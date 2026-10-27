function res = mvst_stabilize(sys, opts)
%MVST_STABILIZE  Exact state-feedback stabilization for the book's system classes.
%   RES = MVST_STABILIZE(SYS) designs a mean-square stabilizing state feedback
%   u = K x for the classes that admit a closed-form / eigen-solver route in the
%   monograph, routing each class to its book Riccati engine, and then CERTIFIES
%   the closed loop by an exact call to MVST_STABILITY.  SYS is accepted in the
%   same forms as MVST_STABILITY (an MVST.SYSTEM object or a thin/rich struct).
%
%   Routing (the engine each class uses, mirroring the *_gains example scripts):
%     'ct'       MVST_CARE_CT            single stochastic algebraic Riccati;
%     'mjls'     MVST_SOLVE_INFINITE     coupled continuous-time MJLS Riccati;
%     'switched' MVST_CARE_CT per mode   per-mode infinite-horizon gains (valid
%                                        for a sufficiently large dwell-time);
%     'sampled'  MVST_IGRE_PERIODIC_SD   periodic sampled-data impulsive I-GRE;
%     'poisson'  MVST_CARE_CT_POISSON    GBM + Poissonian-jump algebraic Riccati.
%
%   The DISCRETE-TIME ('dt','dtmjls') and the DWELL-TIME-CONSTRAINED impulsive
%   designs have no closed-form Riccati route in the monograph: their stabilizing
%   state feedback is the timer-dependent sum-of-squares program of
%   MVST_STABZ_SOS (timer-dependent gains, minimum / range dwell-time).  This
%   dispatcher DELEGATES those classes to MVST_STABZ_SOS and flags the result
%   with RES.method = 'sos' (the closed-form classes leave RES.method unset).
%   Supply an explicit dwell specification through OPTS.dwell to override the
%   inferred one ('discrete' for dt/dtmjls, 'minimum' at the system period, or
%   'arbitrary', for impulsive).
%
%   RES = MVST_STABILIZE(SYS,OPTS) accepts:
%     OPTS.class  force the class (same tags as MVST_STABILITY).
%     OPTS.T      dwell-time / period used to certify the closed loop for the
%                 dwell classes ('switched','sampled'); defaults to the system's
%                 own period.  Forwarded to MVST_STABILITY.
%     other fields are forwarded to the underlying Riccati engine.
%
%   RES has fields
%     .K           the gain(s): an m x n matrix for the single-mode classes, a
%                  1 x M cell of m x n gains for the multi-mode classes.
%     .P           the Riccati solution(s) (matrix or 1 x M cell), when available.
%     .closedloop  the autonomous closed-loop system struct fed to the check
%                  (drift A already including B*K, channels folded likewise).
%     .check       the MVST_STABILITY result on .closedloop; .check.stable must
%                  be true for a successful design.
%     .class       resolved canonical class tag.
%     .info        engine diagnostics (residual / iterations) where reported.
%
%   Examples:
%     S = mvst.System.ct([0 1;-1 0.2], [0;1], 'Q', eye(2), 'R', 1, ...
%                        'Ac', 0.3*[0 0;1 0], 'Bc', [0;0.2]);
%     r = S.stabilize();  r.K, r.check.stable
%
%     Pi = [-1 1; 2 -2];
%     S  = mvst.System.mjls({[0 1;-1 0.2],[0 1;-2 -0.1]}, {[0;1],[0;1.2]}, Pi, ...
%                           'Q', eye(2), 'R', 1, 'Ac', 0.25*[0 0;1 0], 'Bc', [0;0.15]);
%     r = S.stabilize();  r.check.stable
%
%   See also MVST_STABILITY, MVST_STABZ_SOS, MVST_CARE_CT, MVST_SOLVE_INFINITE,
%   MVST_IGRE_PERIODIC_SD, MVST_CARE_CT_POISSON, MVST.SYSTEM.

    arguments
        sys
        opts struct = struct()
    end

    [T, modes] = normalize(sys);
    kind = resolve_kind(T, modes, opts);
    res  = struct('K', [], 'P', [], 'closedloop', [], 'check', [], ...
                  'class', kind, 'info', struct());

    switch kind
        case 'ct'
            res = ct_route(res, modes{1}, opts);
        case 'mjls'
            res = mjls_route(res, T, modes, opts);
        case 'switched'
            res = switched_route(res, T, modes, opts);
        case 'sampled'
            res = sampled_route(res, modes{1}, opts);
        case {'poisson', 'levy'}
            res = poisson_route(res, modes{1});
        case {'dt', 'dtmjls', 'impulsive'}
            res = sos_route(res, T, modes, kind, opts);
        otherwise
            error('mvst:stabilize:class', 'Unknown stabilization class "%s".', kind);
    end
end

% ======================================================================== %
%  Per-class routes                                                          %
% ======================================================================== %
function res = ct_route(res, md, opts)
%CT_ROUTE  Single continuous-time stochastic ARE via MVST_CARE_CT.
    [P, K, info] = mvst_care_ct(md, opts);
    res.K = K; res.P = P; res.info = info;
    Acl = md.A + md.B * K;
    Cc  = fold_channels(md, K);
    res.closedloop = auto_ct(Acl, Cc);
    res.check = mvst_stability(res.closedloop, struct('class', 'ct'));
end

function res = mjls_route(res, T, modes, opts)
%MJLS_ROUTE  Coupled continuous-time MJLS gains via MVST_SOLVE_INFINITE.
    sys = rebuild_system(T, modes);
    out = mvst_solve_infinite(sys, opts);
    res.K = out.K; res.P = out.P;
    res.info = struct('residual', out.residual, 'iters', out.iters, 'cost', out.cost);
    cl = cell(1, numel(modes));
    for i = 1:numel(modes)
        Acl = modes{i}.A + modes{i}.B * out.K{i};
        cl{i} = auto_ct(Acl, fold_channels(modes{i}, out.K{i}));
    end
    res.closedloop = mvst_system(cell_modes(cl), 'Pi', T.Pi);
    res.check = mvst_stability(res.closedloop, struct('class', 'mjls'));
end

function res = switched_route(res, T, modes, opts)
%SWITCHED_ROUTE  Per-mode infinite-horizon gains via MVST_CARE_CT.
    M = numel(modes); K = cell(1, M); P = cell(1, M); cl = cell(1, M);
    for i = 1:M
        [P{i}, K{i}] = mvst_care_ct(modes{i}, opts);
        Acl = modes{i}.A + modes{i}.B * K{i};
        cl{i} = auto_ct(Acl, fold_channels(modes{i}, K{i}));
    end
    res.K = K; res.P = P;
    Jc = []; if isfield(T, 'J'); Jc = T.J; end
    args = {}; if ~isempty(Jc); args = {'J', Jc}; end
    res.closedloop = mvst_system(cell_modes(cl), args{:});
    copts = struct('class', 'switched');
    Tp = getf(opts, 'T', default_period(modes));
    if ~isempty(Tp); copts.T = Tp; end
    res.check = mvst_stability(res.closedloop, copts);
end

function res = sampled_route(res, md, opts)
%SAMPLED_ROUTE  Periodic sampled-data impulsive I-GRE via MVST_IGRE_PERIODIC_SD.
    n = size(md.A, 1); m = size(md.B, 2); Tbar = md.Tbar;
    Ac = getf(md, 'Ac', zeros(n, n, 0));
    Bc = getf(md, 'Bc', zeros(n, m, 0));
    out = mvst_igre_periodic_sd(md.A, md.B, Ac, Bc, md.Q, md.R, Tbar, opts);
    K = out.K; res.K = K; res.P = out.P0;
    res.info = struct('residual', out.residual, 'iters', out.iters);
    % augmented autonomous closed loop z = [x; u] (flow holds u; jump sets u=Kx)
    na  = n + m;
    Ab0 = [md.A md.B; zeros(m, na)];
    Nc  = size(Ac, 3); Abl = zeros(na, na, max(Nc, 1));
    if Nc == 0
        Abl = zeros(na, na, 0);
    else
        for l = 1:Nc; Abl(:, :, l) = [Ac(:, :, l) Bc(:, :, l); zeros(m, na)]; end
    end
    Jb = [eye(n) zeros(n, m); K zeros(m, m)];
    cl.A = Ab0; cl.Ac = Abl; cl.J = Jb; cl.Tbar = Tbar;
    res.closedloop = cl;
    res.check = mvst_stability(cl, struct('class', 'sampled', 'T', Tbar));
end

function res = poisson_route(res, md)
%POISSON_ROUTE  GBM + Poissonian-jump algebraic Riccati via MVST_CARE_CT_POISSON.
    Cc = channels_cell(getf(md, 'Ac', []));
    Dc = channels_cell(getf(md, 'Bc', []));
    if numel(Dc) < numel(Cc)                     % pad input-noise channels
        Dc(end+1:numel(Cc)) = {zeros(size(md.B))};
    end
    [Jf, lamf] = pjumps(md);
    Jc = channels_cell(Jf);
    [P, K, info] = mvst_care_ct_poisson(md.A, md.B, Cc, Dc, Jc, lamf, md.Q, md.R);
    res.K = K; res.P = P; res.info = info;
    Acl = md.A + md.B * K;
    Cccl = cell(1, numel(Cc));
    for l = 1:numel(Cc); Cccl{l} = Cc{l} + Dc{l} * K; end
    cl.A  = Acl;
    cl.Ac = cat3(Cccl, size(Acl, 1));
    cl.Jf = Jf; cl.lamf = lamf;
    res.closedloop = cl;
    res.check = mvst_stability(cl, struct('class', 'poisson'));
end

% ======================================================================== %
%  Delegated SOS route (discrete / dwell-time-constrained synthesis)         %
% ======================================================================== %
function res = sos_route(res, T, modes, kind, opts)
%SOS_ROUTE  Delegate to MVST_STABZ_SOS for the classes with no Riccati route.
%   The discrete-time ('dt','dtmjls') and impulsive (dwell-time-constrained)
%   classes are stabilized by the timer-dependent sum-of-squares program of
%   MVST_STABZ_SOS.  Builds the rich struct and a dwell specification, runs that
%   engine, and reports RES.method = 'sos' alongside the timer-gain fields.
    if isfield(T, 'modes'); rich = mvst_torich(T); else; rich = mvst_torich(modes{1}); end
    dwell = getf(opts, 'dwell', []);
    if isempty(dwell)
        switch kind
            case {'dt', 'dtmjls'}
                dwell = struct('type', 'discrete');
            otherwise                                   % impulsive
                Tp = getf(opts, 'T', default_period(modes));
                if ~isempty(Tp)
                    dwell = struct('type', 'minimum', 'Tbar', Tp);
                else
                    dwell = struct('type', 'arbitrary');
                end
        end
    end
    out = mvst_stabz_sos(rich, dwell, opts);
    res.method = 'sos';
    res.feas   = out.feas;
    res.Kfun   = out.Kfun;
    if isfield(out, 'Kf'); res.Kf = out.Kf; end
    if isfield(out, 'KJ'); res.KJ = out.KJ; end
    res.sos    = out;
    if isfield(out, 'info'); res.info = out.info; end
    try; res.K = out.Kfun(0); catch; res.K = []; end   % timer gain at tau = 0
end

% ======================================================================== %
%  Closed-loop assembly helpers                                             %
% ======================================================================== %
function Cc = fold_channels(md, K)
%FOLD_CHANNELS  Closed-loop state-noise channels Ac_l + Bc_l*K (as an array).
    n = size(md.A, 1);
    Ac = getf(md, 'Ac', zeros(n, n, 0));
    Bc = getf(md, 'Bc', []);
    q = size(Ac, 3);
    Cc = zeros(n, n, q);
    for l = 1:q
        Cc(:, :, l) = Ac(:, :, l);
        if ~isempty(Bc) && size(Bc, 3) >= l
            Cc(:, :, l) = Cc(:, :, l) + Bc(:, :, l) * K;
        end
    end
end

function md = auto_ct(A, Cc)
%AUTO_CT  Autonomous thin mode with drift A and state-noise channels Cc.
%   A phantom single input (all zeros, drift and noise alike) is attached so the
%   MVST_MODE channel-count check is satisfied for any number of channels; the
%   input is ignored by the autonomous stability reading, which uses only .A and
%   the state-noise channels .Ac.
    n = size(A, 1); q = size(Cc, 3);
    md = mvst_mode(A, zeros(n, 1), zeros(n), 1, ...
                   'Ac', Cc, 'Bc', zeros(n, 1, q));
end

function modes = cell_modes(cl)
%CELL_MODES  Pass through a 1xM cell of mode structs (identity; for clarity).
    modes = cl;
end

function A3 = cat3(C, n)
%CAT3  Stack a 1xq cell of n x n matrices into an n x n x q array.
    q = numel(C); A3 = zeros(n, n, q);
    for l = 1:q; A3(:, :, l) = C{l}; end
end

function C = channels_cell(Ac)
%CHANNELS_CELL  n x n x q array -> 1xq cell (empty for empty input).
    C = {};
    if isempty(Ac); return; end
    for l = 1:size(Ac, 3); C{l} = Ac(:, :, l); end %#ok<AGROW>
end

function sys = rebuild_system(T, modes)
%REBUILD_SYSTEM  Reassemble a thin multi-mode system from its parts.
    args = {};
    if isfield(T, 'Pi');  args = [args, {'Pi',  T.Pi}];  end
    if isfield(T, 'pi0'); args = [args, {'pi0', T.pi0}]; end
    if isfield(T, 'ZT');  args = [args, {'ZT',  T.ZT}];  end
    sys = mvst_system(modes, args{:});
end

% ======================================================================== %
%  Input normalisation and class resolution (shared with MVST_STABILITY)    %
% ======================================================================== %
function [T, modes] = normalize(sys)
    if isa(sys, 'mvst.System')
        T = sys.thin();
    elseif isstruct(sys)
        if exist('mvst.System', 'class') == 8
            T = mvst.System.fromstruct(sys).thin();
        else
            T = sys;
        end
    else
        error('mvst:stabilize:input', ...
            'SYS must be an mvst.System or a thin/rich struct.');
    end
    if isfield(T, 'modes'); modes = T.modes; else; modes = {T}; end
end

function kind = resolve_kind(T, modes, opts)
    allowed = {'ct','dt','mjls','dtmjls','switched','sampled', ...
               'impulsive','poisson','levy'};
    if isfield(opts, 'class') && ~isempty(opts.class)
        kind = lower(char(opts.class));
        if ~ismember(kind, allowed)
            error('mvst:stabilize:class', 'Unknown OPTS.class "%s".', kind);
        end
        return
    end
    m1 = modes{1};
    M  = numel(modes);
    hasP    = isfield(T, 'P')  && ~isempty(T.P);
    hasPi   = isfield(T, 'Pi') && any(abs(T.Pi(:)) > 0);
    hasAc   = isfield(m1, 'Ac') && size(m1.Ac, 3) > 0;
    [~, lamf] = pjumps(m1);
    hasJf   = ~isempty(lamf);
    hasTbar = isfield(m1, 'Tbar') && ~isempty(m1.Tbar);
    hasJ    = isfield(m1, 'J') && ~isempty(m1.J);
    if M > 1
        if     hasP;  kind = 'dtmjls';
        elseif hasPi; kind = 'mjls';
        else;         kind = 'switched';
        end
    else
        if     hasJf && hasAc; kind = 'levy';
        elseif hasJf;          kind = 'poisson';
        elseif hasTbar;        kind = 'sampled';
        elseif hasJ;           kind = 'impulsive';
        else;                  kind = 'ct';
        end
    end
end

function [Jf, lamf] = pjumps(md)
    if isfield(md, 'Jf') && ~isempty(md.Jf) && ...
       isfield(md, 'lamf') && ~isempty(md.lamf)
        Jf = md.Jf; lamf = md.lamf;
    else
        n = size(md.A, 1); Jf = zeros(n, n, 0); lamf = [];
    end
end

function Tp = default_period(modes)
    m1 = modes{1};
    if isfield(m1, 'Tbar') && ~isempty(m1.Tbar); Tp = m1.Tbar;
    elseif isfield(m1, 'Td') && ~isempty(m1.Td); Tp = m1.Td;
    else; Tp = [];
    end
end

function v = getf(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
