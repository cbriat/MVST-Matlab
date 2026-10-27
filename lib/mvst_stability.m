function res = mvst_stability(sys, opts)
%MVST_STABILITY  Exact mean-square stability analysis for the book's system classes.
%   RES = MVST_STABILITY(SYS) decides mean-square (exponential) stability of an
%   autonomous matrix-valued system SYS by an EXACT spectral reading of the
%   relevant lifted second-moment operator (the reduced representation on the
%   symmetric subspace, MVST_LIFTSYM / MVST_SVECBASIS).  No sum-of-squares and
%   no optimisation are used; the routine reproduces, class by class, the exact
%   closed-form tests that appear inline in the example scripts (EX00_STABILITY,
%   EX_STAB_LEVY, ...).
%
%   SYS may be
%     * an MVST.SYSTEM object (normalised through its THIN() bridge),
%     * a thin system struct (MVST_SYSTEM output, has .modes),
%     * a thin mode struct  (MVST_MODE output, has .A), or
%     * a rich struct (SOS-layer, has .A0);
%   structs are normalised through MVST.SYSTEM.FROMSTRUCT(...).THIN() when the
%   +mvst package is on the path, and consumed directly otherwise.  For the
%   analysis SYS is treated as AUTONOMOUS: its drift .A (or .A0) is the flow /
%   closed-loop generator and any input matrix is ignored.  Feed a closed loop
%   (A already including B*K) to certify a design; MVST_STABILIZE does this.
%
%   RES = MVST_STABILITY(SYS,OPTS) accepts an options struct:
%     OPTS.class  force the criterion, one of
%                 'ct' | 'dt' | 'mjls' | 'dtmjls' | 'switched' |
%                 'sampled' | 'impulsive' | 'poisson' | 'levy'.
%                 When omitted the class is inferred from the data (see below).
%                 Discrete-time single-mode systems ('dt') have no dedicated
%                 factory or tag and MUST be requested with OPTS.class='dt'.
%     OPTS.T      dwell-time / sampling period for the dwell-dependent classes
%                 ('switched','sampled','impulsive').  A SCALAR gives a single
%                 decision; a VECTOR produces a stability-vs-dwell curve in
%                 RES.curve.  Defaults to the system's own period (.Tbar / .Td).
%                 For 'switched' the decision is for the SINGLE PERIODIC SCHEDULE
%                 that visits the modes in INDEX ORDER (1,2,...,M) with the SAME
%                 dwell-time T in each mode (the monodromy over one such cycle),
%                 NOT a certificate for arbitrary / worst-case switching.
%
%   CONVENTION.  All margins are SECOND-MOMENT (mean-square) rates: the analysis
%   is performed on the LIFTED second-moment operator, so the reported rate is
%   TWICE the corresponding deterministic (first-moment) exponent.  A negative
%   margin therefore certifies mean-square exponential stability.
%     OPTS.tol    stability margin threshold (default 0): the system is declared
%                 stable when RES.margin < OPTS.tol.
%
%   RES has fields
%     .stable    logical decision (RES.margin < OPTS.tol).  Vector when OPTS.T
%                is a grid, aligned to that grid.
%     .margin    the stability margin: the spectral ABSCISSA of the lifted
%                generator for the continuous-flow classes ('ct','mjls','poisson',
%                'levy'), or the LOG SPECTRAL RADIUS for the discrete / monodromy
%                classes ('dt','dtmjls','sampled','impulsive','switched').  In
%                either convention a NEGATIVE margin means (mean-square) stable.
%     .spectrum  eigenvalues of the relevant lifted operator (the generator for
%                the abscissa classes, the one-period monodromy or the lifted
%                transition for the radius classes).  Empty for a dwell grid.
%     .kind      human-readable criterion used, e.g. 'ct:spectral-abscissa'.
%     .class     resolved canonical class tag (char).
%     .n         state dimension.
%     .T         dwell-time / period used (present for the dwell classes; a
%                vector when OPTS.T is a grid).
%     .curve     present only for a dwell grid: struct with fields .T, .margin,
%                .radius, .stable (all vectors aligned to the grid).
%
%   Examples:
%     % Continuous-time diffusion: mean-square stability from the lifted abscissa
%     S = mvst.System.ct([0 1;-2 -1], zeros(2,0), 'Ac', 0.15*eye(2));
%     r = S.stability();  r.stable, r.margin
%
%     % Sampled-data monodromy versus the sampling period (stability curve)
%     r = mvst_stability(cl, struct('class','sampled','T',linspace(0.05,2,120)));
%     plot(r.curve.T, r.curve.radius)
%
%   See also MVST_STABILIZE, MVST_LIFTSYM, MVST_SVECBASIS, MVST.SYSTEM.

    arguments
        sys
        opts struct = struct()
    end

    tol = getf(opts, 'tol', 0);

    % ---- normalise the input to a thin model (modes + optional coupling) ----
    [T, modes] = normalize(sys);
    kind = resolve_kind(T, modes, opts);
    n = dimof(modes{1});

    res = struct('stable', false, 'margin', [], 'spectrum', [], ...
                 'kind', '', 'class', kind, 'n', n);

    % ---- Perron-Frobenius irreducibility (finite Burnside test) ------------
    % OPTS.pf = 'irreducible' tests whether the LTI operator (eq:stab:lindblad)
    % is irreducible on S^n_{>=0}: by prop:stab:PFirr the operator is irreducible
    % iff its defining matrix family {A_0,A_1,...,A_N} (flow) / {J_0,...,J_N}
    % (jump) has no common invariant subspace, decided by the finite Burnside
    % span test of prop:stab:PFburnside (I + products of length <= n^2 span
    % C^{n x n}).
    if ~isempty(getf(opts, 'pf', ''))
        fam = [{modes{1}.A}, mode_channels(modes{1})];
        [irr, algdim] = pf_irreducible(fam, n);
        res.irreducible = irr;
        res.pf = struct('irreducible', irr, 'algdim', algdim, 'full', n^2, ...
                        'family', {fam});
        res.kind = [kind ':pf-irreducible-burnside'];
        return
    end

    switch kind
        % ---------------- continuous-flow classes (abscissa) ----------------
        case {'ct', 'poisson', 'levy'}
            m1 = modes{1};
            G  = flowgen(m1);
            res.spectrum = eig(G);
            res.margin   = max(real(res.spectrum));
            res.kind     = [kind ':spectral-abscissa'];

        case 'mjls'                       % continuous-time coupled generator
            G = mjls_generator(modes, coupling_pi(T));
            res.spectrum = eig(G);
            res.margin   = max(real(res.spectrum));
            res.kind     = 'mjls:coupled-abscissa';

        % ------------------ discrete-time classes (log-radius) --------------
        case 'dt'
            m1 = modes{1};
            Tl = mvst_liftsym(zeros(n), [{m1.A}, mode_channels(m1)]);
            res.spectrum = eig(Tl);
            res.margin   = log(max(abs(res.spectrum)));
            res.kind     = 'dt:log-spectral-radius';

        case 'dtmjls'
            if ~isfield(T, 'P') || isempty(T.P)
                error('mvst:stability:coupling', ...
                    'A discrete-time MJLS requires the DT transition matrix .P.');
            end
            G = dtmjls_transition(modes, T.P);
            res.spectrum = eig(G);
            res.margin   = log(max(abs(res.spectrum)));
            res.kind     = 'dtmjls:coupled-log-radius';

        % --------------- dwell-dependent monodromy classes ------------------
        case {'sampled', 'impulsive', 'switched'}
            res = dwell_analysis(res, kind, T, modes, opts, tol);
            return

        otherwise
            error('mvst:stability:class', 'Unknown stability class "%s".', kind);
    end

    res.stable = res.margin < tol;
end

% ----------------------------------------------------------------------- %
% Finite Burnside irreducibility test (prop:stab:PFirr + prop:stab:PFburnside).
% The unital algebra generated by the family FAM (span of I and all products of
% the generators) equals C^{n x n} iff the family has no common nontrivial
% invariant subspace, i.e. iff the operator is irreducible.  The span is built
% by closing {I} under left-multiplication by the generators; the dimension is
% bounded by n^2, so the closure terminates.
function [irr, d] = pf_irreducible(fam, n)
    fam = fam(~cellfun(@(M) isempty(M) || nnz(M)==0, fam));   % drop empty/zero
    tol = 1e-9;
    basis = {eye(n)};  V = reshape(eye(n), [], 1);
    changed = true;
    while changed
        changed = false;
        snap = basis;                       % snapshot of the current basis
        for a = 1:numel(snap)
            for i = 1:numel(fam)
                Pnew = fam{i} * snap{a};
                v = reshape(Pnew, [], 1);
                if rank([V, v], tol) > size(V, 2)
                    V(:, end+1) = v;  basis{end+1} = Pnew;  changed = true; %#ok<AGROW>
                    if size(V, 2) >= n^2; changed = false; break; end
                end
            end
            if size(V, 2) >= n^2; break; end
        end
    end
    d = size(V, 2);  irr = (d == n^2);
end

% ======================================================================== %
%  Dwell-dependent monodromy classes                                        %
% ======================================================================== %
function res = dwell_analysis(res, kind, T, modes, opts, tol)
%DWELL_ANALYSIS  Monodromy radius over one period for the dwell classes.
    Tgrid = getf(opts, 'T', []);
    if isempty(Tgrid); Tgrid = default_period(kind, modes); end
    if isempty(Tgrid)
        error('mvst:stability:dwell', ...
            'Class "%s" needs a dwell-time / period via OPTS.T or the system.', kind);
    end

    mono = @(Td) monodromy(kind, T, modes, Td);
    res.kind = [kind ':monodromy-log-radius'];

    if isscalar(Tgrid)
        Md = mono(Tgrid);
        res.spectrum = eig(Md);
        res.margin   = log(max(abs(res.spectrum)));
        res.stable   = res.margin < tol;
        res.T        = Tgrid;
    else
        rad = zeros(size(Tgrid));
        for k = 1:numel(Tgrid)
            rad(k) = max(abs(eig(mono(Tgrid(k)))));
        end
        margin = log(rad);
        res.margin   = margin;
        res.stable   = margin < tol;
        res.spectrum = [];               % per-point spectra are not stored
        res.T        = Tgrid;
        res.curve    = struct('T', Tgrid, 'margin', margin, ...
                              'radius', rad, 'stable', margin < tol);
    end
end

function Md = monodromy(kind, T, modes, Td)
%MONODROMY  One-period lifted monodromy matrix at dwell-time TD.
    if strcmp(kind, 'switched')
        Jsw = coupling_reset(T, modes{1});          % common switch reset (or [])
        Md = eye(lift_dim(modes{1}));
        for i = 1:numel(modes)
            Md = expm(Td * flowgen(modes{i})) * Md;
            if ~isempty(Jsw); Md = mvst_liftsym(zeros(dimof(modes{i})), {Jsw}) * Md; end
        end
    else                                            % sampled / impulsive
        m1 = modes{1};
        J  = mode_reset(m1);
        if isempty(J)
            error('mvst:stability:dwell', ...
                'A %s system requires a reset map .J.', kind);
        end
        Jl = mvst_liftsym(zeros(dimof(m1)), {J});
        Md = Jl * expm(Td * flowgen(m1));
    end
end

% ======================================================================== %
%  Lifted-operator builders                                                 %
% ======================================================================== %
function G = flowgen(md)
%FLOWGEN  Continuous lifted generator He(A X) + sum C X C' + Poisson term.
%   Covers the continuous-flow classes uniformly: 'ct' (no jumps), 'poisson'
%   (spontaneous jumps only), 'levy' (diffusion channels and jumps together).
    A0 = drift(md);
    n  = size(A0, 1);
    G  = mvst_liftsym(A0, mode_channels(md));
    nr = size(G, 1);
    [Jf, lamf] = mode_pjumps(md);
    for r = 1:numel(lamf)
        G = G + lamf(r) * (mvst_liftsym(zeros(n), {Jf(:, :, r)}) - eye(nr));
    end
end

function G = mjls_generator(modes, Pi)
%MJLS_GENERATOR  Continuous-time coupled MJLS second-moment generator.
%   Block (i,j) = Pi(j,i)*I (coefficient of X_j in Xdot_i), plus the per-mode
%   lifted generator on the diagonal (see EX00_STABILITY, panel (a)).
    M  = numel(modes);
    nr = lift_dim(modes{1});
    G  = zeros(M * nr);
    for i = 1:M
        Li = flowgen(modes{i});
        for j = 1:M
            blk = Pi(j, i) * eye(nr);
            if i == j; blk = blk + Li; end
            G((i-1)*nr+1:i*nr, (j-1)*nr+1:j*nr) = blk;
        end
    end
end

function G = dtmjls_transition(modes, P)
%DTMJLS_TRANSITION  Discrete-time coupled MJLS second-moment transition.
%   X_i^+ = sum_j p_{ji} T_j(X_j), so block (i,j) = P(j,i)*T_j with T_j the
%   lifted per-mode transition sum kron over the mode's channels (drift and
%   noise alike; see EX_STABZ_DTMJLS).
    M  = numel(modes);
    n  = dimof(modes{1});
    nr = lift_dim(modes{1});
    Tj = cell(1, M);
    for j = 1:M
        Tj{j} = mvst_liftsym(zeros(n), [{modes{j}.A}, mode_channels(modes{j})]);
    end
    G = zeros(M * nr);
    for i = 1:M
        for j = 1:M
            G((i-1)*nr+1:i*nr, (j-1)*nr+1:j*nr) = P(j, i) * Tj{j};
        end
    end
end

% ======================================================================== %
%  Input normalisation and class resolution                                 %
% ======================================================================== %
function [T, modes] = normalize(sys)
%NORMALIZE  Reduce any accepted input to a thin model (system/mode) struct.
    if isa(sys, 'mvst.System')
        T = sys.thin();
    elseif isstruct(sys)
        if exist('mvst.System', 'class') == 8
            T = mvst.System.fromstruct(sys).thin();
        else
            T = sys;
        end
    else
        error('mvst:stability:input', ...
            'SYS must be an mvst.System or a thin/rich struct.');
    end
    if isfield(T, 'modes'); modes = T.modes; else; modes = {T}; end
end

function kind = resolve_kind(T, modes, opts)
%RESOLVE_KIND  Canonical class tag from OPTS.class or the model data.
    allowed = {'ct','dt','mjls','dtmjls','switched','sampled', ...
               'impulsive','poisson','levy'};
    if isfield(opts, 'class') && ~isempty(opts.class)
        kind = lower(char(opts.class));
        if ~ismember(kind, allowed)
            error('mvst:stability:class', 'Unknown OPTS.class "%s".', kind);
        end
        return
    end
    m1 = modes{1};
    M  = numel(modes);
    hasP    = isfield(T, 'P')  && ~isempty(T.P);
    hasPi   = isfield(T, 'Pi') && any(abs(T.Pi(:)) > 0);
    hasAc   = ~isempty(mode_channels(m1));
    [~, lamf] = mode_pjumps(m1);
    hasJf   = ~isempty(lamf);
    hasTbar = isfield(m1, 'Tbar') && ~isempty(m1.Tbar);
    hasJ    = ~isempty(mode_reset(m1));
    if M > 1
        if     hasP;  kind = 'dtmjls';
        elseif hasPi; kind = 'mjls';
        else;         kind = 'switched';
        end
    else
        % A single-mode system may still carry a COMMON reset at the SYSTEM
        % level (T.J, e.g. a one-mode switched-tagged system).  MODE_RESET only
        % reads the mode-level .J, so such a reset would be silently dropped and
        % the system analyzed as 'ct'.  Refuse honestly: the reset defines an
        % impulsive/dwell behaviour that needs a dwell-time or an explicit class.
        hasSysReset = isfield(T, 'J') && ~isempty(T.J);
        if     hasJf && hasAc; kind = 'levy';
        elseif hasJf;          kind = 'poisson';
        elseif hasTbar;        kind = 'sampled';
        elseif hasJ;           kind = 'impulsive';
        elseif hasSysReset
            error('mvst:stability:resetDropped', ...
                ['This single-mode system carries a common reset .J but no ', ...
                 'dwell-time / period, so the reset would be silently dropped ', ...
                 'and the flow analyzed as continuous-time.  Supply the ', ...
                 'dwell-time via OPTS.T (the system is then analyzed as ', ...
                 '''impulsive''), or set OPTS.class explicitly.']);
        else;                  kind = 'ct';
        end
    end
end

% ======================================================================== %
%  Small data accessors (thin or rich mode struct)                          %
% ======================================================================== %
function A = drift(md)
    if     isfield(md, 'A');  A = md.A;
    elseif isfield(md, 'A0'); A = md.A0;
    else;  error('mvst:stability:input', 'Mode carries no drift (.A / .A0).');
    end
end

function n = dimof(md)
    n = size(drift(md), 1);
end

function d = lift_dim(md)
    n = dimof(md); d = n * (n + 1) / 2;
end

function C = mode_channels(md)
%MODE_CHANNELS  State-noise channels of a mode as a 1xq cell (possibly empty).
%   A bare struct WITHOUT an .Ac field (or with an empty .Ac) has NO channels:
%   guard on isfield && ~isempty so we return {} rather than fabricating a
%   spurious empty channel (size([],3) is 1, which would inject a {[]} channel
%   and mislead the class resolution / lift builders).
    if isfield(md, 'Ac') && ~isempty(md.Ac)
        Ac = md.Ac;
    else
        Ac = zeros(dimof(md), dimof(md), 0);
    end
    C = cell(1, size(Ac, 3));
    for l = 1:size(Ac, 3)
        C{l} = Ac(:, :, l);
    end
end

function [Jf, lamf] = mode_pjumps(md)
%MODE_PJUMPS  Spontaneous (Poisson) jump maps and rates of a mode.
    if isfield(md, 'Jf') && ~isempty(md.Jf) && ...
       isfield(md, 'lamf') && ~isempty(md.lamf)
        Jf = md.Jf; lamf = md.lamf;
    else
        Jf = zeros(dimof(md), dimof(md), 0); lamf = [];
    end
end

function J = mode_reset(md)
%MODE_RESET  Scheduled state reset of a mode (first slice), or [].
    if isfield(md, 'J') && ~isempty(md.J); J = md.J(:, :, 1); else; J = []; end
end

function Pi = coupling_pi(T)
    if isfield(T, 'Pi') && ~isempty(T.Pi); Pi = T.Pi; else
        error('mvst:stability:coupling', 'A continuous MJLS requires the generator .Pi.');
    end
end

function J = coupling_reset(T, m1)
%COUPLING_RESET  Common switch reset for a switched system (system or mode), or [].
    if isfield(T, 'J') && ~isempty(T.J); J = T.J(:, :, 1);
    else; J = mode_reset(m1);
    end
end

function Tp = default_period(~, modes)
%DEFAULT_PERIOD  System-carried dwell / period for a dwell class, or [].
    m1 = modes{1};
    if isfield(m1, 'Tbar') && ~isempty(m1.Tbar); Tp = m1.Tbar;
    elseif isfield(m1, 'Td') && ~isempty(m1.Td); Tp = m1.Td;
    else; Tp = [];
    end
end

function v = getf(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
