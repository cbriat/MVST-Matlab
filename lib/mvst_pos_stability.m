function res = mvst_pos_stability(sys, opts)
%MVST_POS_STABILITY  Copositive linear-programming stability analysis for
%   linear positive systems (Chapter "Linear Positive Systems" of the
%   monograph).  A single dispatcher over the copositive-LP stability
%   certificates: a linear copositive Lyapunov function V(x)=p'x whose
%   existence is a linear program in the nonnegative weight(s) p.
%
%   RES = MVST_POS_STABILITY(SYS, OPTS) selects a certificate family with
%   OPTS.class and returns a uniform result.
%
%   SUPPORTED CLASSES (OPTS.class):
%     'ct'        Metzler LTI drift A.  Implements the linear copositive
%                 Lyapunov certificate of Lemma lem:pos:metzstab(d):
%                 A is Hurwitz iff there is p>0 with A'p<0.
%                 SYS.A (n x n, Metzler).
%     'dt'        nonnegative LTI transition matrix M (Schur).  Discrete-time
%                 analogue of Lemma lem:pos:metzstab: M>=0 is Schur iff there
%                 is p>0 with M'p<p.  SYS.A (n x n, entrywise >= 0).
%     'mjls'      positive Markov jump system.  Implements the coupled
%                 copositive program eq:pos:mjlscop of
%                 Prop. prop:pos:classmjls:
%                 A_i'p_i + sum_j pi_ij p_j < 0, p_i>0, i=1..N.
%                 SYS.A (1xN cell of Metzler drifts), SYS.Pi (N x N generator).
%     'switched'  positive switched system, x' = A_sigma x.
%                 OPTS.dwell = 'arbitrary' (default): common linear copositive
%                   Lyapunov function, A_i'p<0 for all i, p>0 (arbitrary
%                   switching; a common p exists only when no obstruction).
%                 OPTS.dwell = 'min': the mode-dependent minimum-dwell-time
%                   program of Prop. prop:pos:classswsd (the switched slice of
%                   Cor. cor:pos:dwellLPsuff(b)): weights p_i>0 with
%                   A_i'p_i<0 and, for every admissible transition i->j,
%                   e^{A_i' Tmin} p_j < p_i.  SYS.A (1xN cell of Metzler
%                   drifts), OPTS.Tmin (>0).  OPTS.trans (K x 2 list of
%                   admissible ordered transitions [i j]; default all i~=j).
%     'impulsive' positive impulsive system with constant dwell-time.  When
%                 the Metzler flow A is UNSTABLE and the reset J>=0 is
%                 contractive, implements the maximum-dwell-time copositive LP
%                 of Ex. ex:pos:dwellnum (the instance of Briat:17NAHS):
%                 p>0, A'p>0, (J e^{A Tbar})'p < p.  SYS.A, SYS.J, SYS.Tbar.
%                 When A is Hurwitz, implements instead the constant-dwell-time
%                 program eq:pos:dwellLPper of Cor. cor:pos:dwellLP(b):
%                 A'p<0, J'p - p < 0 (constant weight slice), certifying
%                 stability under constant dwell-time Tbar.
%
%   OPTS FIELDS (all optional unless noted):
%     .class      one of the class tags above (required; inferred when
%                 possible: SYS.Pi -> 'mjls', SYS.J -> 'impulsive', cell
%                 SYS.A -> 'switched', matrix SYS.A -> 'ct').
%     .dwell      'arbitrary' | 'min' (class 'switched').
%     .Tmin       minimum dwell-time (class 'switched', dwell 'min').
%     .Tbar       constant dwell-time (class 'impulsive'; may also be SYS.Tbar).
%     .trans      admissible ordered transitions (class 'switched', 'min').
%     .margin     strict-inequality margin (default 1e-6).
%     .pmin       lower bound on the weight entries (default 1).
%     .threshold  logical; for 'switched'/'min' and 'impulsive' also bisect
%                 the dwell-time at which the copositive LP loses / gains
%                 feasibility and return it in RES.Tstar (default false).
%     .Trange     [lo hi] bisection interval (default [1e-2 10]).
%     .Titers     number of bisection steps (default 60).
%     plus any field consumed by MVST_POS_LP (.solver, .lpalgo, .lpopts).
%
%   RES FIELDS:
%     .stable   logical: the copositive LP is feasible (a certificate exists).
%     .p        the copositive weight: an n x 1 vector, or a 1xN cell of
%               n x 1 vectors for the multi-weight classes ('mjls',
%               'switched'); [] if infeasible.
%     .lp       the LP data actually solved (LINPROG naming; see MVST_POS_LP).
%     .margin   the achieved certificate margin: the smallest slack
%               rhs - lhs over all inequality rows at RES.p (>0 certifies the
%               strict inequalities); [] if infeasible.
%     .Tstar    dwell-time feasibility threshold when OPTS.threshold=true,
%               else [].
%     .class    the resolved class tag.
%
%   Example (Metzler-Hurwitz drift):
%     res = mvst_pos_stability(struct('A',[-1 0;3 -1]), struct('class','ct'));
%     % res.stable == true, res.p > 0, A'*res.p < 0.
%
%   See also MVST_POS_GAIN, MVST_POS_STABILIZE, MVST_POS_SUPPLY, MVST_POS_LP,
%   MVST_POS.

    arguments
        sys  (1,1) struct
        opts (1,1) struct = struct()
    end

    margin = mvst_getf(opts, 'margin', 1e-6);
    pmin   = mvst_getf(opts, 'pmin',   1);
    cls    = resolve_class(sys, opts);

    res = struct('stable', false, 'p', [], 'lp', [], 'margin', [], ...
                 'Tstar', [], 'class', cls);

    switch cls
        case 'ct'
            A = getA(sys);
            [res, lp] = solve_lclf(A.', margin, pmin, opts, res);
            res.lp = lp;

        case 'dt'
            M = getA(sys);
            % M'p < p  <=>  (M' - I) p <= -margin
            [res, lp] = solve_lclf(M.' - eye(size(M,1)), margin, pmin, opts, res);
            res.lp = lp;

        case 'mjls'
            Ac = ascell(getA(sys));
            Pi = mvst_getf(sys, 'Pi', mvst_getf(sys, 'Pig', []));
            if isempty(Pi)
                error('mvst:pos:mjls', ...
                      'Class ''mjls'' requires the generator SYS.Pi (N x N).');
            end
            [res, lp] = solve_mjls(Ac, Pi, margin, pmin, opts, res);
            res.lp = lp;

        case 'switched'
            Ac    = ascell(getA(sys));
            dwell = lower(mvst_getf(opts, 'dwell', 'arbitrary'));
            switch dwell
                case 'arbitrary'
                    [res, lp] = solve_common(Ac, margin, pmin, opts, res);
                    res.lp = lp;
                case 'min'
                    Tmin = mvst_getf(opts, 'Tmin', mvst_getf(sys, 'Tmin', []));
                    if isempty(Tmin)
                        error('mvst:pos:dwell', ...
                              ['Class ''switched'' with dwell ''min'' ', ...
                               'requires OPTS.Tmin (> 0).']);
                    end
                    tr = mvst_getf(opts, 'trans', all_trans(numel(Ac)));
                    [res, lp] = solve_switched_min(Ac, Tmin, tr, ...
                                                   margin, pmin, opts, res);
                    res.lp = lp;
                    if mvst_getf(opts, 'threshold', false)
                        res.Tstar = bisect_min(Ac, tr, margin, pmin, opts);
                    end
                otherwise
                    error('mvst:pos:dwell', ...
                          'OPTS.dwell must be ''arbitrary'' or ''min''.');
            end

        case 'impulsive'
            A = getA(sys);
            J = mvst_getf(sys, 'J', []);
            if isempty(J)
                error('mvst:pos:impulsive', ...
                      'Class ''impulsive'' requires the reset SYS.J (>= 0).');
            end
            Tbar = mvst_getf(opts, 'Tbar', mvst_getf(sys, 'Tbar', []));
            if isempty(Tbar)
                error('mvst:pos:impulsive', ...
                      'Class ''impulsive'' requires the dwell-time Tbar.');
            end
            if max(real(eig(A))) > 0
                [res, lp] = solve_imp_max(A, J, Tbar, margin, pmin, opts, res);
            else
                [res, lp] = solve_imp_const(A, J, Tbar, margin, pmin, opts, res);
            end
            res.lp = lp;
            if mvst_getf(opts, 'threshold', false)
                res.Tstar = bisect_imp(A, J, margin, pmin, opts);
            end

        case 'cone'
            % Stability of dot x = A x on the polyhedral proper cone
            % Omega = cone(V) (V is n x m, m >= n generators): find a Metzler
            % generator H with A V = V H (invariance, Vandergraft:68 /
            % BermanPlemmons:94), then run the copositive Lyapunov LP on H.
            A = getA(sys);
            V = mvst_getf(sys, 'V', []);
            if isempty(V)
                error('mvst:pos:cone', ...
                      'Class ''cone'' requires the generator matrix SYS.V (n x m).');
            end
            [res, lp] = solve_cone(A, V, margin, pmin, opts, res);
            res.lp = lp;

        otherwise
            error('mvst:pos:class', ...
                  ['OPTS.class must be ''ct'', ''dt'', ''mjls'', ', ...
                   '''switched'', ''impulsive'' or ''cone'' (got ''%s'').'], cls);
    end
end

% ===================================================================== %
%  Class inference                                                        %
% ===================================================================== %
function cls = resolve_class(sys, opts)
    cls = mvst_getf(opts, 'class', '');
    if ~isempty(cls); cls = lower(cls); return; end
    if ~isempty(mvst_getf(sys, 'V', []))
        cls = 'cone';
    elseif ~isempty(mvst_getf(sys, 'Pi', mvst_getf(sys, 'Pig', [])))
        cls = 'mjls';
    elseif ~isempty(mvst_getf(sys, 'J', []))
        cls = 'impulsive';
    elseif iscell(mvst_getf(sys, 'A', []))
        cls = 'switched';
    else
        cls = 'ct';
    end
end

% ===================================================================== %
%  LP builders                                                           %
% ===================================================================== %
function [res, lpdata] = solve_lclf(Aub, margin, pmin, opts, res)
% Single-weight copositive LP:  Aub p <= -margin,  p >= pmin.
    n  = size(Aub, 2);
    lp = struct('nvar', n, 'Aub', Aub, 'bub', -margin*ones(size(Aub,1),1), ...
                'lb', pmin*ones(n,1));
    o  = mvst_pos_lp(lp, opts);
    lpdata = o.lp;
    if o.feas
        res.stable = true;
        res.p      = o.x;
        res.margin = min(-Aub*o.x);
    end
end

% ===================================================================== %
%  Cone class: Metzler generator factorization + copositive LP           %
% ===================================================================== %
function [res, lpdata] = solve_cone(A, V, margin, pmin, opts, res)
% Stability of dot x = A x on Omega = cone(V), V in R^{n x m}: first find a
% Metzler H in R^{m x m} with A V = V H (the invariance factorization), then
% certify Hurwitz stability of H by the copositive LP exists p>0 with H'p<0
% (Lemma lem:pos:metzstab(d)).  The pulled-back certificate is V(x)=<p,xi> on
% xi>=0 with x=V xi.  With V=I_n this collapses to the 'ct' branch exactly.
    [ok, H, lpf] = factor_metzler(A, V, opts);
    res.coneInvariant = ok;
    res.V = V;  res.H = [];
    if ~ok
        res.stable = false;  res.margin = [];  lpdata = lpf;
        return;
    end
    res.H = H;
    [res, lpdata] = solve_lclf(H.', margin, pmin, opts, res);   % exists p>0, H'p<0
end

function [ok, H, lpdata] = factor_metzler(A, V, opts)
% Feasibility LP for a Metzler H (m x m) with A V = V H.  Variables are the m^2
% entries of H (column-major).  Equality kron(I_m,V) vec(H) = vec(A V); Metzler
% is H_{ij} >= 0 for i ~= j (the diagonal is free).
    m = size(V, 2);
    Aeq = kron(eye(m), V);
    beq = reshape(A*V, [], 1);
    idx = reshape(1:m*m, m, m);                 % idx(i,j) = position of H(i,j)
    off = idx(~eye(m));                          % off-diagonal positions
    Aub = zeros(numel(off), m*m);
    for r = 1:numel(off); Aub(r, off(r)) = -1; end   % -H_{ij} <= 0
    lp = struct('f', zeros(m*m,1), 'Aub', Aub, 'bub', zeros(numel(off),1), ...
                'Aeq', Aeq, 'beq', beq, 'lb', -inf(m*m,1), 'ub', inf(m*m,1));
    o = mvst_pos_lp(lp, opts);
    lpdata = o.lp;  ok = o.feas;  H = [];
    if ok; H = reshape(o.x, m, m); end
end

function [res, lpdata] = solve_common(Ac, margin, pmin, opts, res)
% Common linear copositive Lyapunov function:  A_i' p <= -margin (all i).
    M = numel(Ac); n = size(Ac{1}, 1);
    Aub = zeros(M*n, n);
    for i = 1:M
        Aub((i-1)*n+(1:n), :) = Ac{i}.';
    end
    lp = struct('nvar', n, 'Aub', Aub, 'bub', -margin*ones(M*n,1), ...
                'lb', pmin*ones(n,1));
    o  = mvst_pos_lp(lp, opts);
    lpdata = o.lp;
    if o.feas
        res.stable = true;
        res.p      = o.x;
        res.margin = min(-Aub*o.x);
    end
end

function [res, lpdata] = solve_mjls(Ac, Pi, margin, pmin, opts, res)
% Coupled copositive program eq:pos:mjlscop:
%   A_i' p_i + sum_j pi_ij p_j <= -margin,  p_i >= pmin,  i = 1..N.
    N = numel(Ac); n = size(Ac{1}, 1); nv = N*n;
    if ~isequal(size(Pi), [N N])
        error('mvst:pos:mjls', ...
              'SYS.Pi must be %d x %d for %d modes.', N, N, N);
    end
    Aub = zeros(N*n, nv);
    for i = 1:N
        rows = (i-1)*n + (1:n);
        Aub(rows, (i-1)*n+(1:n)) = Ac{i}.';           % A_i' p_i
        for j = 1:N
            Aub(rows, (j-1)*n+(1:n)) = Aub(rows, (j-1)*n+(1:n)) ...
                                       + Pi(i,j)*eye(n); % + pi_ij p_j
        end
    end
    lp = struct('nvar', nv, 'Aub', Aub, 'bub', -margin*ones(N*n,1), ...
                'lb', pmin*ones(nv,1));
    o  = mvst_pos_lp(lp, opts);
    lpdata = o.lp;
    if o.feas
        res.stable = true;
        res.p      = splitp(o.x, N, n);
        res.margin = min(-Aub*o.x);
    end
end

function [res, lpdata] = solve_switched_min(Ac, Tmin, tr, margin, pmin, opts, res)
% Mode-dependent minimum-dwell-time program prop:pos:classswsd:
%   A_i' p_i <= -margin  (all i),
%   e^{A_i' Tmin} p_j - p_i <= -margin  for each transition i->j in TR.
    N = numel(Ac); n = size(Ac{1}, 1); nv = N*n;
    K = size(tr, 1);
    Aub = zeros(N*n + K*n, nv);
    for i = 1:N
        Aub((i-1)*n+(1:n), (i-1)*n+(1:n)) = Ac{i}.';   % flow rows
    end
    base = N*n;
    for k = 1:K
        i = tr(k,1); j = tr(k,2);
        rows = base + (k-1)*n + (1:n);
        Aub(rows, (j-1)*n+(1:n)) = expm(Ac{i}.'*Tmin);  % e^{A_i' Tmin} p_j
        Aub(rows, (i-1)*n+(1:n)) = Aub(rows, (i-1)*n+(1:n)) - eye(n); % - p_i
    end
    lp = struct('nvar', nv, 'Aub', Aub, 'bub', -margin*ones(N*n+K*n,1), ...
                'lb', pmin*ones(nv,1));
    o  = mvst_pos_lp(lp, opts);
    lpdata = o.lp;
    if o.feas
        res.stable = true;
        res.p      = splitp(o.x, N, n);
        res.margin = min(-Aub*o.x);
    end
end

function [res, lpdata] = solve_imp_max(A, J, Tbar, margin, pmin, opts, res)
% Maximum-dwell-time copositive LP ex:pos:dwellnum (unstable flow, reset):
%   p >= pmin,  A' p >= margin,  (J e^{A Tbar})' p <= p - margin.
    n = size(A,1);
    Je = (J*expm(A*Tbar)).';
    Aub = [ -A.' ;            % -A'p <= -margin  <=>  A'p >= margin
             Je - eye(n) ];   % (Je^{AT})'p - p <= -margin
    bub = [ -margin*ones(n,1) ; -margin*ones(n,1) ];
    lp = struct('nvar', n, 'Aub', Aub, 'bub', bub, 'lb', pmin*ones(n,1));
    o  = mvst_pos_lp(lp, opts);
    lpdata = o.lp;
    if o.feas
        res.stable = true;
        res.p      = o.x;
        res.margin = min(bub - Aub*o.x);
    end
end

function [res, lpdata] = solve_imp_const(A, J, Tbar, margin, pmin, opts, res)
% Constant-dwell-time slice eq:pos:dwellLPper (Hurwitz flow, constant weight):
%   A' p <= -margin,  (J e^{A Tbar})' p - p <= -margin,  p >= pmin.
    n = size(A,1);
    Je = (J*expm(A*Tbar)).';
    Aub = [ A.' ; Je - eye(n) ];
    bub = -margin*ones(2*n,1);
    lp = struct('nvar', n, 'Aub', Aub, 'bub', bub, 'lb', pmin*ones(n,1));
    o  = mvst_pos_lp(lp, opts);
    lpdata = o.lp;
    if o.feas
        res.stable = true;
        res.p      = o.x;
        res.margin = min(bub - Aub*o.x);
    end
end

% ===================================================================== %
%  Dwell-time feasibility thresholds (bisection on the LP)               %
% ===================================================================== %
function Tstar = bisect_min(Ac, tr, margin, pmin, opts)
% Smallest T at which the switched minimum-dwell-time LP is feasible.
    Tr = mvst_getf(opts, 'Trange', [1e-2 10]);
    Ni = mvst_getf(opts, 'Titers', 60);
    feas = @(T) getfield(solve_switched_min(Ac, T, tr, margin, pmin, opts, ...
                 base_res()), 'stable');                                %#ok<GFLD>
    lo = Tr(1); hi = Tr(2);
    if ~feas(hi); Tstar = NaN; return; end
    for it = 1:Ni
        m = 0.5*(lo+hi);
        if feas(m); hi = m; else; lo = m; end
    end
    Tstar = hi;
end

function Tstar = bisect_imp(A, J, margin, pmin, opts)
% Largest T at which the impulsive maximum-dwell-time LP stays feasible.
    Tr = mvst_getf(opts, 'Trange', [1e-2 10]);
    Ni = mvst_getf(opts, 'Titers', 60);
    feas = @(T) getfield(solve_imp_max(A, J, T, margin, pmin, opts, ...
                 base_res()), 'stable');                                %#ok<GFLD>
    lo = Tr(1); hi = Tr(2);
    if ~feas(lo); Tstar = NaN; return; end
    for it = 1:Ni
        m = 0.5*(lo+hi);
        if feas(m); lo = m; else; hi = m; end
    end
    Tstar = lo;
end

% ===================================================================== %
%  Small utilities                                                       %
% ===================================================================== %
function r = base_res()
    r = struct('stable', false, 'p', [], 'lp', [], 'margin', [], ...
               'Tstar', [], 'class', '');
end

function A = getA(sys)
    A = mvst_getf(sys, 'A', []);
    if isempty(A)
        error('mvst:pos:noA', 'SYS must carry a drift .A.');
    end
end

function P = splitp(x, N, n)
    P = cell(1, N);
    for i = 1:N; P{i} = x((i-1)*n + (1:n)); end
end

function tr = all_trans(N)
    tr = zeros(N*(N-1), 2); k = 0;
    for i = 1:N
        for j = 1:N
            if i ~= j; k = k+1; tr(k,:) = [i j]; end
        end
    end
end

function Ac = ascell(A)
    if iscell(A); Ac = A; else; Ac = {A}; end
end
