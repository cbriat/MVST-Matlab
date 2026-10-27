function res = mvst_stab_sos_modes(sys, dwell, opts)
%MVST_STAB_SOS_MODES  Generic SOSTOOLS stability-ANALYSIS test on a PRODUCT CONE
%  (S^n)^M for mode-indexed linear Hermitian matrix-valued systems: the coupled
%  certificates {P_i} of the Stability chapter's Markov-jump and switched
%  results.  This is the product-cone counterpart of MVST_STAB_SOS (single
%  certificate): here every mode carries its own Lyapunov certificate P_i and
%  the certificates are COUPLED through the mode transition data (Markov
%  generator / transition matrix, or the switch conditions).  As in
%  MVST_STAB_SOS there are NO gain variables and NO inputs: this is a
%  feasibility test, decided by common-slack maximization.
%
%  SYS fields (any omitted field defaults to "none"; NO input/B fields):
%    .A0  n x n x M          per-mode flow drift (CT classes / switched) or
%                            per-mode transition (DT class)
%    .Ac  n x n x q x M      per-mode flow noise channels (Brownian: pass the
%                            state-multiplicative maps; DT class: the DT noise
%                            channels M_l = Ac_l).  q may be 0.
%    .Pi  M x M              Markov GENERATOR (CT-MJLS, MSI); rows sum to 0,
%                            Pi(i,j) = rate i -> j.
%    .P   M x M              Markov TRANSITION matrix (DT-MJLS); rows sum to 1.
%    .J   n x n              impulse reset applied at every switch (MSI) or the
%                            switch reset (switched; default identity).
%
%  DWELL:
%    .type  'ctmjls' | 'dtmjls' | 'msi' | 'switched'
%       ('ctmjls'  = CT Markov jump diffusion, constant {P_i}, generator coupling
%                    sum_j Pi(i,j) P_j (row-i coupling), certified feasible iff
%                    the coupled lifted generator has spectral abscissa < 0.
%        'dtmjls'  = DT Markov jump, constant {P_i}, coupled transition
%                    X_i^+ = sum_j P(j,i) T_j(X_j); the certificate side is
%                    sum_i P(j,i) T_j^*(P_i) - rho P_j <= -slack I (row-j
%                    coupling), feasible iff the coupled lifted radius < 1.
%        'msi'     = CT Markov switching with an impulse RESET J at every switch:
%                    the off-diagonal coupling carries J^* (.) J.
%        'switched' = minimum-dwell-time SWITCHED analysis, mode-indexed timer
%                    certificates P_i(theta), theta = tau/Tbar in [0,1], with the
%                    all-pairs switch conditions P_j(0) - rho P_i(1) <= -slack I
%                    (reset J on the switch; identity by default).  This is the
%                    mode-indexed timer test that a single-certificate program
%                    cannot express; feasible iff Tbar exceeds the minimum
%                    dwell-time, cross-checked by the exact switching monodromy.)
%    .Tbar  scalar (switched)                 dwell interval length.
%    .deg   timer polynomial degree of P_i (switched; default 4).
%    .alpha,.rho  guaranteed hybrid rate (defaults 0, 1), at the COVARIANCE
%           level: the certified state decays as X ~ e^{-alpha t} rho^kappa, so
%           the flow shift is +alpha P and the jump/switch contraction is rho
%           (NOT 2 alpha / rho^2); reported rates are the covariance rates
%           directly, with no halving or rooting.
%
%  Conditions certified (P-side, no congruence), per mode:
%    ctmjls  A_i^*P_i + P_i A_i + sum_l Ac_{il}^*P_i Ac_{il} + aP_i
%               + sum_j Pi(i,j) P_j <= -slack I,  0 <= P_i <= I.
%    dtmjls  sum_j P(i,j) ( A_i^*P_j A_i + sum_l Ac_{il}^*P_j Ac_{il} )
%               - rho P_i <= -slack I,  0 <= P_i <= I.
%               (row-i coupling: P(i,j) is the outgoing transition i -> j; this
%                is the destination sum for the standard DT coupled Lyapunov
%                inequality, equivalent to the source-sum primal above.)
%    msi     A_i^*P_i + P_i A_i + sum_l Ac_{il}^*P_i Ac_{il} + aP_i
%               + Pi(i,i) P_i + sum_{j!=i} Pi(i,j) J^*P_j J <= -slack I.
%    switched flow (per mode i, on [0,1] via theta(1-theta) multiplier):
%               (1/Tbar)P_i' + A_i^*P_i + P_i A_i + sum_l Ac_{il}^*P_i Ac_{il}
%               + aP_i <= -slack I;   bounds 0 <= P_i <= I on [0,1];
%             switch (all ordered pairs i -> j != i):
%               J^*P_j(0)J - rho P_i(1) <= -slack I.
%
%  RES: .feas, .slack, .type, .Pi (cell of certificates), .Pfun (switched:
%       cell of @(tau)->P_i), .info.
%
%  OPTS (struct or MVST.OPTIONS, optional): .solver selects the SDP backend
%  ('sedumi' (default) | 'sdpt3' | 'mosek' | 'nosdp'), normalized by MVST_OPTS.
%
%  SYS may also be a THIN multi-mode system (MVST_SYSTEM): if it lacks the
%  .A0 field it is routed through MVST_TORICH, which stacks the modes into
%  A0 (n x n x M) / Ac (n x n x q x M) and exposes the generator .Pi.

if nargin < 3; opts = struct(); end
oo = mvst_opts(opts); solver = oo;   % normalized options travel to the solve (reltol/verbosity)
if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
n = size(sys.A0,1);  M = size(sys.A0,3);
sys = fill_defaults(sys, n, M);
q = size(sys.Ac,3);
alpha = getf(dwell,'alpha',0);  rho = getf(dwell,'rho',1);  tol = 1e-7;

switch dwell.type
    case {'ctmjls','dtmjls','msi'}
        res = solve_constant(sys, dwell, n, M, q, alpha, rho, tol, solver);
    case 'switched'
        res = solve_switched(sys, dwell, n, M, q, alpha, rho, tol, solver);
    otherwise
        error('mvst:stab_sos_modes:type','unknown dwell type "%s"', dwell.type);
end
end

% ==========================================================================
function res = solve_constant(sys, dwell, n, M, q, alpha, rho, tol, solver)
% Constant per-mode certificates {P_i} for the CT-MJLS, DT-MJLS and MSI classes.
pvar th; sl = dpvar('sl');            %#ok<NASGU> th is a formal placeholder
prog = sosprogram(th);  prog = sosdecvar(prog, sl);
P = cell(1,M);
for i = 1:M
    [prog, P{i}] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
end

for i = 1:M
    L = mode_flow(sys, P, i, alpha, q, rho, dwell.type, M);
    prog = matsos(prog, -L - sl*eye(n));
    prog = matsos(prog, P{i});                 % 0 <= P_i
    prog = matsos(prog, eye(n) - P{i});        % P_i <= I (scale fix)
end

prog = sossetobj(prog, -sl);
o = mvst_sosopts(solver);
prog = mvst_sossolve(prog, o);

sval = double(sosgetsol(prog, sl));
Pv = cell(1,M);
for i = 1:M; Pv{i} = double(sosgetsol(prog, P{i})); end
res.type = dwell.type;  res.slack = sval;  res.feas = sval > tol;
res.Pi = Pv;  res.info = struct('P',{Pv},'alpha',alpha,'rho',rho);
end

% --------------------------------------------------------------------------
function L = mode_flow(sys, P, i, alpha, q, rho, type, M)
% Coupled certificate residual for mode i (constant certificates).
A = sys.A0(:,:,i);
switch type
    case 'ctmjls'
        L = A'*P{i} + P{i}*A + alpha*P{i};
        for l = 1:q; C = sys.Ac(:,:,l,i); L = L + C'*P{i}*C; end
        for j = 1:M; L = L + sys.Pi(i,j)*P{j}; end            % row-i coupling
    case 'msi'
        L = A'*P{i} + P{i}*A + alpha*P{i} + sys.Pi(i,i)*P{i};
        for l = 1:q; C = sys.Ac(:,:,l,i); L = L + C'*P{i}*C; end
        for j = 1:M
            if j ~= i; L = L + sys.Pi(i,j)*(sys.J'*P{j}*sys.J); end
        end
    case 'dtmjls'
        % sum_j P(i,j) ( A_i^* P_j A_i + sum_l Ac_il^* P_j Ac_il ) - rho P_i
        L = -rho*P{i};
        for j = 1:M
            Tstar = A'*P{j}*A;
            for l = 1:q; C = sys.Ac(:,:,l,i); Tstar = Tstar + C'*P{j}*C; end
            L = L + sys.P(i,j)*Tstar;
        end
end
end

% ==========================================================================
function res = solve_switched(sys, dwell, n, M, q, alpha, rho, tol, solver)
% Minimum-dwell-time SWITCHED analysis: mode-indexed timer certificates
% P_i(theta) on [0,1] (theta = tau/Tbar) with all-pairs switch conditions.
T = dwell.Tbar;  deg = getf(dwell,'deg',4);
pvar th; sl = dpvar('sl');
prog = sosprogram(th);  prog = sosdecvar(prog, sl);
gt = th*(1-th);
P = cell(1,M);
for i = 1:M
    [prog, P{i}] = sospolymatrixvar(prog, monomials(th,0:deg), [n n], 'symmetric');
end

% flow + positivity + normalization, per mode, on [0,1]
for i = 1:M
    A = sys.A0(:,:,i);
    Phi = (1/T)*diff(P{i},th) + A'*P{i} + P{i}*A + alpha*P{i};
    for l = 1:q; C = sys.Ac(:,:,l,i); Phi = Phi + C'*P{i}*C; end
    [prog,S] = sospolymatrixvar(prog, monomials(th,0:deg), [n n], 'symmetric');
    prog = matsos(prog, -Phi - sl*eye(n) - gt*S);   prog = matsos(prog, S);
    [prog,Lo] = sospolymatrixvar(prog, monomials(th,0:deg), [n n], 'symmetric');
    prog = matsos(prog, P{i} - 1e-3*eye(n) - gt*Lo); prog = matsos(prog, Lo);
    [prog,Hi] = sospolymatrixvar(prog, monomials(th,0:deg), [n n], 'symmetric');
    prog = matsos(prog, eye(n) - P{i} - gt*Hi);      prog = matsos(prog, Hi);
end

% switch conditions: all ordered pairs i -> j != i, reset J
for i = 1:M
    Pi1 = subs(P{i},th,1);
    for j = 1:M
        if j ~= i
            Pj0 = subs(P{j},th,0);
            prog = matsos(prog, rho*Pi1 - sys.J'*Pj0*sys.J - sl*eye(n));
        end
    end
end

prog = sossetobj(prog, -sl);
o = mvst_sosopts(solver);
prog = mvst_sossolve(prog, o);

sval = double(sosgetsol(prog, sl));
Pv = cell(1,M);  Pfun = cell(1,M);
for i = 1:M
    Pv{i} = sosgetsol(prog, P{i});
    Pfun{i} = @(tau) evalmat(Pv{i}, th, tau/T, n);
end
res.type = dwell.type;  res.slack = sval;  res.feas = sval > tol;
res.Pi = Pv;  res.Pfun = Pfun;
res.info = struct('P',{Pv},'T',T,'deg',deg,'alpha',alpha,'rho',rho);
end

% ==========================================================================
function prog = matsos(prog, EXPR); prog = sosmatrixineq(prog, EXPR); end

function Mout = evalmat(P, th, theta, n)
if isnumeric(P); Mout = P; else; Mout = double(subs(P, th, theta)); end
Mout = reshape(Mout, n, n);
end

function sys = fill_defaults(sys, n, M)
if ~isfield(sys,'Ac')||isempty(sys.Ac); sys.Ac = zeros(n,n,0,M); end
if ndims(sys.Ac) < 4 && size(sys.Ac,3) >= 1 && M == 1
    sys.Ac = reshape(sys.Ac, n, n, size(sys.Ac,3), 1);
end
if ~isfield(sys,'J')||isempty(sys.J);   sys.J = eye(n); end
end

function v = getf(s, f, d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end
