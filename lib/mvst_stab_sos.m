function res = mvst_stab_sos(sys, dwell, opts)
%MVST_STAB_SOS  Generic SOSTOOLS stability-ANALYSIS test for linear impulsive
%  Hermitian matrix-valued systems, covering every dwell-time type of the
%  Stability chapter.  This is the ANALYSIS counterpart of MVST_STABZ_SOS
%  (the state-feedback synthesis routine): here there are NO gain variables and
%  NO inputs, only a Lyapunov certificate P (constant, or polynomial in the
%  timer tau).  The routine is a feasibility test for the stability conditions
%  of the Stability chapter (Theorems th:stab:main, th:stab:dwell,
%  th:stab:MDT, th:stab:RDT and the Poisson / discrete-time corollaries).
%
%  Because there is no gain to eliminate, the conditions are written directly on
%  the certificate P (the adjoint / observability side of the KYP inequality),
%  as plain n x n matrix inequalities.  No congruence, no Schur complement:
%  the second-moment noise terms sum_l Ac_l^* P Ac_l and the spontaneous-jump
%  terms sum_r lamf_r (Jf_r^* P Jf_r - P) enter P linearly.
%
%  SYS fields (any omitted field defaults to "none"; NO input/B fields):
%    .A0   n x n         flow drift
%    .Ac   n x n x q     flow noise channels (Brownian OR Levy: pass
%                        sqrt(c_l) A_l, and fold the Levy mean drift
%                        sum_l mu_l A_l into A0)
%    .Jf   n x n x pf    SPONTANEOUS (Poisson) jump maps, rate-weighted into the
%                        flow through .lamf
%    .lamf 1 x pf        Poisson rates of the spontaneous jumps
%    .J    n x n x p     SCHEDULED (dwell-time) jump maps (p>1 => noisy jumps,
%                        summed as sum_r J_r^* P(0) J_r)
%
%  DWELL:
%    .type  'arbitrary' | 'discrete' | 'periodic' | 'minimum' | 'range'
%           ('arbitrary'  = constant certificate (CT persistent flowing); a
%                           purely spontaneous-Poisson system is 'arbitrary'
%                           with .Jf/.lamf and no scheduled .J.
%            'discrete'    = discrete-time / persistent-jumping test, condition
%                           sum_l M_l^* P M_l - rho P <= -slack I with the
%                           transition M_0 = A0 and channels M_l = Ac_l.
%            'periodic'    = fixed dwell-time Tbar, single jump at tau = Tbar.
%            'minimum'     = minimum dwell-time Tbar (ramp on [0,Tbar] plus the
%                           frozen tail at tau = Tbar).
%            'range'       = ranged dwell-time [Tmin,Tmax], jump allowed for any
%                           tau in [Tmin,Tmax] via the dwell multiplier.)
%    .Tbar  scalar (periodic, minimum) ; .Tmin,.Tmax (range)
%    .deg   timer polynomial degree of P (default 2; forced 0 for
%           arbitrary/discrete)
%    .alpha,.rho  guaranteed hybrid rate (defaults 0, 1), at the COVARIANCE
%           level: the certified state decays as X ~ e^{-alpha t} rho^kappa
%           (kappa = number of jumps), so the flow shift is +alpha P and the
%           jump contraction is rho (NOT 2 alpha / rho^2); the reported rates
%           are the covariance rates directly, with no halving or rooting.
%    .persist 'flowing' | 'jumping' | 'none' (default 'none'): switches the
%           strict slack OFF the non-expansive side (persistent flowing => the
%           jump is only required non-expansive; persistent jumping => the flow
%           is only required non-expansive), all decay carried by the other side.
%
%  Conditions certified (P-side, no congruence):
%    flow    dP/dtau + A0^*P + P A0 + sum_l Ac_l^* P Ac_l + alpha P
%               + sum_r lamf_r (Jf_r^* P Jf_r - P) <= -slack I,  tau in [0,T]
%            (interval Positivstellensatz multiplier theta(1-theta) on [0,1] for
%             the timer types).
%    jump    sum_r J_r^* P(0) J_r - rho P(tau_jump) <= -slack I
%            (periodic/minimum: single end point tau = Tbar, i.e. theta = 1;
%             range: for all tau in [Tmin,Tmax] via the dwell multiplier;
%             arbitrary with scheduled J: constant P at both ends).
%    tail    minimum dwell only: frozen certificate P(Tbar) with dP/dtau = 0 must
%            still satisfy the flow inequality (persistence of flowing).
%    discrete   sum_l M_l^* P M_l - rho P <= -slack I, M_0 = A0, M_l = Ac_l.
%    positive/normalize   0 <= P <= I (lower bound forces a genuine, nonsingular
%            Lyapunov certificate; upper bound fixes the scale so the slack
%            maximization is bounded).  Both are enforced on [0,1] for the timer
%            types through the theta(1-theta) multiplier.
%
%  The problem maximizes the common slack; the system is certified stable iff the
%  optimal slack exceeds a tolerance.
%
%  RES: .feas, .slack, .type, .Pfun (@(tau)->P), .info.
%
%  OPTS (struct or MVST.OPTIONS, optional): .solver selects the SDP backend
%  ('sedumi' (default) | 'sdpt3' | 'mosek' | 'nosdp'), normalized by MVST_OPTS, which
%  errors ('mvst:opts:solverMissing') if the requested solver is not installed.
%
%  SYS may also be a THIN mode/system struct (MVST_MODE / MVST_SYSTEM): if it
%  lacks the .A0 field it is routed through MVST_TORICH to the rich naming.

if nargin < 3; opts = struct(); end
oo = mvst_opts(opts);
if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
n = size(sys.A0,1);
sys = fill_defaults(sys, n);
q = size(sys.Ac,3);  pf = size(sys.Jf,3);  p = size(sys.J,3);
alpha = getf(dwell,'alpha',0);  rho = getf(dwell,'rho',1);  tol = 1e-7;
% Persistent-flowing / persistent-jumping switches: the "persistent" side is only
% required non-expansive (no strict slack), so all decay is carried by the other
% side.  (Mirrors the synthesis routine th:stab:relax / th:stab:persjump.)
persist = getf(dwell,'persist','none');
flowSlack = double(~strcmp(persist,'jumping'));    % flow strict unless persistent jumping
jumpSlack = double(~strcmp(persist,'flowing'));    % jump strict unless persistent flowing

isTimer = any(strcmp(dwell.type,{'periodic','minimum','range'}));
deg = getf(dwell,'deg', 2*isTimer);
switch dwell.type
    case {'periodic','minimum'}; T = dwell.Tbar;
    case 'range';                T = dwell.Tmax;
    otherwise;                   T = 1;
end

pvar th; sl = dpvar('sl');
prog = sosprogram(th); prog = sosdecvar(prog, sl);
mons = monomials(th, 0:deg);
[prog, P] = sospolymatrixvar(prog, mons, [n n], 'symmetric');
gt = th*(1-th);

% ------------------------------------------------------------- flow condition
if strcmp(dwell.type,'discrete')
    Ld = discrete_cond(sys, P, rho, q);
    prog = matsos(prog, -Ld - (flowSlack*sl)*eye(n));
else
    Pdot = diff(P, th)/T;
    Lf = flow_cond(sys, P, Pdot, alpha, q, pf);
    fsl = flowSlack*sl;
    switch dwell.type
        case 'arbitrary'
            prog = matsos(prog, -Lf - fsl*eye(n));
        case {'periodic','minimum','range'}
            [prog, S1] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
            prog = matsos(prog, -Lf - fsl*eye(n) - gt*S1);
            prog = matsos(prog, S1);
    end
    if strcmp(dwell.type,'minimum')                    % frozen tail at theta=1
        LfT = flow_cond(sys, subs(P,th,1), zeros(n), alpha, q, pf);
        prog = matsos(prog, -LfT - fsl*eye(n));
    end
end

% -------------------------------------------------- scheduled-jump condition
if p > 0 && ~strcmp(dwell.type,'discrete')
    jsl = jumpSlack*sl;
    switch dwell.type
        case 'arbitrary'
            Lj = jump_cond(sys, P, P, rho, p);              % constant P both ends
            prog = matsos(prog, -Lj - jsl*eye(n));
        case {'periodic','minimum'}
            Lj = jump_cond(sys, subs(P,th,0), subs(P,th,1), rho, p);
            prog = matsos(prog, -Lj - jsl*eye(n));
        case 'range'
            Lj = jump_cond(sys, subs(P,th,0), P, rho, p);   % pre-jump P(tau) varies
            gm = (th - dwell.Tmin/T)*(dwell.Tmax/T - th);
            [prog, Sj] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
            prog = matsos(prog, -Lj - jsl*eye(n) - gm*Sj);
            prog = matsos(prog, Sj);
    end
end

% --------------------------------- positivity 0 <= P and normalization P <= I
if isTimer
    [prog, Xlo] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, P - gt*Xlo);            prog = matsos(prog, Xlo);
    [prog, Xhi] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, eye(n) - P - gt*Xhi);   prog = matsos(prog, Xhi);
else
    prog = matsos(prog, P);
    prog = matsos(prog, eye(n) - P);
end

prog = sossetobj(prog, -sl);
o = mvst_sosopts(oo);
prog = mvst_sossolve(prog, o);

sval = double(sosgetsol(prog, sl));
Pv   = sosgetsol(prog, P);

res.type  = dwell.type;
res.slack = sval;
res.feas  = sval > tol;
res.Pfun  = @(tau) evalmat(Pv, th, tau/T, n);
res.info  = struct('P',Pv,'T',T,'deg',deg,'alpha',alpha,'rho',rho);
end

% ==========================================================================
function L = flow_cond(sys, P, Pdot, alpha, q, pf)
% dP/dtau + A0^*P + P A0 + sum_l Ac_l^* P Ac_l + alpha P
%    + sum_r lamf_r (Jf_r^* P Jf_r - P)   (covariance-level shift +alpha P)
L = Pdot + sys.A0'*P + P*sys.A0 + alpha*P;
for l = 1:q; L = L + sys.Ac(:,:,l)'*P*sys.Ac(:,:,l); end
for r = 1:pf
    L = L + sys.lamf(r)*(sys.Jf(:,:,r)'*P*sys.Jf(:,:,r) - P);
end
end

function L = jump_cond(sys, P0, Pend, rho, p)
% sum_r J_r^* P0 J_r - rho Pend   (P0 = post-jump certificate at timer 0,
% Pend = pre-jump certificate at the end of the flow interval; covariance-level
% jump scaling rho)
L = -rho*Pend;
for r = 1:p; L = L + sys.J(:,:,r)'*P0*sys.J(:,:,r); end
end

function L = discrete_cond(sys, P, rho, q)
% sum_{l=0..q} M_l^* P M_l - rho P, with M_0 = A0 and M_l = Ac_l
% (covariance-level jump scaling rho)
L = sys.A0'*P*sys.A0 - rho*P;
for l = 1:q; L = L + sys.Ac(:,:,l)'*P*sys.Ac(:,:,l); end
end

function prog = matsos(prog, EXPR); prog = sosmatrixineq(prog, EXPR); end

function M = evalmat(P, th, theta, n)
if isnumeric(P); M = P; else; M = double(subs(P, th, theta)); end
M = reshape(M, n, n);
end

function sys = fill_defaults(sys, n)
if ~isfield(sys,'Ac')  ||isempty(sys.Ac);   sys.Ac  = zeros(n,n,0); end
if ~isfield(sys,'Jf')  ||isempty(sys.Jf);    sys.Jf  = zeros(n,n,0); end
if ~isfield(sys,'lamf')||isempty(sys.lamf);  sys.lamf= ones(1,size(sys.Jf,3)); end
if ~isfield(sys,'J')   ||isempty(sys.J);     sys.J   = zeros(n,n,0); end
end

function v = getf(s, f, d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end
