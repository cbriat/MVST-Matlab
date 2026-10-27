function res = mvst_dissyn_sos(sys, supply, dwell, mode, opts)
%MVST_DISSYN_SOS  Generic SOSTOOLS dissipative state-feedback SYNTHESIS program
%  for linear impulsive Hermitian matrix-valued systems: the synthesis
%  counterpart of the dissipativity ANALYSIS routine MVST_DISS_SOS and the
%  state-feedback sibling of the output-feedback routine MVST_DOF_SOS.  It
%  searches (SOSTOOLS + SeDuMi, matrix-SOS in the timer) for a state-feedback
%  gain K(tau) (and a jump gain K_J) rendering the closed loop dissipative with
%  respect to the quadratic supply priced by a port weight, over every
%  dwell-time regime of the chapter "Dissipativity-Based Synthesis".
%
%  ---------------------------------------------------------------------------
%  WHAT IS SOLVED (the book's convexified conditions).
%  ---------------------------------------------------------------------------
%  The closed-loop dissipation DLMIs (eq:disssyn:clDLMI) are bilinear in the
%  storage matrix S = -P and the gain K.  Following the same chapter
%  (Section~sec:disssyn:convex, Lemmas lem:disssyn:flow / lem:disssyn:jump) the
%  Schur STACKING of the completely positive channels, the CONGRUENCE by the
%  inverse storage matrix, and the CHANGE OF VARIABLES
%
%        W := S^{-1},   Y := K W,   Y_J := K_J W^-        (eq:disssyn:cov:flow,
%                                                          eq:disssyn:cov:jump)
%
%  turn them into the affine "synthesis blocks".  With the shorthands
%  (eq:disssyn:short:flow)
%        G0    = C W + Du Y,          V_l  = A_l W + B_l Y,
%        Phi11 = -Wdot + He(A0 W + B0 Y),
%        Phi12 = Bw - G0^* Pi2,       Phi22 = -Pi3 - He(Dw^* Pi2),
%  the FLOW synthesis block (eq:disssyn:Phif) is
%
%     Phi_f = [ Phi11     Phi12     row_l[V_l^*]    (N G0)^*  ]
%             [ *         Phi22     row_l[Bwl^*]    (N Dw)^*  ]  <= 0,
%             [ *         *         -I_q (x) W       0        ]
%             [ *         *         *               -I_nu     ]
%
%  and, with G_J = C_J W^- + DJu Y_J, Z_r = J_r W^- + BJ_r Y_J, the JUMP
%  synthesis block (eq:disssyn:Phij) is
%
%     Phi_j = [ -W^-      -G_J^* PiJ2    row_r[Z_r^*]      (N_J G_J)^*  ]
%             [ *          PhiJ22        row_r[Jwr^*]      (N_J DJw)^*  ]  <= 0,
%             [ *          *            -I_{p+1}(x) W^+     0           ]
%             [ *          *             *                -I_nuJ        ]
%
%  PhiJ22 = -PiJ3 - He(DJw^* PiJ2).  The port weight
%  Pi = [Pi1 Pi2; Pi2^* Pi3] prices (z;w); it is CONVEXIFIABLE (Pi1 <= 0,
%  Definition def:disssyn:convexifiable) with full-row-rank factor -Pi1 = N^*N,
%  N in C^{nu x pz}.  The four dwell-time theorems assembled here are:
%    'arbitrary' fixed jump sequence / CT LTI, constant W  (th:disssyn:seq);
%    'periodic'  Delta_k = Tp, flow on [0,Tp], jump Phi_j(W(Tp),W(0))
%                                                          (th:disssyn:periodic);
%    'minimum'   Delta_k >= Tmin, flow on [0,Tmin] + frozen tail
%                Phi_f^0(Tmin) (Wdot=0) + jump Phi_j(W(Tmin),W(0))
%                                                          (th:disssyn:MDT);
%    'range'     Delta_k in [Tmin,Tmax], flow on [0,Tmax] + jump Phi_j(theta)
%                for theta in [Tmin,Tmax]                  (th:disssyn:RDT).
%  The gains are recovered by eq:disssyn:gains, K = Y W^{-1}, K_J = Y_J (W^-)^{-1}.
%
%  ---------------------------------------------------------------------------
%  MODEL (rich naming; any omitted field defaults away).
%  ---------------------------------------------------------------------------
%  SYS fields:
%    .A0   n x n         flow drift A0
%    .B0   n x mu        flow control input  B0        (state feedback u = K x)
%    .Bw   n x m         flow disturbance input Bw
%    .Ac   n x n x q     flow state-multiplicative (Brownian) channels A_l
%    .Bc   n x mu x q    control selectors B_l  (paired with .Ac; default 0)
%    .Bwc  n x m x q     disturbance-multiplicative channels B_{w,l}
%                        (paired with .Ac; default 0)
%    .C    pz x n        flow performance output C
%    .Du   pz x mu       control feedthrough  Du       (default 0)
%    .Dw   pz x m        disturbance feedthrough Dw    (default 0)
%    .J    n x n [x P]   scheduled-jump maps J_r (P = p+1 channels; default I
%                        when a jump condition is required)
%    .BJ   n x mu [x P]  jump control selectors B_{J,r} (default 0)
%    .Jw   n x m  [x P]  jump-disturbance channels J_{w,r} (default 0)
%    .CJ   pzj x n       jump performance output C_J   (default none)
%    .DJu  pzj x mu      jump control feedthrough      (default 0)
%    .DJw  pzj x m       jump disturbance feedthrough  (default 0)
%  A thin mode/system struct (.A drift, .B input) or an mvst.System is routed
%  through MVST_TORICH first (A->A0, B->B0, .Ac/.Bc/.J/.BJ passed through).
%
%  SUPPLY:
%    .type 'l2gain' | 'passivity' | 'qsr'
%          ('qsr' = general quadratic supply [z;w]'*[Q S;S' R]*[z;w]: fields
%                         .Q (pz x pz), .S (pz x m), .R (m x m), optional jump
%                         weights .Qk .Sk .Rk; port weight Pi1 = Q, Pi2 = S,
%                         Pi3 = R; synthesis needs Q <= 0 (convexifiable,
%                         -Q = N'N) and mode 'feas'; validated by MVST_QSR_SUPPLY;
%           'l2gain'    = port weight Pi = blkdiag(-I_pz, gamma^2 I_m),
%                         Pi2 = 0, Pi3 = gamma^2 I, N = I_pz, nu = pz
%                         (eq:disssyn:M, Prop prop:disssyn:supplyid(a));
%           'passivity' = Pi = [-rho I, I; I, -nu I] (requires pz = m),
%                         Pi2 = I_m, Pi3 = -nu I, N = sqrt(rho) I (OSP), with
%                         .nu (ISP index) and .rho (OSP index), default 0
%                         (Prop prop:disssyn:supplyid(b)).)
%    .gamma  L2-gain level (used in 'feas' mode; free variable in 'gamma' mode).
%    .nu,.rho  passivity indices (default 0).
%
%  DWELL:
%    .type 'arbitrary' | 'periodic' | 'minimum' | 'range'
%    .Tbar scalar (periodic, minimum) ; .Tmin,.Tmax (range)
%    .deg  timer polynomial degree of W (default 2; forced 0 for 'arbitrary')
%    .alpha  optional flow shift +2 alpha W (default 0)
%
%  MODE: 'gamma' (default; bisection on gamma for the L2-gain supply) |
%        'feas'  (feasibility at the fixed supply.gamma).
%
%  OPTS (5th arg, all optional; see FILL_OPTS for the documented defaults):
%    .feastol  slack threshold for feasibility (default 1e-6, coherent with the
%              achieved SeDuMi residual so a sub-residual slack is not trusted).
%    .restol   relative LMI residual tolerance for the post-solve gate (1e-6).
%    .gaptol   allowed gap of the exact re-analysis over the certified level
%              (default 5e-2: achieved <= gamma*(1+gaptol)).
%    .epsW,.Rcap  certificate floor / cap (defaults 1e-6 / 1e4).
%    .stabtol  relative tolerance on the closed-loop lifted abscissa (1e-7).
%    .ngrid    timer-grid points for the numeric checks (default 11).
%    .validate master switch for the post-solve validation gate (default true).
%
%  POST-SOLVE VALIDATION GATE.  feas=true is NEVER returned on the strength of
%  the SOS slack alone.  Before returning, the actual certificate and LMI
%  blocks are recomputed at the RETURNED (W,Y,Y_J) and certified: (a) W is
%  positive definite on the timer grid; (b) the flow / frozen-tail / jump
%  blocks satisfy their inequalities to the residual tolerance; (c) the
%  closed loop is (mean-square) stable (lifted spectral abscissa, arbitrary
%  case) and (d) for the exact CT l2gain case re-analyses to a gain within
%  gaptol of the certified level.  A gate failure demotes the level to
%  infeasible, so the gamma-bisection converges to a genuinely certified
%  level; if nothing passes, res.feas=false and res.info.gate names the gate.
%
%  ---------------------------------------------------------------------------
%  RETURNS
%  ---------------------------------------------------------------------------
%  RES:
%    .feas       feasibility flag (common-slack maximization > tol).
%    .gamma      the returned (or fixed) supply level; aliased .gammastar.
%    .slack      certified slack margin.
%    .K          @(tau) -> flow state-feedback gain K(tau) = Y(tau) W(tau)^{-1}.
%    .KJ         jump gain K_J (empty when no scheduled jump); a function handle
%                @(theta) in the 'range' regime (dwell-dependent).
%    .S          @(tau) -> storage matrix S(tau) = W(tau)^{-1} = -P(tau) (>= 0).
%    .W          @(tau) -> the synthesis certificate W(tau) (>0).
%    .info       raw solved variables and solver info.
%
%  VALIDATION (mvst_dissyn_sos_selftest, run by mvst_dissyn_sos('selftest')):
%    the synthesized closed loop is re-analysed by the exact backend
%    MVST_L2GAIN_CT and by the analysis routine MVST_DISS_SOS; the achieved
%    dissipativity level agrees with the certified one to solver tolerance.
%
%  See also MVST_DISS_SOS, MVST_GAME_SOS, MVST_DOF_SOS, MVST_STABZ_SOS.

arguments
    sys
    supply struct = struct('type','l2gain')
    dwell  struct = struct('type','arbitrary')
    mode   char   = 'gamma'
    opts   struct = struct()
end

if ischar(sys) && strcmp(sys,'selftest')
    res = mvst_dissyn_sos_selftest(); return;
end

opts = fill_opts(opts);

% Thin struct / mvst.System -> rich naming (A->A0, B->B0, channels passed on).
if ~isstruct(sys) || ~isfield(sys,'A0')
    if isobject(sys); sys = obj2struct(sys); end
    if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
end
if ~isfield(sys,'B0') || isempty(sys.B0)
    error('mvst:dissyn:noB0', ...
        'SYS must carry a flow control input .B0 (state-feedback synthesis).');
end

n  = size(sys.A0,1);
mu = size(sys.B0,2);
m  = infer_m(sys, n);
sys = fill_defaults(sys, n, mu, m);
supply = check_supply(supply, sys, n, m);

switch mode
    case 'gamma'
        res = bisect_gamma(sys, supply, dwell, n, mu, m, opts);
    case 'feas'
        g = getf(supply,'gamma',1);
        res = solve_once(sys, supply, dwell, n, mu, m, g, true, opts);
        res = finish_gains(res, sys, mu, n, dwell);
        res.gammastar = g;
    otherwise
        error('mvst:dissyn:mode','unknown mode "%s"', mode);
end
end

% =========================================================================
function res = bisect_gamma(sys, supply, dwell, n, mu, m, opts)
%BISECT_GAMMA  Smallest gamma for which the synthesis DLMIs are SOS-feasible
%   AND pass the post-solve validation gate.  Feasibility is decided by slack
%   maximization (robust near the boundary) FOLLOWED by numeric certification
%   (validate_solution): a level is accepted only if the returned certificate
%   is genuinely positive definite, satisfies the LMIs at the returned point,
%   yields a stable closed loop, and (for the exact CT case) re-analyses to a
%   gain no larger than the certified level.  A gate failure at a level marks
%   that level infeasible, so the bisection converges to a CERTIFIED gamma.
if ~strcmp(supply.type,'l2gain')
    error('mvst:dissyn:gammamode','gamma mode requires the l2gain supply.');
end
lo = norm(sys.Dw) + 1e-4;
hi = max(2*lo, 1);
Rbest = []; grow = 0;
while true
    r = solve_once(sys, supply, dwell, n, mu, m, hi, true, opts);
    if r.feas; Rbest = r; break; end
    hi = 2*hi; grow = grow + 1;
    if grow > 30; break; end
end
for it = 1:26
    mid = 0.5*(lo + hi);
    r = solve_once(sys, supply, dwell, n, mu, m, mid, true, opts);
    if r.feas; hi = mid; Rbest = r; else; lo = mid; end
end
res = Rbest;
if isempty(res)
    res = solve_once(sys, supply, dwell, n, mu, m, hi, false, opts);
end
res.gamma = hi; res.gammastar = hi;
res = finish_gains(res, sys, mu, n, dwell);
end

% =========================================================================
function res = solve_once(sys, supply, dwell, n, mu, m, gam, useSlackObj, opts)
%SOLVE_ONCE  Build and solve the synthesis program at a fixed supply level,
%   THEN certify the returned solution numerically (validate_solution).  The
%   SOS slack decision variable alone is NOT trusted: on ill-conditioned
%   instances SeDuMi can report a positive slack while returning an indefinite
%   certificate and an LMI residual well above the slack, which previously let
%   a destabilizing gain slip through.  A returned feas=true now REQUIRES both
%   a positive slack (> opts.feastol) and a clean validation gate.
tol    = opts.feastol;  % slack threshold, coherent with the residual tolerance
Rcap   = opts.Rcap;     % certificate cap 0 < epsW I <= W <= Rcap I
margin = opts.feastol;
epsW   = opts.epsW;
alpha  = getf(dwell,'alpha',0);
[Pf, Pj] = port_weights(supply, sys, gam, m);   % flow/jump port-weight data

isTimer = any(strcmp(dwell.type,{'periodic','minimum','range'}));
deg = getf(dwell,'deg', 2*isTimer);
switch dwell.type
    case {'periodic','minimum'}; T = dwell.Tbar;
    case 'range';                T = dwell.Tmax;
    otherwise;                   T = 1;
end

needJump = jump_present(sys, dwell);

pvar th; sl = dpvar('sl');
prog = sosprogram(th); prog = sosdecvar(prog, sl);
mons = monomials(th, 0:deg);
[prog, W] = sospolymatrixvar(prog, mons, [n n], 'symmetric');
[prog, Y] = sospolymatrixvar(prog, mons, [mu n]);
YJ = zeros(mu,n);
if needJump; [prog, YJ] = sospolymatrixvar(prog, mons, [mu n]); end
gt = th*(1-th);
if useSlackObj; slk = sl; else; slk = margin; end

% -------------------------------------------------------------- flow condition
Wdot = diff(W, th)/T;
Phi  = flow_block(sys, W, Wdot, Y, Pf, alpha, n, mu, m);
switch dwell.type
    case 'arbitrary'
        prog = matle(prog, Phi, slk);
    case {'periodic','minimum','range'}
        nP = size(Phi,1);
        [prog, S1] = sospolymatrixvar(prog, monomials(th,0:2*deg), [nP nP], 'symmetric');
        prog = matsos(prog, -Phi - slk*eye(nP) - gt*S1);
        prog = matsos(prog, S1);
    otherwise
        error('mvst:dissyn:dwell','unknown dwell type "%s"', dwell.type);
end
if strcmp(dwell.type,'minimum')                % frozen tail Phi_f^0(Tmin)
    Phi0 = flow_block(sys, subs(W,th,1), zeros(n), subs(Y,th,1), Pf, alpha, n, mu, m);
    prog = matle(prog, Phi0, slk);
end

% -------------------------------------------------------------- scheduled jump
if needJump
    switch dwell.type
        case 'arbitrary'
            Dj = jump_block(sys, W, W, YJ, Pj, n, mu);
            prog = matle(prog, Dj, slk);
        case {'periodic','minimum'}
            Dj = jump_block(sys, subs(W,th,1), subs(W,th,0), subs(YJ,th,1), Pj, n, mu);
            prog = matle(prog, Dj, slk);
        case 'range'
            Dj = jump_block(sys, W, subs(W,th,0), YJ, Pj, n, mu);
            gm = (th - dwell.Tmin/T)*(dwell.Tmax/T - th);
            nj = size(Dj,1);
            [prog, Sj] = sospolymatrixvar(prog, monomials(th,0:2*deg), [nj nj], 'symmetric');
            prog = matsos(prog, -Dj - slk*eye(nj) - gm*Sj);
            prog = matsos(prog, Sj);
    end
end

% ---------------------------- certificate scaling 0 < epsW I <= W <= Rcap I
if isTimer
    [prog, Xlo] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, W - epsW*eye(n) - gt*Xlo);   prog = matsos(prog, Xlo);
    [prog, Xhi] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, Rcap*eye(n) - W - gt*Xhi);   prog = matsos(prog, Xhi);
else
    prog = matsos(prog, W - epsW*eye(n));
    prog = matsos(prog, Rcap*eye(n) - W);
end

if useSlackObj; prog = sossetobj(prog, -sl); end
o = mvst_sosopts(mvst_opts(opts, true, 'sos'));
[prog, info] = mvst_sossolve(prog, o);

Wv = sosgetsol(prog, W);  Yv = sosgetsol(prog, Y);
YJv = zeros(mu,n); if needJump; YJv = sosgetsol(prog, YJ); end
if useSlackObj
    sval = double(sosgetsol(prog, sl));
    feas = sval > tol;
else
    sval = margin;
    feas = (info.pinf == 0 && info.dinf == 0 && info.numerr < 2);
end

res.type  = dwell.type;
res.slack = sval;
res.feas  = feas;
res.gamma = gam;
res.info  = struct('W',Wv,'Y',Yv,'YJ',YJv,'T',T,'deg',deg, ...
                   'alpha',alpha,'needJump',needJump,'solverinfo',info);

% ---------------- POST-SOLVE VALIDATION GATE -----------------------------
% Never return feas=true on the strength of the SOS slack alone.  Recompute
% the actual certificate and LMI blocks at the RETURNED (W,Y,Y_J) and certify
% them; on any failure demote the level to infeasible (the bisection then
% converges to a genuinely certified gamma) and record which gate failed.
res.info.gate = struct('ok',true,'reason','not-run');
if feas && opts.validate
    [gok, ginfo] = validate_solution(res, sys, supply, dwell, n, mu, m, gam, opts);
    res.info.gate = ginfo;
    if ~gok
        res.feas = false;
    end
end
end

% =========================================================================
function opts = fill_opts(opts)
%FILL_OPTS  Documented defaults for the solver / validation tolerances.
%   .feastol  slack threshold deciding feasibility (default 1e-6).  Coherent
%             with the achieved SeDuMi accuracy: a claimed slack below the
%             solver residual is noise, so the threshold is not set tighter.
%   .restol   LMI residual tolerance for the numeric gate (default 1e-6),
%             applied RELATIVE to the block norm at each grid point.
%   .gaptol   allowed gap of the exact re-analysis over the certified level
%             (default 5e-2): achieved <= gamma*(1+gaptol).
%   .epsW     certificate floor 0 < epsW I <= W (default 1e-6); the PD gate
%             requires min eig(W) >= epsW/2 on the timer grid.
%   .Rcap     certificate cap W <= Rcap I (default 1e4).
%   .stabtol  relative tolerance on the closed-loop lifted spectral abscissa
%             (default 1e-7): abscissa must be < stabtol*max(1,||M||).
%   .ngrid    number of timer-grid points for the numeric checks (default 11).
%   .validate master switch for the post-solve gate (default true).
d = struct('feastol',1e-6,'restol',1e-6,'gaptol',5e-2,'epsW',1e-6, ...
           'Rcap',1e4,'stabtol',1e-7,'ngrid',11,'validate',true);
fn = fieldnames(d);
for k = 1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k}))
        opts.(fn{k}) = d.(fn{k});
    end
end
end

% =========================================================================
function [ok, gi] = validate_solution(res, sys, supply, dwell, n, mu, m, gam, opts)
%VALIDATE_SOLUTION  Numeric certification of a returned synthesis solution.
%   Recomputes, at the RETURNED (W,Y,Y_J), the four gates of the fix:
%     (a) certificate positivity   min eig(W(tau)) >= epsW/2 on the grid;
%     (b) LMI residual             max eig(block(tau)) <= restol*||block||
%                                  for the flow block (grid), the frozen tail
%                                  (minimum dwell) and the jump block(s);
%     (c) closed-loop stability    for the CT LTI / SDE 'arbitrary' case the
%                                  lifted operator He(Acl X)+sum Ccl_l X Ccl_l'
%                                  (mvst_liftsym) must be Hurwitz.  For the
%                                  timer-varying dwell regimes the DLMI gates
%                                  (a)+(b), holding at the returned certificate
%                                  across the grid, ARE the stability/
%                                  dissipativity certificate (an independent
%                                  monodromy is not recomputed);
%     (d) exact re-analysis        for l2gain + pure-CT 'arbitrary' the closed
%                                  loop is re-analysed by mvst_l2gain_ct and
%                                  the achieved gain must be <= gam*(1+gaptol).
%   Returns ok=false and a populated GI (which gate, the offending value) on
%   the FIRST failing gate.
gi = struct('ok',true,'reason','pass','gate','', ...
            'mineigW',Inf,'maxeigLMI',-Inf,'abscissa',NaN,'achieved',NaN, ...
            'condW',NaN,'gam',gam);
ok = true;
pvar th;   % symbolic timer used by evalmat/diff on the returned polynomials
Wv = res.info.W; Yv = res.info.Y; YJv = res.info.YJ; T = res.info.T;
alpha = res.info.alpha; needJump = res.info.needJump;
epsW = opts.epsW; restol = opts.restol;
[Pf, Pj] = port_weights(supply, sys, gam, m);

isTimer = any(strcmp(dwell.type,{'periodic','minimum','range'}));
if isTimer; thgrid = linspace(0,1,opts.ngrid); else; thgrid = 0; end
if isnumeric(Wv); Wdotpoly = []; else; Wdotpoly = diff(Wv, th); end

% ---- (a) certificate positivity + (b) flow-block residual on the grid ----
for thk = thgrid
    Wk = symm(evalmat(Wv, th, thk, n, n));
    ew = min(eig(Wk));
    gi.mineigW = min(gi.mineigW, ew);
    gi.condW   = max(gi.condW, cond(Wk));
    if ew < epsW/2
        [ok,gi] = fail(gi,'certificate-PD', ...
            sprintf('min eig W(tau=%.3g)=%.3e < epsW/2=%.1e', thk*T, ew, epsW/2));
        return;
    end
    if isempty(Wdotpoly); Wdotk = zeros(n);
    else; Wdotk = evalmat(Wdotpoly, th, thk, n, n)/T; end
    Yk  = evalmat(Yv, th, thk, mu, n);
    Phi = symm(flow_block(sys, Wk, Wdotk, Yk, Pf, alpha, n, mu, m));
    [ok,gi] = check_le(Phi, restol, 'flow-LMI', thk*T, gi);
    if ~ok; return; end
end

% ---- (b) frozen-tail block for the minimum-dwell regime -----------------
if strcmp(dwell.type,'minimum')
    W1 = symm(evalmat(Wv, th, 1, n, n));
    Y1 = evalmat(Yv, th, 1, mu, n);
    Phi0 = symm(flow_block(sys, W1, zeros(n), Y1, Pf, alpha, n, mu, m));
    [ok,gi] = check_le(Phi0, restol, 'frozen-tail-LMI', T, gi);
    if ~ok; return; end
end

% ---- (b) jump-block residual --------------------------------------------
if needJump
    switch dwell.type
        case 'arbitrary'
            Wm = symm(evalmat(Wv,th,0,n,n)); YJk = evalmat(YJv,th,0,mu,n);
            Dj = symm(jump_block(sys, Wm, Wm, YJk, Pj, n, mu));
            [ok,gi] = check_le(Dj, restol, 'jump-LMI', 0, gi); if ~ok; return; end
        case {'periodic','minimum'}
            Wm = symm(evalmat(Wv,th,1,n,n)); Wp = symm(evalmat(Wv,th,0,n,n));
            YJk = evalmat(YJv,th,1,mu,n);
            Dj = symm(jump_block(sys, Wm, Wp, YJk, Pj, n, mu));
            [ok,gi] = check_le(Dj, restol, 'jump-LMI', T, gi); if ~ok; return; end
        case 'range'
            Wp = symm(evalmat(Wv,th,0,n,n));
            for tj = linspace(dwell.Tmin/T, dwell.Tmax/T, opts.ngrid)
                Wm = symm(evalmat(Wv,th,tj,n,n)); YJk = evalmat(YJv,th,tj,mu,n);
                Dj = symm(jump_block(sys, Wm, Wp, YJk, Pj, n, mu));
                [ok,gi] = check_le(Dj, restol, 'jump-LMI', tj*T, gi);
                if ~ok; return; end
            end
    end
end

% ---- (c) closed-loop mean-square stability (CT LTI / SDE 'arbitrary') ---
if strcmp(dwell.type,'arbitrary') && ~needJump
    K0  = evalmat(Yv,th,0,mu,n) / symm(evalmat(Wv,th,0,n,n));
    Acl = sys.A0 + sys.B0*K0;
    q   = size(sys.Ac,3);
    Cc  = cell(1,q);
    for l = 1:q; Cc{l} = sys.Ac(:,:,l) + sys.Bc(:,:,l)*K0; end
    M   = mvst_liftsym(Acl, Cc);
    absc = max(real(eig(M)));
    gi.abscissa = absc;
    if absc > opts.stabtol*max(1,norm(M,'fro'))
        [ok,gi] = fail(gi,'closed-loop-stability', ...
            sprintf('lifted spectral abscissa %.4g >= 0 (unstable closed loop)', absc));
        return;
    end
end

% ---- (d) exact re-analysis for the l2gain pure-CT 'arbitrary' case ------
q = size(sys.Ac,3);
if strcmp(supply.type,'l2gain') && strcmp(dwell.type,'arbitrary') ...
        && ~needJump && q == 0
    K0  = evalmat(Yv,th,0,mu,n) / symm(evalmat(Wv,th,0,n,n));
    Acl = sys.A0 + sys.B0*K0;   Ccl = sys.C + sys.Du*K0;
    try
        achieved = mvst_l2gain_ct(Acl, sys.Bw, Ccl, sys.Dw);
    catch ME   % non-Hurwitz closed loop etc. -> gate failure
        [ok,gi] = fail(gi,'exact-reanalysis', ...
            sprintf('mvst_l2gain_ct rejected the closed loop (%s)', ME.message));
        return;
    end
    gi.achieved = achieved;
    if achieved > gam*(1+opts.gaptol)
        [ok,gi] = fail(gi,'exact-reanalysis', ...
            sprintf('achieved gain %.4g > gamma*(1+gaptol)=%.4g', ...
                    achieved, gam*(1+opts.gaptol)));
        return;
    end
end
end

% -------------------------------------------------------------------------
function [ok,gi] = check_le(B, restol, gatename, tau, gi)
%CHECK_LE  Certify B <= 0 numerically: max eig(B) <= restol*max(1,||B||).
lam = max(eig(B));
gi.maxeigLMI = max(gi.maxeigLMI, lam);
tolb = restol*max(1, norm(B,'fro'));
if lam > tolb
    [ok,gi] = fail(gi, gatename, ...
        sprintf('max eig block(tau=%.3g)=%.3e > tol=%.1e', tau, lam, tolb));
else
    ok = true;
end
end

function [ok,gi] = fail(gi, gatename, reason)
ok = false; gi.ok = false; gi.gate = gatename; gi.reason = reason;
end

function S = symm(A)
S = (A + A')/2;
end

% =========================================================================
function C = flow_block(sys, W, Wdot, Y, Pf, alpha, n, mu, m)
%FLOW_BLOCK  Flow synthesis block Phi_f (require <= 0), eq:disssyn:Phif.
q  = size(sys.Ac,3);
AW = sys.A0*W + sys.B0*Y;                     % A0 W + B0 Y
Phi11 = -Wdot + (AW + AW.') + 2*alpha*W;
G0  = sys.C*W + sys.Du*Y;                     % pz x n
Phi12 = sys.Bw - G0.'*Pf.Pi2;                 % n x m
Phi22 = -Pf.Pi3 - (sys.Dw.'*Pf.Pi2 + Pf.Pi2.'*sys.Dw);
core = [ Phi11,      Phi12 ; ...
         Phi12.',    Phi22 ];
% ---- Brownian / disturbance-multiplicative Schur channels -> -I_q (x) W ----
Vrow = []; Bwrow = []; Wpiv = [];
for l = 1:q
    V = sys.Ac(:,:,l)*W + sys.Bc(:,:,l)*Y;    % V_l
    Vrow  = [Vrow,  V.'];                      %#ok<AGROW>
    Bwrow = [Bwrow, sys.Bwc(:,:,l).'];         %#ok<AGROW>
    Wpiv  = blkdiagcat(Wpiv, -W);
end
if q == 0; Vrow = zeros(n,0); Bwrow = zeros(m,0); end
% ---- output channel -> -I_nu  (only when nu > 0) ----
if Pf.nu > 0
    Orow = [ (Pf.N*G0).' ; (Pf.N*sys.Dw).' ];   % (n+m) x nu
else
    Orow = zeros(n+m, 0);
end
NoiseTop = [Vrow ; Bwrow];                       % (n+m) x qn
C = [ core,          NoiseTop,               Orow ; ...
      NoiseTop.',    Wpiv,                   zeros(q*n, Pf.nu) ; ...
      Orow.',        zeros(Pf.nu, q*n),     -eye(Pf.nu) ];
C = (C + C.')/2;
end

% -------------------------------------------------------------------------
function D = jump_block(sys, Wm, Wp, YJ, Pj, n, mu) %#ok<INUSD>
%JUMP_BLOCK  Jump synthesis block Phi_j(W^-,W^+) (require <= 0), eq:disssyn:Phij.
%  Wm = W^- (pre-jump), Wp = W^+ (post-jump).  Channels r = 0..p from .J.  The
%  JUMP disturbance dimension md (columns of .Jw) is separate from the flow
%  disturbance m; when md = 0 (a pure state reset, no discrete disturbance) the
%  w-row/column is absent, exactly as the analysis block reduces (jump_dM).
P  = size(sys.J,3);
md = size(sys.Jw,2);
GJ = sys.CJ*Wm + sys.DJu*YJ;                  % pzj x n  (empty pzj -> 0 rows)
PhiJ12 = -GJ.'*Pj.Pi2;                         % n x md
PhiJ22 = -Pj.Pi3 - (sys.DJw.'*Pj.Pi2 + Pj.Pi2.'*sys.DJw);   % md x md
core = [ -Wm,         PhiJ12 ; ...
         PhiJ12.',    PhiJ22 ];
Zrow = []; Jwrow = []; Wpiv = [];
for r = 1:P
    Z = sys.J(:,:,r)*Wm + sys.BJ(:,:,r)*YJ;   % Z_r
    Zrow  = [Zrow,  Z.'];                      %#ok<AGROW>
    Jwrow = [Jwrow, sys.Jw(:,:,r).'];          %#ok<AGROW>
    Wpiv  = blkdiagcat(Wpiv, -Wp);
end
if P == 0; Zrow = zeros(n,0); Jwrow = zeros(md,0); end
if Pj.nu > 0
    Orow = [ (Pj.N*GJ).' ; (Pj.N*sys.DJw).' ];   % (n+md) x nuJ
else
    Orow = zeros(n+md, 0);
end
JumpTop = [Zrow ; Jwrow];                        % (n+md) x Pn
D = [ core,         JumpTop,               Orow ; ...
      JumpTop.',    Wpiv,                  zeros(P*n, Pj.nu) ; ...
      Orow.',       zeros(Pj.nu, P*n),    -eye(Pj.nu) ];
D = (D + D.')/2;
end

% =========================================================================
function [Pf, Pj] = port_weights(supply, sys, gam, m)
%PORT_WEIGHTS  Flow and jump port-weight data (Pi2, Pi3, factor N, rank nu),
%  eq:disssyn:port / Def def:disssyn:convexifiable / Prop prop:disssyn:supplyid.
pz  = size(sys.C,1);
Pf  = one_port(supply, gam, m, pz);
% jump supply on md (jump-disturbance dim) and pzj (jump-output dim); active iff
% a jump output or a jump disturbance is priced, else a pure Lyapunov reset.
pzj = size(sys.CJ,1);
md  = size(sys.Jw,2);
if pzj == 0 && md == 0
    Pj = struct('Pi2', zeros(0,0), 'Pi3', zeros(0,0), 'N', zeros(0,0), 'nu', 0);
else
    Pj = one_port(supply, gam, md, pzj, 'jump');
end
end

function Pp = one_port(supply, gam, m, pz, qwhere)
if nargin < 5; qwhere = 'flow'; end
%ONE_PORT  Port-weight blocks for one condition (flow or jump) with pz outputs.
nu_isp = getf(supply,'nu',0);  rho = getf(supply,'rho',0);
switch supply.type
    case 'l2gain'
        Pp.Pi2 = zeros(pz, m);
        Pp.Pi3 = gam^2*eye(m);
        Pp.N   = eye(pz);        Pp.nu = pz;      % -Pi1 = I_pz = N^*N
    case 'passivity'
        if pz ~= m
            error('mvst:dissyn:passdim', ...
                'passivity requires pz = m (got pz=%d, m=%d).', pz, m);
        end
        Pp.Pi2 = eye(m);
        Pp.Pi3 = -nu_isp*eye(m);
        if rho > 0                                % OSP: -Pi1 = rho I = N^*N
            Pp.N = sqrt(rho)*eye(pz); Pp.nu = pz;
        else
            Pp.N = zeros(0,pz);       Pp.nu = 0;
        end
    case 'qsr'
        % Pi1 = Q (must be <= 0: convexifiable), Pi2 = S, Pi3 = R.
        [Qs,Ss,Rs] = mvst_qsr_supply(supply, pz, m, qwhere);
        Qs = (Qs+Qs')/2;
        [V,D] = eig(Qs); dd = diag(D);
        tolq = 1e-10*max(1,norm(Qs,'fro'));
        if any(dd > tolq)
            error('mvst:dissyn:qsrconvex', ...
                ['qsr synthesis needs Q <= 0 (convexifiable supply, -Q = N''N); ' ...
                 'got lambda_max(Q) = %g.'], max(dd));
        end
        act = dd < -tolq;
        Pp.Pi2 = Ss;  Pp.Pi3 = Rs;
        Pp.N   = diag(sqrt(-dd(act)))*V(:,act)';  Pp.nu = nnz(act);
    otherwise
        error('mvst:dissyn:supply','unknown supply "%s"', supply.type);
end
end

% =========================================================================
function res = finish_gains(res, sys, mu, n, dwell)
%FINISH_GAINS  Attach the recovered gains K = Y W^{-1}, K_J = Y_J (W^-)^{-1}.
if ~isfield(res,'info') || ~isfield(res.info,'W')
    res.K = []; res.KJ = []; res.S = []; res.W = []; return;
end
Wv = res.info.W; Yv = res.info.Y; YJv = res.info.YJ; T = res.info.T;
pvar th;
res.W = @(tau) evalmat(Wv, th, tau/T, n, n);
res.S = @(tau) inv(evalmat(Wv, th, tau/T, n, n));
res.K = @(tau) evalmat(Yv, th, tau/T, mu, n) / evalmat(Wv, th, tau/T, n, n);
res.KJ = [];
if res.info.needJump
    switch dwell.type
        case 'arbitrary'
            res.KJ = evalmat(YJv,th,0,mu,n) / evalmat(Wv,th,0,n,n);
        case {'periodic','minimum'}
            res.KJ = evalmat(YJv,th,1,mu,n) / evalmat(Wv,th,1,n,n);
        case 'range'
            % dwell-dependent jump gain K_J(Delta) = Y_J(Delta) W(Delta)^{-1}
            res.KJ = @(theta) evalmat(YJv,th,theta/T,mu,n) / evalmat(Wv,th,theta/T,n,n);
    end
end
end

% =========================================================================
function prog = matle(prog, EXPR, slk)
% EXPR <= -slk I  as a matrix SOS (constant multiplier: plain LMI).
nn = size(EXPR,1);
prog = sosmatrixineq(prog, -EXPR - slk*eye(nn));
end
function prog = matsos(prog, EXPR); prog = sosmatrixineq(prog, EXPR); end

function Mout = evalmat(P, th, theta, r, c)
if isnumeric(P); Mout = P; else; Mout = double(subs(P, th, theta)); end
Mout = reshape(Mout, r, c);
end

function B = blkdiagcat(B, X)
% Append X on the block diagonal of B (both square; sizes may differ).
if isempty(B); B = X; return; end
[a,~] = size(B); [c,~] = size(X);
B = [ B, zeros(a,c) ; zeros(c,a), X ];
end

% =========================================================================
function tf = jump_present(sys, dwell)
%JUMP_PRESENT  Whether a scheduled-jump condition must be enforced.
tf = any(strcmp(dwell.type,{'periodic','minimum','range'})) || ...
     (~isempty(sys.J) && ~isequal_id(sys.J)) || ...
     size(sys.CJ,1) > 0 || any(sys.Jw(:) ~= 0);
if strcmp(dwell.type,'arbitrary') && isequal_id(sys.J) ...
        && size(sys.CJ,1)==0 && ~any(sys.Jw(:) ~= 0)
    tf = false;                                   % pure CT: no jump condition
end
end

function tf = isequal_id(J)
% true iff J is a single identity channel
tf = (size(J,3) == 1) && isequal(J(:,:,1), eye(size(J,1)));
end

function m = infer_m(sys, n)
if isfield(sys,'Bw') && ~isempty(sys.Bw); m = size(sys.Bw,2);
elseif isfield(sys,'Dw') && ~isempty(sys.Dw); m = size(sys.Dw,2);
else; m = 1; end
end

function supply = check_supply(supply, sys, n, m) %#ok<INUSD>
if ~isfield(supply,'type') || isempty(supply.type); supply.type = 'l2gain'; end
end

function sys = fill_defaults(sys, n, mu, m)
z = @(a,b,c) zeros(a,b,c);
if ~isfield(sys,'Bw') ||isempty(sys.Bw);  sys.Bw  = zeros(n,m);  end
if ~isfield(sys,'Ac') ||isempty(sys.Ac);  sys.Ac  = z(n,n,0);    end
% Accept the ANALYSIS-routine field name .Dwc (SHARED disturbance-multiplicative
% Brownian channels D_{w,l}, n x m x qd, cell or 3-D array) as an alias for the
% per-channel .Bwc, mirroring MVST_DISS_SOS_MODES: the l-th Brownian increment
% drives both the state channel .Ac(:,:,l) and the disturbance channel
% .Dwc(:,:,l), so they load the flow pivot Theta jointly (the D_{w,l}^*SD_{w,l}
% loading is the completely positive Schur channel against -I_q(x)W).  When more
% disturbance channels than state channels are supplied (a PURE
% disturbance-multiplicative channel, Ac_l = 0), .Ac is zero-padded so the two
% stay index-paired.  A directly supplied .Bwc takes precedence (no override).
Dwc0 = getf(sys,'Dwc',[]);
if ~isempty(Dwc0) && (~isfield(sys,'Bwc') || isempty(sys.Bwc))
    if iscell(Dwc0); Dwc0 = cat(3, Dwc0{:}); end
    qd = size(Dwc0,3); q0 = size(sys.Ac,3); qmax = max(q0,qd);
    if q0 < qmax; sys.Ac(:,:,q0+1:qmax) = 0; end        % pad state channels
    if qd < qmax; Dwc0(:,:,qd+1:qmax)   = 0; end        % pad disturbance channels
    sys.Bwc = Dwc0;
end
q = size(sys.Ac,3);
if ~isfield(sys,'Bc') ||isempty(sys.Bc)||size(sys.Bc,3)~=q;  sys.Bc  = z(n,mu,q); end
if ~isfield(sys,'Bwc')||isempty(sys.Bwc)||size(sys.Bwc,3)~=q; sys.Bwc = z(n,m,q);  end
if ~isfield(sys,'C')  ;                   sys.C   = zeros(0,n);  end
if ~isfield(sys,'Du') ||isempty(sys.Du);  sys.Du  = zeros(size(sys.C,1),mu); end
if ~isfield(sys,'Dw') ||isempty(sys.Dw);  sys.Dw  = zeros(size(sys.C,1),m);  end
% scheduled jump
if ~isfield(sys,'J')  ||isempty(sys.J);   sys.J   = eye(n);      end
P = size(sys.J,3);
if ~isfield(sys,'BJ') ||isempty(sys.BJ)||size(sys.BJ,3)~=P; sys.BJ = z(n,mu,P); end
if ~isfield(sys,'CJ') ;                   sys.CJ  = zeros(0,n);  end
pzj = size(sys.CJ,1);
% jump-disturbance dimension md is inferred from .Jw/.DJw (default 0: pure
% state reset).  It is INDEPENDENT of the flow disturbance dimension m.
Jw0  = getf(sys,'Jw',[]);  DJw0 = getf(sys,'DJw',[]);
if     ~isempty(Jw0);  md = size(Jw0,2);
elseif ~isempty(DJw0); md = size(DJw0,2);
else;                  md = 0; end
if ~isfield(sys,'Jw') ||isempty(sys.Jw)||size(sys.Jw,2)~=md||size(sys.Jw,3)~=P
    sys.Jw = z(n,md,P); end
if ~isfield(sys,'DJu')||isempty(sys.DJu); sys.DJu = zeros(pzj,mu); end
if ~isfield(sys,'DJw')||isempty(sys.DJw)||size(sys.DJw,2)~=md; sys.DJw = zeros(pzj,md); end
end

function s = obj2struct(obj)
% Best-effort mvst.System -> struct (field copy of public properties).
s = struct();
p = properties(obj);
for k = 1:numel(p); s.(p{k}) = obj.(p{k}); end
end

function v = getf(s, f, d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end
