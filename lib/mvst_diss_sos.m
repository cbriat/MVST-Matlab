function res = mvst_diss_sos(sys, supply, dwell, mode, opts)
%MVST_DISS_SOS  Generic SOSTOOLS dissipativity-ANALYSIS program for linear
%  impulsive Hermitian matrix-valued systems, covering the L2-gain and
%  passivity supplies and every dwell-time type of the Dissipativity chapter.
%  This is the dissipativity counterpart of MVST_STAB_SOS (stability analysis)
%  and MVST_STABZ_SOS (stabilization synthesis): the storage certificate P is
%  the dual variable of the Dissipativity chapter (Def.~def:IMP:stor, storage
%  matrix S = -P >= 0), and the routine solves the dual differential LMIs
%  C_M(t,P) >= 0 on flow and D_M(k,P^+,P^-) >= 0 at jumps (eq:IMP:LM,
%  Theorem th:IMP:diss:LMI and the dwell-time DLMIs th:diss:periodic:LMI,
%  th:diss:MDT:LMI, th:diss:RDT:LMI) as a matrix sum-of-squares program in the
%  timer with SOSTOOLS + SeDuMi.
%
%  The certificate P is n x n and sign-free (no congruence, no Schur
%  complement): the dissipation operators enter P linearly, and the
%  disturbance/output blocks are carried EXPLICITLY as the (n+m) augmented
%  matrix inequality, so the state-multiplicative second-moment terms
%  sum_l Ac_l^* P Ac_l, the disturbance-multiplicative pivot loading
%  sum_l Dwc_l^* P Dwc_l, the Poisson terms sum_r lamf_r(Jf_r^* P Jf_r - P),
%  and the Poisson jump-disturbance terms sum_r lamf_r Jf_r^* P Bdf_r (with the
%  pivot loading sum_r lamf_r Bdf_r^* P Bdf_r) all enter without eliminating
%  the free block.  Concretely the flow operator built here is
%    C_M(P) = [ Pdot + He(A0^*P) + sum_l Ac_l^*P Ac_l
%                     + sum_r lamf_r(Jf_r^*P Jf_r - P) + M11 ,
%               P Bw + sum_l Ac_l^*P Dwc_l + sum_r lamf_r Jf_r^*P Bdf_r + M12 ;
%               *  ,  sum_l Dwc_l^*P Dwc_l + sum_r lamf_r Bdf_r^*P Bdf_r + M22 ]  >=  0,
%  which is the general (impulsive, GBM, Poisson) bounded-real flow LMI
%  eq:diss:br:gen / eq:diss:bri:gbm / eq:diss:bri:poisson BEFORE the Schur
%  reduction of Remark rem:diss:br:cases; the storage matrix S = -P >= 0 makes
%  it the negated bounded-real inequality of that remark.
%
%  SYS fields (any omitted field defaults to "none"):
%    .A0   n x n         flow drift
%    .Ac   n x n x q     flow state-multiplicative (Brownian) channels C_l
%    .Dwc  n x m x q     flow disturbance-multiplicative channels D_{w,l}
%                        (paired with .Ac; default 0: load the pivot Theta_S)
%    .Bw   n x m         flow disturbance input
%    .C    p x n         flow performance output
%    .Dw   p x m         flow output feedthrough
%    .Jf   n x n x pf    SPONTANEOUS (Poisson) jump maps J_r
%    .lamf 1 x pf        Poisson rates
%    .Bdf  n x m x pf    Poisson jump-disturbance channels D_r (paired with
%                        .Jf/.lamf; default 0), so the spontaneous reset is
%                        x^+ = J_r x + D_r w_c, driven by the flow disturbance
%    .J    n x n         SCHEDULED (dwell-time / sampling) jump map (default I
%                        when a jump condition is required)
%    .Bd   n x md        discrete-disturbance channel at the jump
%    .CJ   pd x n        jump performance output
%    .DJ   pd x md       jump output feedthrough
%
%  SUPPLY:
%    .type 'l2gain' | 'passivity' | 'qsr'
%          ('qsr' = the general quadratic supply [z;w]'*[Q S;S' R]*[z;w] with
%                        z = C x + Dw w; fields .Q (p x p), .S (p x m), .R (m x m),
%                        optional jump weights .Qk .Sk .Rk; blocks
%                        (M11,M12,M22) = (C'QC, C'(Q Dw + S), Dw'Q Dw + He(S'Dw) + R);
%                        l2gain is Q=-I,S=0,R=gamma^2 I, passivity is Q=-rho I,
%                        S=I, R=-nu I.  Validated by MVST_QSR_SUPPLY.  Feasibility
%                        mode only ('gamma' mode needs l2gain).)
%          ('l2gain'   = M = diag(0,gamma^2 I) - Ctilde^*Ctilde, the hybrid
%                        L2-gain / bounded-real supply eq:diss:M:l2, blocks
%                        (M11,M12,M22) = (-C^*C, -C^*Dw, gamma^2 I - Dw^*Dw);
%           'passivity'= M = [0 C^*; C He(Dw)] eq:diss:M:pass, with strict
%                        refinements .nu (ISP: M22 -= nu I) and .rho (OSP:
%                        M -= rho Ctilde^*Ctilde).)
%    .gamma  L2-gain level (used in 'feas' mode; the free variable in 'gamma'
%            mode).  .nu,.rho  passivity indices (default 0).
%
%  DWELL:
%    .type 'arbitrary' | 'discrete' | 'periodic' | 'minimum' | 'range'
%          ('arbitrary' = constant storage (CT LTI, GBM, Poisson: no scheduled
%                         jump, the Poisson jumps live in the flow via Jf/lamf);
%           'discrete'  = pure DT bounded-real: jump map A0, input Bw, output
%                         (C,Dw), constant P, only D_M(P,P) >= 0;
%           'periodic'  = fixed dwell-time Tbar, single jump at tau = Tbar (the
%                         sampled-data / sampled-output case: J = I, output CJ);
%           'minimum'   = minimum dwell-time Tbar, ramp on [0,Tbar] plus the
%                         frozen-tail flow condition C_M^0(Tbar) >= 0;
%           'range'     = ranged dwell-time [Tmin,Tmax], jump for any tau in
%                         [Tmin,Tmax] via the dwell multiplier.)
%    .Tbar scalar (periodic, minimum) ; .Tmin,.Tmax (range)
%    .deg  timer polynomial degree of P (default 2; forced 0 for
%          arbitrary/discrete)
%    .alpha  optional flow shift +2 alpha P (default 0)
%
%  MODE: 'feas'  (dissipativity feasibility with slack maximization) |
%        'gamma' (bisection on gamma for the L2-gain supply).
%
%  RES: .feas, .slack, .gamma, .Pfun (@(tau)->P), .Sfun (@(tau)-> S = -P),
%       .info.
%
%  OPTS (struct or MVST.OPTIONS, optional): .solver selects the SDP backend
%  ('sedumi' (default) | 'sdpt3' | 'mosek' | 'nosdp'), normalized by MVST_OPTS.
%
%  VALIDATION (see mvst_diss_sos_selftest, run by mvst_diss_sos('selftest')):
%    case                          expected                 tolerance
%    (1a) CT L2-gain, A=[-3 0;0 -3], Bw=[0;1], C=[.5 1], Dw=.5
%                                  norm(ss,'inf') via mvst_l2gain_ct   ~1e-2
%    (1b) CT L2-gain, A=[0 1;-2 -1], Bw=[0;1], C=[1 0], Dw=.5
%                                  mvst_l2gain_ct                      ~1e-2
%    (2a) passivity of a positive-real system (A=[0 1;-1 -1],Bw=[0;1],
%         C=[1 0],Dw=.5)          FEASIBLE
%    (2b) non-passive system (Dw=0, C=[0 1])   INFEASIBLE
%    (3)  periodic sampled-output gain, A=[0 1;-1 -1],Bw=[0;1],C=[1 0],Tbar=.3
%                                  mvst_l2gain_sd_hybrid               ~3e-2
%         (the SOS value is the converged infinite-horizon periodic gain; the
%          reference undercounts at a finite horizon: Nper=60 gives 2.053,
%          Nper=200 gives 2.102, SOS gives 2.108.)
%    (4)  Poisson mean-square gain (ex_diss_poisson data)
%                                  mvst_l2gain_ct_poisson              ~3e-2

if nargin == 1 && ischar(sys) && strcmp(sys,'selftest')
    res = mvst_diss_sos_selftest(); return;
end
if nargin < 4 || isempty(mode); mode = 'feas'; end
if nargin < 5; opts = struct(); end
oo = mvst_opts(opts); solver = oo;   % normalized options travel to the solve

% SYS may be a THIN mode/system struct (MVST_MODE / MVST_SYSTEM): if it lacks
% the .A0 field it is routed through MVST_TORICH to the rich naming.  (A thin
% mode carries no disturbance/output data, so the supply blocks default to
% zero unless .Bw/.C/.Dw were attached to the struct before the call.)
if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
n = size(sys.A0,1);
m = infer_m(sys, n);
sys = fill_defaults(sys, n, m);

switch mode
    case 'gamma'
        res = bisect_gamma(sys, supply, dwell, n, m, solver);
    case 'feas'
        res = solve_once(sys, supply, dwell, n, m, supplygamma(supply), true, solver);
    otherwise
        error('mvst:diss_sos:mode','unknown mode "%s"', mode);
end
end

% =========================================================================
function res = bisect_gamma(sys, supply, dwell, n, m, solver)
%BISECT_GAMMA  Smallest gamma for which the L2-gain DLMIs are SOS-feasible.
if ~strcmp(supply.type,'l2gain')
    error('mvst:diss_sos:gammamode','gamma mode requires the l2gain supply.');
end
Dw = sys.Dw;
lo = norm(Dw) + 1e-4;
hi = max(2*lo, 1);
Rbest = [];
grow = 0;
% Feasibility is decided by slack maximization (robust; the pinf/dinf flags
% give false positives near the boundary), exactly as in MVST_STAB_SOS.
while true
    r = solve_once(sys, supply, dwell, n, m, hi, true, solver);
    if r.feas; Rbest = r; break; end
    hi = 2*hi; grow = grow + 1;
    if grow > 30; break; end
end
for it = 1:24
    mid = 0.5*(lo + hi);
    r = solve_once(sys, supply, dwell, n, m, mid, true, solver);
    if r.feas; hi = mid; Rbest = r; else; lo = mid; end
end
res = Rbest;
if isempty(res)
    res = solve_once(sys, supply, dwell, n, m, hi, false, solver);
end
res.gamma = hi;
end

% =========================================================================
function res = solve_once(sys, supply, dwell, n, m, gam, useSlackObj, solver)
%SOLVE_ONCE  Build and solve the dual DLMI program at a fixed supply level.
%   useSlackObj = true  : maximize a common slack (feasibility margin / 'feas')
%   useSlackObj = false : fixed small margin, read solver feasibility ('gamma')
tol   = 1e-7;
Rcap  = 1e4;                                        % storage-magnitude cap
margin = 1e-7;                                       % strict-feasibility margin
alpha = getf(dwell,'alpha',0);
[M, Mk] = supply_blocks(supply, sys, gam, m);        % (M11,M12,M22),(jump)

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

if useSlackObj; slk = sl; else; slk = margin; end

% -------------------------------------------------------------- flow C_M >= 0
if strcmp(dwell.type,'discrete')
    Ld = discrete_cM(sys, P, M);
    prog = matge(prog, Ld, slk);
else
    Pdot = diff(P, th)/T;
    Cf = flow_cM(sys, P, Pdot, alpha, M, m);
    switch dwell.type
        case 'arbitrary'
            prog = matge(prog, Cf, slk);
        case {'periodic','minimum','range'}
            [prog, S1] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n+m n+m], 'symmetric');
            prog = matsos(prog, Cf - slk*eye(n+m) - gt*S1);
            prog = matsos(prog, S1);
    end
    if strcmp(dwell.type,'minimum')                  % frozen tail C_M^0(Tbar)>=0
        Cf0 = flow_cM(sys, subs(P,th,1), zeros(n), alpha, M, m);
        prog = matge(prog, Cf0, slk);
    end
end

% --------------------------------------------------- scheduled-jump D_M >= 0
needJump = ~strcmp(dwell.type,'discrete') && jump_present(sys, dwell);
if needJump
    switch dwell.type
        case 'arbitrary'
            Dj = jump_dM(sys, P, P, Mk);
            prog = matge(prog, Dj, slk);
        case {'periodic','minimum'}
            Dj = jump_dM(sys, subs(P,th,0), subs(P,th,1), Mk);
            prog = matge(prog, Dj, slk);
        case 'range'
            Dj = jump_dM(sys, subs(P,th,0), P, Mk);
            gm = (th - dwell.Tmin/T)*(dwell.Tmax/T - th);
            nj = size(Dj,1);
            [prog, Sj] = sospolymatrixvar(prog, monomials(th,0:2*deg), [nj nj], 'symmetric');
            prog = matsos(prog, Dj - slk*eye(nj) - gm*Sj);
            prog = matsos(prog, Sj);
    end
end

% ---------------- storage nonnegativity 0 <= -P <= Rcap I (scale-bounded)
if isTimer
    [prog, Xlo] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, -P - gt*Xlo);              prog = matsos(prog, Xlo);
    [prog, Xhi] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, Rcap*eye(n) + P - gt*Xhi); prog = matsos(prog, Xhi);
else
    prog = matsos(prog, -P);
    prog = matsos(prog, Rcap*eye(n) + P);
end

if useSlackObj; prog = sossetobj(prog, -sl); end
o = mvst_sosopts(solver);
[prog, info] = mvst_sossolve(prog, o);

Pv = sosgetsol(prog, P);
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
res.Pfun  = @(tau) evalmat(Pv, th, tau/T, n);
res.Sfun  = @(tau) -evalmat(Pv, th, tau/T, n);
res.info  = struct('P',Pv,'T',T,'deg',deg,'alpha',alpha,'solverinfo',info);
end

% =========================================================================
function C = flow_cM(sys, P, Pdot, alpha, M, m)
% C_M(P) = [ Pdot + He(A0^*P) + sum Ac_l^*P Ac_l + Poisson + 2aP + M11 ,
%            P Bw + sum Ac_l^*P Dwc_l + M12 ;
%            *  ,  sum Dwc_l^*P Dwc_l + M22 ]
n  = size(sys.A0,1);
q  = size(sys.Ac,3);  pf = size(sys.Jf,3);
G11 = Pdot + sys.A0'*P + P*sys.A0 + 2*alpha*P + M.M11;
G12 = P*sys.Bw + M.M12;
G22 = M.M22;
for l = 1:q
    G11 = G11 + sys.Ac(:,:,l)'*P*sys.Ac(:,:,l);
    G12 = G12 + sys.Ac(:,:,l)'*P*sys.Dwc(:,:,l);
    G22 = G22 + sys.Dwc(:,:,l)'*P*sys.Dwc(:,:,l);
end
for r = 1:pf
    G11 = G11 + sys.lamf(r)*(sys.Jf(:,:,r)'*P*sys.Jf(:,:,r) - P);
    G12 = G12 + sys.lamf(r)*(sys.Jf(:,:,r)'*P*sys.Bdf(:,:,r));
    G22 = G22 + sys.lamf(r)*(sys.Bdf(:,:,r)'*P*sys.Bdf(:,:,r));
end
C = [ G11, G12 ; G12', G22 ];
C = (C + C')/2;
end

% -------------------------------------------------------------------------
function D = jump_dM(sys, Ppost, Ppre, Mk)
% D_M(P^+,P^-) = [ J^*P^+J - P^- + M11k , J^*P^+ Bd + M12k ;
%                  *  ,  Bd^*P^+ Bd + M22k ]
% If no discrete-disturbance channel (md = 0) the block reduces to the state
% part J^*P^+J - P^- + M11k.
J  = sys.J;  Bd = sys.Bd;  md = size(Bd,2);
G11 = J'*Ppost*J - Ppre + Mk.M11;
if md == 0
    D = (G11 + G11')/2; return;
end
G12 = J'*Ppost*Bd + Mk.M12;
G22 = Bd'*Ppost*Bd + Mk.M22;
D = [ G11, G12 ; G12', G22 ];
D = (D + D')/2;
end

% -------------------------------------------------------------------------
function D = discrete_cM(sys, P, M)
% Pure DT bounded-real: transition A0, noise Ac, input Bw, output (C,Dw).
% D_M(P,P) = [ A0^*P A0 + sum Ac_l^*P Ac_l - P + M11 ,
%              A0^*P Bw + sum Ac_l^*P Dwc_l + M12 ;
%              *  ,  Bw^*P Bw + sum Dwc_l^*P Dwc_l + M22 ]
q = size(sys.Ac,3);
G11 = sys.A0'*P*sys.A0 - P + M.M11;
G12 = sys.A0'*P*sys.Bw + M.M12;
G22 = sys.Bw'*P*sys.Bw + M.M22;
for l = 1:q
    G11 = G11 + sys.Ac(:,:,l)'*P*sys.Ac(:,:,l);
    G12 = G12 + sys.Ac(:,:,l)'*P*sys.Dwc(:,:,l);
    G22 = G22 + sys.Dwc(:,:,l)'*P*sys.Dwc(:,:,l);
end
D = [ G11, G12 ; G12', G22 ];
D = (D + D')/2;
end

% =========================================================================
function [M, Mk] = supply_blocks(supply, sys, gam, m)
%SUPPLY_BLOCKS  Flow (M11,M12,M22) and jump (M11k,M12k,M22k) supply blocks.
n = size(sys.A0,1);
C = sys.C; Dw = sys.Dw; p = size(C,1);
if isempty(C);  C  = zeros(0,n); end
if isempty(Dw); Dw = zeros(p,m); end
nu  = getf(supply,'nu',0);  rho = getf(supply,'rho',0);
switch supply.type
    case 'l2gain'
        M.M11 = -(C'*C);
        M.M12 = -(C'*Dw);
        M.M22 = gam^2*eye(m) - Dw'*Dw;
    case 'passivity'
        % M = [0 C^*; C He(Dw)]  (requires p == m); ISP/OSP refinements.
        M.M11 = zeros(n) - rho*(C'*C);
        M.M12 = C' - rho*(C'*Dw);
        M.M22 = (Dw + Dw') - nu*eye(m) - rho*(Dw'*Dw);
    case 'qsr'
        [Qs,Ss,Rs] = mvst_qsr_supply(supply, p, m, 'flow');
        M.M11 = C'*Qs*C;
        M.M12 = C'*(Qs*Dw + Ss);
        M.M22 = Dw'*Qs*Dw + Dw'*Ss + Ss'*Dw + Rs;
    otherwise
        error('mvst:diss_sos:supply','unknown supply "%s"', supply.type);
end
% Jump supply from the jump output (CJ,DJ) and discrete disturbance (Bd).
CJ = sys.CJ; DJ = sys.DJ; md = size(sys.Bd,2);
if isempty(CJ) && md == 0
    Mk = struct('M11',zeros(n),'M12',zeros(n,0),'M22',zeros(0));
    return;
end
pj = size(CJ,1);
if isempty(CJ); CJ = zeros(0,n); pj = 0; end
if isempty(DJ) || size(DJ,1) ~= pj || size(DJ,2) ~= md; DJ = zeros(pj,md); end
switch supply.type
    case 'l2gain'
        Mk.M11 = -(CJ'*CJ);
        Mk.M12 = -(CJ'*DJ);
        Mk.M22 = gam^2*eye(md) - DJ'*DJ;
    case 'passivity'
        Mk.M11 = zeros(n);
        Mk.M12 = CJ';
        Mk.M22 = (DJ + DJ');
    case 'qsr'
        [Qk,Sk,Rk] = mvst_qsr_supply(supply, pj, md, 'jump');
        Mk.M11 = CJ'*Qk*CJ;
        Mk.M12 = CJ'*(Qk*DJ + Sk);
        Mk.M22 = DJ'*Qk*DJ + DJ'*Sk + Sk'*DJ + Rk;
end
end

% =========================================================================
function tf = jump_present(sys, dwell)
tf = any(strcmp(dwell.type,{'periodic','minimum','range'})) || ...
     (~isempty(sys.J) && ~isequal(sys.J, eye(size(sys.A0,1)))) || ...
     ~isempty(sys.CJ) || size(sys.Bd,2) > 0;
% arbitrary CT with a genuine scheduled jump: keep the constant-P jump test.
if strcmp(dwell.type,'arbitrary') && isequal(sys.J, eye(size(sys.A0,1))) ...
        && isempty(sys.CJ) && size(sys.Bd,2)==0
    tf = false;
end
end

function g = supplygamma(supply)
g = getf(supply,'gamma',1);
end

function m = infer_m(sys, n)
if isfield(sys,'Bw') && ~isempty(sys.Bw); m = size(sys.Bw,2);
elseif isfield(sys,'Dw') && ~isempty(sys.Dw); m = size(sys.Dw,2);
elseif isfield(sys,'C')  && ~isempty(sys.C);  m = 0;
else; m = 1; end
end

% =========================================================================
function prog = matge(prog, EXPR, slk)
% EXPR >= slk*I  as a matrix SOS (constant multiplier: plain LMI).
nn = size(EXPR,1);
prog = sosmatrixineq(prog, EXPR - slk*eye(nn));
end
function prog = matsos(prog, EXPR); prog = sosmatrixineq(prog, EXPR); end

function Mout = evalmat(P, th, theta, n)
if isnumeric(P); Mout = P; else; Mout = double(subs(P, th, theta)); end
Mout = reshape(Mout, n, n);
end

function sys = fill_defaults(sys, n, m)
z = @(a,b,c) zeros(a,b,c);
if ~isfield(sys,'Ac') ||isempty(sys.Ac);  sys.Ac  = z(n,n,0); end
q = size(sys.Ac,3);
if ~isfield(sys,'Dwc')||isempty(sys.Dwc)||size(sys.Dwc,3)~=q; sys.Dwc = z(n,m,q); end
if ~isfield(sys,'Bw') ||isempty(sys.Bw);  sys.Bw  = zeros(n,m); end
if ~isfield(sys,'C')  ;                   sys.C   = zeros(0,n); end
if ~isfield(sys,'Dw') ||isempty(sys.Dw);  sys.Dw  = zeros(size(sys.C,1),m); end
if ~isfield(sys,'Jf') ||isempty(sys.Jf);  sys.Jf  = z(n,n,0); end
if ~isfield(sys,'lamf')||isempty(sys.lamf);sys.lamf= ones(1,size(sys.Jf,3)); end
if ~isfield(sys,'Bdf')||isempty(sys.Bdf)||size(sys.Bdf,3)~=size(sys.Jf,3)
    sys.Bdf = z(n,m,size(sys.Jf,3)); end
if ~isfield(sys,'J')  ||isempty(sys.J);   sys.J   = eye(n); end
if ~isfield(sys,'Bd') ||isempty(sys.Bd);  sys.Bd  = zeros(n,0); end
if ~isfield(sys,'CJ') ;                   sys.CJ  = []; end
if ~isfield(sys,'DJ') ;                   sys.DJ  = []; end
end

function v = getf(s, f, d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end
