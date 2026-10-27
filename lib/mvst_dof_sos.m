function res = mvst_dof_sos(sys, dwell, opts)
%MVST_DOF_SOS  Generic SOSTOOLS timer-dependent output-feedback bracket
%  SYNTHESIS for the scalar dwell-time output-feedback problem of the Dynamic
%  Output Feedback chapter, with controller reconstruction.  Each DOF example
%  formats the plant into SYS, picks a dwell-time regime in DWELL, and this
%  function searches (SOSTOOLS + SeDuMi, polynomial pvar mode) for the
%  polynomial certificate corners X(tau), Y(tau) and the transformed controller
%  data (flow: Chat,Ahat,Bhat,Dhat; jump: ChatJ,AhatJ,BhatJ,DhatJ), enforcing
%  on the WHOLE dwell interval the convex sufficient bracket of Theorem
%  th:dof:dwell and Proposition prop:dof:bracket:
%
%    flow:  [ Psi_f(tau)   L(tau)^*                 ]  <= -alpha I,
%           [ L(tau)      -(2 mu I - mu^2 Y(tau))   ]
%    jump:  [ -Omega(X^-,Y^-)   Xi_J^*              ]  <= -alpha I,
%           [  Xi_J            -Omega(X^+,Y^+)      ]
%    coupling: Omega(X(tau),Y(tau)) = [X I; I Y] >= eps I,
%
%  for three dwell-time regimes (Theorem th:dof:dwell (b),(c),(d)):
%    periodic  Delta_k = T       : flow on [0,T],   jump at tau=T;
%    minimum   Delta_k >= T       : + frozen-tail flow Psi_f(T)|_{Xdot=Ydot=0};
%    range     Delta_k in [Tm,TM] : flow on [0,TM], jump for theta in [Tm,TM].
%
%  Everything is polynomial in the timer (no inverse Y(tau)^{-1}); the interval
%  constraints are matrix Positivstellensatz certificates (an SOS-matrix
%  multiplier times tau(T-tau), no gridding).  After a feasible solve the actual
%  controller K(tau)=(A_c,B_c,C_c,D_c) is reconstructed by the change of
%  variables (Lemma lem:dof:recovery), with the derivative correction
%  Ahat -> Ahat + Y*Xdot applied first.  SOSTOOLS runs in pvar mode: the
%  symbolic engine crashes on the larger jump matrix, the pvar backend does not.
%
%  SYS fields (scalar plant; one state-multiplicative flow channel, one reset):
%    .a0,.a1   flow drift A0(tau) = a0 + a1 tau (timer-affine, may be unstable)
%    .sig      multiplicative-noise intensity
%    .mu       fixed scalar bracket parameter
%    .B0,.Cy   flow control input / measurement output
%    .J0       reset map
%    .BJ,.CJ   jump control input / jump measurement output
%
%  DWELL:
%    .type  'periodic' | 'minimum' | 'range'
%    .T     dwell length (periodic, minimum; also the reset instant)
%    .Tmin,.Tmax  dwell bounds (range)
%    .deg   timer polynomial degree (default 2)
%
%  OPTS (struct or MVST.OPTIONS, optional; normalized by MVST_OPTS):
%    .alpha  strict bracket slack   (default 1e-3)
%    .eps    coupling slack         (default 1e-3)
%    .solver SDP backend 'sedumi' (default) | 'sdpt3' | 'mosek'
%
%  RES: .feas, .info (SeDuMi info), .maxM, .minOm (interval verification of the
%       solved certificate over [0,Tflow]), .Mf/@(tau) flow bracket, .Omf/@(tau)
%       coupling, .Xf/@(tau) X, .Yf/@(tau) Y, .ctrl (solved transformed
%       controller polynomials Ahat,Bhat,Chat,Dhat + jump AhatJ,BhatJ,ChatJ,
%       DhatJ), .sym (raw solved pvars), and .K (reconstructed controller
%       A_c,B_c,C_c,D_c + jump gains AcJ,BcJ,CcJ,DcJ; empty when infeasible).

if nargin < 3 || isempty(opts); opts = struct(); end
opts = mvst_opts(opts);   % normalize + validate .solver ('sedumi' default)

% ---- assemble the internal parameter struct (plant + regime + slacks) ----
P      = sys;
P.T    = getf(dwell,'T',    1.0);
P.Tmin = getf(dwell,'Tmin', 0.5);
P.Tmax = getf(dwell,'Tmax', 1.5);
P.d    = getf(dwell,'deg',  2);
P.alpha= getf(opts,'alpha', 1e-3);
P.eps  = getf(opts,'eps',   1e-3);
solver = opts;   % normalized options (solver, reltol, verbosity) travel to the solve
regime = dwell.type;

out = solveDwell(regime, P, solver);

res.feas  = out.solved;
res.info  = out.info;
res.maxM  = out.maxM;
res.minOm = out.minOm;
res.Mf    = out.Mf;
res.Omf   = out.Omf;
res.Xf    = out.Xf;
res.Yf    = out.Yf;
res.sym   = out.sym;
res.ctrl  = struct('Ahat',out.sym.Ah,'Bhat',out.sym.Bh, ...
                   'Chat',out.sym.Ch,'Dhat',out.sym.Dh, ...
                   'AhatJ',out.sym.AhJ,'BhatJ',out.sym.BhJ, ...
                   'ChatJ',out.sym.ChJ,'DhatJ',out.sym.DhJ);

% ---- controller reconstruction (only meaningful once feasible) ----
if out.solved
    res.K = reconstruct(out, P);
else
    res.K = [];
end
end

% =========================================================================
function out = solveDwell(regime, P, solver)
%SOLVEDWELL  Solve the bracket synthesis for one dwell-time regime by SOSTOOLS
%  (polynomial pvar mode; boundary values are linear combinations of the
%  polynomial coefficients, never subs on the program variable).
pvar tau
Tflow = P.T; if strcmp(regime,'range'), Tflow = P.Tmax; end
prog = sosprogram(tau);
Z0 = monomials(tau, 0);
[prog, X,  cX ] = polyvar(prog, tau, P.d, Z0);
[prog, Y,  cY ] = polyvar(prog, tau, P.d, Z0);
[prog, Ch, cCh] = polyvar(prog, tau, P.d, Z0);
[prog, Ah, cAh] = polyvar(prog, tau, P.d, Z0);
[prog, Bh] = sospolymatrixvar(prog, Z0, [1 1]);
[prog, Dh] = sospolymatrixvar(prog, Z0, [1 1]);
[prog, ChJ] = sospolymatrixvar(prog, Z0, [1 1]);
[prog, AhJ] = sospolymatrixvar(prog, Z0, [1 1]);
[prog, BhJ] = sospolymatrixvar(prog, Z0, [1 1]);
[prog, DhJ] = sospolymatrixvar(prog, Z0, [1 1]);

a  = P.a0 + P.a1*tau;
Xd = polyval_der(cX, tau); Yd = polyval_der(cY, tau);

% ---- flow bracket M(tau) <= -alpha I and coupling on [0,Tflow] ----
Psi11 = -Xd + 2*(a*X + P.B0*Ch);
Psi12 =  a + P.B0*Dh*P.Cy + Ah;
Psi22 =  Yd + 2*(Y*a + Bh*P.Cy);
M  = [Psi11 Psi12 P.sig*X; Psi12 Psi22 P.sig; P.sig*X P.sig -(2*P.mu-P.mu^2*Y)];
Om = [X 1; 1 Y];
[prog, Sm] = sossosmatrixvar(prog, monomials(tau,0:1), [3 3]);
[prog, So] = sossosmatrixvar(prog, monomials(tau,0:1), [2 2]);
prog = sosmatrixineq(prog, -M  - P.alpha*eye(3) - Sm*tau*(Tflow-tau));
prog = sosmatrixineq(prog,  Om - P.eps*eye(2)   - So*tau*(Tflow-tau));

% ---- jump condition (regime-dependent) ----
X0 = polyval_at(cX,0); Y0 = polyval_at(cY,0);
switch regime
  case {'periodic','minimum'}
    Xm = polyval_at(cX,P.T); Ym = polyval_at(cY,P.T);
    XiJ = [P.J0*Xm+P.BJ*ChJ, P.J0+P.BJ*DhJ*P.CJ; AhJ, Y0*P.J0+BhJ*P.CJ];
    Jmp = [-[Xm 1;1 Ym], XiJ.'; XiJ, -[X0 1;1 Y0]];
    [prog, Sj] = sossosmatrixvar(prog, monomials(tau,0:1), [4 4]);
    prog = sosmatrixineq(prog, -Jmp - P.alpha*eye(4) - Sj*tau*(P.T-tau));
    if strcmp(regime,'minimum')
      ChT=polyval_at(cCh,P.T); AhT=polyval_at(cAh,P.T); aT=P.a0+P.a1*P.T;
      Mf = [2*aT*Xm+2*P.B0*ChT, aT+P.B0*Dh*P.Cy+AhT, P.sig*Xm; ...
            aT+P.B0*Dh*P.Cy+AhT, 2*Ym*aT+2*Bh*P.Cy,  P.sig; ...
            P.sig*Xm,            P.sig,               -(2*P.mu-P.mu^2*Ym)];
      [prog, Sf] = sossosmatrixvar(prog, monomials(tau,0:1), [3 3]);
      prog = sosmatrixineq(prog, -Mf - P.alpha*eye(3) - Sf*tau*(P.T-tau));
    end
  case 'range'
    XiJ = [P.J0*X+P.BJ*ChJ, P.J0+P.BJ*DhJ*P.CJ; AhJ, Y0*P.J0+BhJ*P.CJ];
    Jmp = [-[X 1;1 Y], XiJ.'; XiJ, -[X0 1;1 Y0]];
    [prog, Sj] = sossosmatrixvar(prog, monomials(tau,0:1), [4 4]);
    prog = sosmatrixineq(prog, -Jmp - P.alpha*eye(4) - Sj*(tau-P.Tmin)*(P.Tmax-tau));
end

opt = mvst_sosopts(solver);
[prog, info] = mvst_sossolve(prog, opt);
out.solved = (info.pinf==0 && info.dinf==0 && info.numerr<2);
out.info = info;

% ---- recover solved polynomials + interval verification ----
S = @(v) sosgetsol(prog, v);
Xs=S(X); Ys=S(Y); Chs=S(Ch); Ahs=S(Ah); Bhs=S(Bh); Dhs=S(Dh);
out.sym = struct('tau',tau,'X',Xs,'Y',Ys,'Ch',Chs,'Ah',Ahs,'Bh',Bhs,'Dh',Dhs, ...
                 'ChJ',S(ChJ),'AhJ',S(AhJ),'BhJ',S(BhJ),'DhJ',S(DhJ));
Xds=polyval_der(cellfun(@(c) S(c), cX, 'uni',0), tau);
Yds=polyval_der(cellfun(@(c) S(c), cY, 'uni',0), tau);
P11=-Xds+2*(a*Xs+P.B0*Chs); P12=a+P.B0*Dhs*P.Cy+Ahs; P22=Yds+2*(Ys*a+Bhs*P.Cy);
Msol=[P11 P12 P.sig*Xs; P12 P22 P.sig; P.sig*Xs P.sig -(2*P.mu-P.mu^2*Ys)];
Omsol=[Xs 1;1 Ys];
out.Mf  = @(t) evalp(Msol , tau, t);
out.Omf = @(t) evalp(Omsol, tau, t);
out.Xf  = @(t) evalp(Xs, tau, t);
out.Yf  = @(t) evalp(Ys, tau, t);
tt = linspace(0, Tflow, 150);
out.maxM  = max(arrayfun(@(t) max(real(eig(out.Mf(t)))),  tt));
out.minOm = min(arrayfun(@(t) min(real(eig(out.Omf(t)))), tt));
end

% =========================================================================
function K = reconstruct(out, P)
%RECONSTRUCT  Controller K(tau)=(A_c,B_c,C_c,D_c) from the solved (X,Y,hats).
%  B_c, A_c divide by N(tau)=I-Y X, hence rational: returned as numeric-eval
%  function handles of tau.
s = out.sym; tau = s.tau;
ev = @(p,t) evalp(p, tau, t);
Dc = ev(s.Dh,0);
Xdp = diff(s.X, tau);
Xf=@(t) ev(s.X,t); Yf=@(t) ev(s.Y,t); Chf=@(t) ev(s.Ch,t); Ahf=@(t) ev(s.Ah,t);
Xdf=@(t) ev(Xdp,t); Bhc = ev(s.Bh,0);
N  = @(t) 1 - Yf(t).*Xf(t);
Cc = @(t) Chf(t) - Dc.*P.Cy.*Xf(t);
Bc = @(t) (Bhc - Yf(t).*P.B0.*Dc)./N(t);
Ac = @(t) (Ahf(t) + Yf(t).*Xdf(t) - N(t).*Bc(t).*P.Cy.*Xf(t) ...
           - Yf(t).*P.B0.*Cc(t) - Yf(t).*(P.a0+P.a1.*t+P.B0.*Dc.*P.Cy).*Xf(t))./N(t);
K = struct('tau',tau,'Dc',Dc,'Cc',Cc,'Bc',Bc,'Ac',Ac);
% jump gains at the reset (F,G,H)=(J0,BJ,CJ), X^-=X(T), X^+=X(0)
X0=Xf(0); Y0=Yf(0); XT=Xf(P.T); Np=1-Y0*X0;
K.DcJ = ev(s.DhJ,0);
K.CcJ = ev(s.ChJ,0) - K.DcJ*P.CJ*XT;
K.BcJ = (ev(s.BhJ,0) - Y0*P.BJ*K.DcJ)/Np;
K.AcJ = (ev(s.AhJ,0) - Np*K.BcJ*P.CJ*XT - Y0*P.BJ*K.CcJ ...
         - Y0*(P.J0+P.BJ*K.DcJ*P.CJ)*XT)/Np;
end

% =========================================================================
function [prog, p, c] = polyvar(prog, tau, d, Z0)
%POLYVAR  Degree-d scalar polynomial in tau from explicit coefficient vars.
c = cell(1,d+1); p = 0;
for k = 0:d
    [prog, ck] = sospolymatrixvar(prog, Z0, [1 1]);
    c{k+1} = ck; p = p + ck*tau^k;
end
end

function pd = polyval_der(c, tau)
%POLYVAL_DER  Derivative of the polynomial with coefficient cell c.
pd = 0;
for k = 1:numel(c)-1, pd = pd + k*c{k+1}*tau^(k-1); end
end

function pv = polyval_at(c, t)
%POLYVAL_AT  Value of the polynomial with coefficient cell c at tau=t.
pv = 0;
for k = 0:numel(c)-1, pv = pv + c{k+1}*t^k; end
end

function A = evalp(M, tau, t)
%EVALP  Numeric value of a pvar (matrix) polynomial M at tau=t.
A = double(subs(M, tau, t));
end

function v = getf(s, f, d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end
