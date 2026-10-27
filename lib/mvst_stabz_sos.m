function res = mvst_stabz_sos(sys, dwell, opts)
%MVST_STABZ_SOS  Generic SOSTOOLS state-feedback stabilization synthesis for
%  linear impulsive Hermitian matrix-valued systems, covering every dwell-time
%  type of the Stabilization chapter.  Each stabilization example formats its
%  data into SYS, picks a dwell-time type in DWELL, and this function solves the
%  synthesis condition (the flow / jump / tail linear matrix inequalities of
%  Theorems th:stab:stabz, th:stab:stabzdwell, th:stab:stabzMDT,
%  th:stab:stabzRDT, th:stab:stabzpoisson) as a matrix sum-of-squares program
%  in the timer with SOSTOOLS + SeDuMi.
%
%  Single flow control u = K x and, where present, single spontaneous-jump
%  control uf = Kf x and single scheduled-jump control uJ = KJ x; the input
%  matrices act as selectors.
%
%  SYS fields (any omitted field defaults to "none"):
%    .A0  n x n         flow drift
%    .B0  n x m         flow input selector
%    .Ac  n x n x q     flow noise channels (Brownian OR Levy: pass sqrt(c_l) A_l,
%                       and fold the Levy mean drift sum_l mu_l A_l into A0)
%    .Bc  n x m x q     input-multiplicative noise selectors
%    .Jf  n x n x pf    SPONTANEOUS (Poisson) jump maps, rate-weighted into flow
%    .BJf n x m x pf    spontaneous-jump input selectors
%    .lamf 1 x pf       Poisson rates of the spontaneous jumps
%    .J   n x n x p     SCHEDULED (dwell-time) jump maps (p>1 => noisy jumps)
%    .BJ  n x m x p     scheduled-jump input selectors
%
%  DWELL:
%    .type  'arbitrary' | 'discrete' | 'periodic' | 'minimum' | 'range'
%           ('arbitrary' = constant certificate / CT persistent flowing;
%            a purely spontaneous-Poisson system is 'arbitrary' with .Jf/.lamf
%            and no scheduled .J; 'discrete' = persistent-jumping DT synthesis.)
%    .Tbar  scalar (periodic, minimum) ;  .Tmin,.Tmax (range)
%    .deg   timer polynomial degree (default 2; forced 0 for arbitrary/discrete)
%    .alpha,.rho  guaranteed hybrid rate (defaults 0, 1), at the COVARIANCE
%           level: the certified state decays as X ~ e^{-alpha t} rho^kappa,
%           so the flow (1,1) block carries +alpha Y and the jump (1,1) block
%           is scaled by rho (NOT 2 alpha / rho^2); reported rates are the
%           covariance rates directly, with no halving or rooting.
%
%  RES: .feas, .slack, .type, .Kfun (@(tau)->K), .Kf (spontaneous-jump gain),
%       .KJ (scheduled-jump gain), .info.
%
%  OPTS (struct or MVST.OPTIONS, optional): .solver selects the SDP backend
%  ('sedumi' (default) | 'sdpt3' | 'mosek' | 'nosdp'), normalized by MVST_OPTS.
%
%  SYS may also be a THIN mode/system struct (MVST_MODE / MVST_SYSTEM): if it
%  lacks the .A0 field it is routed through MVST_TORICH to the rich naming.

if nargin < 3; opts = struct(); end
oo = mvst_opts(opts);
if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
n = size(sys.A0,1);
if ~isfield(sys,'B0')||isempty(sys.B0); sys.B0 = zeros(n,0); end
m = size(sys.B0,2);
sys = fill_defaults(sys, n, m);
q = size(sys.Ac,3);  pf = size(sys.Jf,3);  p = size(sys.J,3);
alpha = getf(dwell,'alpha',0);  rho = getf(dwell,'rho',1);  tol = 1e-7;
% Persistent-flowing / persistent-jumping switches (th:stab:relax /
% th:stab:persjump): the "persistent" side is only required non-expansive
% (no strict slack), so all decay is carried by the other side.
persist = getf(dwell,'persist','none');
if strcmp(persist,'flowing') && ~isfield(dwell,'rho');   rho   = 1; end   % jumps non-expansive
if strcmp(persist,'jumping') && ~isfield(dwell,'alpha'); alpha = 0; end   % flow  non-expansive
flowSlack = double(~strcmp(persist,'jumping'));    % flow strict unless persistent jumping
jumpSlack = double(~strcmp(persist,'flowing'));    % jump strict unless persistent flowing
% Levy flow noise enters the second moment only through its variance rate:
% the caller passes Ac = sqrt(c_l) A_l (Brownian is c=1) and folds the Levy
% mean-rate drift sum_l mu_l A_l into A0.  Levy/Poisson/noisy JUMPS are the
% Jf/lamf (spontaneous) and multi-channel J (scheduled) data.

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
[prog, Y] = sospolymatrixvar(prog, mons, [n n], 'symmetric');
[prog, U] = sospolymatrixvar(prog, mons, [m n]);
Uf = zeros(m,n); UJ = zeros(m,n);
if pf > 0; [prog, Uf] = sospolymatrixvar(prog, mons, [m n]); end
if p  > 0; [prog, UJ] = sospolymatrixvar(prog, mons, [m n]); end
gt = th*(1-th);

% ------------------------------------------------------------- flow condition
if strcmp(dwell.type,'discrete')
    Phi = discrete_block(sys, Y, U, rho, n, q);
    prog = matsos(prog, -Phi - (flowSlack*sl)*eye(size(Phi,1)));
else
    Yp = diff(Y, th);
    Phif = flow_block(sys, Y, Yp/T, U, Uf, alpha, n, q, pf);
    fsl = flowSlack*sl;
    switch dwell.type
        case 'arbitrary'
            prog = matsos(prog, -Phif - fsl*eye(size(Phif,1)));
        case {'periodic','minimum','range'}
            [prog, S1] = sospolymatrixvar(prog, monomials(th,0:2*deg), size(Phif), 'symmetric');
            prog = matsos(prog, -Phif - fsl*eye(size(Phif,1)) - gt*S1);
            prog = matsos(prog, S1);
    end
    if strcmp(dwell.type,'minimum')                    % frozen tail at theta=1
        PhifT = flow_block(sys, subs(Y,th,1), zeros(n), subs(U,th,1), subs(Uf,th,1), alpha, n, q, pf);
        prog = matsos(prog, -PhifT - fsl*eye(size(PhifT,1)));
    end
end

% -------------------------------------------------- scheduled-jump condition
if p > 0 && ~strcmp(dwell.type,'discrete')
    jsl = jumpSlack*sl;
    switch dwell.type
        case 'arbitrary'
            Phij = jump_block(sys, Y, Y, UJ, rho, n, p);
            prog = matsos(prog, -Phij - jsl*eye(size(Phij,1)));
        case {'periodic','minimum'}
            Phij = jump_block(sys, subs(Y,th,1), subs(Y,th,0), subs(UJ,th,1), rho, n, p);
            prog = matsos(prog, -Phij - jsl*eye(size(Phij,1)));
        case 'range'
            Phij = jump_block(sys, Y, subs(Y,th,0), UJ, rho, n, p);
            gm = (th - dwell.Tmin/T)*(dwell.Tmax/T - th);
            [prog, Sj] = sospolymatrixvar(prog, monomials(th,0:2*deg), size(Phij), 'symmetric');
            prog = matsos(prog, -Phij - jsl*eye(size(Phij,1)) - gm*Sj);
            prog = matsos(prog, Sj);
    end
end

% ------------------------------------------- normalization Y <= I (scale-fix)
if isTimer
    [prog, Xi] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, eye(n) - Y - gt*Xi);  prog = matsos(prog, Xi);
else
    prog = matsos(prog, eye(n) - Y);
end

prog = sossetobj(prog, -sl);
o = mvst_sosopts(oo);
prog = mvst_sossolve(prog, o);

sval = double(sosgetsol(prog, sl));
Yv  = sosgetsol(prog, Y);  Uv = sosgetsol(prog, U);
Ufv = zeros(m,n); if pf>0; Ufv = sosgetsol(prog, Uf); end
UJv = zeros(m,n); if p >0; UJv = sosgetsol(prog, UJ); end

res.type  = dwell.type;
res.slack = sval;
res.feas  = sval > tol;
res.Kfun  = @(tau) evalmat(Uv,th,tau/T,m,n) / evalmat(Yv,th,tau/T,n,n);
res.Kf    = zeros(m,n); if pf>0; res.Kf = evalmat(Ufv,th,0,m,n)/evalmat(Yv,th,0,n,n); end
res.KJ    = zeros(m,n); if p >0; res.KJ = evalmat(UJv,th,1,m,n)/evalmat(Yv,th,1,n,n); end
res.info  = struct('Y',Yv,'U',Uv,'Uf',Ufv,'UJ',UJv,'T',T,'deg',deg);
end

% ==========================================================================
function P = flow_block(sys, Y, Ydot, U, Uf, alpha, n, q, pf)
% [ -Ydot+He(A0Y+B0U)+aY-(sum lamf)Y , row_l[V_l'] , row_r[lamf_r Zf_r'] ;
%   col_l[V_l] , -I_q kron Y , 0 ;
%   col_r[lamf_r Zf_r] , 0 , -diag_r(lamf_r Y) ]   (covariance-level shift +aY)
AY = sys.A0*Y + sys.B0*U; lamf = sys.lamf;
P11 = -Ydot + (AY+AY.') + alpha*Y - sum(lamf)*Y;
Vrow=[];Vcol=[]; for l=1:q; V=sys.Ac(:,:,l)*Y+sys.Bc(:,:,l)*U; Vrow=[Vrow,V.']; Vcol=[Vcol;V]; end %#ok<AGROW>
Zrow=[];Zcol=[]; for r=1:pf; Z=sys.Jf(:,:,r)*Y+sys.BJf(:,:,r)*Uf; Zrow=[Zrow,lamf(r)*Z.']; Zcol=[Zcol;lamf(r)*Z]; end %#ok<AGROW>
P = [ P11, Vrow ; Vcol, -blkdiagrep(Y,q) ];
if pf > 0
    D = blkdiaglam(Y, lamf);
    zt = zeros(q*n, pf*n);
    P = [ P11, Vrow, Zrow ; Vcol, -blkdiagrep(Y,q), zt ; Zcol, zt.', -D ];
end
end

function P = jump_block(sys, Yend, Y0, UJ, rho, n, p)
% [ -rho Yend , row_r[Z_r'] ; col_r[Z_r] , -I_p kron Y0 ], Z_r=J_r Yend+BJ_r UJ
% (covariance-level jump scaling rho)
Zrow=[];Zcol=[];
for r=1:p; Zr=sys.J(:,:,r)*Yend+sys.BJ(:,:,r)*UJ; Zrow=[Zrow,Zr.']; Zcol=[Zcol;Zr]; end %#ok<AGROW>
P = [ -rho*Yend, Zrow ; Zcol, -blkdiagrep(Y0,p) ];
end

function P = discrete_block(sys, Y, U, rho, n, q)
% DT persistent-jumping synthesis: [ -rho Y , row_{l=0..q}[M_l'] ; col , -I kron Y ]
% (covariance-level jump scaling rho)
M0 = sys.A0*Y + sys.B0*U; Mrow = M0.'; Mcol = M0;
for l=1:q; Ml=sys.Ac(:,:,l)*Y+sys.Bc(:,:,l)*U; Mrow=[Mrow, Ml.']; Mcol=[Mcol; Ml]; end %#ok<AGROW>
P = [ -rho*Y, Mrow ; Mcol, -blkdiagrep(Y, 1+q) ];
end

function B = blkdiagrep(Y, k)
n = size(Y,1); if k==0; B = zeros(0); return; end
B = [];
for i = 1:k
    row = [];
    for j = 1:k; if i==j; row=[row,Y]; else; row=[row,zeros(n)]; end; end %#ok<AGROW>
    B = [B; row]; %#ok<AGROW>
end
end

function B = blkdiaglam(Y, lam)
n = size(Y,1); k = numel(lam); B = [];
for i = 1:k
    row = [];
    for j = 1:k; if i==j; row=[row,lam(i)*Y]; else; row=[row,zeros(n)]; end; end %#ok<AGROW>
    B = [B; row]; %#ok<AGROW>
end
end

function prog = matsos(prog, EXPR); prog = sosmatrixineq(prog, EXPR); end

function M = evalmat(P, th, theta, r, c)
if isnumeric(P); M = P; else; M = double(subs(P, th, theta)); end
M = reshape(M, r, c);
end

function sys = fill_defaults(sys, n, m)
if ~isfield(sys,'Ac')||isempty(sys.Ac);  sys.Ac =zeros(n,n,0); end
if ~isfield(sys,'Bc')||isempty(sys.Bc)||size(sys.Bc,3)~=size(sys.Ac,3); sys.Bc=zeros(n,m,size(sys.Ac,3)); end
if ~isfield(sys,'Jf')||isempty(sys.Jf);   sys.Jf =zeros(n,n,0); end
if ~isfield(sys,'BJf')||isempty(sys.BJf)||size(sys.BJf,3)~=size(sys.Jf,3); sys.BJf=zeros(n,m,size(sys.Jf,3)); end
if ~isfield(sys,'lamf')||isempty(sys.lamf); sys.lamf=ones(1,size(sys.Jf,3)); end
if ~isfield(sys,'J') ||isempty(sys.J);    sys.J  =zeros(n,n,0); end
if ~isfield(sys,'BJ')||isempty(sys.BJ)||size(sys.BJ,3)~=size(sys.J,3);  sys.BJ =zeros(n,m,size(sys.J,3)); end
end

function v = getf(s, f, d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end
