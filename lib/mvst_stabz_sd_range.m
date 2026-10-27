function res = mvst_stabz_sd_range(sys, opts)
%MVST_STABZ_SD_RANGE  Aperiodic (range) sampled-data refresh-gain synthesis for
%  a held-input matrix-valued plant, as a matrix sum-of-squares program in the
%  timer.  Implements the two range-sampling designs of the Stabilization
%  chapter (Section sec:stabz:sd:range) at the COVARIANCE level:
%
%    'scheduled'  Proposition prop:stabz:sdrange (period-scheduled refresh gain
%                 K_J(theta), the general range dwell-time synthesis).
%    'constant'   Corollary cor:stabz:sdrange:const (schedule-free constant
%                 refresh gain K_J; the flat-certificate convex restriction).
%    'slack'      Finsler/slack constant-gain synthesis with a PERIOD-DEPENDENT
%                 certificate (the LPV-paper dilation device): a schedule-free
%                 constant refresh gain K_J = U G^{-1} carried by a single
%                 constant slack G and a constant U, while the timer certificate
%                 Y(theta) stays FREE (period-dependent) across the window.  A
%                 strict convex tightening of the constant design that lifts its
%                 frozen-certificate conservatism (see (C) below).
%
%  HELD-INPUT LIFT (Section sec:stabz:sd:lift).  The plant carries a
%  zero-order-held command u.  Adjoining it to the plant state x gives the
%  lifted state z = col(x,u) of dimension N = n + m, with flow data
%     tilde A   = [ A   B_u ; 0 0 ],   tilde A_l = [ A_l  B_{u,l} ; 0 0 ],
%  (the held command does not move between samples, so the lifted flow carries
%  NO flow gain), and the open controlled reset
%     z^+ = J_0 z + B_J u_J,   J_0 = [ I_n 0 ; 0 0 ],   B_J = [ 0 ; I_m ],
%  with jump feedback u_J = K_J z.  All synthesis freedom sits in the single
%  jump channel (p = 0), pair (J_0, B_J).  This backend builds the lift from the
%  plant data internally.
%
%  SLACK CONVENTION.  The held command has ZERO flow dynamics (its lifted block
%  is zero), so the flow cannot contract between samples: the flow (and frozen-
%  tail) LMIs are imposed only as NON-EXPANSIVE (no strict margin) and the whole
%  decay is carried by the STRICT refresh (jump) LMI.  Feasibility is decided by
%  maximizing the jump slack s; feasible iff s > 0.  A strict margin on the flow
%  would be structurally infeasible (the frozen-tail (u,u) block is identically
%  zero), so this jump-carried decay is mandatory, not a modeling choice.
%
%  THE TWO PROGRAMS (normalized timer th in [0,1]; real timer tau = T th,
%  real-time derivative Ydot = (dY/dth)/T):
%
%   (A) 'scheduled', T = Tmax.  Certificate Y(tau) on [0,Tmax], refresh matrix
%       U_J(theta) polynomial in theta; the flow LMI holds on (0,Tmax) and the
%       jump LMI holds for EVERY admissible period theta in [Tmin,Tmax]:
%         flow:  [ -Ydot+He(tilde A Y) , row_l (tilde A_l Y)^* ;
%                   *                   , -I_q kron Y            ] <= 0,
%         jump:  [ -Y(theta) , (J_0 Y(theta)+B_J U_J(theta))^* ;
%                   *         , -Y(0)                            ] <= -s I.
%       Recovered gain  K_J(theta) = U_J(theta) Y(theta)^{-1}.
%
%   (B) 'constant', T = Tmin.  Certificate Y(tau) on [0,Tmin], FROZEN (Ydot=0)
%       on the window [Tmin,Tmax] so Y(theta) = Y(Tmin), and a CONSTANT refresh
%       matrix U_J.  The flow LMI carries the derivative on (0,Tmin); the
%       frozen-tail flow LMI (derivative deleted) holds on [Tmin,Tmax]; the
%       single jump LMI closes the window:
%         flow  (0,Tmin):     [ -Ydot+He(tilde A Y) , row_l (tilde A_l Y)^* ;
%                                *                   , -I_q kron Y ] <= 0,
%         frozen [Tmin,Tmax]: [  He(tilde A Y(Tmin)) , row_l (tilde A_l Y(Tmin))^* ;
%                                *                    , -I_q kron Y(Tmin) ] <= 0,
%         jump:               [ -Y(Tmin) , (J_0 Y(Tmin)+B_J U_J)^* ;
%                                *        , -Y(0) ] <= -s I.
%       Recovered gain  K_J = U_J Y(Tmin)^{-1}  (independent of theta).
%       (The frozen-tail block is time-invariant for the LTI held-input plant,
%        so a single LMI covers the whole window.)
%
%   (C) 'slack', T = Tmax.  Certificate Y(tau) FREE on [0,Tmax] (period-
%       dependent, exactly the scheduled certificate), plus ONE CONSTANT square
%       slack Ys (N x N, NOT symmetric) and a CONSTANT U (m x N).  The flow LMI
%       is the scheduled one on (0,Tmax); the jump LMI is the author's DILATED
%       (Finsler) constant-gain family, imposed for EVERY period theta in
%       [Tmin,Tmax].  With He[W] := W + W^* and Z := J_0 Ys + B_J U (= (J_0 +
%       B_J K_J) Ys):
%         jump:  [ -rho Y(0) + He[Z] ,   Z - Ys^*           ;
%                  (Z - Ys^*)^*       ,   Y(theta) - He[Ys] ] <= -s I.
%       Recovered CONSTANT gain  K_J = U Ys^{-1}  (independent of theta).
%
%       PROVENANCE.  Range/covariance adaptation of the dilated timer-independent
%       synthesis of Briat, Nonlinear Analysis: Hybrid Systems 41 (2021) 101040,
%       Theorem 18 (built on the dilated Theorem 7): the two Finsler slacks of
%       the analysis dilation are collapsed to ONE constant matrix, the whole
%       block is congruence-transformed by diag(Ys,Ys) (Ys = X^{-1}), and the
%       gain is recovered through the change of variables U = K_J Ys.
%
%       SUFFICIENCY (implication chain; the book states it).  The scheduled jump
%       block is equivalent, via a Schur complement on its (1,1) block, to the
%       COVARIANCE jump contraction
%         (J_0 + B_J K_J) Y(theta) (J_0 + B_J K_J)^* <= rho Y(0)   (all theta).
%       Pre- and post-multiply the dilated block by [ I , J_cl ] and its adjoint
%       [ I ; J_cl^* ], where J_cl = J_0 + B_J K_J.  Every slack term telescopes:
%       the He[Z], the +-Z, and the -+Ys^* / He[Ys] contributions cancel in
%       pairs (Z = J_cl Ys, so J_cl Ys and Ys^* J_cl^* meet their negatives),
%       leaving exactly
%         J_cl Y(theta) J_cl^* - rho Y(0) <= 0,
%       the covariance contraction, with K_J = U Ys^{-1} read off from Z.  This
%       is the Finsler DIRECTION (dilated feasibility => original jump LMI).  The
%       strict margin survives the congruence:
%         [ I , J_cl ] (block) [ I ; J_cl^* ] <= -s (I + J_cl J_cl^*) <= -s I.
%       At a SINGLE theta the dilation is LOSSLESS (Finsler equivalence, slack
%       Ys free), so a degenerate window is no more conservative than mode
%       'constant'; a POSITIVE window is certified with a period-dependent
%       Y(theta) yet a genuinely CONSTANT gain, removing the frozen-certificate
%       conservatism of 'constant'.
%       NOTE.  The slack enters as He[Ys] in the (2,2) corner and pairs the DATA
%       map J_0 Ys + B_J U with the transformed certificate.  This exact
%       arrangement (NOT a naive [ -Y(theta) , * ; * , -G-G^*+Y(0) ] block) is
%       what makes the slack cancel; the alternative placement does not.
%
%  SYS fields (the PLANT data, before lifting):
%    .A   n x n        plant drift A                       (required)
%    .Bu  n x m        plant input selector B_u            (required)
%    .Ac  n x n x q    plant diffusion drift channels A_l  (optional, default none)
%    .Bc  n x m x q    plant diffusion input selectors B_{u,l} (optional, zeros)
%
%  OPTS fields:
%    .Tmin, .Tmax  sampling-period range, 0 < Tmin <= Tmax   (required)
%    .mode   'scheduled' | 'constant' | 'slack'               (default 'scheduled')
%    .deg    timer polynomial degree of Y (and of U_J in the scheduled mode;
%            U_J is constant in the constant mode; in the slack mode Ys and U
%            are constant and Y stays degree deg)             (default 4)
%    .slack  requested strict margin s subtracted on every block (default 0:
%            the objective then maximizes the achieved slack)
%    .alpha  covariance-level flow rate: the flow (1,1) block carries +alpha Y
%            (default 0; base feasibility program needs no rate)
%    .rho    covariance-level jump contraction: the jump (1,1) block is scaled
%            by rho (default 1)
%    .solver SDP backend for SOSTOOLS ('sedumi' (default) | 'sdpt3' | 'mosek' | 'nosdp')
%
%  RES fields:
%    .feas   logical, optimal slack > tol
%    .slack  optimal slack s (the certified decay margin proxy)
%    .mode, .Tmin, .Tmax
%    .KJfun  @(theta) -> refresh gain K_J(theta) (constant/slack ignore theta)
%    .KJ     representative gain: K_J(Tmin) here for all modes
%    .info   struct('Y',Y,'UJ',UJ,'T',T,'deg',deg,'N',N,'n',n,'m',m,'th',th);
%            in the slack mode UJ holds the constant U and .G holds the slack Ys
%
%  See also MVST_STABZ_SOS (generic dwell-time synthesis), MVST_SVECBASIS,
%  MVST_LIFTSYM (used by the callers to certify the recovered gains by the exact
%  one-period monodromy on the lifted symmetric space).

if nargin < 2; opts = struct(); end
oo = mvst_opts(opts);
if ~isfield(sys,'A')  || isempty(sys.A);  error('mvst:sdrange:A','SYS.A required.');  end
if ~isfield(sys,'Bu') || isempty(sys.Bu); error('mvst:sdrange:Bu','SYS.Bu required.'); end
A = sys.A; Bu = sys.Bu; n = size(A,1); m = size(Bu,2);
Ac = getf(sys,'Ac',zeros(n,n,0)); q = size(Ac,3);
Bc = getf(sys,'Bc',zeros(n,m,q));
if size(Bc,3) ~= q; Bc = zeros(n,m,q); end

Tmin = getf(opts,'Tmin',[]); Tmax = getf(opts,'Tmax',[]);
if isempty(Tmin) || isempty(Tmax); error('mvst:sdrange:T','OPTS.Tmin and OPTS.Tmax required.'); end
if ~(Tmin > 0 && Tmax >= Tmin); error('mvst:sdrange:Trange','Need 0 < Tmin <= Tmax.'); end
mode  = lower(getf(opts,'mode','scheduled'));
deg   = getf(opts,'deg',4);
sreq  = getf(opts,'slack',0);
alpha = getf(opts,'alpha',0);
rho   = getf(opts,'rho',1);
tol   = 1e-7;

% ---------------------------------------------------- held-input lift
N  = n + m;
tA = [A Bu; zeros(m,N)];
tAl = cell(1,q);
for l = 1:q; tAl{l} = [Ac(:,:,l) Bc(:,:,l); zeros(m,N)]; end
J0 = [eye(n) zeros(n,m); zeros(m,N)];
BJ = [zeros(n,m); eye(m)];

isConst = strcmp(mode,'constant');
isSlack = strcmp(mode,'slack');
if ~(isConst || isSlack || strcmp(mode,'scheduled'))
    error('mvst:sdrange:mode','OPTS.mode must be scheduled | constant | slack.');
end
if isConst; T = Tmin; else; T = Tmax; end   % slack keeps Y free on [0,Tmax]

pvar th; sl = dpvar('sl');
prog = sosprogram(th); prog = sosdecvar(prog, sl);
mons = monomials(th, 0:deg);
[prog, Y]  = sospolymatrixvar(prog, mons, [N N], 'symmetric');
if isConst
    [prog, UJ] = sospolymatrixvar(prog, monomials(th,0), [m N]);   % constant refresh
elseif isSlack
    [prog, Ys] = sospolymatrixvar(prog, monomials(th,0), [N N]);   % constant slack (not symmetric)
    [prog, U]  = sospolymatrixvar(prog, monomials(th,0), [m N]);   % constant U, K_J = U Ys^{-1}
else
    [prog, UJ] = sospolymatrixvar(prog, mons, [m N]);              % scheduled refresh
end
gt = th*(1-th);
sfix = sl + sreq;                     % strict jump margin = optimized slack + request

% Held-input SD convention: the held command has ZERO flow dynamics (its lifted
% block is zero), so the flow cannot contract; the flow LMI is imposed only as
% NON-EXPANSIVE (no strict margin) and ALL the decay is carried by the strict
% refresh (jump) LMI.  This is the persist='jumping' relaxation, mandatory here
% because the frozen-tail flow block has a zero (u,u) diagonal and cannot be
% made strictly negative definite.

% ------------------------------------------------------------- flow LMI
Yp = diff(Y, th);
Phif = flow_block(tA, tAl, Y, Yp/T, alpha, N, q);
[prog, S1] = sospolymatrixvar(prog, monomials(th,0:2*deg), size(Phif), 'symmetric');
prog = sosmatrixineq(prog, -Phif - gt*S1);   % non-expansive flow, no strict margin
prog = sosmatrixineq(prog, S1);

% -------------------------------------- frozen-tail flow LMI (constant mode)
if isConst && Tmax > Tmin
    Ymin  = subs(Y, th, 1);           % Y(Tmin) (flat certificate on the window)
    PhifT = flow_block(tA, tAl, Ymin, zeros(N), alpha, N, q);   % Ydot = 0
    prog  = sosmatrixineq(prog, -PhifT);                         % non-expansive
end

% ------------------------------------------------------------- jump LMI (strict)
Y0 = subs(Y, th, 0);
if isConst
    Yend = subs(Y, th, 1);            % Y(Tmin)
    Phij = jump_block(J0, BJ, Yend, Y0, UJ, rho);
    prog = sosmatrixineq(prog, -Phij - sfix*eye(size(Phij,1)));
elseif isSlack
    % Author's dilated (Finsler) constant-gain jump family over [Tmin,Tmax]:
    % Y(theta) FREE, constant slack Ys and constant U.  Phij is the <=0 block.
    if Tmax > Tmin
        Phij = slack_block(J0, BJ, Y, Y0, Ys, U, rho);   % Y = Y(theta) over window
        gm = (th - Tmin/T)*(Tmax/T - th);
        [prog, Sj] = sospolymatrixvar(prog, monomials(th,0:2*deg), size(Phij), 'symmetric');
        prog = sosmatrixineq(prog, -Phij - sfix*eye(size(Phij,1)) - gm*Sj);
        prog = sosmatrixineq(prog, Sj);
    else                               % degenerate window Tmin==Tmax: single LMI
        Yend = subs(Y, th, 1);         % th=1 maps to tau = T = Tmax
        Phij = slack_block(J0, BJ, Yend, Y0, Ys, U, rho);
        prog = sosmatrixineq(prog, -Phij - sfix*eye(size(Phij,1)));
    end
elseif Tmax > Tmin
    Phij = jump_block(J0, BJ, Y, Y0, UJ, rho);          % over [Tmin,Tmax]
    gm = (th - Tmin/T)*(Tmax/T - th);
    [prog, Sj] = sospolymatrixvar(prog, monomials(th,0:2*deg), size(Phij), 'symmetric');
    prog = sosmatrixineq(prog, -Phij - sfix*eye(size(Phij,1)) - gm*Sj);
    prog = sosmatrixineq(prog, Sj);
else                                   % degenerate window Tmin==Tmax: single LMI
    Yend = subs(Y, th, 1);             % th=1 maps to tau = T = Tmax
    Phij = jump_block(J0, BJ, Yend, Y0, subs(UJ,th,1), rho);
    prog = sosmatrixineq(prog, -Phij - sfix*eye(size(Phij,1)));
end

% ------------------------------------------- normalization Y <= I (scale-fix)
[prog, Xi] = sospolymatrixvar(prog, monomials(th,0:2*deg), [N N], 'symmetric');
prog = sosmatrixineq(prog, eye(N) - Y - gt*Xi);
prog = sosmatrixineq(prog, Xi);

prog = sossetobj(prog, -sl);
o = mvst_sosopts(oo);
prog = mvst_sossolve(prog, o);

sval = double(sosgetsol(prog, sl));
Yv   = sosgetsol(prog, Y);
if isSlack
    Ysv = sosgetsol(prog, Ys);
    Uv  = sosgetsol(prog, U);
else
    UJv = sosgetsol(prog, UJ);
end

res.mode  = mode;
res.Tmin  = Tmin; res.Tmax = Tmax;
res.slack = sval;
res.feas  = sval > tol;
if isConst
    Yend = evalmat(Yv, th, 1, N, N);      % Y(Tmin)
    UJc  = evalmat(UJv, th, 0, m, N);     % constant
    res.KJfun = @(theta) UJc / Yend;
elseif isSlack
    Ysc = evalmat(Ysv, th, 0, N, N);      % constant slack Ys
    Uc  = evalmat(Uv,  th, 0, m, N);      % constant U
    Kc  = Uc / Ysc;                       % K_J = U Ys^{-1} (theta-independent)
    res.KJfun = @(theta) Kc;
else
    res.KJfun = @(theta) evalmat(UJv, th, theta/T, m, N) / evalmat(Yv, th, theta/T, N, N);
end
res.KJ   = res.KJfun(Tmin);
if isSlack
    res.info = struct('Y',Yv,'UJ',Uv,'G',Ysv,'T',T,'deg',deg,'N',N,'n',n,'m',m,'th',th);
else
    res.info = struct('Y',Yv,'UJ',UJv,'T',T,'deg',deg,'N',N,'n',n,'m',m,'th',th);
end
end

% ==========================================================================
function P = flow_block(tA, tAl, Y, Ydot, alpha, N, q)
% [ -Ydot + He(tilde A Y) + alpha Y , row_l (tilde A_l Y)^* ;
%    *                              , -I_q kron Y ]   (covariance-level +alpha Y)
AY = tA*Y;
P11 = -Ydot + (AY + AY.') + alpha*Y;
Vrow = []; Vcol = [];
for l = 1:q
    V = tAl{l}*Y; Vrow = [Vrow, V.']; Vcol = [Vcol; V]; %#ok<AGROW>
end
P = [ P11, Vrow ; Vcol, -blkdiagrep(Y, q) ];
end

function P = jump_block(J0, BJ, Yend, Y0, UJ, rho)
% [ -rho Yend , (J_0 Yend + B_J U_J)^* ; * , -Y0 ]   (covariance-level rho)
Z = J0*Yend + BJ*UJ;
P = [ -rho*Yend, Z.' ; Z, -Y0 ];
end

function P = slack_block(J0, BJ, Ytheta, Y0, Ys, U, rho)
% Dilated (single-slack Finsler) constant-gain jump block, covariance domain
% (<=0 convention: this returns P with -P the PSD dilated block).  With
% Z := J_0 Ys + B_J U = (J_0 + B_J K_J) Ys, K_J = U Ys^{-1}, He[W] := W + W^*:
%   -P = [ rho Y(0) ,   Z              ;
%          Z^*       ,   He[Ys] - Y(theta) ] >= 0.
% The (2,2) He-corner He[Ys] - Y(theta) is the slack corner (a -(Ys+Ys^*)-type
% term relative to +Y(theta)).  By He[Ys] - Y(theta) <= Ys^* Y(theta)^{-1} Ys
% (Young: (Ys - Y(theta))^* Y(theta)^{-1}(Ys - Y(theta)) >= 0) and a Schur
% complement on that corner, the slack CANCELS exactly (Ys Ys^{-1} = I) and the
% block implies the covariance jump contraction
%   (J_0 + B_J K_J) Y(theta) (J_0 + B_J K_J)^* <= rho Y(0)   (all theta),
% the very inequality certified by the scheduled block, but now with Y(theta)
% period-dependent (free SOS) and Ys, U (hence K_J) CONSTANT.
Z  = J0*Ys + BJ*U;
He = @(W) W + W.';
P  = [ -rho*Y0,   -Z            ;
       -Z.',       Ytheta - He(Ys) ];
end

function B = blkdiagrep(Y, k)
n = size(Y,1); if k == 0; B = zeros(0); return; end
B = [];
for i = 1:k
    row = [];
    for j = 1:k; if i==j; row=[row,Y]; else; row=[row,zeros(n)]; end; end %#ok<AGROW>
    B = [B; row]; %#ok<AGROW>
end
end

function M = evalmat(P, th, theta, r, c)
if isnumeric(P); M = P; else; M = double(subs(P, th, theta)); end
M = reshape(M, r, c);
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
