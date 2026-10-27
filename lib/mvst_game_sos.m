function res = mvst_game_sos(sys, game, dwell, mode, opts)
%MVST_GAME_SOS  Generic SOSTOOLS state-feedback zero-sum GAME / minimax program
%  for linear impulsive Hermitian matrix-valued systems, the game analogue of
%  MVST_STAB_SOS (stability), MVST_STABZ_SOS (stabilization synthesis) and
%  MVST_DISS_SOS (dissipativity analysis).  It returns the critical
%  disturbance-attenuation level GAMMASTAR of the H-infinity / minimax game
%  (minimizer control u, maximizer disturbance w) and the associated central
%  gains, over every dwell-time type of the zero-sum games chapter.
%
%  ---------------------------------------------------------------------------
%  WHAT IS COMPUTED, AND THE HONESTY CAVEAT (read this).
%  ---------------------------------------------------------------------------
%  The book's dual game certificate (Theorem th:game:LMI, Corollary
%  cor:CT:game:LMI, Theorem th:game:inf:LMI and the causal dwell-time
%  Theorem th:game:periodic) is the game Riccati INEQUALITY
%
%        Isaacs cond.~(eq:Isaacs:gen) for I, I^j ;
%        I / G  >= 0 ,   I^j / G^j >= 0 ,
%
%  where G = [I_uu I_uw ; I_uw^* I_ww] is the INDEFINITE Isaacs block
%  (eq:Ghat), I_uu = R_u >= 0 (minimizer), I_ww = -R_w <= 0 (maximizer).  As
%  the book stresses (the remark after th:game:LMI, lines "In contrast to the
%  optimal control dual ... this is NOT a plain LMI"), the generalized Schur
%  complement of the indefinite G is rational in P and paired with the
%  saddle-inertia condition, so it is NOT jointly convex in the dual variable.
%  It therefore cannot be handed to an SDP/SOS solver verbatim.
%
%  This routine solves instead the SUFFICIENT convex relaxation obtained by
%  the classical bounded-real change of variables Y = P^{-1}, W = K Y (the
%  same congruence MVST_STABZ_SOS uses for stabilization).  In the separated
%  H-infinity setting I_uw = 0 (the game / bounded-real case singled out in
%  the same remark, eq. with I/G = I11 - I1u pinv(I_uu) I1u^* - I1w pinv(I_ww)
%  I1w^*), this relaxation is EXACT: minimizing gamma over the synthesis LMI
%  reproduces the stabilizing solution of the continuous-time HJI-GRE
%  (cor:CT:game:fh, eq:CT:game:GRE) and hence the true gamma*.  In general
%  (timer-dependent / dwell-time certificates, or nonzero cross term) the SOS
%  value is an UPPER BOUND on gamma*, converging as the certificate degree
%  grows.  This is the mirror of the honesty note in MVST_DISS_SOS: the exact
%  per-class solvers MVST_GAME_GAMMA_CT* integrate the Riccati and return the
%  exact gamma*; this routine returns a certified upper bound.
%
%  ---------------------------------------------------------------------------
%  MODEL  (mirrors the OTHER SOS routines; any omitted field defaults away).
%  ---------------------------------------------------------------------------
%  SYS fields:
%    .A0   n x n         flow drift A
%    .Bu   n x mu        flow control input  B_u   (minimizer)
%    .Bw   n x mw        flow disturbance input B_w (maximizer)
%    .Ac   n x n x q     flow state-multiplicative (Brownian/Levy) channels;
%                        pass sqrt(c_l) A_l, fold the Levy mean drift into A0
%    .Bc   n x mu x q    input-multiplicative noise selectors (paired with .Ac)
%    .Jf   n x n x pf    SPONTANEOUS (Poisson) jump maps J_r
%    .BJf  n x mu x pf   spontaneous-jump control selectors (paired with .Jf)
%    .lamf 1 x pf        Poisson rates
%    .J    n x n [x p]   SCHEDULED (dwell-time / sampling) jump map(s) (default I
%                        when a jump condition is required)
%    .BJ   n x mu [x p]  scheduled-jump control selectors
%
%  GAME weights (the minimax cost is |z|^2 - gamma^2 |w|^2 with
%    z = [Q^{1/2} x ; R_u^{1/2} u], exactly the weights of MVST_GAME_GAMMA_CT):
%    .Q    n x n  >= 0   state weight
%    .Ru   mu x mu > 0   control weight
%    .gamma  fixed level (used only in 'feas' mode).
%
%  DWELL  (identical taxonomy to MVST_STABZ_SOS / MVST_DISS_SOS):
%    .type 'arbitrary' | 'discrete' | 'periodic' | 'minimum' | 'range'
%          ('arbitrary' = constant certificate: CT LTI / GBM / Poisson game,
%                         Poisson jumps live in the flow via Jf/lamf, no
%                         scheduled jump; EXACT for the CT LTI separated game;
%           'discrete'  = pure DT bounded-real game (transition A0, inputs
%                         Bu,Bw, output (Q,Ru), constant certificate);
%           'periodic'  = fixed dwell-time Tbar, single scheduled jump at Tbar
%                         (sampled-data game), th:game:periodic;
%           'minimum'   = minimum dwell-time Tbar, ramp on [0,Tbar] + frozen
%                         tail;
%           'range'     = ranged dwell-time [Tmin,Tmax] via the dwell mult.)
%    .Tbar (periodic,minimum) ; .Tmin,.Tmax (range)
%    .deg  timer polynomial degree of the certificate (default 2; forced 0 for
%          arbitrary/discrete)
%    .alpha optional flow shift (default 0)
%
%  MODE: 'gamma' (default; bisection on gamma for the critical level, mirroring
%        MVST_GAME_GAMMA_CT) | 'feas' (feasibility at the fixed game.gamma).
%
%  ---------------------------------------------------------------------------
%  RETURNS
%  ---------------------------------------------------------------------------
%  RES:
%    .gammastar  critical level (upper bound on the true gamma*; == gamma* for
%                the CT LTI separated game).  Also aliased as .gamma.
%    .feas       feasibility flag at the returned level.
%    .P          @(tau) -> P(tau) = Y(tau)^{-1}, the game-value certificate
%                (dual variable of th:game:inf:LMI).
%    .Ku         @(tau) -> central MINIMIZER (control) gain K_u = W Y^{-1}.
%                Cleanly recoverable from the synthesis variables.
%    .Kw         @(tau) -> central MAXIMIZER (disturbance) gain
%                K_w = gamma^{-2} Bw^* Y^{-1} (eq:game:Ks reading Bw^* P).
%                CAVEAT: this equals the true saddle maximizer gain only at the
%                converged critical level (where Y^{-1} is the HJI-GRE
%                solution); at a strictly feasible slack level it is the
%                worst-case gain w.r.t. the certificate, not the game value.
%    .Kf,.KJ     spontaneous- / scheduled-jump control gains (empty if absent).
%    .slack, .hist, .info.
%
%  OPTS (struct or MVST.OPTIONS, optional): .solver selects the SDP backend
%  ('sedumi' (default) | 'sdpt3' | 'mosek' | 'nosdp'), normalized by MVST_OPTS.
%
%  VALIDATION (mvst_game_sos_selftest, run by mvst_game_sos('selftest')):
%    the CT arbitrary game against MVST_GAME_GAMMA_CT on small instances
%    (SOS bound >= exact gamma*, agreeing to solver tolerance in the separated
%    case), plus a degree-refinement check.

if nargin == 1 && ischar(sys) && strcmp(sys,'selftest')
    res = mvst_game_sos_selftest(); return;
end
if nargin < 4 || isempty(mode); mode = 'gamma'; end
if nargin < 3 || isempty(dwell); dwell = struct('type','arbitrary'); end
if nargin < 5; opts = struct(); end
if ~isstruct(sys)
    error('mvst:game_sos:input', ...
        'SYS must be a rich struct with fields .A0, .Bu, .Bw.');
end
for f = {'A0','Bu','Bw'}
    if ~isfield(sys,f{1}) || isempty(sys.(f{1}))
        error('mvst:game_sos:input', 'SYS.%s is required.', f{1});
    end
end
if nargin < 2 || ~isstruct(game)
    error('mvst:game_sos:input', ...
        'GAME must be a struct carrying the weights .Q and .Ru.');
end
oo = mvst_opts(opts); solver = oo;   % normalized options travel to the solve

n  = size(sys.A0,1);
mu = size(sys.Bu,2);
mw = size(sys.Bw,2);
sys = fill_defaults(sys, n, mu);
if ~isfield(game,'Q')  || isempty(game.Q);  game.Q  = zeros(n);  end
if ~isfield(game,'Ru') || isempty(game.Ru); game.Ru = eye(mu);   end

switch mode
    case 'gamma'
        res = bisect_gamma(sys, game, dwell, n, mu, mw, solver);
    case 'feas'
        g = getf(game,'gamma',1);
        res = solve_once(sys, game, dwell, n, mu, mw, g, true, solver);
        res.gammastar = g; res.gamma = g;
        res = finish_gains(res, sys, mu, mw, dwell);
    otherwise
        error('mvst:game_sos:mode','unknown mode "%s"', mode);
end
end

% =========================================================================
function res = bisect_gamma(sys, game, dwell, n, mu, mw, solver)
%BISECT_GAMMA  Smallest gamma for which the synthesis game LMIs are SOS-feasible.
%   Mirrors MVST_GAME_GAMMA_CT's bracket-then-bisect, but feasibility is
%   decided by slack maximization (robust near the boundary), as in the rest
%   of the SOS family.
lo = 1e-3; hi = 1.0; Rbest = []; grow = 0; hist = [];
while true
    r = solve_once(sys, game, dwell, n, mu, mw, hi, true, solver);
    hist(end+1,:) = [hi, r.slack]; %#ok<AGROW>
    if r.feas; Rbest = r; break; end
    hi = 2*hi; grow = grow + 1;
    if grow > 30; break; end
end
for it = 1:26
    mid = 0.5*(lo + hi);
    r = solve_once(sys, game, dwell, n, mu, mw, mid, true, solver);
    hist(end+1,:) = [mid, r.slack]; %#ok<AGROW>
    if r.feas; hi = mid; Rbest = r; else; lo = mid; end
end
res = Rbest;
if isempty(res)
    res = solve_once(sys, game, dwell, n, mu, mw, hi, false, solver);
end
res.gammastar = hi; res.gamma = hi;
res.hist = sortrows(hist,1);
res = finish_gains(res, sys, mu, mw, dwell);
end

% =========================================================================
function res = solve_once(sys, game, dwell, n, mu, mw, gam, useSlackObj, solver)
%SOLVE_ONCE  Build and solve the bounded-real GAME synthesis program at a
%   fixed disturbance level gam.  Variables Y = P^{-1} (>0), W = K Y, the
%   change of variables that convexifies the indefinite dual DLMI of
%   th:game:LMI into the sufficient certificate documented in the header.
tol    = 1e-7;
Rcap   = 1e4;                      % certificate-magnitude cap 0 < Y <= Rcap I
margin = 1e-7;
epsY   = 1e-6;                     % Y >= epsY I  (so P = Y^{-1} exists)
alpha  = getf(dwell,'alpha',0);
Cq = realsqrt_psd(game.Q);        % Q = Cq^*Cq
Cr = realsqrt_psd(game.Ru);       % Ru = Cr^*Cr

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
[prog, W] = sospolymatrixvar(prog, mons, [mu n]);
q  = size(sys.Ac,3);  pf = size(sys.Jf,3);  p = size(sys.J,3);
Wf = zeros(mu,n); WJ = zeros(mu,n);
if pf > 0; [prog, Wf] = sospolymatrixvar(prog, mons, [mu n]); end
if p  > 0; [prog, WJ] = sospolymatrixvar(prog, mons, [mu n]); end
gt = th*(1-th);

if useSlackObj; slk = sl; else; slk = margin; end

% ------------------------------------------------ flow (or DT) game condition
if strcmp(dwell.type,'discrete')
    Phi = discrete_game(sys, Y, W, game, gam, Cq, Cr, n, mu, mw, q);
    prog = matle(prog, Phi, slk);              % Phi <= -slk I
else
    Ydot = diff(Y, th)/T;
    Phi = flow_game(sys, Y, Ydot, W, Wf, alpha, gam, Cq, Cr, n, mu, mw, q, pf);
    switch dwell.type
        case 'arbitrary'
            prog = matle(prog, Phi, slk);
        case {'periodic','minimum','range'}
            nP = size(Phi,1);
            [prog, S1] = sospolymatrixvar(prog, monomials(th,0:2*deg), [nP nP], 'symmetric');
            prog = matsos(prog, -Phi - slk*eye(nP) - gt*S1);
            prog = matsos(prog, S1);
    end
    if strcmp(dwell.type,'minimum')            % frozen tail at theta = 1
        Phi0 = flow_game(sys, subs(Y,th,1), zeros(n), subs(W,th,1), subs(Wf,th,1), ...
                         alpha, gam, Cq, Cr, n, mu, mw, q, pf);
        prog = matle(prog, Phi0, slk);
    end
end

% ------------------------------------------------- scheduled-jump condition
if p > 0 && ~strcmp(dwell.type,'discrete')
    switch dwell.type
        case 'arbitrary'
            Dj = jump_game(sys, Y, Y, WJ, n, p);
            prog = matle(prog, Dj, slk);
        case {'periodic','minimum'}
            Dj = jump_game(sys, subs(Y,th,1), subs(Y,th,0), subs(WJ,th,1), n, p);
            prog = matle(prog, Dj, slk);
        case 'range'
            Dj = jump_game(sys, Y, subs(Y,th,0), WJ, n, p);
            gm = (th - dwell.Tmin/T)*(dwell.Tmax/T - th);
            nj = size(Dj,1);
            [prog, Sj] = sospolymatrixvar(prog, monomials(th,0:2*deg), [nj nj], 'symmetric');
            prog = matsos(prog, -Dj - slk*eye(nj) - gm*Sj);
            prog = matsos(prog, Sj);
    end
end

% ------------------------- certificate scaling 0 < epsY I <= Y <= Rcap I
if isTimer
    [prog, Xlo] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, Y - epsY*eye(n) - gt*Xlo);   prog = matsos(prog, Xlo);
    [prog, Xhi] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n n], 'symmetric');
    prog = matsos(prog, Rcap*eye(n) - Y - gt*Xhi);   prog = matsos(prog, Xhi);
else
    prog = matsos(prog, Y - epsY*eye(n));
    prog = matsos(prog, Rcap*eye(n) - Y);
end

if useSlackObj; prog = sossetobj(prog, -sl); end
o = mvst_sosopts(solver);
[prog, info] = mvst_sossolve(prog, o);

Yv = sosgetsol(prog, Y);  Wv = sosgetsol(prog, W);
Wfv = zeros(mu,n); if pf>0; Wfv = sosgetsol(prog, Wf); end
WJv = zeros(mu,n); if p >0; WJv = sosgetsol(prog, WJ); end
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
res.P     = @(tau) inv(evalmat(Yv, th, tau/T, n, n));
res.info  = struct('Y',Yv,'W',Wv,'Wf',Wfv,'WJ',WJv,'T',T,'deg',deg, ...
                   'alpha',alpha,'solverinfo',info);
end

% =========================================================================
function C = flow_game(sys, Y, Ydot, W, Wf, alpha, gam, Cq, Cr, n, mu, mw, q, pf)
% Bounded-real GAME flow synthesis block (require C <= 0).  Core state block
%   G11 = -Ydot + He(A0 Y + Bu W) + 2 alpha Y - (sum lamf) Y ,
% augmented with:
%   - disturbance w:  off-diag Bw at (state,w),  diag -gam^2 I
%       (Schur of the negative w-block adds + gam^{-2} Bw Bw^* to G11);
%   - state performance:  off-diag (Cq Y)^*, diag -I  (adds Y Q Y);
%   - control performance: off-diag (Cr W)^*, diag -I (adds W^* Ru W);
%   - Brownian second moment (l): off-diag (Ac_l Y + Bc_l W)^*, diag -Y;
%   - Poisson (r): off-diag lamf_r (Jf_r Y + BJf_r Wf)^*, diag -lamf_r Y.
% Schur-eliminating every trailing block reproduces the game Riccati
% inequality  He(A0 Y + Bu W) + Y Q Y + W^* Ru W + gam^{-2} Bw Bw^*
%             + (noise second moments)  <= 0,  the Y = P^{-1}, W = K Y image
% of I/G >= 0 (cor:CT:game:LMI); the +gam^{-2} Bw Bw^* term entering with the
% opposite sign to the control term, exactly as flagged after th:game:LMI.
lamf = sys.lamf;
AY  = sys.A0*Y + sys.Bu*W;
G11 = -Ydot + (AY + AY.') + 2*alpha*Y;
if pf > 0; G11 = G11 - sum(lamf)*Y; end
g2  = gam^2;
% ---- gain core: state / disturbance / state-perf / control-perf ----
CqY = Cq*Y;   nq = size(Cq,1);
CrW = Cr*W;
core = [ G11,        sys.Bw,      CqY.',        CrW.'      ; ...
         sys.Bw.',  -g2*eye(mw),  zeros(mw,nq), zeros(mw,mu); ...
         CqY,        zeros(nq,mw),-eye(nq),      zeros(nq,mu); ...
         CrW,        zeros(mu,mw), zeros(mu,nq),-eye(mu)     ];
Ncore = size(core,1);
% ---- Brownian + Poisson second-moment Schur channels (couple to state) ----
Rrow = []; Dblk = [];
for l = 1:q
    V = sys.Ac(:,:,l)*Y + sys.Bc(:,:,l)*W;
    Rrow = [Rrow, V.'];               %#ok<AGROW> couples to state only
    Dblk = blkdiagcat(Dblk, -Y);
end
for r = 1:pf
    Z = sys.Jf(:,:,r)*Y + sys.BJf(:,:,r)*Wf;
    Rrow = [Rrow, lamf(r)*Z.'];       %#ok<AGROW>
    Dblk = blkdiagcat(Dblk, -lamf(r)*Y);
end
if isempty(Rrow)
    C = core;
else
    nc = size(Rrow,2);
    Rfull = [Rrow; zeros(Ncore-n, nc)];        % zeros against w / perf blocks
    C = [ core, Rfull ; Rfull.', Dblk ];
end
C = (C + C.')/2;
end

% -------------------------------------------------------------------------
function D = jump_game(sys, Yend, Y0, WJ, n, p)
% Scheduled-jump synthesis block (require D <= 0), state reset only:
%   [ -Yend , row_r[Z_r^*] ; col_r[Z_r] , -I_p (x) Y0 ], Z_r = J_r Yend + BJ_r WJ.
% Schur -> sum_r Z_r Y0^{-1} Z_r^* - Yend <= 0, the Y = P^{-1} image of the
% jump Riccati inequality (dwell-time game, th:game:periodic).  No jump-channel
% L2 term is carried (sampled-data games reset the state and score in flow).
Zrow = []; Zcol = [];
for r = 1:p
    Zr = sys.J(:,:,r)*Yend + sys.BJ(:,:,r)*WJ;
    Zrow = [Zrow, Zr.']; Zcol = [Zcol; Zr];   %#ok<AGROW>
end
D = [ -Yend, Zrow ; Zcol, -blkdiagrep(Y0,p) ];
D = (D + D.')/2;
end

% -------------------------------------------------------------------------
function D = discrete_game(sys, Y, W, game, gam, Cq, Cr, n, mu, mw, q)
% Pure DT bounded-real GAME synthesis (require D <= 0).  Closed-loop transition
% M0 = A0 Y + Bu W, disturbance Bw, performance (Cq Y ; Cr W):
%   [ -Y      M0        Bw          0     ;
%     M0^*    -Y        0           Cz^*  ;
%     Bw^*    0        -gam^2 I     0     ;
%     0       Cz       0           -I     ]  <= 0,   Cz = [Cq Y ; Cr W].
% Also carries DT state-multiplicative noise M_l = Ac_l Y + Bc_l W by extending
% the transition row/col with additional -Y blocks (second-moment sum).
% NOTE: no exact DT game solver exists in this tree to cross-check against;
% this member is provided by the standard DT bounded-real synthesis form and
% is validated only through its CT reduction.
M0 = sys.A0*Y + sys.Bu*W;
Cz = [Cq*Y ; Cr*W];  nz = size(Cz,1);
g2 = gam^2;
% transition channels (drift + Brownian), each Schur'd against a -Y block
Mrow = M0.';  nch = 1;
for l = 1:q
    Ml = sys.Ac(:,:,l)*Y + sys.Bc(:,:,l)*W;
    Mrow = [Mrow, Ml.']; nch = nch + 1;        %#ok<AGROW>
end
% assemble: [ -Y (x) I_nch (only first is transition target) ...]
% Simpler: keep the standard single-transition BRL and append noise as extra
% -Y Schur blocks coupling to the state row.
D11 = -Y;
% state row couples to: transition block(s) (dim nch*n), disturbance (mw), 0
top = [ D11, Mrow, sys.Bw, zeros(n,nz) ];
midT = [ Mrow.', -blkdiagrep(Y,nch), zeros(nch*n,mw), [Cz.'; zeros((nch-1)*n, nz)] ];
midW = [ sys.Bw.', zeros(mw,nch*n), -g2*eye(mw), zeros(mw,nz) ];
botZ = [ zeros(nz,n), [Cz, zeros(nz,(nch-1)*n)], zeros(nz,mw), -eye(nz) ];
D = [ top; midT; midW; botZ ];
D = (D + D.')/2;
end

% =========================================================================
function res = finish_gains(res, sys, mu, mw, dwell)
%FINISH_GAINS  Attach the central-gain function handles to RES.
if ~isfield(res,'info') || ~isfield(res.info,'Y')
    res.Ku = []; res.Kw = []; res.Kf = []; res.KJ = []; return;
end
Yv = res.info.Y; Wv = res.info.W; T = res.info.T; n = size(sys.A0,1);
pvar th;
res.Ku = @(tau) evalmat(Wv,th,tau/T,mu,n) / evalmat(Yv,th,tau/T,n,n);
res.Kw = @(tau) (sys.Bw.' / evalmat(Yv,th,tau/T,n,n)) / (res.gamma^2);
res.Kf = []; res.KJ = [];
if isfield(res.info,'Wf') && any(sys.lamf>0) && size(sys.Jf,3)>0
    Wfv = res.info.Wf;
    res.Kf = evalmat(Wfv,th,0,mu,n) / evalmat(Yv,th,0,n,n);
end
if isfield(res.info,'WJ') && size(sys.J,3) > 0
    WJv = res.info.WJ; thv = 1;
    if strcmp(dwell.type,'arbitrary'); thv = 0; end
    res.KJ = evalmat(WJv,th,thv,mu,n) / evalmat(Yv,th,thv,n,n);
end
end

% =========================================================================
function prog = matle(prog, EXPR, slk)
% EXPR <= -slk I  as a matrix SOS (constant multiplier: plain LMI).
nn = size(EXPR,1);
prog = sosmatrixineq(prog, -EXPR - slk*eye(nn));
end
function prog = matsos(prog, EXPR); prog = sosmatrixineq(prog, EXPR); end

function M = evalmat(P, th, theta, r, c)
if isnumeric(P); M = P; else; M = double(subs(P, th, theta)); end
M = reshape(M, r, c);
end

function S = realsqrt_psd(Q)
% Real symmetric PSD square root Q = S^*S (S returned, size n x n).
Q = (Q + Q.')/2;
[U,D] = eig(Q); d = max(real(diag(D)),0);
S = diag(sqrt(d))*U.';
S = real(S);
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

function B = blkdiagcat(B, X)
% Append X on the block diagonal of B (both square; sizes may differ).
if isempty(B); B = X; return; end
[a,~] = size(B); [c,~] = size(X);
B = [ B, zeros(a,c) ; zeros(c,a), X ];
end

function sys = fill_defaults(sys, n, mu)
if ~isfield(sys,'Ac')||isempty(sys.Ac);  sys.Ac = zeros(n,n,0); end
if ~isfield(sys,'Bc')||isempty(sys.Bc)||size(sys.Bc,3)~=size(sys.Ac,3)
    sys.Bc = zeros(n,mu,size(sys.Ac,3)); end
if ~isfield(sys,'Jf')||isempty(sys.Jf);  sys.Jf = zeros(n,n,0); end
if ~isfield(sys,'BJf')||isempty(sys.BJf)||size(sys.BJf,3)~=size(sys.Jf,3)
    sys.BJf = zeros(n,mu,size(sys.Jf,3)); end
if ~isfield(sys,'lamf')||isempty(sys.lamf); sys.lamf = ones(1,size(sys.Jf,3)); end
if ~isfield(sys,'J')||isempty(sys.J);    sys.J = zeros(n,n,0); end
if ~isfield(sys,'BJ')||isempty(sys.BJ)||size(sys.BJ,3)~=size(sys.J,3)
    sys.BJ = zeros(n,mu,size(sys.J,3)); end
end

function v = getf(s, f, d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end
