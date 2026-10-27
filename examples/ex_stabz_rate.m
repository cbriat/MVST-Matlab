function ex_stabz_rate()
%EX_STABZ_RATE  Guaranteed-rate (alpha,rho) minimum-dwell-time state-feedback
%  synthesis for the impulsive system of ex_stabz_sos, Stabilization chapter
%  (Theorem th:stabz:ratedwell).
%
%  Same system:
%     flow:  dx = (A0 x + B0 u) dt + (A1 x + B1 u) dW1,   A1 = sigma*I,
%     jump:  x+ = J x,   J = jfac*I   (expansive, uncontrolled).
%  The (alpha,rho)-shifted synthesis blocks are those of ex_stabz_sos with the
%  flow (1,1) block gaining +alpha*Y(theta) and the jump (1,1) block scaled to
%  -rho*Y(1) (COVARIANCE-level rate: X ~ e^{-alpha t} rho^kappa):
%     Phi_f^a = [ -(1/Tbar) Y' + He(A0 Y + B0 U0) + a Y , (A1 Y+B1 U0)' ;
%                  A1 Y+B1 U0                          , -Y            ] <= -s I,
%     Phi_j^r = [ -rho Y(1) , (J Y(1))' ; J Y(1) , -Y(0) ] <= -s I.
%  alpha is the guaranteed covariance-level flow decay rate of the second moment
%  and rho<1 the guaranteed per-jump contraction factor.  Because the INPUT-multiplicative
%  noise B1 (sigma>0 channel) bounds the achievable second-moment decay, the
%  maximum guaranteed rate alpha*(Tbar) is FINITE; it increases with the
%  dwell-time Tbar (a longer flow interval buys a larger guaranteed rate against
%  the fixed expansive jump).
%
%  METHOD.  PURE SOSTOOLS MATRIX SUM-OF-SQUARES in the normalized timer
%  theta = tau/Tbar in [0,1] (Positivstellensatz multiplier theta(1-theta); NO
%  gridding of the timer).  Feasibility is decided by SLACK MAXIMIZATION.  For
%  fixed rho, the dwell-time Tbar is SWEPT over scalar values (sweeping the
%  scalar parameter to draw the curve is allowed; it is not a timer grid), and
%  at each Tbar the maximum alpha*(Tbar) is found by bisection on the slack
%  sign.  Solver: SOSTOOLS + SeDuMi.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---------------- system data (as in ex_stabz_sos) ----------------
A0 = [0 1; 2 0.1];  B0 = [0;1];  sigma = 0.3;  A1 = sigma*eye(2);
beta = 0.62;  B1 = beta*[0;1];
jfac = 1.1;   J = jfac*eye(2);
% Covariance-level rate convention: rho is squared and aCap doubled relative to
% the former half-rate convention (rho_old=0.95, aCap_old=0.6) so the certified
% PHYSICAL decay is unchanged (alpha printed doubles, rho printed squares).
tol = 1e-7;  rho = 0.9025;  aCap = 1.2;

fprintf('\n==== guaranteed-rate (alpha,rho) dwell-time synthesis (SOSTOOLS matrix-SOS in tau) ====\n');
fprintf('sigma=%g  beta(B1)=%g  J=%g*I  rho=%g\n', sigma, beta, jfac, rho);

% ---------------- rho-feasible minimum dwell-time (alpha = 0) ----------------
lo = 0.1; hi = 1.5;
assert(slackval(A0,B0,A1,B1,J,lo,0,rho) <= tol, 'rho-min lower bracket feasible');
assert(slackval(A0,B0,A1,B1,J,hi,0,rho) >  tol, 'rho-min upper bracket infeasible');
for it = 1:22
    mid = 0.5*(lo+hi);
    if slackval(A0,B0,A1,B1,J,mid,0,rho) > tol, hi = mid; else, lo = mid; end
end
TminRho = hi;
fprintf('rho-feasible minimum dwell-time (alpha=0): Tbar_min(rho) = %.4f\n', TminRho);

% ---------------- alpha*(Tbar) over a sweep of Tbar values -------------------
Tgrid = linspace(1.05*TminRho, 3.2*TminRho, 9);
aStar = zeros(size(Tgrid));
for i = 1:numel(Tgrid)
    aStar(i) = max_alpha(A0,B0,A1,B1,J,Tgrid(i),rho,tol,aCap);
    fprintf('   Tbar=%.4f  ->  alpha*=%.4f\n', Tgrid(i), aStar(i));
end

% highlighted representative operating point (middle of the sweep)
iHi = ceil(numel(Tgrid)/2);
Tbar_hi = Tgrid(iHi);  alpha_hi = aStar(iHi);
fprintf('   highlighted point: (alpha=%.3f, rho=%.2f, Tbar=%.3f)\n', alpha_hi, rho, Tbar_hi);

% ---------------- figure ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, Tgrid, aStar, '-o', 'LineWidth',2.0, 'Color',gc.order(1,:), ...
     'MarkerFaceColor',gc.order(1,:),'MarkerSize',4);
plot(ax, Tbar_hi, alpha_hi, 'o', 'MarkerSize',9, 'MarkerFaceColor',gc.order(3,:), ...
     'MarkerEdgeColor','k','HandleVisibility','off');
xline(ax, TminRho, ':','Color',gc.onyx,'LineWidth',1.5,'HandleVisibility','off');
% label anchored right of and below the highlighted point (offsets scaled to
% the axis ranges) so the text clears the increasing curve and the borders
xl = xlim(ax);  yl = ylim(ax);
dxl = 0.04*diff(xl);  dyl = -0.06*diff(yl);
text(ax, Tbar_hi+dxl, alpha_hi+dyl, sprintf('$(\\bar T,\\alpha^\\star)=(%.2f,%.2f)$',Tbar_hi,alpha_hi), ...
     'Interpreter','latex','FontSize',9, ...
     'HorizontalAlignment','left','VerticalAlignment','top');
xlabel(ax,'dwell-time $\bar T$','Interpreter','latex');
ylabel(ax,'max guaranteed flow rate $\alpha^\star(\bar T)$','Interpreter','latex');
title(ax,sprintf('Guaranteed-rate synthesis, $\\rho=%.2f$',rho), ...
     'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(figdir,'stabz_rate.pdf'));

% ---------------- data macros ----------------
w = mvst_texdata(fullfile(figdir,'data','stabz_rate.tex'));
w.num('StabzRateRho', rho);
w.num('StabzRateTbar', Tbar_hi);
w.num('StabzRateAlpha', alpha_hi);
w.num('StabzRateTminRho', TminRho);
w.num('StabzRateSigma', sigma);
w.num('StabzRateJ', jfac);
w.close();
fprintf('   wrote stabz_rate.pdf and data/stabz_rate.tex\n');
end

% =========================================================================
function a = max_alpha(A0,B0,A1,B1,J,Tbar,rho,tol,aCap)
%MAX_ALPHA  Largest flow decay rate alpha for which the (alpha,rho)-synthesis
%  is feasible at dwell-time Tbar (bisection on the slack sign).
if slackval(A0,B0,A1,B1,J,Tbar,0,rho) <= tol, a = 0; return; end
lo = 0; hi = aCap;
if slackval(A0,B0,A1,B1,J,Tbar,hi,rho) > tol, a = hi; return; end   % saturated
for it = 1:20
    mid = 0.5*(lo+hi);
    if slackval(A0,B0,A1,B1,J,Tbar,mid,rho) > tol, lo = mid; else, hi = mid; end
end
a = lo;
end

% =========================================================================
function s = slackval(A0,B0,A1,B1,J,Tbar,alpha,rho)
%SLACKVAL  Optimal slack of the (alpha,rho)-synthesis as a PURE SOSTOOLS matrix
%  sum-of-squares in the normalized timer theta = tau/Tbar in [0,1] (interval
%  multiplier theta(1-theta)).  Feasible iff s > 0.  Y normalized to Y <= I so
%  the slack sign is scale-fixed.  Solver: SOSTOOLS + SeDuMi.
n = size(A0,1);
pvar th; sl = dpvar('sl');
prog = sosprogram(th);
prog = sosdecvar(prog, sl);
[prog,Y]    = sospolymatrixvar(prog, monomials(th,0:2), [n n], 'symmetric');
[prog,U]    = sospolymatrixvar(prog, monomials(th,0:2), [1 n]);
[prog,Sig1] = sospolymatrixvar(prog, monomials(th,0),   [2*n 2*n], 'symmetric');
[prog,Xi1]  = sospolymatrixvar(prog, monomials(th,0),   [n n], 'symmetric');

Yp = diff(Y,th);
AY = A0*Y + B0*U;   V = A1*Y + B1*U;
% covariance-level rate: flow (1,1) block carries +alpha Y, jump block scaled rho
Phif = [ -(1/Tbar)*Yp + (AY+AY') + alpha*Y , V.' ;
          V                                 , -Y  ];
gt = th*(1-th);
prog = sosmatrixineq(prog, -Phif - sl*eye(2*n) - gt*Sig1);
prog = sosmatrixineq(prog, Sig1);
prog = sosmatrixineq(prog, eye(n) - Y - gt*Xi1);
prog = sosmatrixineq(prog, Xi1);
Y1 = subs(Y,th,1);  Y0 = subs(Y,th,0);
Phij = [ -rho*Y1 , (J*Y1).' ; J*Y1 , -Y0 ];
prog = sosmatrixineq(prog, -Phij - sl*eye(2*n));
prog = sossetobj(prog, -sl);
o.solver = 'sedumi';  o.params.fid = 0;
prog = sossolve(prog, o);
s = double(sosgetsol(prog, sl));
end
