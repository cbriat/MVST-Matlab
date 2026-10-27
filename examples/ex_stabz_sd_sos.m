function ex_stabz_sd_sos()
%EX_STABZ_SD_SOS  SAMPLED-DATA STABILIZATION: design of the sampled state
%  feedback gain for a target sampling period, certified by the SOSTOOLS
%  matrix-SOS sampled-data analysis of ex_stab_sd_sos (Stabilization chapter).
%
%  Continuous plant  dx = A x dt + B u,  zero-order-hold feedback u(t)=K x(t_k),
%  augmented impulsive model  xi = [x; u_hold]:
%     flow      xidot = Atil xi,  Atil = [A B; 0 0],
%     sampling  xi+  = Jtil xi,   Jtil = [I 0; K 0].
%
%  METHOD (honest).  The gain K enters the sampling map Jtil BILINEARLY with the
%  certificate, so there is no clean convex SOS synthesis.  We use the
%  design-by-analysis route: a one-parameter family of stabilizing gains K(omega)
%  (continuous pole placement at -omega, -2 omega) is swept, and for each the
%  MAXIMUM stable sampling period is CERTIFIED by the degree-4 SOSTOOLS matrix
%  sum-of-squares sampled-data analysis in the normalized timer (identical to
%  ex_stab_sd_sos: flow matrix-SOS on [0,1] via the interval multiplier
%  theta(1-theta), sampling LMI Jtil' P(0) Jtil - P(1) <= -s I, slack-max,
%  bisection on the period; NO gridding of the timer).  The trade-off is
%  monotone: faster feedback (larger omega) tolerates a SMALLER maximum sampling
%  period, so "maximize the period over K" is unbounded (K -> 0).  We therefore
%  fix a TARGET sampling period and report the FASTEST gain (largest bandwidth)
%  whose certified maximum period meets the target.  The chosen closed loop is
%  shown decaying under sampled feedback.  Solver: SOSTOOLS + SeDuMi.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---------------- system data ----------------
A = [0 1; 0 0.1];  B = [0;1];
Ttarget = 0.80;                          % target sampling period (design budget)
tol = 1e-6;
omegas = linspace(0.4, 2.4, 9);          % feedback-bandwidth family
fprintf('\n==== sampled-data STABILIZATION design (SOSTOOLS matrix-SOS certification) ====\n');
fprintf('A=[%g %g; %g %g]  B=[%g;%g]  target sampling period = %.2f\n', ...
    A(1,1),A(1,2),A(2,1),A(2,2),B(1),B(2),Ttarget);

% ---------------- sweep bandwidth, certify max sampling period ---------------
maxT = zeros(size(omegas));  Ks = zeros(numel(omegas),2);
for i = 1:numel(omegas)
    K = -place(A,B,[-omegas(i) -2*omegas(i)]);  Ks(i,:) = K;
    maxT(i) = maxperiod_sos(A,B,K,tol);
    fprintf('  omega=%.2f  K=[%7.3f %7.3f]  certified maxT=%.4f\n', omegas(i), K(1), K(2), maxT(i));
end

% ---------------- pick the fastest gain meeting the target period ------------
feas = find(maxT >= Ttarget);
iDes = feas(end);                        % largest bandwidth with maxT >= target
Kdes = Ks(iDes,:);  Tdes = maxT(iDes);   omDes = omegas(iDes);
fprintf('  design: omega=%.2f  K=[%.3f %.3f]  certified maxT=%.4f (>= target %.2f)\n', ...
    omDes, Kdes(1), Kdes(2), Tdes, Ttarget);

% ---------------- closed-loop decay under sampled feedback at the target -----
[tt, trX] = sd_trajectory(A,B,Kdes,Ttarget,6);

% ---------------- figure ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 11 8]);
ax1 = axes(fig); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, omegas, maxT, '-o','LineWidth',2.0,'Color',gc.order(1,:), ...
     'MarkerFaceColor',gc.order(1,:),'MarkerSize',4);
yline(ax1, Ttarget, '--','Color',gc.order(2,:),'LineWidth',1.5,'HandleVisibility','off');
plot(ax1, omDes, Tdes, 'o','MarkerSize',9,'MarkerFaceColor',gc.order(3,:), ...
     'MarkerEdgeColor','k','HandleVisibility','off');
xlabel(ax1,'feedback bandwidth $\omega$','Interpreter','latex');
ylabel(ax1,'certified max sampling period','Interpreter','latex');
title(ax1,'Design trade-off (SOS-certified)','Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(figdir,'sd_stabz_sos_tradeoff.pdf'));
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 11 8]);
ax2 = axes(fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, tt, trX, '-','LineWidth',1.8,'Color',gc.order(1,:));
set(ax2,'YScale','log');
for c = 1:6, xline(ax2, c*Ttarget, ':','Color',[0.6 0.6 0.6],'HandleVisibility','off'); end
xlabel(ax2,'time $t$','Interpreter','latex');
ylabel(ax2,'$\mathrm{tr}\,X(t)$','Interpreter','latex');
title(ax2,sprintf('Sampled feedback, $\\bar T=%.2f$',Ttarget),'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig2, fullfile(figdir,'sd_stabz_sos_decay.pdf'));

% ---------------- data macros ----------------
w = mvst_texdata(fullfile(figdir,'data','sd_stabz_sos.tex'));
w.num('SdStabzTbar', Tdes);
w.num('SdStabzTarget', Ttarget);
w.num('SdStabzOmega', omDes);
w.num('SdStabzGain', Kdes(1));
w.num('SdStabzGainTwo', Kdes(2));
w.close();
fprintf('  wrote sd_stabz_sos.pdf and data/sd_stabz_sos.tex\n');
end

% =========================================================================
function T = maxperiod_sos(A,B,K,tol)
%MAXPERIOD_SOS  Largest sampling period certified by the degree-4 SOSTOOLS
%  matrix-SOS sampled-data analysis (bisection on the slack sign).
Atil = [A B; zeros(1,3)];  Jtil = [eye(2) [0;0]; K 0];
lo = 0.05; hi = 6.0;
if slackval(Atil,Jtil,lo) <= tol, T = 0; return; end
if slackval(Atil,Jtil,hi) >  tol, T = hi; return; end
for it = 1:22
    mid = 0.5*(lo+hi);
    if slackval(Atil,Jtil,mid) > tol, lo = mid; else, hi = mid; end
end
T = lo;
end

% =========================================================================
function s = slackval(Atil,Jtil,Tbar)
%SLACKVAL  Optimal slack of the sampled-data analysis as a PURE SOSTOOLS matrix
%  sum-of-squares in the normalized timer theta = tau/Tbar in [0,1] (degree-4
%  certificate, interval multiplier theta(1-theta)).  Feasible iff s > 0.
n = size(Atil,1);
pvar th; sl = dpvar('sl');
prog = sosprogram(th);  prog = sosdecvar(prog, sl);
[prog,P]  = sospolymatrixvar(prog, monomials(th,0:4), [n n], 'symmetric');
[prog,S1] = sospolymatrixvar(prog, monomials(th,0:2), [n n], 'symmetric');
[prog,L1] = sospolymatrixvar(prog, monomials(th,0:2), [n n], 'symmetric');
[prog,V1] = sospolymatrixvar(prog, monomials(th,0:2), [n n], 'symmetric');
gt = th*(1-th);
Phi = (1/Tbar)*diff(P,th) + Atil.'*P + P*Atil;
prog = sosmatrixineq(prog, -Phi - sl*eye(n) - gt*S1);   prog = sosmatrixineq(prog, S1);
prog = sosmatrixineq(prog, P - 1e-3*eye(n) - gt*L1);    prog = sosmatrixineq(prog, L1);
prog = sosmatrixineq(prog, eye(n) - P - gt*V1);         prog = sosmatrixineq(prog, V1);
P0 = subs(P,th,0);  P1 = subs(P,th,1);
prog = sosmatrixineq(prog, P1 - Jtil.'*P0*Jtil - sl*eye(n));
prog = sossetobj(prog, -sl);
o.solver = 'sedumi';  o.params.fid = 0;
prog = sossolve(prog, o);
s = double(sosgetsol(prog, sl));
end

% =========================================================================
function [tt, trX] = sd_trajectory(A,B,K,Tbar,Ncyc)
%SD_TRAJECTORY  Trace of the plant second moment E[x x'] under sampled feedback
%  u = K x(t_k), propagated exactly on the augmented lifted symmetric space
%  (flow Atil over [0,Tbar] then sampling Jtil), started from x0 = [1;0].
Atil = [A B; zeros(1,3)];  Jtil = [eye(2) [0;0]; K 0];
na = size(Atil,1);
U = mvst_svecbasis(na);
MJ = mvst_liftsym(zeros(na), {Jtil});
LF = mvst_liftsym(Atil, {});
x0 = [1;0];  xi0 = [x0; K*x0];
sx = U'*reshape(xi0*xi0.',[],1);
tt = []; trX = []; t0 = 0; npt = 30;
for c = 1:Ncyc
    ts = linspace(0,Tbar,npt); sxc = sx;
    for j = 1:npt
        if j > 1, sxc = expm(LF*(ts(j)-ts(j-1)))*sxc; end
        Xi = reshape(U*sxc, na, na);
        tt(end+1) = t0+ts(j); trX(end+1) = trace(Xi(1:2,1:2)); %#ok<AGROW>
    end
    sx = MJ*sxc; t0 = t0 + Tbar;
end
end
