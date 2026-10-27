function Tstar = ex_stab_switched_sos()
%EX_STAB_SWITCHED_SOS  Minimum-dwell-time SWITCHED STABILITY ANALYSIS through the
%  generic mode-indexed matrix-SOS analysis program MVST_STAB_SOS_MODES
%  (Stability chapter, Theorems th:stab:periodic / th:stab:MDT).
%
%  Two-mode switched system  dx = A_sigma x,  state continuous at switches
%  (identity jump).  Second-moment flow operator per mode  F_i(X) = He(A_i X),
%  adjoint  F_i*(P) = A_i' P + P A_i.  Both modes are Hurwitz but their average
%  is unstable, so FAST switching (short dwell-time) is unstable while a
%  sufficiently long dwell-time is stable: a genuine minimum dwell-time Tbar*.
%
%  DUAL (Lyapunov) ANALYSIS with mode-indexed certificates.  On the timer,
%  normalized to theta = tau/Tbar in [0,1], each mode carries a symmetric
%  polynomial certificate P_i(theta) (degree 4).  The minimum-dwell-time
%  conditions are
%     flow (per mode i):   (1/Tbar) P_i' + A_i' P_i + P_i A_i <= -s I  on [0,1],
%     switch (i -> j!=i):  P_j(0) - P_i(1) <= -s I.
%  These are exactly the mode-indexed timer conditions of MVST_STAB_SOS_MODES
%  with dwell type 'switched': a single-certificate program cannot express the
%  switch coupling, so the product-cone routine is used.
%
%  METHOD.  The mode-indexed matrix sum-of-squares in the timer is assembled and
%  solved by MVST_STAB_SOS_MODES (no gridding of the timer): each flow condition
%  is a matrix-SOS Positivstellensatz on [0,1] with the interval multiplier
%  theta(1-theta); the bounds 1e-3 I <= P_i <= I (positivity + scale
%  normalization) are imposed the same way; the switch conditions are constant
%  matrix inequalities.  Feasibility is decided by common-slack maximization
%  (feasible iff slack > 0); the minimum dwell-time Tbar* is found by BISECTION
%  ON THE SOS FEASIBILITY.  CROSS-CHECKED against the exact two-mode monodromy on
%  the lifted symmetric space (expm(A1 Tbar) then expm(A2 Tbar), lifted), whose
%  spectral radius crosses 1 at Tbar*; the SOS and exact critical dwell-times
%  must agree to ~1e-2 relative.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---------------- system data ----------------
A1 = [-1 0; 2.5 -1];  A2 = [-1 2.5; 0 -1];
tol = 1e-6;
fprintf('\n==== minimum-dwell-time SWITCHED stability analysis (mode-indexed matrix-SOS in tau) ====\n');
fprintf('A1 eig %s   A2 eig %s   (A1+A2) eig %s\n', ...
    mat2str(eig(A1)',3), mat2str(eig(A2)',3), mat2str(eig(A1+A2)',3));

% ---------------- minimum dwell-time by bisection on the SOS feasibility ------
lo = 0.3; hi = 3.0;
assert(slackval(A1,A2,lo) <= tol, 'lower bracket already feasible');
assert(slackval(A1,A2,hi) >  tol, 'upper bracket infeasible');
for it = 1:22
    mid = 0.5*(lo+hi);
    if slackval(A1,A2,mid) > tol, hi = mid; else, lo = mid; end
end
Tstar = hi;
fprintf('minimum dwell-time Tbar* = %.4f  (slack there = %.2e)\n', Tstar, slackval(A1,A2,Tstar));

% ---------------- cross-check: exact monodromy critical dwell-time -----------
Texact = mono_crit(A1,A2);
fprintf('  cross-check exact monodromy critical dwell-time = %.4f  (SOS %.4f, rel. dev. %.2e)\n', ...
    Texact, Tstar, abs(Tstar-Texact)/Texact);
assert(abs(Tstar-Texact)/Texact < 1e-2, ...
    'SOS min dwell-time disagrees with the exact monodromy beyond 1e-2 relative');

% ---------------- verification: two-mode monodromy radius vs dwell-time ------
Tgrid = linspace(0.5*Tstar, 1.8*Tstar, 50);
rad = arrayfun(@(T) mono_radius(A1,A2,T), Tgrid);
radAt = mono_radius(A1,A2,Tstar);
radBelow = mono_radius(A1,A2,0.8*Tstar);
radAbove = mono_radius(A1,A2,1.2*Tstar);
fprintf('  monodromy radius: at 0.8 Tbar*=%.3f, at Tbar*=%.3f, at 1.2 Tbar*=%.3f\n', ...
    radBelow, radAt, radAbove);

% ---------------- figure ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, Tgrid, rad, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(ax,1,'--','Color',gc.order(2,:),'LineWidth',1.5,'HandleVisibility','off');
xline(ax,Tstar,':','Color',gc.onyx,'LineWidth',1.5,'HandleVisibility','off');
plot(ax,Tstar,radAt,'o','MarkerSize',7,'MarkerFaceColor',gc.order(3,:), ...
     'MarkerEdgeColor','k','HandleVisibility','off');
xlabel(ax,'dwell-time $\bar T$','Interpreter','latex');
ylabel(ax,'two-cycle monodromy spectral radius','Interpreter','latex');
title(ax,sprintf('Switched minimum dwell-time $\\bar T^\\star=%.3f$',Tstar), ...
    'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(figdir,'sw_stab_sos.pdf'));

% ---------------- data macros ----------------
w = mvst_texdata(fullfile(figdir,'data','sw_stab_sos.tex'));
w.num('SwStabTbar', Tstar);
w.num('SwStabRadAt', radAt);
w.num('SwStabRadBelow', radBelow);
w.num('SwStabRadAbove', radAbove);
w.mat('SwStabAone', A1);
w.mat('SwStabAtwo', A2);
w.close();
fprintf('  wrote sw_stab_sos.pdf and data/sw_stab_sos.tex\n');
end

% =========================================================================
function s = slackval(A1,A2,Tbar)
%SLACKVAL  Optimal common slack of the two-mode switched analysis, assembled and
%  solved by the generic mode-indexed matrix-SOS program MVST_STAB_SOS_MODES
%  (dwell type 'switched'): mode-indexed timer certificates P_i(theta) on [0,1]
%  (interval multiplier theta(1-theta)) with the all-pairs switch conditions.
%  Feasible iff s > 0.
sys = struct('A0', cat(3, A1, A2));
res = mvst_stab_sos_modes(sys, struct('type','switched','Tbar',Tbar,'deg',4));
s = res.slack;
end

% =========================================================================
function Tstar = mono_crit(A1,A2)
%MONO_CRIT  Exact minimum dwell-time from the two-mode monodromy: smallest Tbar
%  whose one-cycle lifted spectral radius crosses 1 (bisection).
lo = 0.05; hi = 4.0;
for it = 1:40
    mid = 0.5*(lo+hi);
    if mono_radius(A1,A2,mid) < 1, hi = mid; else, lo = mid; end
end
Tstar = 0.5*(lo+hi);
end

% =========================================================================
function r = mono_radius(A1,A2,Tbar)
%MONO_RADIUS  Exact two-mode monodromy spectral radius on the lifted symmetric
%  space over one switch cycle: flow mode 1 for Tbar, then flow mode 2 for Tbar.
L1 = mvst_liftsym(A1, {});  L2 = mvst_liftsym(A2, {});
r = max(abs(eig(expm(L2*Tbar)*expm(L1*Tbar))));
end
