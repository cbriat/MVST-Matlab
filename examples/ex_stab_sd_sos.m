function Tstar = ex_stab_sd_sos()
%EX_STAB_SD_SOS  Maximum-sampling-period SAMPLED-DATA STABILITY ANALYSIS through
%  the generic matrix-SOS analysis program MVST_STAB_SOS in the timer (Stability
%  chapter, Theorem th:stab:periodic and the ranged-dwell Theorem th:stab:RDT;
%  the sampled-data lifting, Sec. "Lifting and the sampled stability condition").
%
%  Continuous plant  dx = A x dt + B u  under a zero-order-hold sampled state
%  feedback  u(t) = K x(t_k)  on [t_k, t_{k+1}).  Modelled as an impulsive
%  system on the augmented state  xi = [x; u_hold]:
%     flow:     xidot = Atil xi,   Atil = [A B; 0 0],
%     sampling: xi+  = Jtil xi,    Jtil = [I 0; K 0]   (the hold is refreshed).
%  Second-moment flow operator  F(X) = Atil X + X Atil', adjoint
%     F*(P) = Atil' P + P Atil;   sampling jump adjoint  Jtil' P Jtil.
%  The continuous closed loop A+BK is stable, so a small sampling period is
%  stable and a large one is not: there is a MAXIMUM stable sampling period
%  Tbar*.
%
%  DUAL (Lyapunov) ANALYSIS.  On the timer, normalized to theta = tau/Tbar in
%  [0,1], the certificate P(theta) is a symmetric polynomial matrix (degree 4,
%  the degenerate augmented flow needs degree 4 for a tight bound), and the
%  periodic-sampling conditions of Theorem th:stab:periodic read
%     flow (0,Tbar):  (1/Tbar) P' + Atil' P + P Atil <= -s I  on [0,1],
%     sampling:       Jtil' P(0) Jtil - P(1) <= -s I.
%
%  METHOD.  The matrix sum-of-squares in the timer is assembled and solved by the
%  generic analysis routine MVST_STAB_SOS (no gridding of the timer).  The flow
%  condition is a matrix-SOS Positivstellensatz on [0,1] with the interval
%  multiplier theta(1-theta); the certificate bounds 0 <= P <= I are imposed the
%  same way; the sampling condition is a constant matrix inequality.  Feasibility
%  is decided by common-slack maximization (feasible iff slack > 0).  The MAXIMUM
%  stable sampling period Tbar* (stable for Tbar <= Tbar*) is found by BISECTION
%  ON THE SOS FEASIBILITY (dwell type 'periodic', the main sweep).  The ranged
%  dwell test (type 'range', Theorem th:stab:RDT) additionally certifies a whole
%  admissible sampling interval [Tmin,Tmax] with a single certificate.
%  CROSS-CHECKED against the exact augmented monodromy on the lifted symmetric
%  space (expm(Atil Tbar) then Jtil(.)Jtil', lifted), whose spectral radius
%  crosses 1 at Tbar*; SOS and exact critical period must agree to ~1e-2 rel.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---------------- system data ----------------
A = [0 1; 0 0.1];  B = [0;1];            % unstable plant (eigenvalues 0, 0.1)
K = place(A,B,[-1 -2]);  K = -K;         % continuous stabilizing gain (poles -1,-2)
% place returns K with A-BK stable; our convention u=Kx needs A+BK stable
if max(real(eig(A+B*K))) >= -1e-9, K = -place(A,B,[-1 -2]); end
Atil = [A B; zeros(1,3)];
Jtil = [eye(2) [0;0]; K 0];
tol = 1e-6;
fprintf('\n==== maximum-sampling-period SAMPLED-DATA stability analysis (generic matrix-SOS) ====\n');
fprintf('A=[%g %g; %g %g]  B=[%g;%g]  K=[%.4f %.4f]  (A+BK eig %s)\n', ...
    A(1,1),A(1,2),A(2,1),A(2,2),B(1),B(2),K(1),K(2),mat2str(eig(A+B*K)',3));

% ---------------- maximum stable sampling period by SOS bisection (periodic) --
lo = 0.05; hi = 3.0;                       % lo stable (feasible), hi unstable
assert(slackval(Atil,Jtil,lo) >  tol, 'lower bracket infeasible (retune)');
assert(slackval(Atil,Jtil,hi) <= tol, 'upper bracket feasible (raise hi)');
for it = 1:24
    mid = 0.5*(lo+hi);
    if slackval(Atil,Jtil,mid) > tol, lo = mid; else, hi = mid; end
end
Tstar = lo;
fprintf('maximum stable sampling period Tbar* = %.4f  (slack there = %.2e)\n', Tstar, slackval(Atil,Jtil,Tstar));

% ---------------- cross-check: exact augmented monodromy critical period ------
Texact = mono_crit(Atil,Jtil);
fprintf('  cross-check exact monodromy critical period = %.4f  (SOS %.4f, rel. dev. %.2e)\n', ...
    Texact, Tstar, abs(Tstar-Texact)/Texact);
assert(abs(Tstar-Texact)/Texact < 1e-2, ...
    'SOS max sampling period disagrees with the exact monodromy beyond 1e-2 relative');

% ---------------- ranged-dwell certification of an admissible interval --------
Tmin = 0.05; Tmax = 0.85*Tstar;            % aperiodic sampling anywhere in here
sRange = slackval_range(Atil,Jtil,Tmin,Tmax);
fprintf('  ranged-dwell (th:stab:RDT): admissible sampling interval [%.4f, %.4f] certified (slack %.2e)\n', ...
    Tmin, Tmax, sRange);
assert(sRange > tol, 'ranged-dwell certificate infeasible on the reported interval');

% ---------------- verification: augmented monodromy radius vs period ---------
Tgrid = linspace(0.4*Tstar, 1.6*Tstar, 50);
rad = arrayfun(@(T) mono_radius(Atil,Jtil,T), Tgrid);
radAt = mono_radius(Atil,Jtil,Tstar);
radBelow = mono_radius(Atil,Jtil,0.8*Tstar);
radAbove = mono_radius(Atil,Jtil,1.2*Tstar);
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
xlabel(ax,'sampling period $\bar T$','Interpreter','latex');
ylabel(ax,'augmented monodromy spectral radius','Interpreter','latex');
title(ax,sprintf('Sampled-data maximum period $\\bar T^\\star=%.3f$',Tstar), ...
    'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(figdir,'sd_stab_sos.pdf'));

% ---------------- data macros ----------------
w = mvst_texdata(fullfile(figdir,'data','sd_stab_sos.tex'));
w.num('SdStabTbar', Tstar);
w.num('SdStabRadAt', radAt);
w.num('SdStabRadBelow', radBelow);
w.num('SdStabRadAbove', radAbove);
w.num('SdStabGainK', K(1));
w.num('SdStabGainKtwo', K(2));
w.close();
fprintf('  wrote sd_stab_sos.pdf and data/sd_stab_sos.tex\n');
end

% =========================================================================
function s = slackval(Atil,Jtil,Tbar)
%SLACKVAL  Optimal slack of the sampled-data analysis, assembled and solved by
%  the generic matrix-SOS routine MVST_STAB_SOS (dwell type 'periodic', degree-4
%  certificate P(theta) in the normalized timer theta = tau/Tbar in [0,1],
%  interval multiplier theta(1-theta)).  Feasible iff s > 0.  The augmented flow
%  Atil is the drift .A0 and the sampling map Jtil is the scheduled jump .J.
sys = struct('A0', Atil, 'J', Jtil);
res = mvst_stab_sos(sys, struct('type','periodic','Tbar',Tbar,'deg',4));
s = res.slack;
end

% =========================================================================
function s = slackval_range(Atil,Jtil,Tmin,Tmax)
%SLACKVAL_RANGE  Ranged-dwell certificate (dwell type 'range'): a single
%  certificate that certifies the sampling jump for EVERY interval length in
%  [Tmin,Tmax] (aperiodic sampling), the analysis form of Theorem th:stab:RDT.
sys = struct('A0', Atil, 'J', Jtil);
res = mvst_stab_sos(sys, struct('type','range','Tmin',Tmin,'Tmax',Tmax,'deg',4));
s = res.slack;
end

% =========================================================================
function Tstar = mono_crit(Atil,Jtil)
%MONO_CRIT  Exact maximum sampling period from the augmented monodromy: largest
%  Tbar whose one-period lifted spectral radius is still below 1 (bisection).
lo = 0.02; hi = 4.0;
for it = 1:40
    mid = 0.5*(lo+hi);
    if mono_radius(Atil,Jtil,mid) < 1, lo = mid; else, hi = mid; end
end
Tstar = 0.5*(lo+hi);
end

% =========================================================================
function r = mono_radius(Atil,Jtil,Tbar)
%MONO_RADIUS  Exact augmented one-period monodromy spectral radius on the lifted
%  symmetric space: flow expm(Atil Tbar) then the sampling jump Jtil(.)Jtil'.
LF = mvst_liftsym(Atil, {});
MJ = mvst_liftsym(zeros(size(Atil)), {Jtil});
r = max(abs(eig(MJ*expm(LF*Tbar))));
end
