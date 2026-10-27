function ex_stabz_sdrange()
%EX_STABZ_SDRANGE  Aperiodic (range) sampled-data control: conservatism of the
%  CONSTANT refresh gain versus the PERIOD-SCHEDULED refresh gain, Stabilization
%  chapter (Proposition prop:stabz:sdrange vs Corollary cor:stabz:sdrange:const;
%  the author's "How conservative is that?" on the corollary).
%
%  Master data (the double integrator of
%  ex_stabz_sd, reused VERBATIM): plant
%     dx = (A x + B u) dt + (Ac x + B1 u) dW1,   u = K_J x(t_k) held,
%  with A = [0 1; 2 0.1], B = [0;1], Ac = 0.2 I2, B1 = 0.1[0;1], center period
%  Tc = 0.25.  The command is refreshed aperiodically with inter-sample period
%  in a window [Tc - h, Tc + h]; the multiplicative noise is kept as plant data,
%  there is no exogenous disturbance to switch off.
%
%  THREE SYNTHESIS PROGRAMS (mvst_stabz_sd_range, matrix-SOS in the timer):
%    SCHEDULED (Proposition prop:stabz:sdrange) carries a period-dependent
%      refresh matrix U_J(theta) and a timer certificate Y(tau) on [0,Tc+h],
%      the jump LMI holding over the whole window: the recovered gain
%      K_J(theta) = U_J(theta) Y(theta)^{-1} is scheduled on the realized period.
%    CONSTANT (Corollary cor:stabz:sdrange:const) freezes the certificate on the
%      window (Ydot = 0 on [Tmin,Tmax]) and holds U_J constant, so the jump-LMI
%      family collapses to one and the gain is schedule-free.
%    SLACK (dilated/Finsler constant-gain, the author's LPV-paper device, Briat
%      NAHS 41 (2021) Thm 18) keeps the certificate Y(theta) FREE (period-
%      dependent) across the window but recovers a CONSTANT gain K_J = U Ys^{-1}
%      through a single constant slack Ys and constant U in a dilated jump LMI.
%      This lifts the frozen-certificate conservatism of CONSTANT: a genuinely
%      schedule-free gain feasible over a POSITIVE window.
%
%  FINDINGS.  For the held-input sampled-data lift the flow carries NO control
%  (B0 = 0) and the held command is marginal (zero u-block); the SCHEDULED
%  program absorbs the flow expansion in a GROWING timer certificate and stays
%  SOS-feasible across essentially the whole sampling range, while the CONSTANT
%  program, whose FROZEN certificate cannot certify the marginal held-input flow
%  over a nonzero window (a frozen certificate cannot dominate an expanding flow
%  without a flow gain), is SOS-infeasible for any h > 0: its critical half-width
%  is ~0.  That is the conservatism gap, and it lives entirely in the CERTIFICATE
%  (the corollary's frozen-certificate convex restriction), NOT in the constant
%  gain itself: certified by the exact one-period monodromy on the lifted
%  symmetric space (mvst_liftsym/mvst_svecbasis), the fixed refresh gain of the
%  periodic (h=0) design keeps the spectral radius below 1 across the whole
%  tested range, essentially as well as the scheduled gain.  The frozen-
%  certificate corollary is thus a very loose certificate for a gain that in fact
%  works range-wide; the period-scheduled program certifies that gain tightly.
%
%  DELIVERABLES: (i) SOS feasibility of the three programs over a width sweep;
%  (ii) the SOS critical half-width of each (scheduled >= slack >= constant);
%  (iii) at a representative window the scheduled gains K_J(Tmin), K_J(Tc),
%  K_J(Tmax) vs the constant and slack gains, the scheduled and slack certified
%  slacks, and the exact worst monodromy radius of each; plus the exact
%  robustness half-width of the fixed constant gain.  Solver: SOSTOOLS + SeDuMi.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---------------- master data (verbatim from ex_stabz_sd) ----------------
A = [0 1; 2 0.1]; B = [0; 1]; Ac = 0.2*eye(2);
B1 = 0.1*[0; 1];                 % input selector actuating the diffusion channel
Tc = 0.25;                       % center period (the master's bar T)
deg = 4;                         % timer polynomial degree
tol = 1e-7;
plant = struct('A',A,'Bu',B,'Ac',Ac,'Bc',B1);

fprintf('\n==== range sampled-data synthesis: scheduled vs constant refresh gain ====\n');
fprintf('A=[%g %g;%g %g] B=[%g;%g] Ac=0.2 I2 B1=0.1[0;1] Tc=%g deg=%d\n', ...
    A(1,1),A(1,2),A(2,1),A(2,2),B(1),B(2),Tc,deg);

solve = @(h,mode) mvst_stabz_sd_range(plant, ...
    struct('Tmin',Tc-h,'Tmax',Tc+h,'mode',mode,'deg',deg));

% ---------------- (i) width sweep: SOS feasibility of both designs -------------
wfac  = [0 0.2 0.4 0.6 0.8];     % width w = wfac*Tc, half-width h = w/2
halfs = 0.5*wfac*Tc;
nsw   = numel(halfs);
slkS = zeros(1,nsw); slkC = zeros(1,nsw); slkX = zeros(1,nsw);
for i = 1:nsw
    slkS(i) = solve(halfs(i),'scheduled').slack;
    slkC(i) = solve(halfs(i),'constant').slack;
    slkX(i) = solve(halfs(i),'slack').slack;
    fprintf(['(i)  h=%.4f (w=%.2f Tc): scheduled slack=%+.3e feas=%d | constant slack=%+.3e feas=%d' ...
        ' | slack slack=%+.3e feas=%d\n'], ...
        halfs(i), wfac(i), slkS(i), slkS(i)>tol, slkC(i), slkC(i)>tol, slkX(i), slkX(i)>tol);
end
% sanity: at h=0 (periodic case) scheduled and constant must be feasible and
% identical; the slack (dilated constant-gain) design must be at least as good.
assert(slkS(1) > tol, 'scheduled design infeasible at h=0 (periodic case); check assembly.');
assert(slkC(1) > tol, 'constant design infeasible at h=0 (periodic case); check assembly.');
assert(slkX(1) > tol, 'slack design infeasible at h=0 (periodic case); check assembly.');
assert(slkX(1) >= slkC(1) - 1e-6, ...
    'slack (dilated) design must be at least as good as constant at h=0.');

% ---------------- (ii) SOS critical half-width by bisection --------------------
hcap  = 0.90*Tc;                 % keep Tmin > 0 strictly
critS = crit_halfwidth(solve, 'scheduled', hcap, tol);
critC = crit_halfwidth(solve, 'constant',  hcap, tol);
critX = crit_halfwidth(solve, 'slack',     hcap, tol);
gapSos = critS - critC;
fprintf('(ii) SOS critical half-width: scheduled h*=%.4f | constant h*=%.4f | slack h*=%.4f | gap=%.4f\n', ...
    critS, critC, critX, gapSos);
assert(critS >= critC - 1e-6, 'scheduled SOS critical half-width must be >= constant.');
assert(critX >= critC - 1e-6, 'slack SOS critical half-width must be >= constant.');

% ---------------- (iii) representative window: slacks, gains, certification ----
hrep = 0.2*Tc;                   % representative half-width (w = 0.4 Tc)
Tmin = Tc - hrep; Tmax = Tc + hrep;
rS = solve(hrep,'scheduled');    % scheduled feasible here
r0 = solve(0,'constant');        % constant design at its ONLY feasible window (h=0)
rX = solve(hrep,'slack');        % slack (dilated constant-gain) feasible here
fprintf(['(iii) representative window [%.4f, %.4f]: scheduled slack=%.3e | slack slack=%.3e' ...
    ' (constant infeasible for h>0)\n'], Tmin, Tmax, rS.slack, rX.slack);

KsMin = rS.KJfun(Tmin);  KsTc = rS.KJfun(Tc);  KsMax = rS.KJfun(Tmax);
Kconst = r0.KJfun(Tc);           % schedule-free constant gain (h=0 periodic design)
Kslack = rX.KJfun(Tc);           % dilated schedule-free constant gain (feasible over window)
fprintf('   scheduled K_J(Tmin)=[%s]\n', num2str(KsMin,'%8.3f'));
fprintf('   scheduled K_J(Tc)  =[%s]\n', num2str(KsTc, '%8.3f'));
fprintf('   scheduled K_J(Tmax)=[%s]\n', num2str(KsMax,'%8.3f'));
fprintf('   constant  K_J      =[%s]\n', num2str(Kconst,'%8.3f'));
fprintf('   slack     K_J      =[%s]\n', num2str(Kslack,'%8.3f'));

% independent certification: worst monodromy radius over the representative window
radS = window_radius(A,B,Ac,B1, rS.KJfun, Tmin, Tmax, false);   % scheduled
radC = window_radius(A,B,Ac,B1, @(t) Kconst, Tmin, Tmax, true); % fixed constant gain
radX = window_radius(A,B,Ac,B1, @(t) Kslack, Tmin, Tmax, true); % dilated constant gain
fprintf('   exact worst monodromy radius over rep window: scheduled=%.4f constant=%.4f slack=%.4f\n', ...
    radS, radC, radX);
assert(radS < 1, 'recovered scheduled design must certify (monodromy radius < 1) on the window.');
assert(radX < 1, 'recovered slack (dilated constant-gain) design must certify (radius < 1) on the window.');

% exact robustness half-width of the FIXED constant gain (radius < 1)
constExact = crit_radius_halfwidth(A,B,Ac,B1, @(t) Kconst, Tc, hcap);
fprintf('   fixed constant-gain exact robustness half-width h_c=%.4f\n', constExact);

% ---------------- figures ----------------
fig1 = figure('Color','w','Units','centimeters','Position',[2 2 11 8]);
ax1 = axes(fig1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, halfs, slkS, '-o','LineWidth',2.0,'Color',gc.order(1,:),'MarkerFaceColor',gc.order(1,:));
plot(ax1, halfs, slkC, '-s','LineWidth',2.0,'Color',gc.order(2,:),'MarkerFaceColor',gc.order(2,:));
plot(ax1, halfs, slkX, '-^','LineWidth',2.0,'Color',gc.order(3,:),'MarkerFaceColor',gc.order(3,:));
yline(ax1,0,'--','Color',gc.onyx,'LineWidth',1.2,'HandleVisibility','off');
xlabel(ax1,'window half-width $h$','Interpreter','latex');
ylabel(ax1,'certified slack','Interpreter','latex');
legend(ax1,{'scheduled','constant','slack'},'Location','northeast','Interpreter','latex');
title(ax1,'range sampled-data: scheduled vs constant vs slack','Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig1, fullfile(figdir,'stabz_sdrange_slack.pdf'));

fig2 = figure('Color','w','Units','centimeters','Position',[2 2 11 8]);
ax2 = axes(fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
thg = linspace(Tmin,Tmax,40);
Ksched = cell2mat(arrayfun(@(t) rS.KJfun(t), thg', 'UniformOutput', false));
plot(ax2, thg, Ksched(:,1), '-','LineWidth',2.0,'Color',gc.order(1,:));
plot(ax2, thg, Ksched(:,2), '-','LineWidth',2.0,'Color',gc.order(2,:));
plot(ax2, [Tmin Tmax], Kconst(1)*[1 1], '--','LineWidth',1.6,'Color',gc.order(1,:),'HandleVisibility','off');
plot(ax2, [Tmin Tmax], Kconst(2)*[1 1], '--','LineWidth',1.6,'Color',gc.order(2,:),'HandleVisibility','off');
plot(ax2, [Tmin Tmax], Kslack(1)*[1 1], ':','LineWidth',1.6,'Color',gc.order(1,:),'HandleVisibility','off');
plot(ax2, [Tmin Tmax], Kslack(2)*[1 1], ':','LineWidth',1.6,'Color',gc.order(2,:),'HandleVisibility','off');
xlabel(ax2,'realized period $\theta$','Interpreter','latex');
ylabel(ax2,'refresh gain $K_J(\theta)$','Interpreter','latex');
legend(ax2,{'$K_{J,x_1}(\theta)$','$K_{J,x_2}(\theta)$'},'Location','best','Interpreter','latex');
title(ax2,'scheduled (solid), constant (dashed), slack (dotted)','Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig2, fullfile(figdir,'stabz_sdrange_gain.pdf'));

% ---------------- data macros ----------------
w = mvst_texdata(fullfile(figdir,'data','ex_stabz_sdrange.tex'));
w.num('StabzSdrTc',   Tc);
w.num('StabzSdrDeg',  deg);
w.num('StabzSdrCritSched', critS);
w.num('StabzSdrCritConst', critC);
w.num('StabzSdrGap',       gapSos);
w.num('StabzSdrRepHalf',   hrep);
w.num('StabzSdrTmin',      Tmin);
w.num('StabzSdrTmax',      Tmax);
w.num('StabzSdrSlackSched', rS.slack);
w.mat('StabzSdrKsTmin', KsMin);
w.mat('StabzSdrKsTc',   KsTc);
w.mat('StabzSdrKsTmax', KsMax);
w.mat('StabzSdrKconst', Kconst);
w.num('StabzSdrRadSched', radS);
w.num('StabzSdrRadConst', radC);
w.num('StabzSdrConstExact', constExact);
% slack (dilated constant-gain, period-dependent certificate) mode
w.num('StabzSdrCritSlack',  critX);
w.num('StabzSdrSlackSlack', rX.slack);
w.mat('StabzSdrKslack',     Kslack);
w.num('StabzSdrRadSlack',   radX);
w.close();
fprintf('   wrote stabz_sdrange_slack.pdf, stabz_sdrange_gain.pdf and data/ex_stabz_sdrange.tex\n');
end

% =========================================================================
function hstar = crit_halfwidth(solve, mode, hmax, tol)
%CRIT_HALFWIDTH  Largest window half-width at which the MODE SOS program stays
%  feasible, by bisection on the sign of the optimal slack (feasible at h=0).
if solve(hmax,mode).slack > tol; hstar = hmax; return; end     % feasible to the cap
lo = 0; hi = hmax;
for it = 1:16
    mid = 0.5*(lo+hi);
    if solve(mid,mode).slack > tol; lo = mid; else; hi = mid; end
end
hstar = lo;
end

% =========================================================================
function hstar = crit_radius_halfwidth(A,B,Ac,B1, KJfun, Tc, hmax)
%CRIT_RADIUS_HALFWIDTH  Largest half-width h such that the recovered design keeps
%  the exact worst one-period monodromy radius over [Tc-h,Tc+h] below 1, by
%  bisection.  KJfun is applied as a fixed gain (isConst=true) here.
if window_radius(A,B,Ac,B1, KJfun, Tc-hmax, Tc+hmax, true) < 1; hstar = hmax; return; end
lo = 0; hi = hmax;
for it = 1:22
    mid = 0.5*(lo+hi);
    if window_radius(A,B,Ac,B1, KJfun, Tc-mid, Tc+mid, true) < 1; lo = mid; else; hi = mid; end
end
hstar = lo;
end

% =========================================================================
function rad = window_radius(A,B,Ac,B1, KJfun, Tmin, Tmax, isConst)
%WINDOW_RADIUS  Worst one-period monodromy spectral radius over [Tmin,Tmax] for
%  the recovered refresh design, on the lifted symmetric space.  For a realized
%  period T the closed loop flows tilde A over [0,T] then applies the refresh
%  jump G(T) = J_0 + B_J K_J(T); the scheduled gain uses K_J(T), a fixed gain
%  the single K_J.
n = size(A,1); m = size(B,2); N = n+m;
Atil  = [A B; zeros(m,N)];
Actil = [Ac B1; zeros(m,N)];
J0    = [eye(n) zeros(n,m); zeros(m,N)];
BJ    = [zeros(n,m); eye(m)];
Mf = mvst_liftsym(Atil, {Actil});
Ud = mvst_svecbasis(N);
Ts = linspace(Tmin, Tmax, 25);
rad = 0;
for T = Ts
    if isConst; K = KJfun(Tmin); else; K = KJfun(T); end
    G  = J0 + BJ*K;
    MJ = Ud' * kron(G,G) * Ud;
    rad = max(rad, max(abs(eig(MJ*expm(T*Mf)))));
end
end
