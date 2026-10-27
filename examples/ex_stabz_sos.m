function Tstar = ex_stabz_sos()
%EX_STABZ_SOS  Minimum-dwell-time state-feedback synthesis for a
%  deterministic-plus-multiplicative impulsive system, Stabilization chapter
%  (Theorem th:stab:stabzMDT / th:stab:stabzdwell).
%
%  System (state dimension 2):
%     flow:  dx = (A0 x + B0 u) dt + (A1 x + B1 u) dW1,   A1 = sigma*I,
%     jump:  x+ = J x,   J = jfac*I   (expansive, uncontrolled, no jump input).
%  Second-moment closed-loop flow operator with flow gain K0:
%     F_K(X) = He((A0+B0 K0) X) + (A1+B1 K0) X (A1+B1 K0)',   jump  J X J'.
%
%  TUNING / FINDING.  With the CLEAN fully-actuated flow (B1 = 0) the pair
%  (A0,B0) is controllable, the deterministic part is arbitrarily
%  pole-placeable, and the state-multiplicative noise A1 = sigma*I only adds a
%  constant shift sigma^2 to the lifted spectrum: the achievable second-moment
%  decay is UNBOUNDED and no positive minimum dwell-time exists (Tbar* -> 0).
%  A finite Tbar* (and, in ex_stabz_rate, a finite guaranteed rate) requires a
%  bounded achievable decay.  A small INPUT-multiplicative noise B1 = beta*[0;1]
%  supplies it: large gains inject noise proportional to u, so the second-moment
%  decay is bounded.  This channel is exactly the term V_l = A_l Y + B_l U_l of
%  the synthesis blocks eq:stab:stabzdwell, so it stays inside the theorem.
%  A0,B0,A1(sigma),J are as specified, beta = 0.62 added.
%
%  On the timer, normalized to theta = tau/Tbar in [0,1] for conditioning, the
%  synthesis carries POLYNOMIAL matrices Y(theta) in S^2 (symmetric, degree 2)
%  and U0(theta) (degree 2), and imposes
%     flow block (0,Tbar):  Phi_f = [ -(1/Tbar) Y' + He(A0 Y + B0 U0) , (A1 Y+B1 U0)' ;
%                                       A1 Y+B1 U0                     , -Y            ] <= -s I,
%     jump block  theta=1:  Phi_j = [ -Y(1) , (J Y(1))' ; J Y(1) , -Y(0) ] <= -s I,
%  with recovered flow gain K0(theta) = U0(theta) Y(theta)^{-1}.
%
%  METHOD.  PURE SOSTOOLS MATRIX SUM-OF-SQUARES IN THE TIMER (no gridding of the
%  timer anywhere in the synthesis).  The interval condition Phi_f(theta) <= -s I
%  on [0,1] is a matrix-SOS Positivstellensatz certificate
%     -Phi_f(theta) - s I = Sigma0(theta) + theta(1-theta) Sigma1(theta),
%  with Sigma0, Sigma1 matrix-SOS (sosmatrixineq: v'(.)v is SOS in (theta,v));
%  the certificate Y(theta) <= I (a scale normalization) is imposed the same
%  way.  The jump block at theta = 1 is a single constant LMI (also an
%  sosmatrixineq).  Feasibility is decided by SLACK MAXIMIZATION (maximize the
%  scalar margin s; feasible iff s > 0) which keeps SeDuMi well conditioned.
%  The MINIMUM dwell-time Tbar* is found by bisection on the sign of the optimal
%  slack.  The synthesized closed loop is VERIFIED by the exact one-period
%  monodromy on the lifted symmetric space (mvst_liftsym integrated with expm
%  over [0,Tbar*] using the recovered polynomial gain, then the jump J(.)J'):
%  its spectral radius crosses 1 at Tbar*.  Solver: SOSTOOLS + SeDuMi.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---------------- system data ----------------
A0 = [0 1; 2 0.1];  B0 = [0;1];  sigma = 0.3;  A1 = sigma*eye(2);
beta = 0.62;  B1 = beta*[0;1];
jfac = 1.1;   J = jfac*eye(2);
tol = 1e-7;   % slack sign tolerance

fprintf('\n==== minimum-dwell-time stabilization synthesis (SOSTOOLS matrix-SOS in tau) ====\n');
fprintf('A0=[%g %g; %g %g]  B0=[%g;%g]  sigma=%g  beta(B1)=%g  J=%g*I\n', ...
    A0(1,1),A0(1,2),A0(2,1),A0(2,2),B0(1),B0(2),sigma,beta,jfac);

% ---------------- (i) arbitrary-dwell / constant-certificate INFEASIBLE ----
sC = feas_constant(A0,B0,A1,B1,J);
fprintf('(i)  arbitrary-dwell (constant Y) synthesis slack=%.3e -> feasible? %d  (expected 0)\n', ...
    sC, sC>tol);

% ---------------- (ii) flow-alone second-moment stabilizable ----------------
[Kflow, absFlow] = flow_gain(A0,B0,A1,B1);
fprintf('(ii) flow-alone min second-moment abscissa with K0=[%.3f %.3f]: %.4f  (<0 => stabilizable)\n', ...
    Kflow(1),Kflow(2),absFlow);

% ---------------- (iii) minimum dwell-time by bisection on the slack sign ---
lo = 0.1; hi = 1.0;
assert(synth_slack(A0,B0,A1,B1,J,lo,0,1).slack <= tol, 'lower bracket already feasible');
assert(synth_slack(A0,B0,A1,B1,J,hi,0,1).slack >  tol, 'upper bracket infeasible');
for it = 1:22
    mid = 0.5*(lo+hi);
    if synth_slack(A0,B0,A1,B1,J,mid,0,1).slack > tol, hi = mid; else, lo = mid; end
end
Tstar = hi;
sol = synth_slack(A0,B0,A1,B1,J,Tstar,0,1);
sol = build_gain_profile(sol);
fprintf('(iii) minimum dwell-time Tbar* = %.4f  (slack there = %.2e)\n', Tstar, sol.slack);

% ---------------- verification: monodromy spectral radius vs Tbar -----------
Tgrid = linspace(0.6*Tstar, 1.8*Tstar, 40);
rad = arrayfun(@(T) mono_radius(A0,B0,A1,B1,J,sol,T), Tgrid);
radAt = mono_radius(A0,B0,A1,B1,J,sol,Tstar);
radBelow = mono_radius(A0,B0,A1,B1,J,sol,0.85*Tstar);
radAbove = mono_radius(A0,B0,A1,B1,J,sol,1.15*Tstar);
fprintf('   monodromy radius: at 0.85 Tbar*=%.3f, at Tbar*=%.3f, at 1.15 Tbar*=%.3f\n', ...
    radBelow, radAt, radAbove);

% ---------------- closed-loop second-moment trajectory at Tbar* -------------
[tt, trX] = cl_trajectory(A0,B0,A1,B1,J,sol,6);

% representative flow gain at tau = Tbar* (theta = 1)
Kstar = gain_at(sol, 1.0);
fprintf('   representative flow gain K0(Tbar*) = [%.3f  %.3f]\n', Kstar(1), Kstar(2));

% ---------------- figure ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 11 8]);
ax1 = axes(fig); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, Tgrid, rad, '-', 'LineWidth',2.0,'Color',gc.order(1,:));
yline(ax1,1,'--','Color',gc.order(2,:),'LineWidth',1.5,'HandleVisibility','off');
xline(ax1,Tstar,':','Color',gc.onyx,'LineWidth',1.5,'HandleVisibility','off');
plot(ax1, Tstar, radAt, 'o','MarkerSize',7,'MarkerFaceColor',gc.order(3,:), ...
     'MarkerEdgeColor','k','HandleVisibility','off');
xlabel(ax1,'dwell-time $\bar T$','Interpreter','latex');
ylabel(ax1,'monodromy spectral radius','Interpreter','latex');
title(ax1,sprintf('$\\bar T^\\star=%.3f$',Tstar),'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(figdir,'stabz_sos_radius.pdf'));
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 11 8]);
ax2 = axes(fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, tt, trX, '-', 'LineWidth',1.8,'Color',gc.order(1,:));
set(ax2,'YScale','log');
% pad the log-scale y-limits by 3% of the decade range at both ends so the
% trace clears the top and bottom borders (same padding in the combined figure)
ly = log10([min(trX) max(trX)]);  ly = ly + 0.03*diff(ly)*[-1 1];
ylim(ax2, 10.^ly);
xlabel(ax2,'time $t$','Interpreter','latex');
ylabel(ax2,'$\mathrm{tr}\,X(t)$','Interpreter','latex');
title(ax2,'closed loop at $\bar T^\star$ (decay)','Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig2, fullfile(figdir,'stabz_sos_decay.pdf'));

% native combined figure: the two panels above STACKED vertically (same
% content), full text width with readable fonts (no half-width tiling)
fig3 = figure('Color','w','Units','centimeters','Position',[2 2 14 12]);
tl = tiledlayout(fig3,2,1,'TileSpacing','compact','Padding','compact');
ax3 = nexttile(tl); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
plot(ax3, Tgrid, rad, '-', 'LineWidth',2.0,'Color',gc.order(1,:));
yline(ax3,1,'--','Color',gc.order(2,:),'LineWidth',1.5,'HandleVisibility','off');
xline(ax3,Tstar,':','Color',gc.onyx,'LineWidth',1.5,'HandleVisibility','off');
plot(ax3, Tstar, radAt, 'o','MarkerSize',7,'MarkerFaceColor',gc.order(3,:), ...
     'MarkerEdgeColor','k','HandleVisibility','off');
xlabel(ax3,'dwell-time $\bar T$','Interpreter','latex');
ylabel(ax3,'monodromy spectral radius','Interpreter','latex');
title(ax3,sprintf('$\\bar T^\\star=%.3f$',Tstar),'Interpreter','latex','FontWeight','normal');
ax4 = nexttile(tl); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
plot(ax4, tt, trX, '-', 'LineWidth',1.8,'Color',gc.order(1,:));
set(ax4,'YScale','log');
ylim(ax4, 10.^ly);
xlabel(ax4,'time $t$','Interpreter','latex');
ylabel(ax4,'$\mathrm{tr}\,X(t)$','Interpreter','latex');
title(ax4,'closed loop at $\bar T^\star$ (decay)','Interpreter','latex','FontWeight','normal');
set([ax3 ax4],'FontSize',10);
mvst_savepdf(fig3, fullfile(figdir,'stabz_sos_combined.pdf'));

% ---------------- data macros ----------------
w = mvst_texdata(fullfile(figdir,'data','stabz_sos.tex'));
w.num('StabzSosTbar', Tstar);
w.num('StabzSosSigma', sigma);
w.num('StabzSosBeta', beta);
w.num('StabzSosJ', jfac);
w.num('StabzSosFlowAbscissa', absFlow);
w.num('StabzSosRadAt', radAt);
w.num('StabzSosRadBelow', radBelow);
w.num('StabzSosRadAbove', radAbove);
w.num('StabzSosGainK', Kstar(1));
w.num('StabzSosGainKtwo', Kstar(2));
w.close();
fprintf('   wrote stabz_sos.pdf and data/stabz_sos.tex\n');
end

% =========================================================================
function out = synth_slack(A0,B0,A1,B1,J,Tbar,alpha,rho)
%SYNTH_SLACK  Timer-dependent (alpha,rho)-synthesis as a PURE SOSTOOLS matrix
%  sum-of-squares in the normalized timer theta = tau/Tbar in [0,1].  Maximizes
%  the scalar slack s; feasible iff s > 0.  Returns .slack, the solution
%  polynomials .Y, .U (in the pvar .th), and .Tbar for gain recovery.
n = size(A0,1);
pvar th; sl = dpvar('sl');
prog = sosprogram(th);
prog = sosdecvar(prog, sl);
[prog,Y]    = sospolymatrixvar(prog, monomials(th,0:2), [n n], 'symmetric');
[prog,U]    = sospolymatrixvar(prog, monomials(th,0:2), [1 n]);
[prog,Sig1] = sospolymatrixvar(prog, monomials(th,0),   [2*n 2*n], 'symmetric');
[prog,Xi1]  = sospolymatrixvar(prog, monomials(th,0),   [n n], 'symmetric');

Yp = diff(Y,th);                       % dY/dtheta ; real-time Ydot = Yp/Tbar
AY = A0*Y + B0*U;   V = A1*Y + B1*U;
% covariance-level rate: flow (1,1) block carries +alpha Y, jump block scaled rho
Phif = [ -(1/Tbar)*Yp + (AY+AY') + alpha*Y , V.' ;
          V                                 , -Y  ];
gt = th*(1-th);

% flow block matrix-SOS on [0,1] via the interval multiplier theta(1-theta)
prog = sosmatrixineq(prog, -Phif - sl*eye(2*n) - gt*Sig1);
prog = sosmatrixineq(prog, Sig1);
% scale normalization Y(theta) <= I on [0,1] (keeps the slack scale-fixed)
prog = sosmatrixineq(prog, eye(n) - Y - gt*Xi1);
prog = sosmatrixineq(prog, Xi1);
% jump block: single constant LMI at theta = 1
Y1 = subs(Y,th,1);  Y0 = subs(Y,th,0);
Phij = [ -rho*Y1 , (J*Y1).' ; J*Y1 , -Y0 ];
prog = sosmatrixineq(prog, -Phij - sl*eye(2*n));

prog = sossetobj(prog, -sl);
o.solver = 'sedumi';  o.params.fid = 0;
prog = sossolve(prog, o);
out.slack = double(sosgetsol(prog, sl));
out.Y  = sosgetsol(prog, Y);
out.U  = sosgetsol(prog, U);
out.th = th;
out.Tbar = Tbar;
end

% =========================================================================
function s = feas_constant(A0,B0,A1,B1,J)
%FEAS_CONSTANT  Arbitrary-dwell / constant-certificate synthesis: constant Y,
%  U0 (no timer) with the jump block J Y J' <= Y.  Expansive J makes this
%  infeasible; returns the optimal slack (<= 0).  Solved with SOSTOOLS
%  (degree-0 matrices; sosmatrixineq on a constant matrix is the LMI M >= 0).
n = size(A0,1);
pvar th; sl = dpvar('sl');
prog = sosprogram(th);
prog = sosdecvar(prog, sl);
[prog,Y] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
[prog,U] = sospolymatrixvar(prog, monomials(th,0), [1 n]);
AY = A0*Y + B0*U;   V = A1*Y + B1*U;
Phif = [ (AY+AY') , V.' ; V , -Y ];
prog = sosmatrixineq(prog, -Phif - sl*eye(2*n));
prog = sosmatrixineq(prog, [ -Y , (J*Y).' ; J*Y , -Y ] - sl*eye(2*n));
prog = sosmatrixineq(prog, eye(n) - Y);        % Y <= I normalization
prog = sossetobj(prog, -sl);
o.solver = 'sedumi';  o.params.fid = 0;
prog = sossolve(prog, o);
s = double(sosgetsol(prog, sl));
end

% =========================================================================
function [K, absc] = flow_gain(A0,B0,A1,B1)
%FLOW_GAIN  Flow gain minimizing the flow-alone second-moment abscissa on the
%  lifted symmetric space, max real eig( liftsym(A0+B0K,{A1+B1 K}) ), searched
%  over the pole-placement family (both flow poles at -c).  With input noise the
%  abscissa is bounded below (a finite minimum), hence a finite minimum dwell.
best = inf; K = [];
for c = 0.1:0.05:5
    Kc = -acker(A0,B0,[-c -c-1e-3]);
    a = max(real(eig(mvst_liftsym(A0+B0*Kc, {A1+B1*Kc}))));
    if a < best, best = a; K = Kc; end
end
absc = best;
end

% =========================================================================
function sol = build_gain_profile(sol)
%BUILD_GAIN_PROFILE  Precompute the recovered flow gain K0(theta)=U0 Y^{-1} on a
%  dense theta grid (verification side only), for fast lookup in the monodromy.
thg = linspace(0,1,401); Kg = zeros(numel(thg),2);
for i = 1:numel(thg)
    Yv = double(subs(sol.Y, sol.th, thg(i)));
    Uv = double(subs(sol.U, sol.th, thg(i)));
    Kg(i,:) = Uv / Yv;
end
sol.thg = thg;  sol.Kg = Kg;
end

function K = gain_at(sol, theta)
%GAIN_AT  Recovered flow gain at normalized timer theta in [0,1].
th = min(max(theta,0),1);
K = interp1(sol.thg, sol.Kg, th);
end

% =========================================================================
function rho = mono_radius(A0,B0,A1,B1,J,sol,Teval)
%MONO_RADIUS  Spectral radius of the one-period monodromy on the lifted
%  symmetric space: flow of F_{K0(theta)} over [0,Teval] (theta = t/Tstar,
%  frozen at theta = 1 for t > Tstar) then the jump J(.)J'.
Tstar = sol.Tbar;
m = size(mvst_liftsym(A0,{A1}),1);
Psi = eye(m); nstep = 240; ds = Teval/nstep;
for i = 1:nstep
    tmid = (i-0.5)*ds;
    K = gain_at(sol, tmid/Tstar);
    M = mvst_liftsym(A0+B0*K, {A1+B1*K});
    Psi = expm(M*ds)*Psi;
end
MJ = mvst_liftsym(zeros(size(A0)), {J});
rho = max(abs(eig(MJ*Psi)));
end

% =========================================================================
function [tt, trX] = cl_trajectory(A0,B0,A1,B1,J,sol,Ncyc)
%CL_TRAJECTORY  Trace of the closed-loop second moment over Ncyc dwell-time
%  periods at Tbar* (exact lifted flow with the synthesized timer gain, then
%  jump), started from X0 = I.
n = size(A0,1);  Tstar = sol.Tbar;
U = mvst_svecbasis(n);
MJ = mvst_liftsym(zeros(n), {J});
sx = U'*reshape(eye(n),[],1);
tt = []; trX = []; t0 = 0; npt = 40;
for c = 1:Ncyc
    ts = linspace(0,Tstar,npt); sxc = sx;
    for j = 1:npt
        if j > 1
            ds = ts(j)-ts(j-1); tmid = 0.5*(ts(j)+ts(j-1));
            K = gain_at(sol, tmid/Tstar);
            sxc = expm(mvst_liftsym(A0+B0*K,{A1+B1*K})*ds)*sxc;
        end
        X = reshape(U*sxc, n, n);
        tt(end+1) = t0+ts(j); trX(end+1) = trace(X); %#ok<AGROW>
    end
    sx = MJ*sxc; t0 = t0 + Tstar;
end
end
