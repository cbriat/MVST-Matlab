function Tstar = ex_stabz_switched_sos()
%EX_STABZ_SWITCHED_SOS  Minimum-dwell-time SWITCHED STABILIZATION by SOSTOOLS
%  matrix sum-of-squares in the timer (Stabilization chapter, mode-indexed
%  counterpart of Theorem th:stab:stabzMDT).
%
%  Two-mode controlled switched system, mode-dependent state feedback:
%     flow (mode i):  dx = (A_i x + B u) dt + (sigma x + B1 u) dW1,   u = K_i x,
%     switch:         x+ = Jsw x   (identity-like reset).
%  Second-moment closed-loop flow operator per mode with gain K:
%     F_i(X) = He((A_i+B K) X) + (sigma I + B1 K) X (sigma I + B1 K)'.
%
%  TUNING / FINDING.  With PURE identity switches and fully controllable modes,
%  the mode-dependent feedback stabilizes under ARBITRARILY fast switching (a
%  common control-Lyapunov certificate exists), so there is no positive minimum
%  dwell-time (Tbar* -> 0).  A well-posed minimum dwell-time requires (a) the
%  achievable per-mode second-moment decay to be BOUNDED and (b) the switch to
%  cost something.  As in ex_stabz_sos, a small INPUT-multiplicative noise
%  B1 = beta*[0;1] bounds the decay (large gains inject noise), and a small
%  EXPANSIVE switch reset Jsw = jfac*I (a switching transient) makes fast
%  switching lossy.  Both channels stay inside the theorem (the noise is the
%  V_l = A_l Y + B_l U_l term; the reset is the switch map).  Reported data:
%  A1, A2 (both unstable), B, sigma, beta, jfac.
%
%  On the timer, normalized to theta = tau/Tbar in [0,1], each mode carries
%  Y_i(theta) in S^2 (symmetric, degree 2) and U_i(theta) (degree 2).  The
%  minimum-dwell-time SYNTHESIS blocks are
%     flow (mode i):  Phi_f^i = [ -(1/Tbar) Y_i' + He(A_i Y_i + B U_i) , (sigma Y_i + B1 U_i)' ;
%                                  sigma Y_i + B1 U_i                   , -Y_i ] <= -s I,
%     switch (i->j):  [ -Y_i(1) , (Jsw Y_i(1))' ; Jsw Y_i(1) , -Y_j(0) ] <= -s I,
%  with recovered flow gains K_i(theta) = U_i(theta) Y_i(theta)^{-1}.
%
%  METHOD.  PURE SOSTOOLS MATRIX SUM-OF-SQUARES IN THE TIMER (no gridding of the
%  timer): each flow block is a matrix-SOS Positivstellensatz on [0,1] with the
%  interval multiplier theta(1-theta); the bounds 1e-3 I <= Y_i <= I are imposed
%  the same way; the switch blocks are constant LMIs.  Slack-maximized; minimum
%  dwell-time Tbar* by bisection on the slack sign.  Tbar* is a SUFFICIENT
%  certified minimum dwell-time (the SOS relaxation for switched SYNTHESIS is
%  looser than the tight analysis cases, so Tbar* over-estimates the exact
%  boundary).  VERIFIED by the exact switched one-period second-moment monodromy
%  on the lifted symmetric space (flow mode 1 with K1(theta) over [0,Tbar],
%  switch Jsw, flow mode 2 with K2(theta), switch Jsw): its spectral radius at
%  Tbar* is < 1, i.e. the recovered gains stabilize, and the closed-loop trace
%  decays over the switch cycles.  Solver: SOSTOOLS + SeDuMi.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---------------- system data ----------------
A1 = [0 1; 2 0.1];  A2 = [0 1; -3 0.2];  B = [0;1];
sigma = 0.1;  beta = 0.55;  B1 = beta*[0;1];
jfac = 1.1;   Jsw = jfac*eye(2);
tol = 1e-5;
fprintf('\n==== minimum-dwell-time SWITCHED stabilization (SOSTOOLS matrix-SOS in tau) ====\n');
fprintf('A1 eig %s  A2 eig %s  B=[%g;%g]  sigma=%g  beta=%g  Jsw=%g*I\n', ...
    mat2str(eig(A1)',3), mat2str(eig(A2)',3), B(1),B(2), sigma, beta, jfac);

% ---------------- minimum dwell-time by bisection ----------------
lo = 0.4; hi = 1.6;
assert(synth(A1,A2,B,sigma,B1,Jsw,lo).slack <= tol, 'lower bracket already feasible');
assert(synth(A1,A2,B,sigma,B1,Jsw,hi).slack >  tol, 'upper bracket infeasible');
for it = 1:22
    mid = 0.5*(lo+hi);
    if synth(A1,A2,B,sigma,B1,Jsw,mid).slack > tol, hi = mid; else, lo = mid; end
end
Tstar = hi;
sol = synth(A1,A2,B,sigma,B1,Jsw,Tstar);
sol = build_profiles(sol);
fprintf('certified (sufficient) minimum dwell-time Tbar* = %.4f  (slack there = %.2e)\n', Tstar, sol.slack);

% ---------------- verification: synthesized closed loop is stable at Tbar* ---
% The minimum-dwell-time synthesis gives a SUFFICIENT certified dwell-time; the
% exact switched second-moment monodromy with the recovered gains confirms
% stability at Tbar* (spectral radius < 1).  The SOS relaxation for switched
% SYNTHESIS is looser than the tight analysis cases, so Tbar* is a conservative
% sufficient bound, not the exact stability boundary.
radAt = mono_radius(A1,A2,B,sigma,B1,Jsw,sol,Tstar);
K1s = gain_at(sol,1,1.0);  K2s = gain_at(sol,2,1.0);
fprintf('  switched monodromy radius at Tbar* = %.4f  (<1: recovered gains stabilize)\n', radAt);
fprintf('  representative gains K1(Tbar*)=[%.3f %.3f]  K2(Tbar*)=[%.3f %.3f]\n', K1s, K2s);

% ---------------- closed-loop second-moment trajectory at Tbar* -------------
[tt, trX] = cl_trajectory(A1,A2,B,sigma,B1,Jsw,sol,5);

% ---------------- figure ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, tt, trX, '-','LineWidth',1.8,'Color',gc.order(1,:));
set(ax,'YScale','log');
for c = 1:5, xline(ax, c*Tstar, ':','Color',[0.6 0.6 0.6],'HandleVisibility','off'); end
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'$\mathrm{tr}\,X(t)$','Interpreter','latex');
title(ax,sprintf('Switched closed loop at $\\bar T^\\star=%.3f$ (decay)',Tstar), ...
    'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(figdir,'sw_stabz_sos.pdf'));

% ---------------- data macros ----------------
w = mvst_texdata(fullfile(figdir,'data','sw_stabz_sos.tex'));
w.num('SwStabzTbar', Tstar);
w.num('SwStabzSigma', sigma);
w.num('SwStabzBeta', beta);
w.num('SwStabzJ', jfac);
w.num('SwStabzKoneA', K1s(1));  w.num('SwStabzKoneB', K1s(2));
w.num('SwStabzKtwoA', K2s(1));  w.num('SwStabzKtwoB', K2s(2));
w.num('SwStabzRadAt', radAt);
w.close();
fprintf('  wrote sw_stabz_sos.pdf and data/sw_stabz_sos.tex\n');
end

% =========================================================================
function out = synth(A1,A2,B,sigma,B1,Jsw,Tbar)
%SYNTH  Two-mode switched minimum-dwell-time synthesis as a PURE SOSTOOLS matrix
%  sum-of-squares in the normalized timer.  Slack-maximized; feasible iff s > 0.
n = size(A1,1);
pvar th; sl = dpvar('sl');
prog = sosprogram(th);  prog = sosdecvar(prog, sl);
[prog,Y1] = sospolymatrixvar(prog, monomials(th,0:2), [n n], 'symmetric');
[prog,Y2] = sospolymatrixvar(prog, monomials(th,0:2), [n n], 'symmetric');
[prog,U1] = sospolymatrixvar(prog, monomials(th,0:2), [1 n]);
[prog,U2] = sospolymatrixvar(prog, monomials(th,0:2), [1 n]);
[prog,S1] = sospolymatrixvar(prog, monomials(th,0), [2*n 2*n], 'symmetric');
[prog,S2] = sospolymatrixvar(prog, monomials(th,0), [2*n 2*n], 'symmetric');
[prog,L1] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
[prog,L2] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
gt = th*(1-th);
prog = flowblk(prog, A1,B,sigma,B1, Y1,U1,S1,L1, sl,gt,th,Tbar,n);
prog = flowblk(prog, A2,B,sigma,B1, Y2,U2,S2,L2, sl,gt,th,Tbar,n);
Y1_0 = subs(Y1,th,0); Y1_1 = subs(Y1,th,1);
Y2_0 = subs(Y2,th,0); Y2_1 = subs(Y2,th,1);
% switch blocks (expansive reset Jsw), mode 1 -> 2 and 2 -> 1
prog = sosmatrixineq(prog, -[ -Y1_1, (Jsw*Y1_1).'; Jsw*Y1_1, -Y2_0 ] - sl*eye(2*n));
prog = sosmatrixineq(prog, -[ -Y2_1, (Jsw*Y2_1).'; Jsw*Y2_1, -Y1_0 ] - sl*eye(2*n));
prog = sossetobj(prog, -sl);
o.solver = 'sedumi';  o.params.fid = 0;
prog = sossolve(prog, o);
out.slack = double(sosgetsol(prog, sl));
out.Y = {sosgetsol(prog,Y1), sosgetsol(prog,Y2)};
out.U = {sosgetsol(prog,U1), sosgetsol(prog,U2)};
out.th = th;  out.Tbar = Tbar;
end

function prog = flowblk(prog, A,B,sigma,B1, Y,U,S,L, sl,gt,th,Tbar,n)
AY = A*Y + B*U;  V = sigma*Y + B1*U;
Phi = [ -(1/Tbar)*diff(Y,th) + (AY+AY.') , V.' ;
         V                               , -Y  ];
prog = sosmatrixineq(prog, -Phi - sl*eye(2*n) - gt*S);   prog = sosmatrixineq(prog, S);
prog = sosmatrixineq(prog, Y - 1e-3*eye(n) - gt*L);      prog = sosmatrixineq(prog, L);
prog = sosmatrixineq(prog, eye(n) - Y);
end

% =========================================================================
function sol = build_profiles(sol)
thg = linspace(0,1,401); Kg = {zeros(numel(thg),2), zeros(numel(thg),2)};
for m = 1:2
    for i = 1:numel(thg)
        Yv = double(subs(sol.Y{m}, sol.th, thg(i)));
        Uv = double(subs(sol.U{m}, sol.th, thg(i)));
        Kg{m}(i,:) = Uv / Yv;
    end
end
sol.thg = thg;  sol.Kg = Kg;
end

function K = gain_at(sol, m, theta)
K = interp1(sol.thg, sol.Kg{m}, min(max(theta,0),1));
end

% =========================================================================
function r = mono_radius(A1,A2,B,sigma,B1,Jsw,sol,Teval)
%MONO_RADIUS  Switched one-period second-moment monodromy spectral radius on the
%  lifted symmetric space: flow mode 1 (gain K1) over [0,Teval], switch Jsw,
%  flow mode 2 (gain K2) over [0,Teval], switch Jsw.
Am = {A1, A2};
MJ = mvst_liftsym(zeros(size(A1)), {Jsw});
Psi = eye(size(MJ,1));
for m = 1:2
    Psi = MJ * flow_transition(Am{m}, B, sigma, B1, sol, m, Teval) * Psi;
end
r = max(abs(eig(Psi)));
end

function Phi = flow_transition(A,B,sigma,B1,sol,m,Teval)
Tstar = sol.Tbar;  nstep = 200; ds = Teval/nstep;
Phi = eye(size(mvst_liftsym(A,{sigma*eye(size(A,1))}),1));
for i = 1:nstep
    tmid = (i-0.5)*ds;
    K = gain_at(sol, m, tmid/Tstar);
    M = mvst_liftsym(A+B*K, {sigma*eye(size(A,1)) + B1*K});
    Phi = expm(M*ds)*Phi;
end
end

% =========================================================================
function [tt, trX] = cl_trajectory(A1,A2,B,sigma,B1,Jsw,sol,Ncyc)
%CL_TRAJECTORY  Trace of the closed-loop second moment over Ncyc switch cycles
%  at Tbar* (exact lifted flow with the synthesized timer gains, then the
%  expansive switch), started from X0 = I.  One cycle = mode 1 then mode 2.
n = size(A1,1);  Tstar = sol.Tbar;  Am = {A1, A2};
U = mvst_svecbasis(n);
MJ = mvst_liftsym(zeros(n), {Jsw});
sx = U'*reshape(eye(n),[],1);
tt = []; trX = []; t0 = 0; npt = 30;
for c = 1:Ncyc
    for m = 1:2
        ts = linspace(0,Tstar,npt); sxc = sx;
        for j = 1:npt
            if j > 1
                ds = ts(j)-ts(j-1); tmid = 0.5*(ts(j)+ts(j-1));
                K = gain_at(sol, m, tmid/Tstar);
                sxc = expm(mvst_liftsym(Am{m}+B*K,{sigma*eye(n)+B1*K})*ds)*sxc;
            end
            X = reshape(U*sxc, n, n);
            tt(end+1) = t0+ts(j); trX(end+1) = trace(X); %#ok<AGROW>
        end
        sx = MJ*sxc; t0 = t0 + Tstar;
    end
end
end
