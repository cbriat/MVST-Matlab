function ex_dof_sos()
%EX_DOF_SOS  Timer-dependent output-feedback bracket SYNTHESIS across the
%  dwell-time regimes, solved by SOSTOOLS, with controller reconstruction.
%
%  Chapter "Dynamic Output Feedback", Theorem th:dof:dwell and
%  Proposition prop:dof:bracket. Thin caller of the generic engine
%  lib/mvst_dof_sos.m; see that function's header for the bracket, the three
%  dwell-time regimes, and the reconstruction. Solver: SOSTOOLS + SeDuMi.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();

% ---- data (scalar plant, one multiplicative flow channel, controlled reset) ----
sys.a0 = 1.0; sys.a1 = 0.4;   % A0(tau) = a0 + a1 tau (unstable, timer-dependent)
sys.sig = 0.5;                % multiplicative-noise intensity
sys.mu  = 0.5;               % fixed scalar bracket parameter
sys.B0 = 1; sys.Cy = 1;      % flow control in / measurement out
sys.J0 = 0.5; sys.BJ = 1; sys.CJ = 1;   % reset map / jump control / jump measurement
opts.alpha = 1e-3; opts.eps = 1e-3;

T = 1.0; Tmin = 0.5; Tmax = 1.5; deg = 2;
dw.periodic = struct('type','periodic','T',T,'deg',deg);
dw.minimum  = struct('type','minimum', 'T',T,'deg',deg);
dw.range    = struct('type','range','T',T,'Tmin',Tmin,'Tmax',Tmax,'deg',deg);

regimes = {'periodic','minimum','range'};
sol = struct();
for i = 1:numel(regimes)
    r = regimes{i};
    out = mvst_dof_sos(sys, dw.(r), opts);
    sol.(r) = out;
    fprintf('%-9s dwell: solved=%d (pinf=%d dinf=%d numerr=%d) | max lam(M)=%+.3e  min lam(Om)=%+.3e\n', ...
        r, out.feas, out.info.pinf, out.info.dinf, out.info.numerr, out.maxM, out.minOm);
end

% ---- reconstructed controller K(tau) for the periodic solve ----
K = sol.periodic.K;
taus = [0 0.5 1.0];
fprintf('\nReconstructed periodic flow controller K(tau) = (A_c,B_c,C_c,D_c):\n');
fprintf('  D_c = %+.4f (constant)\n', K.Dc);
samplerep('  C_c(tau):', K.Cc, taus);
samplerep('  B_c(tau):', K.Bc, taus);
samplerep('  A_c(tau):', K.Ac, taus);
fprintf('  jump gains: D_cJ=%+.4f  C_cJ=%+.4f  B_cJ=%+.4f  A_cJ=%+.4f\n', ...
        K.DcJ, K.CcJ, K.BcJ, K.AcJ);

% ================== figures ==================
% (1) certificate of the periodic solve over [0,T]
tt = linspace(0, T, 200);
lamM  = arrayfun(@(t) max(real(eig(sol.periodic.Mf(t)))),  tt);
lamOm = arrayfun(@(t) min(real(eig(sol.periodic.Omf(t)))), tt);
fig1 = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig1); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h1 = plot(ax, tt, lamM,  '-','LineWidth',2.0,'Color',gc.order(1,:));
h2 = plot(ax, tt, lamOm, '-','LineWidth',2.0,'Color',gc.order(2,:));
yline(ax,0,'--','Color','k','LineWidth',1.2);
xlabel(ax,'timer $\tau$','Interpreter','latex');
ylabel(ax,'eigenvalue bound','Interpreter','latex');
lg = legend(ax,[h1 h2], ...
   {'$\lambda_{\max}$ of the flow bracket $M(\tau)$ ($<0$)', ...
    '$\lambda_{\min}$ of the coupling $\Omega(\tau)$ ($>0$)'}, ...
   'Interpreter','latex','Location','east');
set(lg,'Box','off'); set(ax,'FontSize',10);
title(ax,'Periodic dwell-time DOF bracket solved by SOS on $[0,T]$', ...
   'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig1, fullfile(mvst_figdir(),'dof_sos.pdf'));

% (2) reconstructed timer-dependent controller gains
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig2); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
p1 = plot(ax, tt, arrayfun(K.Cc, tt), '-','LineWidth',2.0,'Color',gc.order(1,:));
p2 = plot(ax, tt, arrayfun(K.Bc, tt), '-','LineWidth',2.0,'Color',gc.order(2,:));
p3 = plot(ax, tt, arrayfun(K.Ac, tt), '-','LineWidth',2.0,'Color',gc.order(3,:));
yline(ax, K.Dc, ':','Color',[.5 .5 .5],'LineWidth',1.5);
xlabel(ax,'timer $\tau$','Interpreter','latex');
ylabel(ax,'controller gain','Interpreter','latex');
lg = legend(ax,[p1 p2 p3], {'$C_c(\tau)$','$B_c(\tau)$','$A_c(\tau)$'}, ...
   'Interpreter','latex','Location','best');
set(lg,'Box','off'); set(ax,'FontSize',10);
title(ax,'Reconstructed timer-dependent controller $K(\tau)$', ...
   'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig2, fullfile(mvst_figdir(),'dof_sos_gains.pdf'));
end

% =========================================================================
function samplerep(label, fun, taus)
%SAMPLEREP  Print a (rational) gain sampled at a few timer values.
fprintf('%s ', label);
for t = taus, fprintf(' tau=%.2f: %+.4f ', t, fun(t)); end
fprintf('\n');
end
