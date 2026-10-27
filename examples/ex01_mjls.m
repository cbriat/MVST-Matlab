%% EX01_MJLS  Optimal control of a Markov jump system with diffusion.
%  Solves the coupled infinite-horizon / average-cost I-GRE for the two-mode
%  Markov-jump-with-diffusion system, then shows stochastic sample paths of
%  the state and the optimal control.  The long-term average cost rate is
%  reported separately in ex_average (fig:ex:average), so it is not replotted
%  here.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
sys = D.mjls;

ra = mvst_solve_average(sys);
fprintf('Coupled I-GRE: %d iters, residual %.2e\n', ra.iters, ra.residual);
for i = 1:sys.M; fprintf('  K_%d = [% .4f % .4f]\n', i, ra.K{i}); end
fprintf('  rho* = %.4f\n', ra.rho);

% Stochastic sample paths (for state and control).
sp  = mvst_sample_mjls(sys, ra.K, [1;0], 12, 1e-3, 30, 11);
% Downsample for plotting (dense paths bloat vector PDFs).
dsi = max(1, round(numel(sp.t)/1500)); ii = 1:dsi:numel(sp.t);
sp.t = sp.t(ii); sp.x = sp.x(:,ii,:); sp.u = sp.u(:,ii,:);

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 10]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
for r = 1:size(sp.x,3)
    plot(ax1, sp.t, squeeze(sp.x(1,:,r)), '-','LineWidth',0.5,'Color',[gc.order(1,:) 0.20]);
    plot(ax1, sp.t, squeeze(sp.x(2,:,r)), '-','LineWidth',0.5,'Color',[gc.order(2,:) 0.20]);
end
h1 = plot(ax1, sp.t, mean(sp.x(1,:,:),3), '-','LineWidth',2.2,'Color',gc.order(1,:));
h2 = plot(ax1, sp.t, mean(sp.x(2,:,:),3), '-','LineWidth',2.2,'Color',gc.order(2,:));
ylabel(ax1,'states $x_1(t),\ x_2(t)$','Interpreter','latex');lg1=legend(ax1,[h1 h2],{'$x_1$ mean','$x_2$ mean'},'Interpreter','latex','Location','northeast'); set(lg1,'Box','off');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for r = 1:size(sp.u,3)
    plot(ax2, sp.t, squeeze(sp.u(1,:,r)), '-','LineWidth',0.5,'Color',[gc.order(1,:) 0.25]);
end
hu = plot(ax2, sp.t, mean(sp.u(1,:,:),3), '-','LineWidth',2.2,'Color',gc.order(1,:));
xlabel(ax2,'time $t$','Interpreter','latex');
ylabel(ax2,'control $u(t)=K_{\sigma}x$', 'Interpreter', 'latex');
lg2=legend(ax2,hu,{'ensemble mean'},'Location','northeast'); set(lg2,'Box','off');
set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_mjls.pdf'));
