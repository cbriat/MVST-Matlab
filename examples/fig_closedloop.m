%% FIG_CLOSEDLOOP  Closed-loop sample paths and control input.
%  Under the infinite-horizon optimal feedback, simulates several stochastic
%  sample paths of the representative diffusion and shows (top) the state
%  trajectories with their mean and (bottom) the optimal control input
%  u = K x over time.  Saves figures/oc_closedloop.pdf.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();

ri = mvst_solve_infinite(D.ct);  K = ri.K{1};
md = D.ct.modes{1};
x0 = [1; 0];  T = 12;  dt = 1e-3;  nreal = 40;
sp = mvst_sample_ct(md, K, x0, T, dt, nreal, 7);
% Downsample for plotting (dense paths bloat vector PDFs).
ds = max(1, round(numel(sp.t)/1500)); ii = 1:ds:numel(sp.t);
sp.t = sp.t(ii); sp.x = sp.x(:,ii,:); sp.u = sp.u(:,ii,:);
sp.xmean = sp.xmean(:,ii); sp.umean = sp.umean(:,ii);

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 12]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
for r = 1:nreal
    plot(ax1, sp.t, squeeze(sp.x(1,:,r)), '-', 'LineWidth',0.5, 'Color',[gc.order(1,:) 0.20]);
    plot(ax1, sp.t, squeeze(sp.x(2,:,r)), '-', 'LineWidth',0.5, 'Color',[gc.order(2,:) 0.20]);
end
h1 = plot(ax1, sp.t, sp.xmean(1,:), '-', 'LineWidth',2.4, 'Color',gc.order(1,:));
h2 = plot(ax1, sp.t, sp.xmean(2,:), '-', 'LineWidth',2.4, 'Color',gc.order(2,:));
ylabel(ax1, 'states $x_1(t),\ x_2(t)$', 'Interpreter','latex');
lg1 = legend(ax1, [h1 h2], {'$x_1$ mean','$x_2$ mean'}, 'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg1,'Box','off'); set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for r = 1:nreal
    plot(ax2, sp.t, squeeze(sp.u(1,:,r)), '-', 'LineWidth',0.5, 'Color',[gc.order(1,:) 0.25]);
end
humean = plot(ax2, sp.t, sp.umean(1,:), '-', 'LineWidth',2.4, 'Color',gc.order(1,:));
xlabel(ax2, 'time $t$', 'Interpreter','latex');
ylabel(ax2, 'control $u(t)=Kx(t)$', 'Interpreter','latex');
lg2 = legend(ax2, humean, {'ensemble mean'}, 'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg2,'Box','off'); set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'oc_closedloop.pdf'));
