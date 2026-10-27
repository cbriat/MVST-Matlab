%% EX01_MJLS_GAINS  Markov jump system: P(t), gains, and the three regimes.
%  Solves the finite-horizon, infinite-horizon, and long-term average
%  problems for the two-mode Markov-jump-with-diffusion system.  The
%  finite-horizon gains K_i(t) are time-varying and are plotted together
%  with the (constant) infinite-horizon gains; P_i(t) is plotted likewise.
%  All constant data (systems, optimal gains, rho*) are exported to
%  figures/data/ex_mjls.tex for inclusion in the document.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
sys = D.mjls; T = 5.0;

rf = mvst_solve_finite(sys, T);     % coupled finite horizon: P_i(t), K_i(t)
ri = mvst_solve_infinite(sys);      % infinite horizon: constant P_i, K_i
ra = mvst_solve_average(sys);       % long-term average rate rho*

fprintf('K1_inf = [% .4f % .4f],  K2_inf = [% .4f % .4f]\n', ri.K{1}, ri.K{2});
fprintf('rho* = %.4f,  finite-horizon cost = %.4f\n', ra.rho, rf.cost);

% ---- export data for LaTeX (\input) ----
w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_mjls.tex'));
w.mat('AmjlsOne', sys.modes{1}.A); w.mat('AmjlsTwo', sys.modes{2}.A);
w.mat('BmjlsOne', sys.modes{1}.B); w.mat('BmjlsTwo', sys.modes{2}.B);
w.mat('GenMjls',  sys.Pi);
w.mat('QmjlsOne', sys.modes{1}.Q); w.mat('QmjlsTwo', sys.modes{2}.Q);
w.num('RmjlsOne', sys.modes{1}.R); w.num('RmjlsTwo', sys.modes{2}.R);
w.mat('KmjlsInfOne', ri.K{1}); w.mat('KmjlsInfTwo', ri.K{2});
w.num('RhoMjls', ra.rho); w.num('CostMjls', rf.cost);
w.close();

% ---- figure: gains K_i(t) and P_i(t) over the finite horizon ----
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 11]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, rf.t, squeeze(rf.K{1}(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(1,:));
plot(ax1, rf.t, squeeze(rf.K{1}(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(2,:));
plot(ax1, rf.t, squeeze(rf.K{2}(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(3,:));
plot(ax1, rf.t, squeeze(rf.K{2}(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(4,:));
yline(ax1, ri.K{1}(1),':','Color',gc.order(1,:),'LineWidth',1.0);
yline(ax1, ri.K{1}(2),':','Color',gc.order(2,:),'LineWidth',1.0);
yline(ax1, ri.K{2}(1),':','Color',gc.order(3,:),'LineWidth',1.0);
yline(ax1, ri.K{2}(2),':','Color',gc.order(4,:),'LineWidth',1.0);
ylabel(ax1,'gains $K_i(t)$','Interpreter','latex');
lg = legend(ax1, {'$[K_1]_1$','$[K_1]_2$','$[K_2]_1$','$[K_2]_2$'}, ...
      'Interpreter','latex','Location','northeastoutside'); set(lg,'Box','off');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, rf.t, squeeze(rf.P{1}(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(1,:));
plot(ax2, rf.t, squeeze(rf.P{1}(2,2,:)), '-','LineWidth',2.0,'Color',gc.order(2,:));
plot(ax2, rf.t, squeeze(rf.P{1}(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(3,:));
yline(ax2, ri.P{1}(1,1),':','Color',gc.order(1,:),'LineWidth',1.0);
yline(ax2, ri.P{1}(2,2),':','Color',gc.order(2,:),'LineWidth',1.0);
yline(ax2, ri.P{1}(1,2),':','Color',gc.order(3,:),'LineWidth',1.0);
xlabel(ax2,'time $t$','Interpreter','latex');
ylabel(ax2,'mode-1 dual $P_1(t)$','Interpreter','latex');
lg2 = legend(ax2, {'$[P_1]_{11}$','$[P_1]_{22}$','$[P_1]_{12}$'}, ...
      'Interpreter','latex','Location','northeastoutside'); set(lg2,'Box','off');
set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_mjls_gains.pdf'));
