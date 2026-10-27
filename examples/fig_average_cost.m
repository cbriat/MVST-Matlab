%% FIG_AVERAGE_COST  Long-term average cost under the optimal feedback.
%  Computes P_inf and the optimal gain, then simulates the closed-loop
%  covariance under additive noise W and shows the running average cost
%  rate converging to rho* = <P_inf,W>.  Saves figures/oc_average_cost.pdf.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D  = ex_data();  gc = mvst_gemcolors();

ra = mvst_solve_average(D.ct);
K  = ra.K{1};  rho = ra.rho;
fprintf('Optimal gain K = [% .4f % .4f]\n', K);
fprintf('Average-cost rate rho* = <P_inf,W> = %.4f\n', rho);

sim = mvst_simulate_ct(D.ct.modes{1}, K, [0 40]);

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, sim.t, sim.runavg, '-', 'LineWidth',2.2, 'Color',gc.order(1,:));
yline(ax, rho, '--', 'LineWidth',1.6, 'Color',gc.order(2,:));
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'running average cost rate','Interpreter','latex');
lg = legend(ax, {'$\frac{1}{t}\int_0^t \mathrm{tr}((Q+K^\top RK)\Sigma)\,ds$', ...
      sprintf('$\\rho^\\star=%.3f$', rho)}, ...
      'Interpreter','latex','Location','northeast');
set(lg,'Box','off'); set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'oc_average_cost.pdf'));
