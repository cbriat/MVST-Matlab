%% FIG_FINITE_HORIZON  Finite-horizon I-GRE solution and optimal cost.
%  Solves the finite-horizon problem for the single-mode CT system of
%  EX_DATA by integrating the differential I-GRE backward, and saves
%  figures/oc_finite_horizon.pdf.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D  = ex_data();  gc = mvst_gemcolors();

T   = 5.0;
res = mvst_solve_finite(D.ct, T);

P11 = squeeze(res.P(1,1,:));
P12 = squeeze(res.P(1,2,:));
P22 = squeeze(res.P(2,2,:));
fprintf('Finite-horizon optimal cost  J*(X0) = %.4f  (T=%.1f)\n', res.cost, T);

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, res.t, P11, '-' , 'LineWidth',2.2, 'Color',gc.order(1,:));
plot(ax, res.t, P22, '-' , 'LineWidth',2.2, 'Color',gc.order(2,:));
plot(ax, res.t, P12, '-', 'LineWidth',2.2, 'Color',gc.order(3,:));
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'entries of $P(t)$','Interpreter','latex');
lg = legend(ax, {'$P_{11}(t)$','$P_{22}(t)$','$P_{12}(t)$'}, ...
      'Interpreter','latex','Location','northeastoutside');
set(lg,'Box','off'); set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'oc_finite_horizon.pdf'));
