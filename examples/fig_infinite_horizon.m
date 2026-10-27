%% FIG_INFINITE_HORIZON  Monotone convergence P_T -> P_inf.
%  For a grid of horizons T, solves the finite-horizon I-GRE and records
%  P_T = P(0); overlays the infinite-horizon solution P_inf.  Saves
%  figures/oc_infinite_horizon.pdf.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D  = ex_data();  gc = mvst_gemcolors();

Tgrid = linspace(0.25, 10, 40);
P11 = zeros(size(Tgrid)); P12 = P11; P22 = P11;
for k = 1:numel(Tgrid)
    r = mvst_solve_finite(D.ct, Tgrid(k));
    P11(k) = r.P0(1,1); P12(k) = r.P0(1,2); P22(k) = r.P0(2,2);
end
ri = mvst_solve_infinite(D.ct);  Pinf = ri.P{1};
fprintf('P_inf =\n'); disp(Pinf);
fprintf('I-GRE residual = %.2e\n', ri.residual);

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, Tgrid, P11, '-o','LineWidth',2.0,'Color',gc.order(1,:),'MarkerSize',3,'MarkerFaceColor',gc.order(1,:));
plot(ax, Tgrid, P22, '-s','LineWidth',2.0,'Color',gc.order(2,:),    'MarkerSize',3,'MarkerFaceColor',gc.order(2,:));
plot(ax, Tgrid, P12, '-^','LineWidth',2.0,'Color',gc.order(3,:), 'MarkerSize',3,'MarkerFaceColor',gc.order(3,:));
yline(ax, Pinf(1,1),'--','Color',gc.order(1,:),'LineWidth',1.2);
yline(ax, Pinf(2,2),'--','Color',gc.order(2,:),    'LineWidth',1.2);
yline(ax, Pinf(1,2),'--','Color',gc.order(3,:), 'LineWidth',1.2);
xlabel(ax,'horizon $T$','Interpreter','latex');
ylabel(ax,'entries of $P_T=P(0)$','Interpreter','latex');
lg = legend(ax, {'$[P_T]_{11}$','$[P_T]_{22}$','$[P_T]_{12}$', ...
      '$[P_\infty]_{11}$','$[P_\infty]_{22}$','$[P_\infty]_{12}$'}, ...
      'Interpreter','latex','Location','northeastoutside');
set(lg,'Box','off'); set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'oc_infinite_horizon.pdf'));
