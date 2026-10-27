%% EX04_POISSON_GAINS  Finite-horizon gains and cost-vs-rate for the Poisson example.
%  Top: the finite-horizon gains K(t) of the generalized differential Riccati
%  equation settling onto the constant infinite-horizon gain.  Bottom: the
%  optimal average cost rate rho* as a function of the jump rate lambda,
%  increasing from the jump-free (pure-diffusion) value.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

A  = [0 1; 2 -0.5];  B  = [0;1];
C1 = 0.2*eye(2);     D1 = [0;0];
J1 = 1.25*eye(2);    lam = 0.5;
Q  = eye(2);         R  = 1;
W  = 0.05*eye(2);

% ----- finite-horizon gains vs the infinite-horizon gain -----
[Pinf, Kinf] = mvst_care_ct_poisson(A, B, {C1}, {D1}, {J1}, lam, Q, R); %#ok<ASGLU>
T = 8;
[tg, ~, Kt] = mvst_dre_ct_poisson(A, B, {C1}, {D1}, {J1}, lam, Q, R, T, zeros(2), 200);

% ----- average cost rate versus the jump rate -----
lams = linspace(0, 1.2, 25); rhos = zeros(size(lams));
for k = 1:numel(lams)
    Pk = mvst_care_ct_poisson(A, B, {C1}, {D1}, {J1}, lams(k), Q, R);
    rhos(k) = trace(Pk*W);
end

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 12]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
h = gobjects(2,1);
h(1) = plot(ax1, tg, squeeze(Kt(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(1,:));
h(2) = plot(ax1, tg, squeeze(Kt(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(2,:));
yline(ax1, Kinf(1), ':','LineWidth',1.4,'Color',gc.order(1,:));
yline(ax1, Kinf(2), ':','LineWidth',1.4,'Color',gc.order(2,:));
xlabel(ax1,'time $t$'); ylabel(ax1,'gains $K(t)$');
lg1 = legend(ax1, h, {'$K_1(t)$','$K_2(t)$'},'Interpreter','latex', ...
             'Location','southwest'); set(lg1,'Box','off');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, lams, rhos, '-o','LineWidth',2.0,'MarkerSize',3,'Color',gc.order(1,:));
xlabel(ax2,'jump rate $\lambda_1$'); ylabel(ax2,'average cost rate $\rho^\star$');
set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_poisson_gains.pdf'));
