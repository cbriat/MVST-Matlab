%% EX_TV_MJLS  Time-varying optimal control of a Markov jump diffusion.
%  Time-varying counterpart of EX01_MJLS: the drift of each mode decays in
%  time, A_i(t) = A_i + dA*exp(-t), so the system is genuinely time-varying
%  but settles to the LTI system of EX01_MJLS as t grows.  We solve the
%  coupled differential I-GRE on a finite horizon (giving a time-varying dual
%  P_i(t) and gain K_i(t)) and the infinite-horizon problem (giving the
%  time-varying dual P_{inf,i}(t) that settles to the LTI value).  This
%  illustrates the main contribution of the paper: the theory holds for
%  time-varying systems over both finite and infinite horizons.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data(); gc = mvst_gemcolors();
base = D.mjls; n = base.n; m = base.m; M = base.M;
dA = [0 0; 1.5 0];                         % decaying drift perturbation

% ----- build the time-varying modes A_i(t) = A_i + dA*exp(-t) -----
md = cell(1,M);
for i = 1:M
    b = base.modes{i};
    md{i} = mvst_mode_tv(n, m, ...
        'A', @(t) b.A + dA*exp(-t), 'B', b.B, 'Q', b.Q, 'R', b.R, ...
        'Ac', b.Ac, 'Bc', b.Bc, 'W', b.W, 'X0', b.X0);
end
sys = mvst_system(md, 'Pi', base.Pi, 'pi0', base.pi0);

T = 5.0; Tinf = 8.0;
rf = mvst_solve_finite_tv(sys, T);          % time-varying finite-horizon dual
ri = mvst_solve_infinite_tv(sys, Tinf);     % time-varying infinite-horizon dual
fprintf('TV-MJLS: asymptotic gains K1=[% .4f % .4f] K2=[% .4f % .4f]\n', ...
        ri.Klti{1}, ri.Klti{2});

% ----- export data macros -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_mjls.tex'));
w.mat('TVdAmjls', dA);
w.mat('TVKmjlsInfOne', ri.Klti{1}); w.mat('TVKmjlsInfTwo', ri.Klti{2});
w.num('TVTmjls', T); w.close();

% ----- figure: mode-1 dual P_1(t), finite (top) and infinite (bottom) -----
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 12]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, rf.t, squeeze(rf.P{1}(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(1,:));
plot(ax1, rf.t, squeeze(rf.P{1}(2,2,:)), '-','LineWidth',2.0,'Color',gc.order(2,:));
plot(ax1, rf.t, squeeze(rf.P{1}(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(3,:));
ylabel(ax1,'finite-horizon $P_1(t)$','Interpreter','latex');
lg = legend(ax1,{'$[P_1]_{11}$','$[P_1]_{22}$','$[P_1]_{12}$'}, ...
      'Interpreter','latex','Location','northeastoutside'); set(lg,'Box','off');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, ri.t, squeeze(ri.P{1}(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(1,:));
plot(ax2, ri.t, squeeze(ri.P{1}(2,2,:)), '-','LineWidth',2.0,'Color',gc.order(2,:));
plot(ax2, ri.t, squeeze(ri.P{1}(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(3,:));
yline(ax2, ri.Plti{1}(1,1), '--','Color',gc.order(1,:),'LineWidth',1.0);
yline(ax2, ri.Plti{1}(2,2), '--','Color',gc.order(2,:),'LineWidth',1.0);
yline(ax2, ri.Plti{1}(1,2), '--','Color',gc.order(3,:),'LineWidth',1.0);
xlabel(ax2,'time $t$','Interpreter','latex');
ylabel(ax2,'infinite-horizon $P_{\infty,1}(t)$','Interpreter','latex');
lg2 = legend(ax2,{'$[P_{\infty,1}]_{11}$','$[P_{\infty,1}]_{22}$','$[P_{\infty,1}]_{12}$'}, ...
      'Interpreter','latex','Location','northeastoutside'); set(lg2,'Box','off');
set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_tv_mjls.pdf'));
