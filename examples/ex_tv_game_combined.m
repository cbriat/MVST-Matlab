%% EX_TV_GAME_COMBINED  Consolidated time-varying game-value convergence.
%  Native, vertically-stacked replacement for the former three-across
%  minipage tiling (fig:game:tv:mjls).  Reproduces, in one figure with
%  readable fonts, the convergence of the finite-horizon time-varying game
%  value P_T(0) as the horizon grows, for the three classes:
%    row 1: Markov-jump   (EX_TV_GAME_MJLS),
%    row 2: sampled-data  (EX_TV_GAME_SD, held-input game, Tp = 0.3),
%    row 3: switched      (EX_TV_GAME_SWITCHED),
%  each under the genuinely time-varying drift A(t) = A0 + dA*exp(-t) and
%  gamma = 1.2*gamma*.
%  The curves are taken from
%  the per-class scripts themselves (each is run in its own workspace), so the
%  printed figure carries exactly their data: the diffusion channels of each
%  base example and, for the sampled-data row, the held-input periodic game.
%  The per-class scripts
%  rewrite their own data macros and single-panel PDFs (same values).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

dM = grab(fullfile(here,'ex_tv_game_mjls.m'));
dS = grab(fullfile(here,'ex_tv_game_sd.m'));
dW = grab(fullfile(here,'ex_tv_game_switched.m'));
close all

% ---- combined figure: three subsystems stacked vertically ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 20]);
tl  = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); draw_game_ax(ax1, dM.Tg, dM.P0, gc, false); title(ax1,'(a) Markov-jump','Interpreter','latex','FontWeight','normal');
ax2 = nexttile(tl); draw_game_ax(ax2, dS.Tg, dS.P0, gc, false); title(ax2,'(b) Sampled-data','Interpreter','latex','FontWeight','normal');
ax3 = nexttile(tl); h = draw_game_ax(ax3, dW.Tg, dW.P0, gc, true); title(ax3,'(c) Switched','Interpreter','latex','FontWeight','normal');

Tmax = max([dM.Tg(:); dS.Tg(:); dW.Tg(:)]);
set([ax1 ax2 ax3], 'XLim', [0 Tmax]);                % common horizon axis

lg = legend(ax3, h, {'$P_{11}$','$P_{12}$','$P_{22}$'}, ...
            'Interpreter','latex','Orientation','horizontal');
lg.Layout.Tile = 'south'; set(lg,'Box','off');

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_tv_game_combined.pdf'));

% ========================================================================
function d = grab(script)
%GRAB  Run a per-class script in this function's workspace; return its curve.
run(script);
d.Tg = Tg; d.P0 = P0;
end

% ========================================================================
function h = draw_game_ax(ax, Tg, P0, gc, showx)
%DRAW_GAME_AX  Distinct entries of the game value P0(:,:,k) versus horizon.
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h(1) = plot(ax, Tg, reshape(P0(1,1,:),1,[]), '-o','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(1,:));
h(2) = plot(ax, Tg, reshape(P0(1,2,:),1,[]), '-s','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(2,:));
h(3) = plot(ax, Tg, reshape(P0(2,2,:),1,[]), '-^','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(3,:));
ylabel(ax,'$P_T(0)$ entries','Interpreter','latex');
if showx, xlabel(ax,'horizon $T$','Interpreter','latex'); end
set(ax,'FontSize',11);
end
