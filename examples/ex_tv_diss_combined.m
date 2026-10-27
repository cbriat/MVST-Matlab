%% EX_TV_DISS_COMBINED  Consolidated time-varying dissipativity convergence.
%  Native, vertically-stacked replacement for the former three-across
%  minipage tiling (fig:diss:tv:mjls).  Reproduces, in one figure with
%  readable fonts, the convergence of the finite-horizon time-varying
%  available-storage matrix S_T(0) = -P_T(0) >= 0 as the horizon grows, for
%  the three classes:
%    row 1: Markov-jump        (EX_TV_DISS_MJLS, stable mode, both channels),
%    row 2: sampled-data       (EX_TV_DISS_SD, column D_{w,1}),
%    row 3: switched           (EX_TV_DISS_SWITCHED, mode 1, column D_{w,1}),
%  each under the genuinely time-varying drift A(t) = A0 + dA*exp(-t).
%  The curves are taken from
%  the per-class scripts themselves (each run in its own workspace), so the
%  printed figure carries exactly their data: the diffusion channels of the
%  time-invariant examples, the storage level 1.05 max(gamma*, g0) with g0 the
%  gain of the plotted time-varying system from t_0 = 0, and the storage
%  itself (S >= 0).  This script writes no data macros of
%  its own (the per-class scripts rewrite theirs, same values).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

dM = grab(fullfile(here,'ex_tv_diss_mjls.m'));
dS = grab(fullfile(here,'ex_tv_diss_sd.m'));
dW = grab(fullfile(here,'ex_tv_diss_switched.m'));
close all

% ---- combined figure: three subsystems stacked vertically ----------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 20]);
tl  = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); draw_conv_ax(ax1, dM.Tg, dM.S0, gc, false); title(ax1,'(a) Markov-jump','Interpreter','latex','FontWeight','normal');
ax2 = nexttile(tl); draw_conv_ax(ax2, dS.Tg, dS.S0, gc, false); title(ax2,'(b) Sampled-data','Interpreter','latex','FontWeight','normal');
ax3 = nexttile(tl); h = draw_conv_ax(ax3, dW.Tg, dW.S0, gc, true); title(ax3,'(c) Switched','Interpreter','latex','FontWeight','normal');

lg = legend(ax3, h, {'$-P_{11}$','$-P_{12}$','$-P_{22}$'}, ...
            'Interpreter','latex','Orientation','horizontal');
lg.Layout.Tile = 'south'; set(lg,'Box','off');

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_tv_diss_combined.pdf'));

% ========================================================================
function d = grab(script)
%GRAB  Run a per-class script in this function's workspace; return its curve.
run(script);
d.Tg = Tg; d.S0 = S0;
end

% ========================================================================
function h = draw_conv_ax(ax, Tg, S, gc, showx)
%DRAW_CONV_AX  Distinct entries of the storage matrix S(:,:,k) versus horizon.
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h(1) = plot(ax, Tg, reshape(S(1,1,:),1,[]), '-o','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(1,:));
h(2) = plot(ax, Tg, reshape(S(1,2,:),1,[]), '-s','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(2,:));
h(3) = plot(ax, Tg, reshape(S(2,2,:),1,[]), '-^','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(3,:));
ylabel(ax,'$-P_T(0)$ entries','Interpreter','latex');
if showx, xlabel(ax,'horizon $T$','Interpreter','latex'); end
set(ax,'FontSize',11);
end
