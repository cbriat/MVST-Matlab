%% EX_GAME_POISSON  State-feedback H-infinity game with Poissonian jumps.
%  Zero-sum game (minimizer u, maximizer w) for the GBM + Poissonian-jump
%  system.  The flow carries the diffusion term sum_l C_l'P C_l; the jump
%  now carries BOTH players' channels: at each Poisson event the state is
%  reset x -> J_1 x + BuJ u_J + BwJ w_J and the players solve a static
%  ISAACS SADDLE (minimizer u_J with weight Rd, maximizer w_J with weight
%  -gamma^2 I).  This reduces the jump-rate term to
%  lam (J_1' Phat J_1 - P) with Phat = P - P B (B'PB + D)^{-1} B'P,
%  B = [BuJ BwJ], D = blkdiag(Rd, -gamma^2 I).  The critical level gamma*
%  is found by bisection (boundedness of the backward game Riccati), the
%  saddle gains are evaluated at gamma = 1.2 gamma*.  Figures: convergence
%  of the game value P_T(0), and the game value versus gamma.
% The critical
% levels use the LONG-HORIZON criterion of the game library (horizon 1000,
% relative bracket 1e-6, uniform Isaacs margin lambda_min(Theta) >= 1e-6 gamma^2
% along the whole backward solution).
% Here the pivot is the jump maximizer block gamma^2 - BwJ'P BwJ.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

A  = [0 1; 2 -0.5];  Bu = [0;1];  Bw = [1;0];
C1 = 0.2*eye(2);     J1 = 1.25*eye(2);   lam = 0.5;
Q  = eye(2);         Ru = 1;
% ----- jump channels: both control and disturbance enter the reset -----
BuJ = [0;0.5];       BwJ = [0.5;0];      Rd = 1;

[gstar,~,infP] = mvst_game_gamma_ct_poisson(A, Bu, Bw, {C1}, {J1}, lam, Q, Ru, 1000, ...
                                   {BuJ}, {BwJ}, {Rd}, 1e-6);
fprintf('Poisson game bracket (%.6f,%.6f], min eig jump pivot at the feasible end %.3e\n', ...
    infP.lo, infP.hi, infP.thetaMin);
gam   = 1.2*gstar;
[~, Pt, Kut, Kwt] = mvst_game_dre_ct_poisson(A, Bu, Bw, {C1}, {J1}, lam, ...
                                 Q, Ru, gam, [0 12], [], 1e-8, {BuJ}, {BwJ}, {Rd});
Ku = Kut(:,:,1); Kw = Kwt(:,:,1);
Pinf = Pt(:,:,1);
% ----- static Isaacs jump saddle gains at the settled value -----
Br = [BuJ BwJ]; Dr = blkdiag(Rd, -gam^2);
KJ  = -(Br.'*Pinf*Br + Dr)\(Br.'*Pinf*J1);
KuJ = KJ(1,:);   KwJ = KJ(2,:);
fprintf(['Poisson game (jump-io): gstar %.4f, gam %.4f\n' ...
         '  flow Ku = [% .4f % .4f], Kw = [% .4f % .4f]\n' ...
         '  jump KuJ = [% .4f % .4f], KwJ = [% .4f % .4f]\n'], ...
    gstar, gam, Ku, Kw, KuJ, KwJ);
fprintf('Poisson game: lambda_min(Q - gam^2 Kw''Kw) = %.4f\n', min(eig(Q - gam^2*(Kw.'*Kw))));

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','game_poisson.tex'));
w.mat('GApoi',A); w.mat('GJpoi',J1); w.num('GLamPoi',lam);
w.mat('GBuJpoi',BuJ); w.mat('GBwJpoi',BwJ); w.num('GRdPoi',Rd);
w.num('GGamPoi',gstar); w.mat('GKuPoi',Ku); w.mat('GKwPoi',Kw);
w.mat('GKuJPoi',KuJ); w.mat('GKwJPoi',KwJ);
w.close();

% ----- figure 1: convergence of the game value P_T(0) -----
Tg = linspace(0.1, 8, 40); P0 = zeros(2,2,numel(Tg));
for k = 1:numel(Tg)
    [~, Ptk] = mvst_game_dre_ct_poisson(A, Bu, Bw, {C1}, {J1}, lam, ...
                        Q, Ru, gam, [0 Tg(k)], [], 1e-8, {BuJ}, {BwJ}, {Rd});
    P0(:,:,k) = Ptk(:,:,1);
end
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h = gobjects(3,1);
h(1) = plot(ax, Tg, reshape(P0(1,1,:),1,[]), '-o','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(1,:));
h(2) = plot(ax, Tg, reshape(P0(1,2,:),1,[]), '-s','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(2,:));
h(3) = plot(ax, Tg, reshape(P0(2,2,:),1,[]), '-^','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(3,:));
xlabel(ax,'$T$','Interpreter','latex');
ylabel(ax,'game value $P_T(0)$ entries','Interpreter','latex');
lg = legend(ax, h, {'$P_{11}$','$P_{12}$','$P_{22}$'},'Interpreter','latex', ...
            'Location','southeast');
set(lg,'Box','off'); set(ax,'FontSize',11);
mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_game_poisson_conv.pdf'));

% ----- figure 2: infinite-horizon game value versus gamma -----
gg = linspace(1.02*gstar, 3*gstar, 40); vv = zeros(size(gg));
for k = 1:numel(gg)
    [~, Ptk] = mvst_game_dre_ct_poisson(A, Bu, Bw, {C1}, {J1}, lam, ...
                        Q, Ru, gg(k), [0 12], [], 1e-8, {BuJ}, {BwJ}, {Rd});
    vv(k) = trace(Ptk(:,:,1));
end
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
ax2 = axes(fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
h2 = gobjects(2,1);
h2(1) = plot(ax2, gg, vv, '-','LineWidth',2.0,'Color',gc.order(1,:));
h2(2) = xline(ax2, gstar, '--','LineWidth',1.8,'Color',[0.5 0.5 0.5]);
xlim(ax2, [gstar - 0.03*(gg(end)-gstar), gg(end)]);   % curve fills the axis
xlabel(ax2,'disturbance level $\gamma$','Interpreter','latex');
ylabel(ax2,'game value $\mathrm{tr}\,P_\infty$','Interpreter','latex');
lg2 = legend(ax2, h2, {'$\mathrm{tr}\,P_\infty(\gamma)$','$\gamma^\star$'}, ...
             'Interpreter','latex','Location','northeast');
set(lg2,'Box','off'); set(ax2,'FontSize',11);
mvst_savepdf(fig2, fullfile(mvst_figdir(),'ex_game_poisson_val.pdf'));

% ----- combined 2-up (stacked): value convergence (top) + value vs gamma (bottom) -----
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tlc = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');
axc = nexttile(tlc); hold(axc,'on'); box(axc,'on'); grid(axc,'on');
hc = gobjects(3,1);
hc(1) = plot(axc, Tg, reshape(P0(1,1,:),1,[]), '-o','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(1,:));
hc(2) = plot(axc, Tg, reshape(P0(1,2,:),1,[]), '-s','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(2,:));
hc(3) = plot(axc, Tg, reshape(P0(2,2,:),1,[]), '-^','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(3,:));
xlabel(axc,'$T$','Interpreter','latex');
ylabel(axc,'game value $P_T(0)$ entries','Interpreter','latex');
lg = legend(axc, hc, {'$P_{11}$','$P_{12}$','$P_{22}$'},'Interpreter','latex','Location','southeast');
set(lg,'Box','off'); set(axc,'FontSize',11);
title(axc,'(a) game-value convergence','Interpreter','latex','FontWeight','normal');
axc = nexttile(tlc); hold(axc,'on'); box(axc,'on'); grid(axc,'on');
hc2 = gobjects(2,1);
hc2(1) = plot(axc, gg, vv, '-','LineWidth',2.0,'Color',gc.order(1,:));
hc2(2) = xline(axc, gstar, '--','LineWidth',1.8,'Color',[0.5 0.5 0.5]);
xlim(axc, [gstar - 0.03*(gg(end)-gstar), gg(end)]);
xlabel(axc,'disturbance level $\gamma$','Interpreter','latex');
ylabel(axc,'game value $\mathrm{tr}\,P_\infty$','Interpreter','latex');
lg2 = legend(axc, hc2, {'$\mathrm{tr}\,P_\infty(\gamma)$','$\gamma^\star$'}, ...
             'Interpreter','latex','Location','northeast');
set(lg2,'Box','off'); set(axc,'FontSize',11);
title(axc,'(b) value vs.\ disturbance level $\gamma$','Interpreter','latex','FontWeight','normal');
mvst_savepdf(figc, fullfile(mvst_figdir(),'ex_game_poisson_combined.pdf'));
