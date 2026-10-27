%% EX_TV_GAME_MJLS  Time-varying H-infinity game for a Markov-jump system.
%  Time-varying counterpart of EX_GAME_MJLS: the mode-1 flow drift decays in
%  time, A_1(t) = A_1 + dA*exp(-t), so the game is genuinely time-varying but
%  settles to the LTI Markov-jump game as t grows.  The critical disturbance
%  level gamma* is computed at the ASYMPTOTIC drift A0 (worst-mode level
%  gstar = max(gstar_1,gstar_2), a per-mode reference value: the coupled
%  game's level may be larger or smaller); the convergence figure integrates
%  the TIME-VARYING game (HJI) Riccati backward over growing horizons at
%  gamma = 1.2*gstar, showing the finite-horizon game value settling.
% The critical
% levels use the LONG-HORIZON criterion of the game library (horizon 1000,
% relative bracket 1e-6, uniform Isaacs margin lambda_min(Theta) >= 1e-6 gamma^2
% along the whole backward solution).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- two-mode game instance (same base as EX_GAME_MJLS) -----
A1 = [0 1;-2 -1];   A2 = [0 1;-3 -1.5];     % two modes (asymptotic A0=A1)
Bu = [0;1];         Bw = [1;0];
Q  = eye(2);        Ru = 1;
dA = [0 0; 1.5 0];                           % decaying drift perturbation
Afun = @(t) A1 + dA*exp(-t);                 % time-varying flow drift (mode 1)

% ----- w-actuated diffusion channels (same as the base EX_GAME_MJLS):
%       per-mode state-multiplicative C_{i,1}=0.2 I_2 and shared
%       disturbance-multiplicative column D_{w,1}=0.1[1;0] -----
Ci  = 0.2*eye(2);   Dwc = {0.1*[1;0]};

% ----- critical disturbance level at the asymptotic drift A0 -----
gstar1 = mvst_game_gamma_ct(A1,Bu,Bw,Q,Ru,1000,{Ci},Dwc,1e-6);   % long-horizon criterion
gstar2 = mvst_game_gamma_ct(A2,Bu,Bw,Q,Ru,1000,{Ci},Dwc,1e-6);
gstar  = max(gstar1,gstar2);                 % worst-mode level (reference)
gam    = 1.2*gstar;

fprintf('TV-MJLS game (asymptotic A0): gstar1 %.4f, gstar2 %.4f, gstar %.4f, gam %.4f\n', ...
    gstar1, gstar2, gstar, gam);

% ----- maximizer gain at t_0 = 0 (quoted in the (G4) sentence of the text) -----
[~,~,~,Kwt0] = mvst_game_dre_ct(Afun,Bu,Bw,Q,Ru,gam,[0 30],[],[],{Ci},Dwc);
Kw0 = Kwt0(:,:,1);
fprintf('  Kw(0) = [%.4f %.4f], lambda_min(Q - gam^2 Kw''Kw) = %.4f\n', ...
    Kw0, min(eig(Q - gam^2*(Kw0.'*Kw0))));

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_game_mjls.tex'));
w.mat('TVGdAmjls',dA); w.num('TVGGamMjls',gstar);
w.close();

% ----- figure: convergence of the TIME-VARYING game value P_T(0) -----
Tg = linspace(0.1,8,40); P0 = zeros(2,2,numel(Tg));
for k = 1:numel(Tg)
    [~,Ptk] = mvst_game_dre_ct(Afun,Bu,Bw,Q,Ru,gam,[0 Tg(k)],[],[],{Ci},Dwc);
    P0(:,:,k) = Ptk(:,:,1);
end
game_conv_fig(Tg, P0, ...
    'Time-varying Markov-jump game: convergence of the game value', ...
    fullfile(mvst_figdir(),'ex_tv_game_mjls_conv.pdf'), gc);

% ========================================================================
function game_conv_fig(Tg, P0, ttl, fname, gc)
%GAME_CONV_FIG  Plot the distinct entries of the game value P_T(0) vs horizon T.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h(1) = plot(ax, Tg, reshape(P0(1,1,:),1,[]), '-o', 'LineWidth',1.8, ...
            'MarkerSize',3, 'Color',gc.order(1,:));
h(2) = plot(ax, Tg, reshape(P0(1,2,:),1,[]), '-s', 'LineWidth',1.8, ...
            'MarkerSize',3, 'Color',gc.order(2,:));
h(3) = plot(ax, Tg, reshape(P0(2,2,:),1,[]), '-^', 'LineWidth',1.8, ...
            'MarkerSize',3, 'Color',gc.order(3,:));
xlabel(ax,'$T$','Interpreter','latex');
ylabel(ax,'game value matrix $P_T(0)$ entries','Interpreter','latex');
lg = legend(ax, h, {'$P_{11}$','$P_{12}$','$P_{22}$'}, ...
            'Interpreter','latex','Orientation','horizontal', ...
            'Location','southoutside');
set(lg,'Box','off'); set(ax,'FontSize',11);
mvst_savepdf(fig, fname);
end
