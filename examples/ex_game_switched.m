%% EX_GAME_SWITCHED  State-feedback H-infinity game for a two-mode switched system.
%  Zero-sum game with minimizer control u and maximizer disturbance w.  For
%  each mode the critical disturbance level gamma*_i is found by bisection;
%  gstar = max(gstar_1,gstar_2) is the WORST-MODE level, a LOWER bound for
%  the switched game under a minimum dwell-time (staying in one mode is
%  admissible; switching generally raises the level, and forced switching
%  can lower it below every per-mode level).  The saddle-point gains and
%  the infinite-horizon game value are evaluated at gamma = 1.2*gstar over a
%  long horizon (mode 1).  Two single-panel figures are produced: the
%  convergence of the finite-horizon game value P_T(0) and the infinite-horizon
%  game value versus the disturbance level gamma.
% The critical
% levels use the LONG-HORIZON criterion of the game library (horizon 1000,
% relative bracket 1e-6, uniform Isaacs margin lambda_min(Theta) >= 1e-6 gamma^2
% along the whole backward solution).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- two-mode game instance (scalar control u, scalar disturbance w) -----
A1 = [0 1;-2 -1];   A2 = [0 1;-1 -2];       % two modes
Bu = [0;1];         Bw = [1;0];
Q  = eye(2);        Ru = 1;    Td = 0.5;     % minimum dwell-time

% ----- w-actuated diffusion channel (matching the
%       performance switched system of Section sec:diss:ex:switched): a shared
%       disturbance-multiplicative column D_{w,1} = 0.1[1;0] per mode, no
%       state-multiplicative part -----
Dwc = {0.1*[1;0]};

% ----- worst-mode critical level (lower bound for the switched game) -----
Tlong = 1000; rtol = 1e-6;                   % long-horizon criterion (see header)
[gstar1,~,inf1] = mvst_game_gamma_ct(A1,Bu,Bw,Q,Ru,Tlong,{},Dwc,rtol);
[gstar2,~,inf2] = mvst_game_gamma_ct(A2,Bu,Bw,Q,Ru,Tlong,{},Dwc,rtol);
gstar  = max(gstar1,gstar2);                 % worst-mode level (lower bound)

% ----- dwell-time certified UPPER bound on the achievable level: mode-wise
%       central gains + closed-loop MDT dissipativity certificate ---------
% With the long-horizon criterion passed through
% (..., Dwc, Tlong, rtol) the certificate bisection starts from the worst-mode
% level 1.747201 and returns 2.027026 (bracket (2.020201, 2.027026]).  The
% certificate search (MVST_DISS_DWELL_CERT, Nelder-Mead) is not monotone in the level, the
% certificate is FOUND at 2.027, 2.040, 2.050 with either set of central
% gains, so both values are valid certified upper bounds and 2.027 is the
% tighter one.
% The call below passes Tlong and rtol (long-horizon value 2.027, the tighter certified bound).
[gstarDT,~,infDT] = mvst_game_gamma_ct_dwell({A1,A2}, Bu, Bw, Q, Ru, Td, [], Dwc, Tlong, rtol);
gam    = 1.2*gstar;

% ----- mode-1 infinite-horizon dual and saddle-point gains -----
[~,Pt,Kut,Kwt] = mvst_game_dre_ct(A1,Bu,Bw,Q,Ru,gam,[0 12],[],[],{},Dwc);
Pinf = Pt(:,:,1);                            % t=0 of a long horizon ~ inf-horizon
Ku   = Kut(:,:,1);
Kw   = Kwt(:,:,1);

fprintf('Switched game: dwell-time certified achievable level %.4f (Td=%.2f)\n', ...
    gstarDT, Td);
fprintf('Switched game: gstar1 %.4f, gstar2 %.4f, gstar %.4f, gam %.4f\n', ...
    gstar1, gstar2, gstar, gam);
fprintf(['Switched game brackets: mode1 (%.6f,%.6f] mode2 (%.6f,%.6f]; min eig Theta %.4f, %.4f;\n' ...
         '  MDT certificate bracket (%.6f,%.6f]\n'], inf1.lo, inf1.hi, inf2.lo, inf2.hi, ...
    inf1.thetaMin, inf2.thetaMin, infDT.lo, infDT.hi);
% mode-2 maximizer gain at the same level (quoted in the (G4) sentence)
[~,~,~,Kwt2] = mvst_game_dre_ct(A2,Bu,Bw,Q,Ru,gam,[0 12],[],[],{},Dwc);
Kw2 = Kwt2(:,:,1);
fprintf(['Switched game: Ku1 = [%.4f %.4f], Kw1 = [%.4f %.4f], Kw2 = [%.4f %.4f];\n' ...
         '  lambda_min(Q - gam^2 Kw''Kw): mode1 %.4f, mode2 %.4f\n'], Ku, Kw, Kw2, ...
    min(eig(Q - gam^2*(Kw.'*Kw))), min(eig(Q - gam^2*(Kw2.'*Kw2))));

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','game_sw.tex'));
w.mat('GAswOne',A1); w.mat('GAswTwo',A2);
w.mat('GBu',Bu);     w.mat('GBw',Bw);
w.num('GGamSw',gstar);
w.num('GGamSwDT',gstarDT); w.num('GTdSw',Td);
w.mat('GKuSw',Ku);   w.mat('GKwSw',Kw);
w.close();

% ----- figure 1: convergence of the game value P_T(0) (mode 1) -----
Tg = linspace(0.1,8,40); P0 = zeros(2,2,numel(Tg));
for k = 1:numel(Tg)
    [~,Ptk] = mvst_game_dre_ct(A1,Bu,Bw,Q,Ru,gam,[0 Tg(k)],[],[],{},Dwc);
    P0(:,:,k) = Ptk(:,:,1);
end
game_conv_fig(Tg, P0, ...
    'Switched game: convergence of the game value $P_T(0)$', ...
    fullfile(mvst_figdir(),'ex_game_switched_conv.pdf'), gc);

% ----- figure 2: infinite-horizon game value versus gamma (mode 1) -----
gg = linspace(1.02*gstar, 3*gstar, 40); vv = zeros(size(gg));
for k = 1:numel(gg)
    [~,Ptk] = mvst_game_dre_ct(A1,Bu,Bw,Q,Ru,gg(k),[0 12],[],[],{},Dwc);
    vv(k) = trace(Ptk(:,:,1));
end
game_val_fig(gg, vv, gstar, ...
    'Switched game value vs.\ disturbance level $\gamma$', ...
    fullfile(mvst_figdir(),'ex_game_switched_val.pdf'), gc);

% ----- combined 2-up: value convergence (left) + value vs gamma (right) -----
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tlc = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');
game_conv_ax(nexttile(tlc), Tg, P0, '(a) game-value convergence', gc);
game_val_ax(nexttile(tlc), gg, vv, gstar, '(b) value vs.\ disturbance level $\gamma$', gc);
mvst_savepdf(figc, fullfile(mvst_figdir(),'ex_game_switched_combined.pdf'));

% ========================================================================
function game_conv_fig(Tg, P0, ttl, fname, gc)
%GAME_CONV_FIG  Standalone wrapper around GAME_CONV_AX.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
game_conv_ax(axes(fig), Tg, P0, '', gc);
mvst_savepdf(fig, fname);
end

% ========================================================================
function game_conv_ax(ax, Tg, P0, ttl, gc)
%GAME_CONV_AX  Plot the distinct entries of the game value P_T(0) vs horizon T.
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
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
if ~isempty(ttl); title(ax, ttl, 'Interpreter','latex','FontWeight','normal'); end
end

% ========================================================================
function game_val_fig(gg, vv, gstar, ttl, fname, gc)
%GAME_VAL_FIG  Standalone wrapper around GAME_VAL_AX.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
game_val_ax(axes(fig), gg, vv, gstar, '', gc);
mvst_savepdf(fig, fname);
end

% ========================================================================
function game_val_ax(ax, gg, vv, gstar, ttl, gc)
%GAME_VAL_AX  Infinite-horizon game value tr P_inf versus gamma, with gamma*.
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h(1) = plot(ax, gg, vv, '-', 'LineWidth',1.8, 'Color',gc.order(1,:));
h(2) = xline(ax, gstar, '--', 'LineWidth',1.8, 'Color',[0.5 0.5 0.5]);
xlim(ax, [gstar - 0.03*(gg(end)-gstar), gg(end)]);   % curve fills the axis
xlabel(ax,'$\gamma$','Interpreter','latex');
ylabel(ax,'game value $\mathrm{tr}\,P_\infty$','Interpreter','latex');
lg = legend(ax, h, {'$\mathrm{tr}\,P_\infty(\gamma)$','$\gamma^\star$'}, ...
            'Interpreter','latex','Orientation','horizontal', ...
            'Location','southoutside');
set(lg,'Box','off'); set(ax,'FontSize',11);
if ~isempty(ttl); title(ax, ttl, 'Interpreter','latex','FontWeight','normal'); end
end
