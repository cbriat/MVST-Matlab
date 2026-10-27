%% EX_GAME_MJLS  State-feedback H-infinity game for a two-mode Markov-jump system.
%  Zero-sum game with minimizer control u and maximizer disturbance w.  For
%  each mode the critical disturbance level gamma*_i is found by bisection;
%  gstar = max(gstar_1,gstar_2) is the WORST-MODE level, a per-mode reference
%  value: the critical level of the coupled game solves the coupled HJI
%  equations and can be larger or smaller (stochastic switching can degrade
%  or improve the achievable attenuation).  The saddle-point
%  gains and the infinite-horizon game value are evaluated at gamma = 1.2*gstar
%  over a long horizon (mode 1).  Two single-panel figures are produced: the
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
A1 = [0 1;-2 -1];   A2 = [0 1;-3 -1.5];     % two modes
Bu = [0;1];         Bw = [1;0];
Q  = eye(2);        Ru = 1;
Pig = [-1 1; 2 -2];                          % Markov chain generator

% ----- w-actuated diffusion channels (matching the
%       performance MJLS of Section sec:diss:ex:mjls): per-mode
%       state-multiplicative C_{i,1} = 0.2 I_2 and shared
%       disturbance-multiplicative column D_{w,1} = 0.1[1;0] -----
Ci   = 0.2*eye(2);                           % per-mode state channel
Chan = {{Ci},{Ci}};                          % 1 x M cell of 1 x q channel cells
Dwc  = {0.1*[1;0]};                          % shared disturbance-mult. column

% ----- worst-mode critical level (per-mode reference value) -----
Tlong = 1000; rtol = 1e-6;                   % long-horizon criterion (see header)
[gstar1,~,inf1] = mvst_game_gamma_ct(A1,Bu,Bw,Q,Ru,Tlong,{Ci},Dwc,rtol);
[gstar2,~,inf2] = mvst_game_gamma_ct(A2,Bu,Bw,Q,Ru,Tlong,{Ci},Dwc,rtol);
gstar  = max(gstar1,gstar2);                 % worst-mode level (reference)

% ----- EXACT coupled critical level: bisection on the boundedness of the
%       COUPLED game (HJI) Riccati equations (generator coupling) --------
[gstarC,~,infC] = mvst_game_gamma_ct_mjls({A1,A2}, Bu, Bw, Q, Ru, Pig, Tlong, Chan, Dwc, rtol);
gam    = 1.2*gstarC;

% ----- coupled infinite-horizon duals and mode-1 saddle gains -----------
Pc   = mjls_game_dre({A1,A2}, Bu, Bw, Q, Ru, Pig, gam, 30, Chan, Dwc);
Pinf = Pc{1};
Ku   = -(Ru\(Bu.'*Pinf));
% Maximizer saddle gain of the coupled HJI-GRE with the diffusion channels
% (the w-block of eq:game:Ks):
%   Kw = Theta_1^{-1}(Bw'P_1 + sum_l Dwc{l}'P_1 C_{1,l}),
%   Theta_1 = gam^2 I - sum_l Dwc{l}'P_1 Dwc{l}.
Th1  = gam^2 - Dwc{1}.'*Pinf*Dwc{1};
Kw   = Th1\(Bw.'*Pinf + Dwc{1}.'*Pinf*Chan{1}{1});

fprintf('MJLS game: gstar1 %.4f, gstar2 %.4f, worst %.4f\n', ...
    gstar1, gstar2, gstar);
fprintf('MJLS game: exact coupled level %.4f, gam %.4f\n', gstarC, gam);
fprintf(['MJLS game brackets: mode1 (%.6f,%.6f] mode2 (%.6f,%.6f] coupled (%.6f,%.6f];\n' ...
         '  min eig Theta at the feasible ends: %.4f, %.4f, %.4f\n'], ...
    inf1.lo, inf1.hi, inf2.lo, inf2.hi, infC.lo, infC.hi, ...
    inf1.thetaMin, inf2.thetaMin, infC.thetaMin);
fprintf('MJLS game: Ku = [%.4f %.4f], Kw = [%.4f %.4f], lambda_min(Q - gam^2 Kw''Kw) = %.4f\n', ...
    Ku, Kw, min(eig(Q - gam^2*(Kw.'*Kw))));

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','game_mjls.tex'));
w.mat('GAmjlsOne',A1); w.mat('GAmjlsTwo',A2);
w.mat('GBu',Bu);       w.mat('GBw',Bw);
w.num('GGamMjls',gstar);
w.mat('GGenMjls',Pig); w.num('GGamMjlsC',gstarC);
w.mat('GKuMjls',Ku);   w.mat('GKwMjls',Kw);
w.close();

% ----- figure 1: convergence of the game value P_T(0) (mode 1) -----
Tg = linspace(0.1,8,40); P0 = zeros(2,2,numel(Tg));
for k = 1:numel(Tg)
    [~,Ptk] = mvst_game_dre_ct(A1,Bu,Bw,Q,Ru,gam,[0 Tg(k)],[],[],{Ci},Dwc);
    P0(:,:,k) = Ptk(:,:,1);
end
game_conv_fig(Tg, P0, ...
    'Markov-jump game: convergence of the game value $P_T(0)$', ...
    fullfile(mvst_figdir(),'ex_game_mjls_conv.pdf'), gc);

% ----- figure 2: infinite-horizon game value versus gamma (mode 1) -----
gg = linspace(1.02*gstar, 3*gstar, 40); vv = zeros(size(gg));
for k = 1:numel(gg)
    [~,Ptk] = mvst_game_dre_ct(A1,Bu,Bw,Q,Ru,gg(k),[0 12],[],[],{Ci},Dwc);
    vv(k) = trace(Ptk(:,:,1));
end
game_val_fig(gg, vv, gstar, ...
    'Markov-jump game value vs.\ disturbance level $\gamma$', ...
    fullfile(mvst_figdir(),'ex_game_mjls_val.pdf'), gc);

% ----- combined 2-up: value convergence (left) + value vs gamma (right) -----
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tlc = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');
game_conv_ax(nexttile(tlc), Tg, P0, '(a) game-value convergence', gc);
game_val_ax(nexttile(tlc), gg, vv, gstar, '(b) value vs.\ disturbance level $\gamma$', gc);
mvst_savepdf(figc, fullfile(mvst_figdir(),'ex_game_mjls_combined.pdf'));

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

% ========================================================================
function Pc = mjls_game_dre(Ac, Bu, Bw, Q, Ru, Pig, g, Tlong, Chan, Dwc)
%MJLS_GAME_DRE  Settled duals of the coupled game Riccati at level g, folding
%   the per-mode state-multiplicative channels CHAN (1 x M cell of 1 x q cells)
%   and the shared disturbance-multiplicative channels DWC (1 x q cell): each
%   mode gains sum_l Chan{i}{l}'P_i Chan{i}{l} and the maximizer Hessian
%   Theta_i = g^2 I - sum_l Dwc{l}'P_i Dwc{l} (mirrors MVST_GAME_GAMMA_CT_MJLS).
M = numel(Ac); n = size(Q,1); m = size(Bw,2);
if nargin < 9  || isempty(Chan); Chan = repmat({{}},1,M); end
if nargin < 10 || isempty(Dwc);  Dwc  = {}; end
qd = numel(Dwc); for ii = 1:M; qd = max(qd, numel(Chan{ii})); end
if qd > 0
    if numel(Dwc) < qd; Dwc(end+1:qd) = {zeros(n,m)}; end
    for ii = 1:M; if numel(Chan{ii}) < qd; Chan{ii}(end+1:qd) = {zeros(n)}; end; end
end
hasNoise = any(cellfun(@(x) nnz(x)>0, Dwc));
for ii = 1:M; hasNoise = hasNoise || any(cellfun(@(x) nnz(x)>0, Chan{ii})); end
BuRu = Bu*(Ru\Bu.'); BwBw = Bw*Bw.';
opts = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
              'Events',@(s,p) mvst_odecap(p,1e6));
[~, Y] = ode45(@(s,p) rhs(p), [0 Tlong], zeros(M*n*n,1), opts);
Pc = cell(1,M);
for i = 1:M
    Pi_ = reshape(Y(end,(i-1)*n*n+1:i*n*n).',n,n); Pc{i} = (Pi_+Pi_.')/2;
end
    function dp = rhs(p)
        P = cell(1,M);
        for i = 1:M
            Pi_ = reshape(p((i-1)*n*n+1:i*n*n),n,n); P{i} = (Pi_+Pi_.')/2;
        end
        dp = zeros(M*n*n,1);
        for i = 1:M
            cpl = zeros(n);
            for j = 1:M; cpl = cpl + Pig(i,j)*P{j}; end
            if ~hasNoise
                G = Ac{i}.'*P{i} + P{i}*Ac{i} + cpl + Q ...
                    - P{i}*BuRu*P{i} + (1/g^2)*P{i}*BwBw*P{i};
            else
                mult = zeros(n); v = P{i}*Bw; Th = g^2*eye(m);
                for l = 1:qd
                    mult = mult + Chan{i}{l}.'*P{i}*Chan{i}{l};
                    v    = v    + Chan{i}{l}.'*P{i}*Dwc{l};
                    Th   = Th   - Dwc{l}.'*P{i}*Dwc{l};
                end
                % Escape guard (mirrors the library backend): if the maximizer
                % Hessian Theta loses positivity the game value escapes; return a
                % NaN derivative so ode45 stops (with mvst_odecap) rather than
                % grinding on a near-singular Theta\v solve.
                [~,pc] = chol((Th+Th.')/2);
                if pc ~= 0; dp = nan(M*n*n,1); return; end
                G = Ac{i}.'*P{i} + P{i}*Ac{i} + cpl + Q ...
                    - P{i}*BuRu*P{i} + mult + v*(Th\v.');
            end
            G = (G+G.')/2;
            dp((i-1)*n*n+1:i*n*n) = G(:);
        end
    end
end
