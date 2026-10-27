%% EX_GAME_SD  Held-input sampled-data H-infinity game (Game chapter,
%  Section sec:game:ex:sd).  The game is posed on the held-input
%  sampled-data class of Section sec:pre:cov:sd.  Between samples the held
%  control u is a state, so the flow Riccati on chi = (x,u) eliminates the
%  disturbance w only; the sample value u_d is eliminated by the Schur step
%  at each sampling instant (periodic I-HJI-GRE of Definition def:game:GRE,
%  Theorem th:game:periodic), see MVST_GAME_SD_HELD.  Sampling period
%  Tp = 0.3, the value of the other sampled-data examples of the toolbox
%  (Tbarsd, DTbarsd).
%
%  Computes: the critical level gamma* (bisection on convergence of the
%  backward periodic iteration with Isaacs' condition, relative bracket 1e-6);
%  the saddle gains at gamma = 1.2 gamma* (sample gain Ku on x, maximizer gain
%  Kw(0^+) on chi); the critical level of the game with a continuously acting
%  minimizer (comparison, tolerance 1e-6); cross-check gamma* at Tp = 0.1.
%  Figures: convergence of the finite-horizon value at a sampling instant
%  (x-block, horizon T = N Tp) and the infinite-horizon value tr S vs gamma.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- game instance (scalar control u, scalar disturbance w) -----
A  = [0 1;-1 -1];
Bu = [0;1];         Bw = [1;0];
Q  = eye(2);        Ru = 1;
Tp = 0.3;                                    % sampling period

% ----- w-actuated diffusion channel (matching the performance SD plant of
%       Section sec:diss:ex:sd): D_{w,1} = 0.1[1;0], no state part -----
Dw1 = 0.1*[1;0];

% ----- critical level of the held-input sampled-data game -----
t0 = tic;
[gstar, glo] = mvst_game_gamma_sd_held(A,Bu,Bw,Q,Ru,Tp,Dw1,1e-6,200000);
oLo = mvst_game_sd_held(A,Bu,Bw,Q,Ru,glo,Tp,Dw1,200000);
assert(~oLo.ok && ~strcmp(oLo.reason,'maxiter'), 'lower bracket end must fail by escape/Isaacs');
gam = 1.2*gstar;
o = mvst_game_sd_held(A,Bu,Bw,Q,Ru,gam,Tp,Dw1,200000);
assert(o.ok, 'periodic iteration must converge at 1.2 gamma*');
Ku = o.Ku; Kw = o.Kw0;
fprintf(['SD game (held input, Tp=%.2f): gstar %.6f (infeasible end %.6f, %s), gam %.6f\n', ...
    '  at gam: periods %d, jump residual %.1e, min eig Theta %.4f, Puu(0+) %.4f\n', ...
    '  Ku = %s, Kw(0+) = %s, tr S = %.6f\n'], Tp, gstar, glo, oLo.reason, gam, ...
    o.N, o.jumpRes, o.thetaMin, o.Pplus(3,3), mat2str(Ku,6), mat2str(Kw,6), trace(o.S));

% ----- comparison: continuously acting minimizer (CT game with pivot Theta) -----
ctFeas = @(g) ctfeasible(A,Bu,Bw,Q,Ru,g,Dw1);
lo = 1.0; hi = 2.0; assert(ctFeas(hi) && ~ctFeas(lo));
while hi - lo > 1e-6*hi
    mid = 0.5*(lo+hi); if ctFeas(mid); hi = mid; else; lo = mid; end
end
gstarCT = hi;
fprintf('  continuous minimizer: gstar %.6f (diagnostic value 1.3533)\n', gstarCT);
assert(gstar >= gstarCT - 1e-6, 'held-input level must not be below the continuous one');

% ----- cross-check: Tp = 0.1 (former operating period, gamma* = 1.362739) -----
g01 = mvst_game_gamma_sd_held(A,Bu,Bw,Q,Ru,0.1,Dw1,1e-6,200000);
fprintf('  cross-check Tp=0.10: gstar %.6f (former value 1.362739); Tp=0.30 diagnostic value 1.3867\n', g01);
fprintf('  computation time %.1f s\n', toc(t0));

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','game_sd.tex'));
w.mat('GAsd',A);
w.mat('GBu',Bu);       w.mat('GBw',Bw);
w.num('GTpSd',Tp,1);
w.num('GGamSd',gstar);
w.num('GGamSdCT',gstarCT);
w.mat('GKuSd',Ku);     w.mat('GKwSd',Kw);
w.close();

% ----- figure 1: convergence of the value at a sampling instant -----
Nc = round(8/Tp);                                       % horizons up to T = 8 (27 periods at Tp = 0.3)
oc = mvst_game_sd_held(A,Bu,Bw,Q,Ru,gam,Tp,Dw1,Nc,0);   % no early stop
Tg = (1:oc.N)*Tp; P0 = oc.Shist;
game_conv_fig(Tg, P0, ...
    'Sampled-data game: convergence of the game value $P_T(0)$', ...
    fullfile(mvst_figdir(),'ex_game_sd_conv.pdf'), gc);

% ----- figure 2: infinite-horizon value versus gamma -----
gg = linspace(1.02*gstar, 3*gstar, 40); vv = zeros(size(gg));
for k = 1:numel(gg)
    ok_ = mvst_game_sd_held(A,Bu,Bw,Q,Ru,gg(k),Tp,Dw1,200000);
    assert(ok_.ok, 'value curve: iteration must converge');
    vv(k) = trace(ok_.S);
end
game_val_fig(gg, vv, gstar, ...
    'Sampled-data game value vs.\ disturbance level $\gamma$', ...
    fullfile(mvst_figdir(),'ex_game_sd_val.pdf'), gc);

% ----- combined 2-up: value convergence (top) + value vs gamma (bottom) -----
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tlc = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');
game_conv_ax(nexttile(tlc), Tg, P0, '(a) game-value convergence', gc);
game_val_ax(nexttile(tlc), gg, vv, gstar, '(b) value vs.\ disturbance level $\gamma$', gc);
mvst_savepdf(figc, fullfile(mvst_figdir(),'ex_game_sd_combined.pdf'));

% ========================================================================
function ok = ctfeasible(A,Bu,Bw,Q,Ru,g,Dw1)
%CTFEASIBLE  Continuous-minimizer game: bounded backward Riccati over [0,200].
try
    ws = warning('off','all');
    [t,Pt] = mvst_game_dre_ct(A,Bu,Bw,Q,Ru,g,[0 200],[],1e-9,{},{Dw1});
    warning(ws);
    nrm = arrayfun(@(k) norm(Pt(:,:,k)), 1:size(Pt,3));
    ok = abs(t(1)) < 1e-9 && all(isfinite(nrm)) && max(nrm) < 1e6;
catch
    ok = false;
end
end

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
