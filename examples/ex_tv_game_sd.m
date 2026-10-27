%% EX_TV_GAME_SD  Time-varying held-input sampled-data H-infinity game.
%  Time-varying counterpart of EX_GAME_SD: the flow drift decays in time,
%  A(t) = A0 + dA*exp(-t), so the game is genuinely time-varying but settles to
%  the time-invariant held-input game as t grows.
%  The game is the HELD-INPUT sampled-data game of
%  EX_GAME_SD (Section sec:pre:cov:sd), with the same sampling period Tp = 0.3.
%  On each inter-sample interval the state is chi = (x,u), the flow Riccati
%    dP/dt + Abar(t)'P + P Abar(t) + Qbar + P Bwbar Theta^{-1} Bwbar'P = 0,
%    Theta = gam^2 I - Dwbar'P Dwbar,  Abar(t) = [A(t) Bu; 0 0],
%  eliminates w only, and the sample value u_d is eliminated by the Schur step
%    P(t_k^-) = blkdiag(Pxx - Pxu Puu^{-1} Pux, 0),  P = P(t_k^+),
%  at each sampling instant t_k = k Tp.  The critical level gamma* is that of
%  the ASYMPTOTIC drift A0 (MVST_GAME_GAMMA_SD_HELD, identical call to
%  EX_GAME_SD, so TVGGamSd = GGamSd).  The convergence figure integrates the
%  TIME-VARYING periodic HJI-GRE backward from P(T^-) = 0 over horizons
%  T = N Tp at gamma = 1.2 gamma*, and plots the x-block value at t_0 = 0.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- game instance (same base as EX_GAME_SD) -----
A0 = [0 1;-1 -1];                            % asymptotic LTI drift (t->inf)
Bu = [0;1];         Bw = [1;0];
Q  = eye(2);        Ru = 1;
Tp = 0.3;                                    % sampling period (as EX_GAME_SD)
dA = [0 0; 1.5 0];                           % decaying drift perturbation
Afun = @(t) A0 + dA*exp(-t);                 % time-varying flow drift

% ----- w-actuated diffusion channel (same as the base EX_GAME_SD): a
%       disturbance-multiplicative column D_{w,1}=0.1[1;0], no state part -----
Dw1 = 0.1*[1;0];

% ----- critical disturbance level at the asymptotic drift A0 -----
t0 = tic;
gstar = mvst_game_gamma_sd_held(A0,Bu,Bw,Q,Ru,Tp,Dw1,1e-6,200000);
gam   = 1.2*gstar;
fprintf('TV-SD game (held input, Tp=%.2f, asymptotic A0): gstar %.6f, gam %.6f\n', Tp, gstar, gam);

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_game_sd.tex'));
w.mat('TVGdAsd',dA); w.num('TVGGamSd',gstar);
w.close();

% ----- convergence of the TIME-VARYING game value at t_0 = 0 -----
Nmax = round(8/Tp);                          % horizons T = N Tp up to 8
Tg = (1:Nmax)*Tp; P0 = zeros(2,2,Nmax); okAll = true;
for N = 1:Nmax
    [S0, ok, why] = tv_held_value(Afun,Bu,Bw,Q,Ru,gam,Tp,Dw1,N);
    if ~ok; okAll = false; fprintf('  horizon N=%d: %s\n', N, why); end
    P0(:,:,N) = S0;
end
assert(okAll, 'time-varying held-input iteration must stay feasible at 1.2 gamma*');
% settling of the value at t_0 = 0 with the horizon (it settles to the
% time-varying infinite-horizon value at t_0 = 0, not to the fixed point at A0,
% since A(0) differs from A0)
fprintf('  |P_T(0) - P_{T-Tp}(0)| at T=%.1f: %.2e; tr P_T(0) = %.6f\n', ...
    Tg(end), norm(P0(:,:,end) - P0(:,:,end-1)), trace(P0(:,:,end)));
fprintf('  computation time %.1f s\n', toc(t0));

game_conv_fig(Tg, P0, ...
    'Time-varying sampled-data game: convergence of the game value', ...
    fullfile(mvst_figdir(),'ex_tv_game_sd_conv.pdf'), gc);

% ========================================================================
function [S, ok, why] = tv_held_value(Afun, Bu, Bw, Q, Ru, gam, Tp, Dw1, N)
%TV_HELD_VALUE  x-block game value at t_0 = 0^- of the time-varying held-input
%   game on [0, N Tp], from P(N Tp^-) = 0, by backward flow/jump iteration.
n = size(Bu,1); mu = size(Bu,2); mw = size(Bw,2); N3 = n + mu;
Bwb = [Bw; zeros(mu,mw)]; Dwb = [Dw1; zeros(mu,mw)]; Qb = blkdiag(Q,Ru);
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
Pm = zeros(N3); ok = true; why = ''; S = nan(n);
for j = N:-1:1
    ta = (j-1)*Tp; tb = j*Tp;
    try
        [t, Y] = ode45(@(t,p) rhs(t,p,gam,Afun,Bu,Bwb,Dwb,Qb,n,mu,mw), [tb ta], Pm(:), opts);
    catch
        ok = false; why = 'integrator failure on the flow'; return
    end
    if abs(t(end) - ta) > 1e-9*Tp || any(~isfinite(Y(end,:)))
        ok = false; why = 'escape or Isaacs failure on the flow'; return
    end
    Pp = reshape(Y(end,:), N3, N3); Pp = (Pp+Pp.')/2;          % P(t_{j-1}^+)
    Puu = Pp(n+1:end,n+1:end);
    [~, pc] = chol((Puu+Puu.')/2);
    if pc ~= 0; ok = false; why = 'Puu not positive definite at a sample'; return; end
    S = Pp(1:n,1:n) - Pp(1:n,n+1:end)*(Puu\Pp(n+1:end,1:n)); S = (S+S.')/2;
    Pm = blkdiag(S, zeros(mu));                                  % P(t_{j-1}^-)
end
end

function dp = rhs(t, p, gam, Afun, Bu, Bwb, Dwb, Qb, n, mu, mw)
%RHS  dP/dt = -(Abar(t)'P + P Abar(t) + Qbar + P Bwbar Theta^{-1} Bwbar'P).
N3 = n + mu;
Ab = [Afun(t) Bu; zeros(mu, N3)];
P = reshape(p, N3, N3); P = (P+P.')/2;
Th = gam^2*eye(mw) - Dwb.'*P*Dwb; Th = (Th+Th.')/2;
[~, pc] = chol(Th);
if pc ~= 0; dp = nan(N3*N3,1); return; end        % Isaacs fails on the flow
v = P*Bwb;
G = Ab.'*P + P*Ab + Qb + v*(Th\v.'); G = (G+G.')/2;
dp = -G(:);
end

% ========================================================================
function game_conv_fig(Tg, P0, ttl, fname, gc) %#ok<INUSL>
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
