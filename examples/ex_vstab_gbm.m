function ex_vstab_gbm()
%EX_VSTAB_GBM  Geometric Brownian motion: sample paths and moments, for the
%   two mean-square regimes of the vector-stability appendix (Appendix E,
%   Figures fig:vstab:gbmSPastable and fig:vstab:gbmSPunstable).  
%
%   MODEL.  One-dimensional geometric Brownian motion (eq:vstab:gbm)
%       dx = mu x dt + sigma x dW,   x(0) = x0,
%   integrated by Euler-Maruyama with a nonnegativity floor (max(0,.)):
%       x(t+dt) = max(0, x + mu x dt + sigma x sqrt(dt) N(0,1)).
%   The analytic moments are E[x] = x0 exp(mu t) and
%   E[x^2] = x0^2 exp((2 mu + sigma^2) t), so mu < 0 gives mean stability
%   and 2 mu + sigma^2 < 0 gives mean-square stability.
%
%   TWO REGIMES (faithful to the legacy code and to the book captions).
%     (i)  Asymptotically-stable variance: mu = -2, sigma = sqrt(3), T = 4.
%          2 mu + sigma^2 = -1 < 0, both moments exponentially stable.
%     (ii) Unstable variance:             mu = -2, sigma = sqrt(5), T = 3.
%          2 mu + sigma^2 = +1 > 0, only the first moment is stable; the
%          variance diverges while paths still converge (a.s. stable).
%   Common parameters: dt = 1e-3, x0 = 10, N = 10 paths, seeded rng.
%
%   OUTPUT.  Two vertically stacked combined vector PDFs (figures only, no
%   data macros): figures/vecapp_gbmAStable.pdf and
%   figures/vecapp_gbmUnstable.pdf, each with (a) normalized sample paths
%   and (b) mean and standard deviation.

here   = fileparts(mfilename('fullpath'));
root   = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

local_gbm(figdir, gc, -2, sqrt(3), 4, ...
    'vecapp_gbmAStable.pdf', 42, ...
    'both moments exponentially stable');
local_gbm(figdir, gc, -2, sqrt(5), 3, ...
    'vecapp_gbmUnstable.pdf', 43, ...
    'only the first moment exponentially stable');
end

% ------------------------------------------------------------------------
function local_gbm(figdir, gc, mu, sigma, T, fname, seed, tag)
dt = 1e-3;  x0 = 10;  N = 10;
t  = 0:dt:T;  L = numel(t);
rng(seed,'twister');                     % fixed, documented seed per regime

% ------------------------------------------------ Euler-Maruyama sample paths
x = zeros(N, L);  x(:,1) = x0;
for i = 1:L-1
    dW = sqrt(dt)*randn(N,1);
    x(:,i+1) = max(0, x(:,i) + mu*x(:,i)*dt + sigma*x(:,i).*dW);
end
xn = x/max(x(:));                        % normalized sample paths

% ------------------------------------------------ analytic moments
Ex   = x0*exp(mu*t);
Exx  = x0^2*exp((2*mu+sigma^2)*t);
stdx = sqrt(max(0, Exx - Ex.^2));
fprintf('  GBM mu=%.3g sigma=%.4g (2mu+s^2=%.3g, %s)\n', ...
        mu, sigma, 2*mu+sigma^2, tag);

% --------------------------------------------------------------- combined PDF
fig = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
set(ax,'ColorOrder',gc.order); set(ax,'LineStyleOrder','-');
plot(ax, t, xn, 'LineWidth', 1.0);
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$x(t)$ (normalized)','Interpreter','latex');
title(ax,'(a) normalized sample paths','Interpreter','latex');
xlim(ax,[0 T]); set(ax,'FontSize',13);

ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, t, Ex,   '-', 'Color', gc.sapphire, 'LineWidth', 1.6);
plot(ax, t, stdx, '-', 'Color', gc.ruby,     'LineWidth', 1.6);
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'moments','Interpreter','latex');
title(ax,'(b) mean and standard deviation','Interpreter','latex');
legend(ax,{'Mean','Standard deviation'}, 'Location','best', ...
       'Interpreter','latex');
xlim(ax,[0 T]); set(ax,'FontSize',13);

mvst_savepdf(fig, fullfile(figdir,fname));
fprintf('  wrote %s (T=%g, dt=%g, N=%d paths)\n', fname, T, dt, N);
end
