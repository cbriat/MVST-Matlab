%% EX02_SAMPLEDDATA  Optimal periodic sampled-data control (impulsive I-GRE).
%  Solves the periodic-sampling impulsive I-GRE for the double-integrator
%  diffusion, then shows the (deterministic) state-covariance trace, a state
%  sample path, and the zero-order-hold control input (a staircase, updated
%  every Tbar).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
s = D.sd; n = size(s.A,1); m = size(s.B,2); Tbar = s.Tbar;

out = mvst_igre_periodic_sd(s.A,s.B,s.Ac,s.Bc,s.Q,s.R,Tbar);
K = out.K;
fprintf('Periodic sampled-data gain K = [% .4f % .4f] (iters %d, res %.1e)\n', K, out.iters, out.residual);

% Deterministic augmented covariance (flow + sample-and-hold resets).
Ab0 = [s.A s.B; zeros(m,n+m)];  Abl = [s.Ac s.Bc; zeros(m,n+m)];
seg = struct('Acl',Ab0,'Lc',Abl,'W',zeros(n+m),'L',blkdiag(s.Q,s.R));
Jbar = [eye(n) zeros(n,m); K zeros(m,m)];
tf = 6; tk = Tbar:Tbar:(tf-1e-9);
covsim = mvst_simulate_impulsive(repmat({seg},1,numel(tk)+1), tk, @(X,k) Jbar*X*Jbar.', ...
         blkdiag(s.X0,zeros(m)), tf, struct('npts',40));
trX1 = squeeze(arrayfun(@(k) trace(covsim.X(1:n,1:n,k)), 1:size(covsim.X,3)));

% Stochastic sample paths with the zero-order-hold input.
dt = 1e-3; rng(5); tt = 0:dt:tf; nreal = 20;
xs = zeros(n,numel(tt),nreal); us = zeros(numel(tt),nreal);
for r = 1:nreal
    x = [1;0]; u = K*x;
    for k = 1:numel(tt)
        if any(abs(tt(k)-tk) < dt/2); u = K*x; end       % sample-and-hold update
        xs(:,k,r) = x; us(k,r) = u;
        inc = (s.A*x + s.B*u)*dt + (s.Ac*x + s.Bc*u)*(sqrt(dt)*randn);
        x = x + inc;
    end
end
ip = 1:3:numel(tt);                                       % plot subsample

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 14]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, covsim.t, trX1, '-','LineWidth',2.0,'Color',gc.order(1,:));
ylabel(ax1,'$\mathrm{trace}\ \mathrm{E}[xx^\top]$', 'Interpreter', 'latex');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for r = 1:nreal
    plot(ax2, tt(ip), squeeze(xs(1,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(1,:) 0.25]);
    plot(ax2, tt(ip), squeeze(xs(2,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(2,:) 0.25]);
end
h1 = plot(ax2, tt(ip), mean(xs(1,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(1,:));
h2 = plot(ax2, tt(ip), mean(xs(2,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(2,:));
ylabel(ax2,'states $x_1(t),\ x_2(t)$', 'Interpreter', 'latex');
lg2 = legend(ax2,[h1 h2],{'$x_1$ mean','$x_2$ mean'},'Interpreter','latex','Location','northeast');
set(lg2,'Box','off'); set(ax2,'FontSize',11);

ax3 = nexttile(tl); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
for r = 1:nreal
    stairs(ax3, tt(ip), us(ip,r), '-','LineWidth',0.4,'Color',[gc.order(1,:) 0.30]);
end
stairs(ax3, tt(ip), mean(us(ip,:),2), '-','LineWidth',2.0,'Color',gc.order(1,:));
xlabel(ax3,'time $t$', 'Interpreter', 'latex'); ylabel(ax3,'control $u(t)$ (held)', 'Interpreter', 'latex'); set(ax3,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_sampleddata.pdf'));
