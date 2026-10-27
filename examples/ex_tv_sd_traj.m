%% EX_TV_SD_TRAJ  Sample paths and mean-square trace, time-varying sampled-data.
%  Trajectory companion of EX_TV_SD: the flow drift is time-varying,
%  A(t) = A + dA*exp(-t), and the input is held between samples and updated
%  by u(t_k^+) = K x(t_k) at each sampling instant with the (asymptotic)
%  optimal sampled-data gain K, to which the time-varying design settles.
%  We simulate an ensemble of stochastic sample paths and show (top) the
%  mean-square trace E[x'x] = trace E[x x'], (middle) the state sample paths
%  x_1(t),x_2(t) with their ensemble mean, and (bottom) the zero-order-hold
%  control input.  Companion figure to the dual-variable plot of EX_TV_SD.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
s = D.sd; n = size(s.A,1); m = size(s.B,2); Tbar = s.Tbar;
dA = [0 0; 1.5 0];                          % decaying drift perturbation (as in EX_TV_SD)

% ----- asymptotic (LTI) optimal sampled-data gain the TV design settles to -
out = mvst_igre_periodic_sd(s.A, s.B, s.Ac, s.Bc, s.Q, s.R, Tbar);
K = out.K;
fprintf('TV sampled-data trajectory: settled gain K = [% .4f % .4f]\n', K);

% ----- ensemble of sample paths, time-varying flow + sample-and-hold input -
dt = 1e-3; tf = 6; tt = 0:dt:tf; nt = numel(tt);
tk = Tbar:Tbar:(tf-1e-9);                    % sampling instants
rng(5); nreal = 200; nshow = 20;
xs = zeros(n,nt,nreal); us = zeros(1,nt,nreal);
for r = 1:nreal
    x = [1;0]; u = K*x;
    for k = 1:nt
        if any(abs(tt(k)-tk) < dt/2); u = K*x; end     % sample-and-hold update
        xs(:,k,r) = x; us(1,k,r) = u;
        if k < nt
            A_t = s.A + dA*exp(-tt(k));                 % time-varying flow drift
            inc = (A_t*x + s.B*u)*dt + (s.Ac*x + s.Bc*u)*(sqrt(dt)*randn);
            x = x + inc;
        end
    end
end
trX = squeeze(mean(sum(xs.^2,1),3));         % empirical E[x'x] = trace E[xx']
ip = 1:3:nt;                                 % plot subsample

% ----- figure: trace (top), state sample paths (middle), held input (bottom)
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 14]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, tt(ip), trX(ip), '-','LineWidth',2.0,'Color',gc.order(1,:));
for x0 = tk; xline(ax1, x0, ':','Color',[0.7 0.7 0.7]); end
ylabel(ax1,'$\mathrm{trace}\ \mathrm{E}[xx^\top]$','Interpreter','latex');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for r = 1:nshow
    plot(ax2, tt(ip), squeeze(xs(1,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(1,:) 0.25]);
    plot(ax2, tt(ip), squeeze(xs(2,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(2,:) 0.25]);
end
h1 = plot(ax2, tt(ip), mean(xs(1,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(1,:));
h2 = plot(ax2, tt(ip), mean(xs(2,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(2,:));
ylabel(ax2,'states $x_1(t),\ x_2(t)$','Interpreter','latex');
lg2 = legend(ax2,[h1 h2],{'$x_1$ mean','$x_2$ mean'},'Interpreter','latex','Location','northeast');
set(lg2,'Box','off'); set(ax2,'FontSize',11);

ax3 = nexttile(tl); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
for r = 1:nshow
    stairs(ax3, tt(ip), squeeze(us(1,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(1,:) 0.30]);
end
stairs(ax3, tt(ip), mean(us(1,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(1,:));
xlabel(ax3,'time $t$','Interpreter','latex');
ylabel(ax3,'control $u(t)$ (held)','Interpreter','latex');
set(ax3,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_oc_tv_sd_traj.pdf'));
