%% EX_TV_SWITCHED_TRAJ  Sample paths and mean-square trace, time-varying switched.
%  Trajectory companion of EX_TV_SWITCHED: each mode's flow drift is
%  time-varying, A_i(t) = A_i + dA*exp(-t), and the closed loop uses the
%  mode-dependent (asymptotic) optimal gains K_i to which the time-varying
%  design settles, applied under the dwell-time switching law.  The state is
%  continuous at switches (identity jump); only the active flow changes.  We
%  simulate an ensemble of stochastic sample paths and show (top) the
%  mean-square trace E[x'x] = trace E[x x'], (middle) the state sample paths
%  x_1(t),x_2(t) with their ensemble mean, and (bottom) the continuous
%  control input.  Companion figure to the dual-variable plot of EX_TV_SWITCHED.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
modes = D.sw.modes; Td = D.sw.Tdwell; M = numel(modes); n = modes{1}.n;
dA = [0 0; 1.5 0];                          % decaying drift perturbation (as in EX_TV_SWITCHED)

% ----- mode-dependent (asymptotic LTI) gains the TV design settles to -----
K = cell(1,M);
for i = 1:M
    [~, K{i}] = mvst_care_ct(modes{i});
    fprintf('  K_%d = [% .4f % .4f]\n', i, K{i});
end

% ----- ensemble of switched sample paths, time-varying flow A_i(t) -----
dt = 1e-3; tf = 8; tt = 0:dt:tf; nt = numel(tt);
tk = Td:Td:(tf-1e-9); seq = 1 + mod(0:numel(tk), M);   % active mode per interval
modeof = @(t) seq(min(floor(t/Td)+1, numel(seq)));
rng(3); nreal = 200; nshow = 20;
xs = zeros(n,nt,nreal); us = zeros(1,nt,nreal);
for r = 1:nreal
    x = [1;0];
    for k = 1:nt
        i = modeof(tt(k)); md = modes{i};
        xs(:,k,r) = x; us(1,k,r) = K{i}*x;
        if k < nt
            Acl = (md.A + dA*exp(-tt(k))) + md.B*K{i};
            inc = Acl*x*dt;
            for l = 1:md.Nc
                inc = inc + (md.Ac(:,:,l)+md.Bc(:,:,l)*K{i})*x*(sqrt(dt)*randn);
            end
            x = x + inc;
        end
    end
end
trX = squeeze(mean(sum(xs.^2,1),3));         % empirical E[x'x] = trace E[xx']
ip = 1:3:nt;                                 % plot subsample

% ----- figure: trace (top), state sample paths (middle), control (bottom) -
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
for x0 = tk; xline(ax2, x0, ':','Color',[0.7 0.7 0.7]); end
ylabel(ax2,'states $x_1(t),\ x_2(t)$','Interpreter','latex');
lg2 = legend(ax2,[h1 h2],{'$x_1$ mean','$x_2$ mean'},'Interpreter','latex','Location','northeast');
set(lg2,'Box','off'); set(ax2,'FontSize',11);

ax3 = nexttile(tl); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
for r = 1:nshow
    plot(ax3, tt(ip), squeeze(us(1,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(1,:) 0.30]);
end
plot(ax3, tt(ip), mean(us(1,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(1,:));
for x0 = tk; xline(ax3, x0, ':','Color',[0.7 0.7 0.7]); end
xlabel(ax3,'time $t$','Interpreter','latex');
ylabel(ax3,'control $u(t)=K_{\sigma}x$','Interpreter','latex');
set(ax3,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_oc_tv_switched_traj.pdf'));
