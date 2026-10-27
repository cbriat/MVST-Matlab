%% EX03_SWITCHED  Optimal control of a switched system (dwell-time).
%  Designs mode-dependent gains (per-mode infinite-horizon I-GRE) and shows
%  the (deterministic) state-covariance trace, a state sample path, and the
%  continuous control input under a dwell-time switching law.  The state is
%  continuous at switches (identity jump); only the active flow changes.
%  (The exact minimum/range dwell-time co-design uses the timer-dependent
%  I-GRE; here each mode uses its own stabilizing gain, valid for a
%  sufficiently large dwell-time.)
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
modes = D.sw.modes; Td = D.sw.Tdwell; M = numel(modes); n = modes{1}.n;

K = cell(1,M); L = cell(1,M); Acl = cell(1,M); Lc = cell(1,M);
for i = 1:M
    [~, K{i}] = mvst_care_ct(modes{i});
    md = modes{i}; Acl{i} = md.A + md.B*K{i};
    Lc{i} = zeros(n,n,md.Nc);
    for l = 1:md.Nc; Lc{i}(:,:,l) = md.Ac(:,:,l) + md.Bc(:,:,l)*K{i}; end
    L{i} = md.Q + K{i}.'*md.R*K{i};
    fprintf('  K_%d = [% .4f % .4f]\n', i, K{i});
end

tf = 8; tk = Td:Td:(tf-1e-9);  seq = 1 + mod(0:numel(tk), M);
intervals = cell(1, numel(tk)+1);
for i = 1:numel(intervals)
    mi = seq(i);
    intervals{i} = struct('Acl',Acl{mi},'Lc',Lc{mi},'W',zeros(n),'L',L{mi});
end
covsim = mvst_simulate_impulsive(intervals, tk, @(X,k) X, modes{1}.X0, tf, struct('npts',40));
trX = squeeze(arrayfun(@(k) trace(covsim.X(:,:,k)), 1:size(covsim.X,3)));

% Stochastic sample paths of the switched closed loop (u = K_sigma x).
dt = 1e-3; rng(3); tt = 0:dt:tf; nreal = 20;
xs = zeros(n,numel(tt),nreal); us = zeros(numel(tt),nreal);
modeof = @(t) seq(min(floor(t/Td)+1, numel(seq)));
for r = 1:nreal
    x = [1;0];
    for k = 1:numel(tt)
        i = modeof(tt(k)); md = modes{i};
        xs(:,k,r) = x; us(k,r) = K{i}*x;
        inc = (md.A + md.B*K{i})*x*dt;
        for l = 1:md.Nc; inc = inc + (md.Ac(:,:,l)+md.Bc(:,:,l)*K{i})*x*(sqrt(dt)*randn); end
        x = x + inc;
    end
end
ip = 1:3:numel(tt);

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 14]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, covsim.t, trX, '-','LineWidth',2.0,'Color',gc.order(1,:));
for tkk = tk; xline(ax1, tkk, ':','Color',[0.6 0.6 0.6]); end
ylabel(ax1,'$\mathrm{trace}\ \mathrm{E}[xx^\top]$', 'Interpreter', 'latex');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for r = 1:nreal
    plot(ax2, tt(ip), squeeze(xs(1,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(1,:) 0.25]);
    plot(ax2, tt(ip), squeeze(xs(2,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(2,:) 0.25]);
end
h1 = plot(ax2, tt(ip), mean(xs(1,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(1,:));
h2 = plot(ax2, tt(ip), mean(xs(2,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(2,:));
for tkk = tk; xline(ax2, tkk, ':','Color',[0.6 0.6 0.6]); end
ylabel(ax2,'states $x_1(t),\ x_2(t)$','Interpreter','latex');
lg2 = legend(ax2,[h1 h2],{'$x_1$ mean','$x_2$ mean'},'Interpreter','latex','Location','northeast');
set(lg2,'Box','off'); set(ax2,'FontSize',11);

ax3 = nexttile(tl); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
for r = 1:nreal
    plot(ax3, tt(ip), us(ip,r), '-','LineWidth',0.4,'Color',[gc.order(1,:) 0.30]);
end
plot(ax3, tt(ip), mean(us(ip,:),2), '-','LineWidth',2.0,'Color',gc.order(1,:));
for tkk = tk; xline(ax3, tkk, ':','Color',[0.6 0.6 0.6]); end
xlabel(ax3,'time $t$'); ylabel(ax3,'control $u(t)=K_{\sigma}x$', 'Interpreter', 'latex'); set(ax3,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_switched.pdf'));
