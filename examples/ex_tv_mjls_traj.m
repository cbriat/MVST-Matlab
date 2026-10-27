%% EX_TV_MJLS_TRAJ  Sample paths and mean-square trace, time-varying MJLS.
%  Trajectory companion of EX_TV_MJLS: under the time-varying optimal
%  feedback u = K_{sigma}(t) x of the time-varying Markov jump diffusion
%  A_i(t) = A_i + dA*exp(-t), we simulate an ensemble of stochastic sample
%  paths and show (top) the mean-square trace  E[x'x] = trace E[x x'],
%  (middle) the individual state sample paths x_1(t),x_2(t) with their
%  ensemble mean, and (bottom) the realized control input.  The empirical
%  trace of the second moment is the mean-square quantity the second-moment
%  design controls in expectation, not pathwise.  Companion figure to the
%  dual-variable plot of EX_TV_MJLS.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data(); gc = mvst_gemcolors();
base = D.mjls; n = base.n; m = base.m; M = base.M;
dA = [0 0; 1.5 0];                          % decaying drift perturbation (as in EX_TV_MJLS)

% ----- build the time-varying modes A_i(t) = A_i + dA*exp(-t) -----
md = cell(1,M);
for i = 1:M
    b = base.modes{i};
    md{i} = mvst_mode_tv(n, m, ...
        'A', @(t) b.A + dA*exp(-t), 'B', b.B, 'Q', b.Q, 'R', b.R, ...
        'Ac', b.Ac, 'Bc', b.Bc, 'W', b.W, 'X0', b.X0);
end
sys = mvst_system(md, 'Pi', base.Pi, 'pi0', base.pi0);

% ----- time-varying infinite-horizon dual and gains K_i(t) -----
Tinf = 8.0;
ri = mvst_solve_infinite_tv(sys, Tinf);     % ri.K{i}: m x n x nt on grid ri.t

% ----- gains resampled onto the simulation grid (settle to LTI values) -----
dt = 1e-3; T = Tinf; tt = 0:dt:T; nt = numel(tt);
[tu, iu] = unique(ri.t(:));                 % strictly increasing knots for interp1
Kg = cell(1,M);
for i = 1:M
    Kmat = reshape(ri.K{i}, m*n, []).';     % nt_ri x (m*n), column-major vec of K
    tcl  = min(max(tt, tu(1)), tu(end));
    Kg{i} = interp1(tu, Kmat(iu,:), tcl, 'pchip');   % nt x (m*n)
end

% ----- ensemble of Markov-jump sample paths under u = K_sigma(t) x -----
rng(11); nreal = 200; nshow = 20;
xs = zeros(n,nt,nreal); us = zeros(m,nt,nreal);
for r = 1:nreal
    x = [1;0]; s = pick(base.pi0);
    for k = 1:nt
        Ki = reshape(Kg{s}(k,:), m, n);
        xs(:,k,r) = x; us(:,k,r) = Ki*x;
        if k < nt
            b = base.modes{s};
            Acl = (b.A + dA*exp(-tt(k))) + b.B*Ki;
            inc = Acl*x*dt;
            for l = 1:b.Nc
                Lcl = b.Ac(:,:,l) + b.Bc(:,:,l)*Ki;
                inc = inc + (Lcl*x)*(sqrt(dt)*randn);
            end
            x = x + inc;
            % Markov mode update (Euler on the generator).
            pr = base.Pi(s,:)*dt; pr(s) = 0; pr = max(pr,0);
            if rand < sum(pr); s = pick(pr/sum(pr)); end
        end
    end
end
trX = squeeze(mean(sum(xs.^2,1),3));        % empirical E[x'x] = trace E[xx']
ip = 1:8:nt;                                % plot subsample

% ----- figure: trace (top), state sample paths (middle), control (bottom) -
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 14]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, tt(ip), trX(ip), '-','LineWidth',2.0,'Color',gc.order(1,:));
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
    plot(ax3, tt(ip), squeeze(us(1,ip,r)), '-','LineWidth',0.4,'Color',[gc.order(1,:) 0.30]);
end
plot(ax3, tt(ip), mean(us(1,ip,:),3), '-','LineWidth',2.0,'Color',gc.order(1,:));
xlabel(ax3,'time $t$','Interpreter','latex');
ylabel(ax3,'control $u(t)=K_{\sigma}(t)x$','Interpreter','latex');
set(ax3,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_oc_tv_mjls_traj.pdf'));

% ========================================================================
function i = pick(p)
p = p(:).'/sum(p); i = find(rand <= cumsum(p), 1, 'first');
if isempty(i); i = numel(p); end
end
