%% EX04_POISSON  Optimal control of a stochastic system with Poissonian jumps.
%  Solves the generalized algebraic Riccati equation (diffusion term
%  sum_l C_l'P C_l plus jump-rate term lam (J'PJ - P)) for the GBM +
%  Poissonian-jump system, then shows closed-loop sample paths (with the
%  random resets visible).  The long-term average cost rate is reported
%  separately in ex_average (fig:ex:average), so it is not replotted here.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- instance: unstable drift, geometric noise, expansive resets -----
A  = [0 1; 2 -0.5];  B  = [0;1];
C1 = 0.2*eye(2);     D1 = [0;0];
J1 = 1.25*eye(2);    lam = 0.5;
Q  = eye(2);         R  = 1;
W  = 0.05*eye(2);                       % additive noise (average-cost reading)

[P, K, info] = mvst_care_ct_poisson(A, B, {C1}, {D1}, {J1}, lam, Q, R);
rho = trace(P*W);
fprintf('Poisson-jump I-GRE: residual %.2e, K = [% .4f % .4f], rho* = %.4f\n', ...
    info.residual, K, rho);

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_poisson.tex'));
w.mat('Apoi',A); w.mat('Bpoi',B); w.mat('Cpoi',C1);
w.mat('Jpoi',J1); w.num('LamPoi',lam);
w.mat('Qpoi',Q); w.num('Rpoi',R);
w.mat('KpoiInf',K); w.num('RhoPoi',rho);
w.close();

% ----- stochastic closed-loop sample paths (Euler-Maruyama + thinning) -----
rng(7);
T = 12; dt = 1e-3; tgrid = 0:dt:T; nt = numel(tgrid); nruns = 30;
n = 2; sqW = sqrtm(W);
X = zeros(n,nt,nruns); U = zeros(1,nt,nruns);
Acl = A + B*K; C1cl = C1 + D1*K;
for r = 1:nruns
    x = [1;0];
    for k = 1:nt
        X(:,k,r) = x; U(1,k,r) = K*x;
        dWl = sqrt(dt)*randn;  dBeta = sqrt(dt)*randn(n,1);
        x = x + Acl*x*dt + C1cl*x*dWl + sqW*dBeta;
        if rand < lam*dt; x = J1*x; end          % Poissonian reset
    end
end

% ----- downsample paths for plotting -----
dsi = max(1, round(nt/1500)); ii = 1:dsi:nt;
tp = tgrid(ii); Xp = X(:,ii,:); Up = U(:,ii,:);

fig = figure('Color','w','Units','centimeters','Position',[2 2 16 10]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
for r = 1:nruns
    plot(ax1, tp, squeeze(Xp(1,:,r)), '-','LineWidth',0.5,'Color',[gc.order(1,:) 0.20]);
    plot(ax1, tp, squeeze(Xp(2,:,r)), '-','LineWidth',0.5,'Color',[gc.order(2,:) 0.20]);
end
h1 = plot(ax1, tp, mean(Xp(1,:,:),3), '-','LineWidth',2.2,'Color',gc.order(1,:));
h2 = plot(ax1, tp, mean(Xp(2,:,:),3), '-','LineWidth',2.2,'Color',gc.order(2,:));
ylabel(ax1,'states $x_1(t),\ x_2(t)$','Interpreter','latex');
lg1 = legend(ax1,[h1 h2],{'$x_1$ mean','$x_2$ mean'},'Interpreter','latex', ...
             'Location','northeast'); set(lg1,'Box','off');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
for r = 1:nruns
    plot(ax2, tp, squeeze(Up(1,:,r)), '-','LineWidth',0.5,'Color',[gc.order(1,:) 0.25]);
end
hu = plot(ax2, tp, mean(Up(1,:,:),3), '-','LineWidth',2.2,'Color',gc.order(1,:));
xlabel(ax2,'time $t$','Interpreter','latex');
ylabel(ax2,'control $u(t)=K^{\star}x$','Interpreter','latex');
lg2 = legend(ax2,hu,{'ensemble mean'},'Location','northeast'); set(lg2,'Box','off');
set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_poisson.pdf'));
