%% EX02_SAMPLEDDATA_GAINS  Sampled-data: P(tau) with its jump, gain, regimes.
%  Solves the periodic-sampling impulsive I-GRE.  The optimal gain K is a
%  constant matrix (exported).  The dual variable P(tau) is timer-dependent
%  WITH a jump at each sample; it is plotted over two periods so the jump is
%  visible.  The long-term average rate rho* is obtained from a forced
%  steady-state simulation.  Constants are exported to figures/data/.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
s = D.sd; n = size(s.A,1); m = size(s.B,2); Tbar = s.Tbar;

out = mvst_igre_periodic_sd(s.A,s.B,s.Ac,s.Bc,s.Q,s.R,Tbar);
K = out.K;
fprintf('Sampled-data gain K = [% .4f % .4f]\n', K);

% within-period dual variable P(tau) (timer-dependent), sampled for plotting.
% Here the state is AUGMENTED, z = [x; u] with x = (x_1,x_2) and a scalar
% held input u, so P(tau) is (n+m) x (n+m) = 3 x 3.  We plot every distinct
% entry: the full state block (x_1 and x_2 both shown), the state-input
% coupling block, and the input block.
tau = linspace(0, Tbar, 60);
Px1x1 = zeros(size(tau)); Px2x2 = Px1x1; Px1x2 = Px1x1;
Px1u  = Px1x1;            Px2u  = Px1x1; Puu   = Px1x1;
for k = 1:numel(tau)
    P = out.Pfun(tau(k));
    Px1x1(k) = P(1,1);   Px2x2(k) = P(2,2);   Px1x2(k) = P(1,2);
    Px1u(k)  = P(1,n+1); Px2u(k)  = P(2,n+1); Puu(k)   = P(n+1,n+1);
end

% long-term average rate via a forced steady-state simulation
Ab0 = [s.A s.B; zeros(m,n+m)]; Abl = [s.Ac s.Bc; zeros(m,n+m)];
seg = struct('Acl',Ab0,'Lc',Abl,'W',blkdiag(s.W,zeros(m)),'L',blkdiag(s.Q,s.R));
Jbar = [eye(n) zeros(n,m); K zeros(m,m)];
tf = 30; tk = Tbar:Tbar:(tf-1e-9);
av = mvst_simulate_impulsive(repmat({seg},1,numel(tk)+1), tk, @(X,k) Jbar*X*Jbar.', ...
        blkdiag(s.X0,zeros(m)), tf, struct('npts',24));
rho = av.runavg(end);
fprintf('rho* (forced steady state) = %.4f\n', rho);

% ---- export to LaTeX ----
w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_sampleddata.tex'));
w.mat('Asd', s.A); w.mat('Bsd', s.B); w.mat('Ksd', K);
w.mat('Qsd', s.Q); w.num('Rsd', s.R);
w.num('Tbarsd', Tbar); w.num('RhoSd', rho);
w.close();

% ---- figure: P(tau) over two periods (jump visible) + gain text ----
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
NaNb = NaN;                                   % break to show the jump
tt  = [tau, NaNb, tau+Tbar];
plot(ax, tt, [Px1x1 NaNb Px1x1], '-','LineWidth',2.2,'Color',gc.order(1,:));
plot(ax, tt, [Px2x2 NaNb Px2x2], '-','LineWidth',2.2,'Color',gc.order(2,:));
plot(ax, tt, [Px1x2 NaNb Px1x2], '-','LineWidth',2.2,'Color',gc.order(3,:));
plot(ax, tt, [Px1u  NaNb Px1u ], '-','LineWidth',2.2,'Color',gc.order(4,:));
plot(ax, tt, [Px2u  NaNb Px2u ], '-','LineWidth',2.2,'Color',gc.order(5,:));
plot(ax, tt, [Puu   NaNb Puu  ], '-','LineWidth',2.2,'Color',gc.order(6,:));
xline(ax, Tbar, ':','Color',[0.6 0.6 0.6]);
xlabel(ax,'timer $\tau$ (sampling period $\bar T$)','Interpreter','latex');
ylabel(ax,'entries of $P(\tau)$','Interpreter','latex');
lg = legend(ax, {'$P_{x_1x_1}$','$P_{x_2x_2}$','$P_{x_1x_2}$', ...
      '$P_{x_1u}$','$P_{x_2u}$','$P_{uu}$'}, ...
      'Interpreter','latex','Location','northeastoutside'); set(lg,'Box','off');
set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_sampleddata_gains.pdf'));
