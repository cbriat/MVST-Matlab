%% FIG_PERIODIC_DWELL  Periodic dwell-time I-GRE: P(tau) with its jump.
%  Illustrates the periodic-impulses (timer-based) I-GRE on a generic
%  second-order instance: the dual variable P(tau) flows over one period
%  and is reset by the discrete-time jump Riccati at tau = Tbar.  Plotted
%  over two periods so the jump is explicit.  Saves figures/oc_periodic.pdf.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
gc = mvst_gemcolors();

A = [0 1;0 0]; B = [0;1]; Ac = 0.10*[0 0;1 0]; Bc = [0;0.10];
Q = eye(2); R = 1.0; Tbar = 0.6; n = 2; m = 1;
out = mvst_igre_periodic_sd(A,B,Ac,Bc,Q,R,Tbar);

% Augmented state z = [x; u] with x = (x_1,x_2) and a scalar held input u,
% so P(tau) is (n+m) x (n+m) = 3 x 3.  Plot every distinct entry: the full
% state block (both x_1 and x_2), the state-input coupling, and the input.
tau = linspace(0, Tbar, 60);
Px1x1 = zeros(size(tau)); Px2x2 = Px1x1; Px1x2 = Px1x1;
Px1u  = Px1x1;            Px2u  = Px1x1; Puu   = Px1x1;
for k = 1:numel(tau)
    P = out.Pfun(tau(k));
    Px1x1(k) = P(1,1);   Px2x2(k) = P(2,2);   Px1x2(k) = P(1,2);
    Px1u(k)  = P(1,n+1); Px2u(k)  = P(2,n+1); Puu(k)   = P(n+1,n+1);
end
NaNb = NaN; tt = [tau, NaNb, tau+Tbar];

fig = figure('Color','w','Units','centimeters','Position',[2 2 15 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, tt, [Px1x1 NaNb Px1x1], '-','LineWidth',2.2,'Color',gc.order(1,:));
plot(ax, tt, [Px2x2 NaNb Px2x2], '-','LineWidth',2.2,'Color',gc.order(2,:));
plot(ax, tt, [Px1x2 NaNb Px1x2], '-','LineWidth',2.2,'Color',gc.order(3,:));
plot(ax, tt, [Px1u  NaNb Px1u ], '-','LineWidth',2.2,'Color',gc.order(4,:));
plot(ax, tt, [Px2u  NaNb Px2u ], '-','LineWidth',2.2,'Color',gc.order(5,:));
plot(ax, tt, [Puu   NaNb Puu  ], '-','LineWidth',2.2,'Color',gc.order(6,:));
xline(ax, Tbar, ':','Color',[0.6 0.6 0.6]);
xlabel(ax,'timer $\tau$','Interpreter','latex');
ylabel(ax,'entries of $P(\tau)$','Interpreter','latex');
lg = legend(ax, {'$P_{x_1x_1}$','$P_{x_2x_2}$','$P_{x_1x_2}$', ...
      '$P_{x_1u}$','$P_{x_2u}$','$P_{uu}$'}, ...
      'Interpreter','latex','Location','northeastoutside'); set(lg,'Box','off');
set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'oc_periodic.pdf'));
