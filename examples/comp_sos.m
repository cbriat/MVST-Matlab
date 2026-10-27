%% COMP_SOS  Timer-dependent matrix inequality over one period.
%  Verifies, on the periodic dwell-time sampled-data instance, that the
%  timer-dependent dual variable P(tau) stays positive semidefinite across
%  one period [0,Tbar].  The smallest eigenvalue of the symmetric part of
%  P(tau) is plotted against the timer; it remains nonnegative, certifying
%  the timer-dependent matrix inequality (the continuous analogue of a
%  sum-of-squares feasibility check).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
gc = mvst_gemcolors();

A = [0 1;0 0]; B = [0;1]; Ac = 0.10*[0 0;1 0]; Bc = [0;0.10];
Q = eye(2); R = 1.0; Tbar = 0.6;
out = mvst_igre_periodic_sd(A,B,Ac,Bc,Q,R,Tbar);

% ---- smallest eigenvalue of the symmetric part of P(tau) over a period ----
tau = linspace(0, Tbar, 200);
lam = zeros(size(tau));
for k = 1:numel(tau)
    P = out.Pfun(tau(k));
    lam(k) = min(real(eig((P+P.')/2)));
end
fprintf('min over [0,Tbar] of lambda_min(P(tau)) = %.4e\n', min(lam));

% ---- figure: timer-dependent matrix inequality lambda_min(P(tau)) >= 0 ----
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, tau, lam, '-','LineWidth',2.2,'Color',gc.order(1,:));
yline(ax, 0, '--','Color','k','LineWidth',1.5);
xlabel(ax,'timer $\tau$','Interpreter','latex');
ylabel(ax,'min eigenvalue of $P(\tau)$','Interpreter','latex');
lg = legend(ax, {'$\lambda_{\min}(P(\tau))$','$0$'}, ...
      'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off'); set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'comp_sos.pdf'));
