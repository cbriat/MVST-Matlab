function ex_dof_dilation()
%EX_DOF_DILATION  Conservatism of the scalar restriction of the lossless
%  dilation, for the multiplicative-noise output-feedback synthesis of the
%  chapter "Dynamic Output Feedback".
%
%  The exact synthesis carries the rational pivot Y^{-1}; the convex
%  sufficient bracket replaces it by the affine tangent 2*mu*I - mu^2*Y,
%  which satisfies 2*mu*I - mu^2*Y <= Y^{-1}. In the eigenbasis of Y the
%  pivot gap in the direction of eigenvalue y is
%       g(mu,y) = 1/y - (2*mu - mu^2*y) = (mu*y - 1)^2 / y >= 0,
%  vanishing exactly at mu = 1/y. Hence the scalar restriction is LOSSLESS
%  when Y is a scalar multiple of the identity, and its conservatism is
%  governed by the conditioning of Y. The full matrix-slack dilation
%  (Proposition prop:dof:dilation) removes the gap entirely.
%
%  Panel 1: the best-over-mu pivot gap lambda_max(Y^{-1} - (2 mu I - mu^2 Y))
%           versus the condition number of Y (worst rank-one certificate).
%  Panel 2: the gap versus mu for one ill-conditioned Y, showing the scalar
%           search and the residual.
%
%  Both panels are solver-free (eigenvalues only): they isolate the single
%  source of conservatism of the convex bracket, namely the tangent
%  under-estimate of Y^{-1}. The timer-dependent DLMIs are SOLVED
%  separately with SOSTOOLS (see ex_dof_sos.m).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();

% ---------- Panel 1: conservatism vs conditioning ----------
kappas = linspace(1, 40, 60);
bestgap = zeros(size(kappas));
for i = 1:numel(kappas)
    Y = diag([1, kappas(i)]);                 % eigenvalues 1 and kappa
    bestgap(i) = bestMuGap(Y);
end

fig1 = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
ax = axes(fig1); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, kappas, bestgap, '-', 'LineWidth',2.0, 'Color',gc.order(1,:));
xlabel(ax,'condition number of the certificate $Y$','Interpreter','latex');
ylabel(ax,'best-$\mu$ pivot gap $\lambda_{\max}\!\big(Y^{-1}-(2\mu I-\mu^2Y)\big)$', ...
    'Interpreter','latex');
set(ax,'FontSize',11);
mvst_savepdf(fig1, fullfile(mvst_figdir(),'dof_dilation_kappa.pdf'));

% ---------- Panel 2: gap vs mu for one ill-conditioned Y ----------
% The dilated LMI is feasible only while 2 mu I - mu^2 Y > 0, i.e.
% mu < 2/lambda_max(Y): beyond that wall the (2,2) block loses definiteness
% and the certificate is invalid. Feasibility therefore pins mu near the
% large-eigenvalue tangent 1/lambda_max, far from the 1/lambda_min tangent
% that the small-eigenvalue direction of Y^{-1} would need.
Y = diag([1, 12]); lmax = max(eig(Y));
muWall = 2/lmax;
mus = linspace(0.01, 0.45, 300);
gap = arrayfun(@(mu) max(eig(inv(Y) - (2*mu*eye(2) - mu^2*Y))), mus);
feas = mus < muWall;
% best FEASIBLE mu (constrained minimiser sits at the wall)
[gEdge, iE] = min(gap(feas)); musF = mus(feas);

fig2 = figure('Color','w','Units','centimeters','Position',[2 2 12 7]);
ax = axes(fig2); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, mus(feas),  gap(feas),  '-', 'LineWidth',2.0, 'Color',gc.order(1,:));
plot(ax, mus(~feas), gap(~feas), '--','LineWidth',1.6, 'Color',[.6 .6 .6]);
xline(ax, muWall, ':', 'LineWidth',1.4, 'Color',[.5 .5 .5]);
plot(ax, musF(iE), gEdge, 'o', 'MarkerSize',6, ...
     'MarkerFaceColor',gc.order(1,:), 'MarkerEdgeColor','none');
plot(ax, 1/lmax, 0, 's', 'MarkerSize',6, ...
     'MarkerFaceColor',gc.order(2,:), 'MarkerEdgeColor','none');
text(ax, muWall, 0.15, ' feasibility wall $\mu=2/\lambda_{\max}$', ...
     'Interpreter','latex','FontSize',9,'Color',[.5 .5 .5], ...
     'HorizontalAlignment','left','VerticalAlignment','bottom');
xlabel(ax,'scalar parameter $\mu$','Interpreter','latex');
ylabel(ax,'pivot gap $\lambda_{\max}\!\big(Y^{-1}-(2\mu I-\mu^2Y)\big)$', ...
    'Interpreter','latex');
title(ax,'$Y=\mathrm{diag}(1,12)$','Interpreter','latex','FontWeight','normal');
set(ax,'FontSize',11);
mvst_savepdf(fig2, fullfile(mvst_figdir(),'dof_dilation_muscan.pdf'));

% ---------- console summary ----------
gTen   = bestMuGap(diag([1,10]));
gForty = bestMuGap(diag([1,40]));
fprintf('best-mu pivot gap: kappa=10 -> %.4f,  kappa=40 -> %.4f\n', gTen, gForty);

% ---------- data macros (names must be letters only) ----------
try
    w = mvst_texdata(fullfile(mvst_figdir(),'data','dof_dilation.tex'));
    w.num('DofGapKappaTen',   gTen);
    w.num('DofGapKappaForty', gForty);
    w.close();
catch err
    fprintf('data-macro write skipped: %s\n', err.message);
end
end

% -------------------------------------------------------------------------
function g = bestMuGap(Y)
%BESTMUGAP  min over mu>0 of lambda_max( Y^{-1} - (2 mu I - mu^2 Y) ).
Yi = inv(Y); n = size(Y,1);
f = @(mu) max(eig(Yi - (2*mu*eye(n) - mu^2*Y)));
mus = linspace(1e-3, 2/max(eig(Y)), 400);
g = min(arrayfun(f, mus));
end
