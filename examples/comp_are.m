%% COMP_ARE  Newton--Kleinman iteration on the deterministic LQR ARE.
%  Illustrates the quadratic convergence of the Newton--Kleinman scheme on
%  the deterministic algebraic Riccati equation of the single-mode
%  continuous-time instance (the multiplicative-noise terms Ac,Bc are
%  ignored here).  Starting from a stabilizing gain, each iteration solves a
%  Lyapunov equation (via the Kronecker formula, no Control System Toolbox)
%  and the ARE residual is plotted on a log scale.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
md = D.ct.modes{1};
A = md.A; B = md.B; Q = md.Q; R = md.R; n = size(A,1);

% ---- Newton--Kleinman: start from a stabilizing gain ----
K = [-2 -2];                                 % A+B*K is Hurwitz for D.ct
assert(all(real(eig(A+B*K)) < 0), 'comp_are:stab','Initial gain not stabilizing.');
maxit = 8; res = zeros(1,maxit); nit = 0;
for i = 1:maxit
    Acl = A + B*K;
    P = lyapsolve(Acl, Q + K.'*R*K);          % solves Acl'*P + P*Acl + (Q+K'RK) = 0
    G = A.'*P + P*A + Q - P*B*(R\(B.'*P));     % ARE residual G(P_i)
    nit = nit + 1; res(nit) = norm(G, 'fro');
    K = -(R\(B.'*P));                          % gain update
    if res(nit) < 1e-12; break; end
end
res = res(1:nit);
fprintf('Newton--Kleinman: %d iterations, final residual %.3e\n', nit, res(end));

% ---- figure: ARE residual versus iteration (log scale) ----
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, 1:nit, res, '-o','LineWidth',2.0,'MarkerSize',5,'Color',gc.order(1,:));
set(ax,'YScale','log');
xlabel(ax,'iteration $i$','Interpreter','latex');
ylabel(ax,'ARE residual $\|G(P_i)\|_F$','Interpreter','latex');
lg = legend(ax, {'Newton--Kleinman'}, ...
      'Interpreter','latex','Location','northeast'); set(lg,'Box','off');
set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'comp_are.pdf'));

% ========================================================================
function X = lyapsolve(Acl, Qc)
%LYAPSOLVE  Solve the Lyapunov equation Acl'*X + X*Acl + Qc = 0 (Kronecker).
n = size(Acl,1);
x = -(kron(eye(n),Acl.') + kron(Acl.',eye(n))) \ Qc(:);
X = reshape(x, n, n); X = (X+X.')/2;
end
