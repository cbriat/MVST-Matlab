%% COMP_RIC_COMBINED  Consolidated Riccati-solver convergence figure.
%  Native, vertically-stacked replacement for the former two-panel minipage
%  tiling (fig:comp:are).  Reproduces, in one figure with readable fonts, the
%  Frobenius residual versus iteration on a logarithmic scale for the two
%  numerical Riccati schemes:
%    row 1: Newton--Kleinman on the deterministic infinite-horizon ARE
%           (single-mode CT instance; verbatim core of COMP_ARE),
%    row 2: relaxation (Gauss--Seidel) on the coupled Markov-jump ARE
%           (two-mode instance; verbatim core of COMP_MJLS).
%  This script produces the figure ONLY: it writes no data macros (COMP_ARE
%  and COMP_MJLS write none either, so the golden data is unaffected).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
D = ex_data();

% ---- row 1: Newton--Kleinman on the deterministic ARE (COMP_ARE core) ----
md = D.ct.modes{1};
A = md.A; B = md.B; Qc = md.Q; R = md.R;
K = [-2 -2];                                 % A+B*K Hurwitz for D.ct
assert(all(real(eig(A+B*K)) < 0), 'comp:stab','Initial gain not stabilizing.');
maxit = 8; resNK = zeros(1,maxit); nit = 0;
for i = 1:maxit
    Acl = A + B*K;
    P = lyapsolve(Acl, Qc + K.'*R*K);
    G = A.'*P + P*A + Qc - P*B*(R\(B.'*P));
    nit = nit + 1; resNK(nit) = norm(G,'fro');
    K = -(R\(B.'*P));
    if resNK(nit) < 1e-12, break; end
end
resNK = resNK(1:nit);

% ---- row 2: Gauss--Seidel on the coupled Markov-jump ARE (COMP_MJLS core) -
sys = D.mjls; M = numel(sys.modes); n = sys.modes{1}.n; Pi = sys.Pi;
nsweep = 8; Pc = repmat({zeros(n)},1,M); resGS = zeros(nsweep,M);
for s = 1:nsweep
    for i = 1:M
        mdi = sys.modes{i};
        mdi.A = sys.modes{i}.A + 0.5*Pi(i,i)*eye(n);
        Qadd = zeros(n);
        for j = 1:M, if j ~= i, Qadd = Qadd + Pi(i,j)*Pc{j}; end, end
        mdi.Q = sys.modes{i}.Q + (Qadd+Qadd.')/2;
        Pc{i} = mvst_care_ct(mdi, struct());
    end
    for i = 1:M
        Cpl = zeros(n);
        for j = 1:M, Cpl = Cpl + Pi(i,j)*Pc{j}; end
        resGS(s,i) = norm(mvst_ric_ct(Pc{i}, sys.modes{i}, Cpl),'fro') + eps;
    end
end

% ---- combined figure -----------------------------------------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 15]);
tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on'); set(ax1,'YScale','log');
plot(ax1, 1:nit, resNK, '-o','LineWidth',2.0,'MarkerSize',5,'Color',gc.order(1,:));
ylabel(ax1,'ARE residual $\|G(P_i)\|_F$','Interpreter','latex');
title(ax1,'(a) Newton--Kleinman (deterministic ARE)','Interpreter','latex','FontWeight','normal');
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on'); set(ax2,'YScale','log');
h(1) = plot(ax2, 1:nsweep, resGS(:,1), '-o','LineWidth',1.8,'MarkerSize',4,'Color',gc.order(1,:));
h(2) = plot(ax2, 1:nsweep, resGS(:,2), '-s','LineWidth',1.8,'MarkerSize',4,'Color',gc.order(2,:));
xlabel(ax2,'iteration / sweep','Interpreter','latex');
ylabel(ax2,'coupled ARE residual $\|G_i(P_i)\|_F$','Interpreter','latex');
title(ax2,'(b) Gauss--Seidel relaxation (coupled Markov-jump ARE)','Interpreter','latex','FontWeight','normal');
lg = legend(ax2, h, {'mode 1','mode 2'}, 'Interpreter','latex','Location','northeast');
set(lg,'Box','off'); set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'comp_ric_combined.pdf'));

% ========================================================================
function X = lyapsolve(Acl, Qc)
%LYAPSOLVE  Solve Acl'*X + X*Acl + Qc = 0 (Kronecker; verbatim from COMP_ARE).
n = size(Acl,1);
x = -(kron(eye(n),Acl.') + kron(Acl.',eye(n))) \ Qc(:);
X = reshape(x, n, n); X = (X+X.')/2;
end
