function ex_oc_sos_poisson()
%EX_OC_SOS_POISSON  Optimal control by the SOS/LMI dual conditions for the
%  GBM + Poissonian-jump system class, recovering the Poisson generalized
%  Riccati (I-GRE) solution of the Linear-Cost Optimal Control chapter.
%
%  For a controlled jump-diffusion mode
%     dx = (Ax+Bu)dt + sum_l (Cc_l x + Dc_l u) dW_l
%            + sum_r (Jc_r - I) x  dN_r,   E[dN_r] = lam_r dt,
%  with cost [x;u]'[Q 0;0 R][x;u], the optimal value is <P*,X0> with P* the
%  stabilizing solution of the Poisson I-GRE (mvst_care_ct_poisson):
%     G(P) = He[A'P] + sum_l Cc_l'P Cc_l
%            + sum_r lam_r (Jc_r'P Jc_r - P) + Q - S'Theta^{-1}S = 0,
%     Theta = R + sum_l Dc_l'P Dc_l,   S = B'P + sum_l Dc_l'P Cc_l.
%  The Poisson terms lam_r (Jc_r'P Jc_r - P) enter the accumulation block
%  exactly as the multiplicative-noise terms Cc_l'P Cc_l do; they are LINEAR
%  in P, so the storage LMI stays affine in the decision variable.
%
%  By a Schur complement G(P)>=0 <=> M(P)=[acc S'; S Theta]>=0, and P* is the
%  MAXIMAL P with M(P)>=0. The dual (storage-function) LMI is therefore
%     max <P,X0>   s.t.   M(P) >= 0,
%  posed to SOSTOOLS as the SOS condition [x;u]'M(P)[x;u] in Sigma (pvar mode).
%  Its optimum recovers P* and the gain K*=-Theta^{-1}S to solver accuracy;
%  complementary slackness M(P*)[I;K*]=0 holds exactly. Solver: SOSTOOLS +
%  SeDuMi (shipped under matlab/).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();

% ---- a small stabilizable GBM + Poisson instance (n=2, m=1) ----
D = poisson_data();

% ---- Riccati benchmark ----
[Pr, Kr, infr] = mvst_care_ct_poisson(D.A,D.B,D.Cc,D.Dc,D.Jc,D.lam,D.Q,D.R);

% ---- SOS storage program ----
[Ps, Ks, cost] = sos_poisson(D);

errP = norm(Ps-Pr,'fro'); errK = norm(Ks-Kr,'fro');
[Ms] = kyp_poisson(Ps, D);
slack = norm(Ms*[eye(D.n); Ks],'fro');

fprintf('\n==== Poisson optimal-control SOS-dual recovery of the Riccati solution ====\n');
fprintf('  Riccati residual ||G(Pr)||=%.2e  (horizon %.0f)\n', infr.residual, infr.horizon);
fprintf('  ||P_sos-P_ric|| = %.3e\n', errP);
fprintf('  ||K_sos-K_ric|| = %.3e\n', errK);
fprintf('  cost <P,X0>: SOS primal=%.6f  dual(moment)=%.6f  compl-slack res=%.2e\n', ...
    cost.primal, cost.dual, slack);
fprintf('  Riccati value <Pr,X0> = %.6f\n', trace(D.X0*Pr));

% ================== figure ==================
% Compare the SOS-dual solution against the Poisson Riccati solution entry
% by entry: the bars coincide (recovery error ~1e-7, annotated).
labels = {'$P_{11}$','$P_{12}$','$P_{22}$','$K_1$','$K_2$'};
vs = [Ps(1,1) Ps(1,2) Ps(2,2) Ks(1) Ks(2)];
vr = [Pr(1,1) Pr(1,2) Pr(2,2) Kr(1) Kr(2)];
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
b = bar(ax, [vs; vr].', 'grouped');
b(1).FaceColor = gc.order(1,:); b(2).FaceColor = gc.order(2,:);
set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'FontSize',10);
ylabel(ax,'value','Interpreter','latex');
lg = legend(ax,{'SOS dual','Poisson Riccati'}, ...
    'Interpreter','latex','Location','best'); set(lg,'Box','off');
title(ax, sprintf('Poisson optimal control: SOS dual recovers the Riccati solution (error $%.0e$)', ...
    max(errP,errK)),'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(mvst_figdir(),'oc_sos_poisson.pdf'));
end

% =========================================================================
function [P,K,cost] = sos_poisson(D)
%SOS_POISSON  Storage LMI  max <P,X0> s.t. M(P)>=0  by SOSTOOLS (pvar mode).
n=D.n; m=D.m; V = sos_vars(n,m);
prog = sosprogram(V);
[prog,P] = sospolymatrixvar(prog, monomials(V,0), [n n], 'symmetric');
M = kyp_poisson(P, D);
prog = sosineq(prog, V.'*M*V);
prog = sossetobj(prog, -trace(D.X0*P));
o.solver='sedumi'; o.params.fid=0;
[prog,~] = sossolve(prog,o);
P = double(sosgetsol(prog,P)); P=(P+P.')/2;
[~,K] = kyp_poisson(P, D);
cost = dual_cost(P, K, D);
end

% =========================================================================
function [M,K] = kyp_poisson(P, D)
%KYP_POISSON  M(P) = [acc S'; S Theta] for the GBM + Poisson I-GRE.
%  acc = He[A'P] + Q + sum_l Cc_l'P Cc_l + sum_r lam_r (Jc_r'P Jc_r - P).
n=D.n;
Theta = D.R; S = D.B.'*P; acc = D.A.'*P + P*D.A + D.Q;
for l=1:numel(D.Cc)
    Cc=D.Cc{l}; Dc=D.Dc{l};
    Theta = Theta + Dc.'*P*Dc; S = S + Dc.'*P*Cc; acc = acc + Cc.'*P*Cc;
end
for r=1:numel(D.Jc)
    Jc=D.Jc{r}; acc = acc + D.lam(r)*(Jc.'*P*Jc - P);
end
Theta=(Theta+Theta.')/2;
M=[acc S.'; S Theta]; M=(M+M.')/2;
K=[];
if nargout>=2, K = -(Theta\S); end   % only for numeric P (recovery), not pvar
end

% =========================================================================
function cost = dual_cost(P, K, D)
%DUAL_COST  Primal <P,X0>, and the dual moment cost trace([Q 0;0 R] Z*) with
%  Z* = [I;K] Xocc [I;K]' supported on u=Kx (complementary slackness). Xocc is
%  the mean-square occupation second moment: the closed-loop generator of the
%  jump-diffusion applied to the second moment, including the Poisson jumps.
n=D.n;
Acl = D.A + D.B*K;
Lyap = kron(eye(n),Acl) + kron(Acl,eye(n));        % vec(Acl X + X Acl')
for l=1:numel(D.Cc)
    Ccl = D.Cc{l} + D.Dc{l}*K;
    Lyap = Lyap + kron(Ccl,Ccl);                   % + sum vec(Ccl X Ccl')
end
for r=1:numel(D.Jc)
    Jc = D.Jc{r};
    Lyap = Lyap + D.lam(r)*(kron(Jc,Jc) - eye(n*n));   % Poisson jump moment
end
xocc = -Lyap\D.X0(:); Xocc = reshape(xocc,n,n); Xocc=(Xocc+Xocc.')/2;
QSR  = blkdiag(D.Q, D.R);
Z    = [eye(n); K]*Xocc*[eye(n); K].';
cost.primal = trace(D.X0*P);
cost.dual   = trace(QSR*Z);
end

% =========================================================================
function D = poisson_data()
%POISSON_DATA  A small mean-square stabilizable GBM + Poisson LQ instance.
D.A  = [ 0 1 ; -1 0.2 ];         % unstable drift
D.B  = [ 0 ; 1 ];
D.Cc = { 0.30*[0 0;1 0] };        % one state-multiplicative (GBM) channel
D.Dc = { [0;0.20] };              % input-multiplicative channel
D.Jc = { [0.9 0.1; 0.0 0.8] };    % one Poisson jump map
D.lam = 0.5;                      % Poisson intensity
D.Q  = diag([2.0 1.0]);
D.R  = 1.0;
D.X0 = eye(2);
D.n = size(D.A,1); D.m = size(D.B,2);
end

% =========================================================================
function v = sos_vars(n,m)
%SOS_VARS  Column of pvar indeterminates [x(1:n); u(1:m)] (n=2,m=1 here).
persistent V
if isempty(V)
    pvar x1 x2 u
    V = [x1;x2;u];
end
v = V;
end
