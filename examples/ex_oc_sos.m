function ex_oc_sos()
%EX_OC_SOS  Optimal control by the SOS/LMI dual conditions, recovering the
%  Riccati (I-GRE) solution, for the constant-data examples of the
%  Linear-Cost Optimal Control chapter.
%
%  For a controlled diffusion mode  dx = (Ax+Bu)dt + sum_l(Ac_l x+Bc_l u)dW_l
%  with cost [x;u]'[Q S;S' R][x;u], the optimal value is <P*,X0> with P* the
%  stabilizing solution of the I-GRE (mvst_care_ct):
%     G(P) = He[A'P] + sum Ac_l'P Ac_l + Q - S'Theta^{-1}S = 0,
%     Theta = R + sum Bc_l'P Bc_l,   S = B'P + Scr' + sum Bc_l'P Ac_l.
%  By a Schur complement G(P)>=0 <=> M(P)=[acc S'; S Theta]>=0, and P* is the
%  MAXIMAL P with M(P)>=0. The DUAL (storage-function) LMI is therefore
%     max <P,X0>   s.t.   M(P) >= 0,
%  which SOSTOOLS solves as an SOS program in (x,u): [x;u]'M(P)[x;u] is SOS.
%  Its optimum recovers P* (and the gain K*=-Theta^{-1}S) exactly.
%
%  Dual conditions / complementary slackness: at the optimum M(P*)[I;K*]'=0,
%  so the dual (moment) measure is supported on the optimal law u=K*x, and
%  the dual cost trace([Q S;S' R] Z*) with Z*=[I;K*]Xocc[I;K*]' equals the
%  primal cost <P*,X0>. Both are printed against the Riccati value.
%
%  For the Markov-jump system the same conditions couple across modes through
%  the generator term sum_j pi_ij P_j, and one SOS program over {P_i} solves
%  all modes at once. Solver: SOSTOOLS + SeDuMi (shipped under matlab/).
%
%  The SOS engine now lives in the library: single-mode in MVST_OC_SOS, the
%  coupled Markov-jump case in MVST_OC_SOS_MODES (with the Gauss-Seidel coupled
%  Riccati benchmark returned under OPTS.riccati).  This script is a thin
%  caller that formats the chapter data and draws the recovery figure.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
D = ex_data();

names = {}; errP = []; errK = []; dgap = [];

% ================= single-mode continuous-time =================
md = D.ct.modes{1};
[Pr,Kr] = mvst_care_ct(md);
r1 = mvst_oc_sos(md);
Ps = r1.P; Ks = r1.K; cost = r1.cost;
[names,errP,errK,dgap] = logrow(names,errP,errK,dgap,'single', ...
    norm(Ps-Pr,'fro'), norm(Ks-Kr,'fro'), cost.primal-cost.dual);
report('single-mode', Pr,Ps,Kr,Ks,cost,mvst_ric_ct(Pr,md));

% ================= Markov-jump (coupled) =================
sys = D.mjls; N = numel(sys.modes);
rc = mvst_oc_sos_modes(sys, [], [], struct('riccati',true));
Prc = rc.Pric; Krc = rc.Kric; Psc = rc.P; Ksc = rc.K; costc = rc.cost;
eP = max(arrayfun(@(i) norm(Psc{i}-Prc{i},'fro'), 1:N));
eK = max(arrayfun(@(i) norm(Ksc{i}-Krc{i},'fro'), 1:N));
[names,errP,errK,dgap] = logrow(names,errP,errK,dgap,'mjls',eP,eK,costc.primal-costc.dual);
fprintf('\n[MJLS coupled] max_i ||Ps_i-Pr_i||=%.3e  max_i ||Ks_i-Kr_i||=%.3e\n', eP, eK);
fprintf('   cost sum_i pi0_i <P_i,X0>: Riccati=%.6f  SOS primal=%.6f  dual=%.6f\n', ...
    costc.riccati, costc.primal, costc.dual);

% ================= summary table + figure =================
fprintf('\n==== optimal-control SOS-dual recovery of the Riccati solution ====\n');
for i=1:numel(names)
    fprintf('  %-8s  ||P_sos-P_ric||=%.2e   ||K_sos-K_ric||=%.2e   primal-dual gap=%.2e\n', ...
        names{i}, errP(i), errK(i), dgap(i));
end

fig = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
b = bar(ax, [errP(:) errK(:)], 'grouped');
b(1).FaceColor = gc.order(1,:); b(2).FaceColor = gc.order(2,:);
set(ax,'YScale','log','XTick',1:numel(names),'XTickLabel',names,'FontSize',10);
ylabel(ax,'recovery error (Frobenius)','Interpreter','latex');
lg = legend(ax,{'$\|P_{\mathrm{SOS}}-P_{\mathrm{Ric}}\|$','$\|K_{\mathrm{SOS}}-K_{\mathrm{Ric}}\|$'}, ...
    'Interpreter','latex','Location','best'); set(lg,'Box','off');
title(ax,'Optimal control: SOS dual recovers the Riccati solution', ...
    'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(mvst_figdir(),'oc_sos.pdf'));
end

% =========================================================================
function report(tag, Pr,Ps,Kr,Ks,cost,Gr)
fprintf('\n[%s] ||Ps-Pr||=%.3e  ||Ks-Kr||=%.3e  ||G(Pr)||=%.2e\n', ...
    tag, norm(Ps-Pr,'fro'), norm(Ks-Kr,'fro'), norm(Gr,'fro'));
fprintf('   cost <P,X0>: SOS primal=%.6f  dual(moment)=%.6f  compl-slack res=%.2e\n', ...
    cost.primal, cost.dual, cost.slackres);
end

function [names,eP,eK,dg]=logrow(names,eP,eK,dg,nm,a,b,c)
names{end+1}=nm; eP(end+1)=a; eK(end+1)=b; dg(end+1)=c;
end
