function ex_diss_sos()
%EX_DISS_SOS  Dissipativity by the SOS/LMI dual conditions, recovering the
%  bounded-real Riccati solution, for the examples of the Dissipativity
%  chapter, plus dissipativity-based state-feedback control.
%
%  ANALYSIS (bounded-real / mean-square L2-gain).  The certified gains are now
%  computed through the GENERIC library routines of the dissipativity family:
%  MVST_DISS_SOS (single mode) and MVST_DISS_SOS_MODES (coupled Markov-jump),
%  used exactly as MVST_STAB_SOS / MVST_STAB_SOS_MODES are used for stability.
%  These assemble the flow dissipation operator C_M(P) >= 0 of the chapter
%  (the dual DLMI th:IMP:diss:LMI, storage matrix S = -P >= 0) and bisect on
%  gamma; for the Markov-jump system the modes couple through the generator
%  sum_j Pi(i,j) P_j and one SOS program returns all {P_i}.  The recovered
%  gamma matches mvst_l2gain_ct_mjls (printed as a cross-check).
%
%  CONTROL (bounded-real state feedback).  This part is a SYNTHESIS LMI, not an
%  analysis condition, so it does not route through the analysis routine: for
%  dx=(A x+B u+Bw w)dt, z=C x+Dz u+Dw w, the change of variables W=S^{-1},
%  L=K W linearizes the synthesis; the dual LMI in (W,L,g) is convex, and
%  K=L W^{-1} achieves the minimized gamma, recovering the H-infinity (game)
%  Riccati controller.  It is kept here as the local SOS_BR_CONTROL helper.
%
%  Solver: SOSTOOLS + SeDuMi (shipped under matlab/).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();

names={}; gsos=[]; gric=[];

% ============ analysis: single bounded-real mode (generic MVST_DISS_SOS) ====
A = [-3 0;0 -3]; Bw=[0;1]; C=[0.5 1]; Dw=0.5;
r1 = mvst_diss_sos(struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                   struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
gR1 = mvst_l2gain_ct_mjls({A},Bw,{C},Dw,0);
names{end+1}='single'; gsos(end+1)=r1.gamma; gric(end+1)=gR1;
fprintf('[single BR] gamma: SOS dual=%.6f  Riccati/l2gain=%.6f\n', r1.gamma, gR1);

% ==== analysis: coupled Markov-jump (generic MVST_DISS_SOS_MODES) ====
A1=[1 2;0 1]; A2=[-3 0;0 -3]; C1=[1 0]; C2=[0.5 1]; Bwm=[0;1]; Dwm=0.5;
Pig=[-4 4;1 -1];
rm = mvst_diss_sos_modes(struct('A0',cat(3,A1,A2),'Bw',Bwm,'C',cat(3,C1,C2), ...
                   'Dw',Dwm,'Pig',Pig), ...
                   struct('type','l2gain'), struct('type','ctmjls'), 'gamma');
gRm = mvst_l2gain_ct_mjls({A1,A2},Bwm,{C1,C2},Dwm,Pig);
names{end+1}='mjls'; gsos(end+1)=rm.gamma; gric(end+1)=gRm;
fprintf('[MJLS  BR] gamma: SOS=%.6f  coupled l2gain=%.6f\n', rm.gamma, gRm);

% ============ control: bounded-real state feedback ============
% performance output z = [x1 ; u]  (state + control penalty), disturbance w
Ac=[0 1;-1 0.2]; Bc=[0;1]; Bwc=[0.5;0];
Cc=[1 0; 0 0]; Dz=[0;1]; Dwc=[0;0];
[K,W,gc2] = sos_br_control(Ac,Bc,Bwc,Cc,Dz,Dwc);
Acl = Ac+Bc*K; Ccl=Cc+Dz*K;
gCL = mvst_l2gain_ct_mjls({Acl},Bwc,{Ccl},Dwc,0);   % achieved closed-loop gain
fprintf('[BR control] synthesized gamma=%.6f  achieved closed-loop gain=%.6f  gain K=[% .4f % .4f]\n', ...
    sqrt(gc2), gCL, K(1), K(2));

% ============ summary + figure ============
fprintf('\n==== dissipativity SOS-dual recovery of the bounded-real gain ====\n');
for i=1:numel(names)
    fprintf('  %-8s  gamma_SOS=%.6f  gamma_Riccati=%.6f  |diff|=%.2e\n', ...
        names{i}, gsos(i), gric(i), abs(gsos(i)-gric(i)));
end

fig=figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax=axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
b=bar(ax,[gsos(:) gric(:)],'grouped');
b(1).FaceColor=gc.order(1,:); b(2).FaceColor=gc.order(2,:);
set(ax,'XTick',1:numel(names),'XTickLabel',names,'FontSize',10);
ylabel(ax,'bounded-real gain $\gamma$','Interpreter','latex');
lg=legend(ax,{'$\gamma_{\mathrm{SOS\ dual}}$','$\gamma_{\mathrm{Riccati}}$'}, ...
    'Interpreter','latex','Location','best'); set(lg,'Box','off');
title(ax,'Dissipativity: SOS dual recovers the bounded-real Riccati gain', ...
    'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(mvst_figdir(),'diss_sos.pdf'));
end


% =========================================================================
function [K,W,g] = sos_br_control(A,B,Bw,C,Dz,Dw)
%SOS_BR_CONTROL  Bounded-real state feedback: dual LMI in (W,L) by SOSTOOLS,
%  bisecting on gamma. K = L W^{-1} achieves the minimized gamma.
n=size(A,1); m=size(B,2); mw=size(Bw,2); nz=size(C,1);
lo=norm(Dw)+1e-4; hi=max(2*lo,1);
while ~feas(hi), hi=2*hi; if hi>1e4, break; end, end
for it=1:26
    mid=0.5*(lo+hi); ok=feas(mid);
    if ok, hi=mid; else, lo=mid; end
end
gsyn = 1.08*hi;              % back off from the singular optimum: well-conditioned K
[~,K,W] = feas(gsyn); g = gsyn^2;
    function [ok,K,W] = feas(gam)
        V=brvars(n,mw+nz); prog=sosprogram(V);
        [prog,W]=sospolymatrixvar(prog,monomials(V,0),[n n],'symmetric');
        [prog,L]=sospolymatrixvar(prog,monomials(V,0),[m n]);
        He=A*W+B*L; He=He+He.';
        M=[ He, Bw, (C*W+Dz*L).'; Bw.', -gam^2*eye(mw), Dw.'; C*W+Dz*L, Dw, -eye(nz) ];
        prog=sosineq(prog, V.'*(-M-1e-8*eye(n+mw+nz))*V);
        prog=sosineq(prog, V(1:n).'*(W-1e-8*eye(n))*V(1:n));
        o.solver='sedumi'; o.params.fid=0;
        [prog,info]=sossolve(prog,o);
        ok=(info.pinf==0 && info.dinf==0 && info.numerr<2);
        K=[];
        if ok, W=double(sosgetsol(prog,W)); W=(W+W.')/2; L=double(sosgetsol(prog,L)); K=L/W; end
    end
end

function v = brvars(n, mw)
%BRVARS  Column of pvar indeterminates [x(1:n); w(1:mw)].
persistent V key
k = sprintf('%d_%d', n, mw);
if isempty(V) || ~strcmp(key,k)
    pvar x1 x2 w1 w2 w3
    all = [x1;x2;w1;w2;w3];
    V = all(1:(n+mw)); key = k;
end
v = V;
end
