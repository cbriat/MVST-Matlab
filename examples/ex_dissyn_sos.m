function ex_dissyn_sos()
%EX_DISSYN_SOS  Dissipativity-based state-feedback SYNTHESIS by the generic
%  SOS/LMI convexified conditions of the Dissipative-Synthesis chapter
%  (chapter "Dissipativity-Based Synthesis"), on the running system classes of the Game
%  chapter for direct comparison.
%
%  For each class the bounded-real (L2-gain) state-feedback gain is DESIGNED
%  through the GENERIC library routine MVST_DISSYN_SOS, used exactly as
%  MVST_GAME_SOS / MVST_DISS_SOS / MVST_OC_SOS are used in their chapters: the
%  change of variables W = S^{-1}, Y = K W linearizes the closed-loop
%  dissipation DLMIs (Section sec:disssyn:convex), and the routine bisects on
%  the certified bounded-real level gamma while a post-solve gate re-analyses
%  the recovered closed loop.  The certified synthesis level gamma_SOS is
%  compared with the EXACT state-feedback H-infinity level of the Game chapter,
%  computed by the same game-Riccati backends MVST_GAME_GAMMA_CT /
%  _CT_POISSON used there (state/control penalty Q = I_2, R_u = 1, so the
%  bounded-real synthesis and the zero-sum game share one optimal level).
%
%  As documented in the routine's header and validated by its self-test, the
%  constant-storage relaxation is EXACT for the continuous-time cases (the
%  certified level equals the game-Riccati level), and is a certified UPPER
%  BOUND for the timer-dependent dwell-time certificates (sampled-data,
%  switched), where the reference is the underlying jump-free / worst-mode
%  continuous-time level: the gap is the price of the causal timer certificate,
%  mirroring the losslessness discussion of the chapter.  This is the synthesis
%  analogue of EX_GAME_SOS / EX_DISS_SOS / EX_OC_SOS.
%
%  Solver: SOSTOOLS + SeDuMi (shipped under matlab/).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();

% Common state/control penalty z = [x1 ; x2 ; u]  (Q = I_2, R_u = 1), matching
% the Game chapter cost x^*Qx + u^*R_u u.
C  = [eye(2); zeros(1,2)];   Du = [0;0;1];
Bu = [0;1];  Q = eye(2);  Ru = 1;

% w-actuated diffusion channel inherited from the dissipativity examples:
% a shared disturbance-multiplicative column D_{w,1}=0.1[1;0]
% on the continuous-time, sampled-data, and switched classes.  The synthesis
% routine folds it into the flow block (field .Dwc) and the game references
% carry it through mvst_game_gamma_ct's Dwc argument, so the memoryless levels
% still coincide.  The Poissonian class keeps its master's state / jump channels.
Dwc = 0.1*[1;0];

names={}; gsos=[]; gexa=[]; kind={};

% ===== (1) continuous-time LTI (arbitrary; relaxation is EXACT) ===========
A=[0 1;-1 -1]; Bw=[1;0];
r1 = mvst_dissyn_sos(struct('A0',A,'B0',Bu,'Bw',Bw,'C',C,'Du',Du,'Dwc',Dwc), ...
                     struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
% The game-Riccati reference levels use the LONG-HORIZON criterion of the game
% library (horizon 1000, relative bracket 1e-6, uniform Isaacs margin 1e-6 gamma^2).
g1 = mvst_game_gamma_ct(A,Bu,Bw,Q,Ru,1000,{},{Dwc},1e-6);
names{end+1}='ct'; gsos(end+1)=r1.gammastar; gexa(end+1)=g1; kind{end+1}='exact';
fprintf('[ct   syn] gamma: SOS cert=%.6f  game-Riccati=%.6f  (exact)\n', r1.gammastar, g1);

% ===== (2) sampled-data (periodic dwell-time; SOS = certified u.bound) =====
Tbar=0.5;
r2 = mvst_dissyn_sos(struct('A0',A,'B0',Bu,'Bw',Bw,'C',C,'Du',Du,'J',eye(2),'Dwc',Dwc), ...
                     struct('type','l2gain'), ...
                     struct('type','periodic','Tbar',Tbar,'deg',2), 'gamma');
g2 = mvst_game_gamma_ct(A,Bu,Bw,Q,Ru,1000,{},{Dwc},1e-6);  % jump-free CT reference
names{end+1}='sd'; gsos(end+1)=r2.gammastar; gexa(end+1)=g2; kind{end+1}='u.bound';
fprintf('[sd   syn] gamma: SOS cert=%.6f  CT reference=%.6f  (upper bound)\n', r2.gammastar, g2);

% ===== (3) switched (minimum dwell-time; worst mode; SOS u.bound) =========
A1=[0 1;-2 -1]; A2=[0 1;-1 -2]; Td=0.5;
gS=[];
for Am = {A1,A2}
    ri = mvst_dissyn_sos(struct('A0',Am{1},'B0',Bu,'Bw',Bw,'C',C,'Du',Du,'Dwc',Dwc), ...
                         struct('type','l2gain'), ...
                         struct('type','minimum','Tbar',Td,'deg',2), 'gamma');
    gS(end+1)=ri.gammastar; %#ok<AGROW>
end
% g3 = max(mvst_game_gamma_ct(A1,Bu,Bw,Q,Ru,40,{},{Dwc}), mvst_game_gamma_ct(A2,Bu,Bw,Q,Ru,40,{},{Dwc}));
g3 = max(mvst_game_gamma_ct(A1,Bu,Bw,Q,Ru,1000,{},{Dwc},1e-6), ...
         mvst_game_gamma_ct(A2,Bu,Bw,Q,Ru,1000,{},{Dwc},1e-6));
names{end+1}='switched'; gsos(end+1)=max(gS); gexa(end+1)=g3; kind{end+1}='u.bound';
fprintf('[sw   syn] gamma: SOS cert(worst mode)=%.6f  worst-mode CT=%.6f  (upper bound)\n', ...
        max(gS), g3);

% ===== (4) Poissonian jump-diffusion (state-only reset; arbitrary) ========
% Second-moment flow of the compound-Poisson jump x+ = J x at rate lam plus the
% Brownian channel A_1 = C1: it enters MVST_DISSYN_SOS as the drift shift
% A0 -> A0 - (lam/2) I and the completely-positive channels {C1, sqrt(lam) J},
% the covariance model of Section sec:pre:cov:ctjd.
Ap=[0 1;2 -0.5]; C1=0.2*eye(2); J1=1.25*eye(2); lam=0.5;
A0e = Ap - (lam/2)*eye(2); Ace = cat(3, C1, sqrt(lam)*J1);
r4 = mvst_dissyn_sos(struct('A0',A0e,'B0',Bu,'Bw',Bw,'Ac',Ace,'C',C,'Du',Du), ...
                     struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
g4 = mvst_game_gamma_ct_poisson(Ap,Bu,Bw,{C1},{J1},lam,Q,Ru,1000,[],[],[],1e-6);
names{end+1}='poisson'; gsos(end+1)=r4.gammastar; gexa(end+1)=g4; kind{end+1}='exact';
fprintf('[poi  syn] gamma: SOS cert=%.6f  game-Riccati(state-only)=%.6f\n', r4.gammastar, g4);

% ===== summary =====
fprintf('\n==== dissipative-synthesis SOS: certified level vs exact game-Riccati level ====\n');
for i=1:numel(names)
    fprintf('  %-9s  gamma_SOS=%.6f  gamma_exact=%.6f  |diff|=%.2e  (%s)\n', ...
        names{i}, gsos(i), gexa(i), abs(gsos(i)-gexa(i)), kind{i});
end

% ===== figure: certified synthesis level vs exact game level per class =====
fig=figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax=axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
b=bar(ax,[gsos(:) gexa(:)],'grouped');
b(1).FaceColor=gc.order(1,:); b(2).FaceColor=gc.order(2,:);
set(ax,'XTick',1:numel(names),'XTickLabel',names,'FontSize',10);
ylabel(ax,'bounded-real level $\gamma^\star$','Interpreter','latex');
lg=legend(ax,{'$\gamma_{\mathrm{SOS}}$ (synthesis)','$\gamma_{\mathrm{Riccati}}$ (game)'}, ...
    'Interpreter','latex','Location','best'); set(lg,'Box','off');
title(ax,'Dissipative synthesis: SOS-certified level vs exact game-Riccati level', ...
    'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(mvst_figdir(),'dissyn_sos.pdf'));

% ===== data macros =====
w = mvst_texdata(fullfile(mvst_figdir(),'data','dissyn_sos.tex'));
w.num('DsynGsosCT',gsos(1));  w.num('DsynGexaCT',gexa(1));
w.num('DsynGsosSD',gsos(2));  w.num('DsynGexaSD',gexa(2));
w.num('DsynGsosSw',gsos(3));  w.num('DsynGexaSw',gexa(3));
w.num('DsynGsosPoi',gsos(4)); w.num('DsynGexaPoi',gexa(4));
w.close();
end
