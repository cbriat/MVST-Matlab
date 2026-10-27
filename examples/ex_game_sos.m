function ex_game_sos()
%EX_GAME_SOS  Zero-sum H-infinity GAME by the generic SOS/LMI dual conditions,
%  for the running system classes of the zero-sum games chapter.
%
%  The certified critical disturbance level gamma* is computed through the
%  GENERIC library routine MVST_GAME_SOS, used exactly as MVST_STAB_SOS /
%  MVST_DISS_SOS are used in their chapters.  MVST_GAME_SOS assembles the
%  bounded-real GAME synthesis DLMI (change of variables Y = P^{-1}, W = K Y)
%  and bisects on gamma.  As documented in that routine's header, the SOS
%  value is EXACT for the separated continuous-time LTI game (dwell type
%  'arbitrary') and a certified UPPER BOUND for the timer-dependent
%  (dwell-time) certificates; the exact per-class Riccati levels
%  MVST_GAME_GAMMA_CT / _DWELL / _POISSON are printed alongside as the ground
%  truth.  This is the game analogue of EX_DISS_SOS / EX_OC_SOS.
%
%  Solver: SOSTOOLS + SeDuMi (shipped under matlab/).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();

names={}; gsos=[]; gexa=[]; kind={};

% ===== (1) continuous-time separated game (arbitrary; SOS is EXACT) =====
A=[0 1;-1 -1]; Bu=[0;1]; Bw=[1;0]; Q=eye(2); Ru=1;
r  = mvst_game_sos(struct('A0',A,'Bu',Bu,'Bw',Bw), ...
                   struct('Q',Q,'Ru',Ru), struct('type','arbitrary'), 'gamma');
ge = mvst_game_gamma_ct(A,Bu,Bw,Q,Ru,40);
names{end+1}='ct'; gsos(end+1)=r.gammastar; gexa(end+1)=ge; kind{end+1}='exact';
fprintf('[ct   game] gamma: SOS=%.6f  Riccati=%.6f  (exact)\n', r.gammastar, ge);

% ===== (2) sampled-data game (periodic dwell-time; SOS upper bound) =====
Asd=[0 1;-1 -1]; Tbar=0.5;
r  = mvst_game_sos(struct('A0',Asd,'Bu',Bu,'Bw',Bw,'J',eye(2)), ...
                   struct('Q',Q,'Ru',Ru), ...
                   struct('type','periodic','Tbar',Tbar,'deg',2), 'gamma');
ge = mvst_game_gamma_ct(Asd,Bu,Bw,Q,Ru,40);            % jump-free reference
names{end+1}='sd'; gsos(end+1)=r.gammastar; gexa(end+1)=ge; kind{end+1}='u.bound';
fprintf('[sd   game] gamma: SOS=%.6f  CT reference=%.6f  (upper bound)\n', r.gammastar, ge);

% ===== (3) switched game (minimum dwell-time; SOS upper bound) =====
A1=[0 1;-2 -1]; A2=[0 1;-1 -2]; Td=0.5;
% generic minimum-dwell certificate on the worst-case mode (frozen tail)
r1 = mvst_game_sos(struct('A0',A1,'Bu',Bu,'Bw',Bw), ...
                   struct('Q',Q,'Ru',Ru), ...
                   struct('type','minimum','Tbar',Td,'deg',2), 'gamma');
r2 = mvst_game_sos(struct('A0',A2,'Bu',Bu,'Bw',Bw), ...
                   struct('Q',Q,'Ru',Ru), ...
                   struct('type','minimum','Tbar',Td,'deg',2), 'gamma');
gS = max(r1.gammastar, r2.gammastar);
% clean worst-mode CT reference (a lower bound on the switched level)
ge = max(mvst_game_gamma_ct(A1,Bu,Bw,Q,Ru,40), mvst_game_gamma_ct(A2,Bu,Bw,Q,Ru,40));
names{end+1}='switched'; gsos(end+1)=gS; gexa(end+1)=ge; kind{end+1}='u.bound';
fprintf('[sw   game] gamma: SOS(min-dwell,worst mode)=%.6f  worst-mode CT=%.6f\n', gS, ge);

% ===== (4) Poissonian-jump game (state-only reset; arbitrary) =====
Ap=[0 1;2 -0.5]; C1=0.2*eye(2); J1=1.25*eye(2); lam=0.5;
r  = mvst_game_sos(struct('A0',Ap,'Bu',Bu,'Bw',Bw, ...
                          'Ac',C1,'Jf',J1,'lamf',lam), ...
                   struct('Q',Q,'Ru',Ru), struct('type','arbitrary'), 'gamma');
ge = mvst_game_gamma_ct_poisson(Ap,Bu,Bw,{C1},{J1},lam,Q,Ru,30);   % state-only
names{end+1}='poisson'; gsos(end+1)=r.gammastar; gexa(end+1)=ge; kind{end+1}='u.bound';
fprintf('[poi  game] gamma: SOS=%.6f  Riccati(state-only)=%.6f\n', r.gammastar, ge);

% ===== summary =====
fprintf('\n==== game SOS dual: critical level vs exact Riccati level ====\n');
for i=1:numel(names)
    fprintf('  %-9s  gamma_SOS=%.6f  gamma_exact=%.6f  |diff|=%.2e  (%s)\n', ...
        names{i}, gsos(i), gexa(i), abs(gsos(i)-gexa(i)), kind{i});
end

% ===== figure: SOS vs exact critical level per class =====
fig=figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax=axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
b=bar(ax,[gsos(:) gexa(:)],'grouped');
b(1).FaceColor=gc.order(1,:); b(2).FaceColor=gc.order(2,:);
set(ax,'XTick',1:numel(names),'XTickLabel',names,'FontSize',10);
ylabel(ax,'critical level $\gamma^\star$','Interpreter','latex');
lg=legend(ax,{'$\gamma_{\mathrm{SOS}}$','$\gamma_{\mathrm{Riccati}}$'}, ...
    'Interpreter','latex','Location','best'); set(lg,'Box','off');
title(ax,'Game: SOS dual level vs exact game-Riccati level', ...
    'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(mvst_figdir(),'game_sos.pdf'));

% ===== data macros =====
w = mvst_texdata(fullfile(mvst_figdir(),'data','game_sos.tex'));
w.num('GsosCT',gsos(1));  w.num('GexaCT',gexa(1));
w.num('GsosSD',gsos(2));  w.num('GsosSw',gsos(3)); w.num('GexaSw',gexa(3));
w.num('GsosPoi',gsos(4)); w.num('GexaPoi',gexa(4));
w.close();
end
