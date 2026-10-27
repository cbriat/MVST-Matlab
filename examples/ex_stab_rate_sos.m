function ex_stab_rate_sos()
%EX_STAB_RATE_SOS  Guaranteed-rate study of the MASTER IMPULSIVE system of the
%  worked example ex:stab:imp (Stochastic systems with deterministic impulses),
%  run through the generic matrix-SOS analysis program MVST_STAB_SOS in the timer
%  (Stability chapter, range dwell-time Theorem th:stab:RDT and the hybrid-rate
%  Definition def:stab:gues).
%
%  System (analysis configuration of ex:stab:imp, state-only, u=w=0):
%     flow:  dx = A0 x dt + A1 x dW,   A0 = [0 1; 2 0.1] (a SADDLE: UNSTABLE),
%                                      A1 = 0.3 I,
%     jump:  x+ = J0 x + J1 x nu(k),   J0 = [0.7 0.3; 0 0.7], J1 = [0.2 0.1; -0.1 0.2]
%            (a deterministic reset plus one stochastic reset channel, CONTRACTING).
%  Second-moment flow operator  F(X) = A0 X + X A0' + A1 X A1',
%     jump operator  J(X) = J0 X J0' + J1 X J1'.
%
%  REGIME (flipped from the design-agnostic Levy/expansive-jump case): the flow
%  EXPANDS the second moment (abscissa muF>0) while the resets CONTRACT it
%  (lifted jump radius rhoJ<1).  Stability is therefore of the PERSISTENT-JUMPING
%  / MAXIMUM dwell-time kind: the loop is mean-square stable iff the resets fire
%  often enough, i.e. the dwell-time stays below a critical bound T*, and a
%  SHORTER admissible dwell buys a STRONGER guaranteed per-jump contraction.
%
%  STUDY.  The guaranteed hybrid rate is (0,rho): no flow decay is claimed
%  (the flow expands), and all decay is carried by the jump factor rho<1
%  (persistent-jumping side, .persist='jumping', .alpha=0).  For a maximum
%  admissible dwell-time Tbar the certificate P(tau) on [0,Tbar] must satisfy,
%  for EVERY closing length theta in (0,Tbar], the range dwell-time conditions of
%  Theorem th:stab:RDT with the pre-jump certificate scaled by rho; this is the
%  generic 'range' test with Tmin -> 0 and Tmax = Tbar.  The SMALLEST feasible
%  rho*(Tbar) (the strongest certified contraction) is found by BISECTION ON THE
%  SOS FEASIBILITY (no timer grid), and SWEPT over Tbar in (0,T*).
%
%  CROSS-CHECK.  The exact guaranteed per-jump contraction over dwells of length
%  at most Tbar is the worst one-period monodromy radius, attained at the longest
%  dwell (the flow expands monotonically):
%     rho_exact(Tbar) = max_{theta<=Tbar} radius( MJ * expm(LF*theta) )
%                     = radius( MJ * expm(LF*Tbar) ),
%  with LF = mvst_liftsym(A0,{A1}) the lifted flow generator and
%  MJ = mvst_liftsym(0,{J0,J1}) the lifted jump operator.  The critical dwell T*
%  solves rho_exact(T*) = 1.  SOS rho*(Tbar) and rho_exact(Tbar) must agree to
%  ~3e-2.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---------------- system data (analysis config of ex:stab:imp) ----------------
A0 = [0 1; 2 0.1];  A1 = 0.3*eye(2);
J0 = [0.7 0.3; 0 0.7];  J1 = [0.2 0.1; -0.1 0.2];
sys = struct('A0',A0,'Ac',A1,'J',cat(3,J0,J1));
tol = 1e-6;

LF = mvst_liftsym(A0,{A1});        muF  = max(real(eig(LF)));   % flow abscissa
MJ = mvst_liftsym(zeros(2),{J0,J1}); rhoJ = max(abs(eig(MJ)));  % jump radius
mono = @(T) max(abs(eig(MJ*expm(LF*T))));                       % monodromy radius

fprintf('\n==== guaranteed-rate study, master impulsive system of ex:stab:imp ====\n');
fprintf('flow abscissa muF=%.4f (UNSTABLE), jump radius rhoJ=%.4f (CONTRACTING)\n',muF,rhoJ);
assert(muF>0 && rhoJ<1, 'expected unstable flow with contracting jump');

% ---------------- critical maximum dwell-time T* (plain stability) ----------
Tstar = fzero(@(T) mono(T)-1, [0.02 0.3]);
fprintf('critical maximum dwell-time T* = %.4f (monodromy radius = 1)\n', Tstar);

% ---------------- guaranteed contraction rho*(Tbar) --------------------------
Tsweep = linspace(0.02, 0.95*Tstar, 10);
rhoStar  = nan(size(Tsweep));
rhoExact = arrayfun(mono, Tsweep);          % exact = monodromy at the max dwell
for i = 1:numel(Tsweep)
    rhoStar(i) = min_rho(sys, Tsweep(i), tol);
    fprintf('  Tbar=%.4f : SOS rho*=%.4f  exact=%.4f  (abs. dev. %.2e)\n', ...
        Tsweep(i), rhoStar(i), rhoExact(i), abs(rhoStar(i)-rhoExact(i)));
    assert(abs(rhoStar(i)-rhoExact(i)) < 3e-2, ...
        'SOS guaranteed contraction disagrees with the exact monodromy beyond 3e-2');
end

% ---------------- highlighted operating point ----------------
iH = 5; TbarH = Tsweep(iH); rhoH = rhoStar(iH);
fprintf('  highlighted point: (Tbar=%.4f, rho*=%.4f)\n', TbarH, rhoH);

% ---------------- figure ----------------
figc = figure('Color','w','Units','centimeters','Position',[2 2 13 9]);
axc = axes(figc); hold(axc,'on'); box(axc,'on'); grid(axc,'on');
Tfine = linspace(Tsweep(1), Tstar, 200);
plot(axc, Tfine, arrayfun(mono,Tfine), '-','LineWidth',1.6,'Color',gc.order(2,:));
plot(axc, Tsweep, rhoStar, 'o','MarkerSize',5,'MarkerFaceColor',gc.order(1,:), ...
     'MarkerEdgeColor','k','LineStyle','none');
plot(axc, TbarH, rhoH, 'o','MarkerSize',9,'MarkerFaceColor',gc.order(3,:), ...
     'MarkerEdgeColor','k','HandleVisibility','off');
yline(axc, 1, ':','Color',[.4 .4 .4],'HandleVisibility','off');
xline(axc, Tstar, '--','Color',[.4 .4 .4],'HandleVisibility','off');
xlabel(axc,'maximum admissible dwell-time $\bar T$','Interpreter','latex');
ylabel(axc,'guaranteed per-jump contraction $\rho^\star(\bar T)$','Interpreter','latex');
title(axc,'Guaranteed rate versus maximum dwell-time', ...
      'Interpreter','latex','FontWeight','normal');
lg = legend(axc, {'exact monodromy radius','matrix-SOS certificate $\rho^\star$'}, ...
            'Interpreter','latex','Location','northwest'); set(lg,'Box','off');
set(axc,'FontSize',11);
mvst_savepdf(figc, fullfile(figdir,'stab_rate_sos_combined.pdf'));

% ---------------- data macros ----------------
w = mvst_texdata(fullfile(figdir,'data','stab_rate_sos.tex'));
w.num('StabRateSosMuF',    muF);
w.num('StabRateSosRhoJ',   rhoJ);
w.num('StabRateSosTstar',  Tstar);
w.num('StabRateSosTbar',   TbarH);
w.num('StabRateSosRho',    rhoH);
w.close();
fprintf('   wrote stab_rate_sos_combined.pdf and data/stab_rate_sos.tex\n');
end

% =========================================================================
function rho = min_rho(sys, Tbar, tol)
%MIN_RHO  Smallest jump factor rho for which the persistent-jumping,
%  maximum-dwell (range [1e-3,Tbar]) analysis is SOS-feasible at that Tbar
%  (bisection on the slack sign; slack increases with rho).
if slackval(sys,Tbar,0.999) <= tol, rho = NaN; return; end   % beyond T* (unstable)
lo = 0.30; hi = 0.999;
if slackval(sys,Tbar,lo) > tol, rho = lo; return; end
for it = 1:22
    mid = 0.5*(lo+hi);
    if slackval(sys,Tbar,mid) > tol, hi = mid; else, lo = mid; end
end
rho = hi;
end

% =========================================================================
function s = slackval(sys, Tbar, rho)
%SLACKVAL  Optimal slack of the persistent-jumping, maximum-dwell analysis
%  (generic matrix-SOS routine MVST_STAB_SOS, range dwell-time [1e-3,Tbar],
%  degree-2 certificate P(theta), flow shift alpha=0, jump factor rho, flow
%  only required non-expansive).  Feasible iff s > 0.
res = mvst_stab_sos(sys, struct('type','range','Tmin',1e-3,'Tmax',Tbar, ...
                        'deg',2,'alpha',0,'rho',rho,'persist','jumping'));
s = res.slack;
end
