function ex_perf()
%EX_PERF  Covariance (Gramian) performance measures of a constant dwell-time
%  impulsive system, for the chapter "Further Performance Measures" (Ch. 11).
%
%  Constant dwell-time impulsive system, flow dx=(A x + Bw w)dt, jump x+=J x
%  at each period Tbar, output z=C x.  The (periodic) controllability Gramian
%  Pi(tau) is the state second moment under unit-intensity disturbance:
%     Pi_dot = A Pi + Pi A' + Bw Bw'   on (0,Tbar),   Pi(0+) = J Pi(Tbar-) J',
%  in steady state.  On it, the chapter's covariance measures read
%     H2^2            = (1/Tbar) int_0^Tbar tr(C Pi(tau) C') dtau   (avg power),
%     energy-to-peak  = sqrt( max_tau lambda_max(C Pi(tau) C') )    (gen. H2),
%     peak            = max_tau lambda_max(C Pi(tau) C')            (peak gain),
%  the average, the peak eigenvalue, and the peak power along Pi.
%
%  EVERYTHING is routed through the generic dispatcher MVST_PERF: the headline
%  measures come from the measure-by-name interface MVST_PERF(SYS,'h2'|'e2p'|
%  'peak'), and the multi-period output-power profile plus the white-noise
%  Monte-Carlo cross-check come from the same function's legacy Gramian batch
%  interface MVST_PERF(SYS,OUT,OPTS).  The two interfaces agree to machine
%  precision (see tests/tperf.m, group B1).

here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---- constant dwell-time impulsive system (dispatcher-style struct) ----
sys = struct('A',[-0.5 1; -1 -0.5], 'J',0.8*eye(2), 'Tbar',1.2, ...
             'Bw',[0; 1], 'C',[1 0.5]);
Tbar = sys.Tbar;

% ---- headline covariance measures through the MVST_PERF dispatcher ----
h2  = mvst_perf(sys, 'h2',   struct('npts',200));   % avg-power H2 cost
e2p = mvst_perf(sys, 'e2p',  struct('npts',200));   % energy-to-peak (gen. H2)
pk  = mvst_perf(sys, 'peak', struct('npts',200));   % hybrid peak gain (= e2p^2)
H2sq = h2.gamma^2;  e2psq = e2p.gamma^2;

% ---- multi-period profile + white-noise Monte-Carlo cross-check ----
%      (the legacy Gramian batch interface of the SAME function)
out  = struct('Bw',sys.Bw, 'C',sys.C);
sysB = struct('A',sys.A, 'J',sys.J, 'Tbar',Tbar);
opts = struct('npts',200, 'nprofper',3, 'profpts',120, 'mc',true);
res  = mvst_perf(sysB, out, opts);

fprintf('H2   = %.4f  (H2^2 = %.4f, average output power)\n', h2.gamma, H2sq);
fprintf('energy-to-peak (generalized H2) = %.4f\n', e2p.gamma);
fprintf('peak gain (= energy-to-peak^2)  = %.4f\n', pk.gamma);
fprintf('H2^2 cross-check: dispatcher = %.4f, batch = %.4f, Monte-Carlo = %.4f\n', ...
        H2sq, res.H2sq, res.mc.h2);

% ================= figure: instantaneous output power over 3 periods =====
Np = opts.nprofper; tt = res.profile.t; pw = res.profile.power;
figure('Color','w','Units','centimeters','Position',[2 2 12 9.5]);
ax = axes; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h1 = plot(ax, tt, pw, '-','LineWidth',1.8,'Color',gc.order(1,:));
h2l = yline(ax, H2sq, '--','LineWidth',1.6,'Color',gc.order(2,:));
h3 = yline(ax, e2psq, '-.','LineWidth',2.4,'Color',gc.order(3,:));
for c = 1:Np-1, xline(ax, c*Tbar, ':','Color',[.6 .6 .6]); end
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'output power $\mathrm{tr}\,C\Pi(t)C^\top$','Interpreter','latex');
% Legend entries are padded on the right (\quad), wrapped on two columns.
set(ax,'FontSize',12);
lg = legend(ax,[h1 h2l h3], ...
   {'instantaneous power\quad','$H_2^2$ (period average)\quad', ...
    'energy-to-peak$^2$ (peak $\lambda$)\quad'}, ...
   'Interpreter','latex','NumColumns',2,'Location','southoutside','FontSize',12);
set(lg,'Box','off');
mvst_savepdf(gcf, fullfile(figdir,'perf_measures.pdf'));

% ---- data macros ----
try
    w = mvst_texdata(fullfile(figdir,'data','perf.tex'));
    w.num('PerfHtwo',h2.gamma); w.num('PerfEtwoP',e2p.gamma);
    w.close();
catch, end
end
