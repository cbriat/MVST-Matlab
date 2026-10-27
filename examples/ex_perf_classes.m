function ex_perf_classes()
%EX_PERF_CLASSES  Induced performance (hybrid L2 x l2-gain and trace gain)
%  across the system classes of the monograph, for the chapter "Further
%  Performance Measures" (Ch. 11).
%
%  The covariance measures (H2, energy-to-peak, peak / intensity peak-to-peak)
%  are Gramian-characterized only for the continuous-time and the constant
%  dwell-time impulsive / sampled-data classes (see EX_PERF).  For the
%  multi-mode (Markov jump, switched) and the stochastic (Poissonian-jump)
%  classes the book states INDUCED-gain conditions instead, so the performance
%  measure available for every class is the hybrid L2 x l2-gain and its square,
%  the TRACE GAIN of th:perf:tracegain.  This script computes both for the five
%  classes, all through the generic dispatcher MVST_PERF:
%     gamma  = MVST_PERF(sys,'l2').gamma          (hybrid L2 x l2-gain)
%     tgain  = MVST_PERF(sys,'tracegain').gamma   (= gamma^2)
%  MVST_PERF routes 'l2'/'tracegain' to the exact per-class backends of
%  MVST_L2GAIN (Riccati / lifted-eigenvalue / shooting).  The
%  backends are called directly with the dissipativity-example settings (horizon 1000,
%  tolerance 1e-6), and the switched entry is the SOS minimum dwell-time bound
%  (SeDuMi/SOSTOOLS); MVST_PERF defaults are printed as a diagnostic.  The system instances are the disturbance instances of the four
%  dissipativity-chapter examples (Section sec:diss:examples), plus a plain
%  continuous-time reference, so the numbers are cross-consistent with that
%  chapter.
%
%  Produces perf_classes.pdf: a grouped bar chart of the L2-gain and the trace
%  gain per class, and writes the values as data macros to
%  figures/data/perf_classes.tex.

here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---- the five class instances (disturbance instances) -------------------
% Canonical performance variants of the Ch08 dissipativity examples: canonical
% triple Bw = [1;0], C = [1 0], Dw = 0.5, and the w-actuated diffusion channel
% D_{w,1} = 0.1[1;0] throughout (plus, for the Markov-jump instance, the
% state-multiplicative channel 0.2 I2 per mode).  All gains route through
% mvst_perf -> mvst_l2gain -> the exact backends, so they equal the Ch08 values.
Dwcm  = 0.1*[1;0];
mjAc  = zeros(2,2,1,2); mjAc(:,:,1,1) = 0.2*eye(2); mjAc(:,:,1,2) = 0.2*eye(2);
% continuous-time reference (single stable mode).  NOTE: the CT reference drift
% was not pinned by the family table; its output/feedthrough are made canonical.
ct  = struct('A',[-0.5 1; -1 -0.5], 'Bw',[1;0], 'C',[1 0], 'Dw',0.5, 'Dwc',Dwcm);
% Markov jump diffusion: mode 1 UNSTABLE, mode 2 stable; coupled gain finite
mj  = struct('A0',cat(3,[1 2;0 1],[-3 0;0 -3]), 'Bw',[1;0], ...
             'C',cat(3,[1 0],[0.5 1]), 'Dw',0.5, 'Pig',[-4 4;1 -1], ...
             'Ac',mjAc, 'Dwc',Dwcm);
% switched system under a minimum dwell-time Td (worst-case over schedules)
sw  = struct('A0',cat(3,[0 1;-2 -1],[0 1;-1 -2]), 'Bw',[1;0], ...
             'C',cat(3,[1 0],[1 0]), 'Dw',0.5, 'Td',0.5, 'Dwc',Dwcm);
% sampled-data system (hybrid sampled-output gain), period Tbar
sd  = struct('A',[0 1;-1 -1], 'Bw',[1;0], 'C',[1 0], 'Dw',0.5, 'Tbar',0.3, 'Dwc',Dwcm);
% Poissonian-jump diffusion
po  = struct('A',[0 1;-2 -1], 'Bw',[1;0], 'C',[1 0], 'Dw',0.5, ...
             'Ac',0.15*eye(2), 'Dwc',Dwcm, 'Jf',1.15*eye(2), ...
             'lamf',0.5, 'Bdf',[0;0.2]);

names = {'continuous-time','Markov jump','switched','sampled-data','Poissonian-jump'};
tags  = {'Ct','Mjls','Sw','Sd','Poi'};
syss  = {ct, mj, sw, sd, po};

% The values are
% computed by the SAME backend calls as the dissipativity scripts ex_diss_poisson,
% ex_diss_mjls, ex_diss_sd, ex_diss_sw (shooting horizon 1000, relative
% bisection tolerance 1e-6; the disturbance-diffusion column is kept, which the
% former mvst_l2gain_ct(...,{},Dwc) path dropped).  The dispatched defaults of
% MVST_PERF (horizon 60, tolerance 1e-3) under-estimated the gains, and the
% switched value of MVST_L2GAIN_CT_DWELL is a Nelder-Mead search, not exact:
% the switched entry is now the degree-4 SOS minimum dwell-time bound, an UPPER
% bound, reported with the worst-mode LOWER bound as in Ch08.
Hz = 1000; rt = 1e-6;
gam = zeros(1,numel(syss));
gam(1) = mvst_l2gain_ct_poisson(ct.A,ct.Bw,ct.C,ct.Dw,{zeros(2)},{Dwcm},{},[],{},Hz,rt);
gam(2) = mvst_l2gain_ct_mjls({mj.A0(:,:,1),mj.A0(:,:,2)}, mj.Bw, ...
            {mj.C(:,:,1),mj.C(:,:,2)}, mj.Dw, mj.Pig, Hz, ...
            {0.2*eye(2),0.2*eye(2)}, {Dwcm}, rt);
sysm  = struct('A0',sw.A0,'Bw',sw.Bw,'C',sw.C,'Dw',sw.Dw,'Dwc',Dwcm);
gam(3) = mvst_diss_sos_modes(sysm, struct('type','l2gain'), ...
            struct('type','switched','Tbar',sw.Td,'deg',4), 'gamma').gamma;   % SOS MDT upper bound
gSwLo  = max(mvst_l2gain_ct_poisson(sw.A0(:,:,1),sw.Bw,sw.C(:,:,1),sw.Dw,{zeros(2)},{Dwcm},{},[],{},Hz,rt), ...
             mvst_l2gain_ct_poisson(sw.A0(:,:,2),sw.Bw,sw.C(:,:,2),sw.Dw,{zeros(2)},{Dwcm},{},[],{},Hz,rt)); % worst-mode lower bound
gam(4) = mvst_l2gain_sd_hybrid(sd.A, sd.Bw, sd.C, sd.Tbar, ceil(Hz/sd.Tbar), {}, {Dwcm}, rt);
gam(5) = mvst_l2gain_ct_poisson(po.A,po.Bw,po.C,po.Dw,{po.Ac},{po.Dwc},{po.Jf}, ...
            po.lamf,{po.Bdf},Hz,rt);
tg  = gam.^2;                       % th:perf:tracegain: trace gain = L2-gain^2
for k = 1:numel(syss)
    fprintf('%-13s  L2-gain = %.6f   trace gain = %.6f\n', names{k}, gam(k), tg(k));
end
fprintf('switched: worst-mode lower bound %.6f, SOS upper bound %.6f\n', gSwLo, gam(3));
% diagnostic: the dispatched defaults (former values of this script)
for k = [1 2 4 5]
    fprintf('  [diagnostic] mvst_perf default %-13s = %.6f\n', names{k}, mvst_perf(syss{k},'l2').gamma);
end

% ---- honest coverage demonstration: a covariance measure on a multi-mode
%      class is NOT Gramian-characterized and MUST raise 'unsupported' -------
try
    mvst_perf(mj, 'h2');
    fprintf('WARNING: expected mvst:perf:unsupported for h2 on the mjls class\n');
catch ME
    fprintf('coverage: h2 on the Markov jump class -> %s (as designed)\n', ME.identifier);
end

% ---- data macros --------------------------------------------------------
w = mvst_texdata(fullfile(figdir,'data','perf_classes.tex'));
for k = 1:numel(syss)
    w.num(['PerfGam' tags{k}], gam(k));
    w.num(['PerfTG'  tags{k}], tg(k));
end
w.num('PerfGamSwLo', gSwLo);  w.num('PerfTGSwLo', gSwLo^2);
w.num('PerfHoriz', Hz, 0);
w.close();

% ================= figure: grouped bar chart per class ===================
figure('Color','w','Units','centimeters','Position',[2 2 12.5 8.5]);
ax = axes; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
B = bar(ax, [gam(:) tg(:)], 'grouped', 'BarWidth',0.9);
B(1).FaceColor = gc.order(1,:); B(1).EdgeColor = 'none';
B(2).FaceColor = gc.order(2,:); B(2).EdgeColor = 'none';
set(ax,'XTick',1:numel(names),'XTickLabel',names,'FontSize',12);
ylabel(ax,'performance value','Interpreter','latex');
lg = legend(ax, {'$L_2\times\ell_2$-gain $\gamma$\hspace{0.5em}', ...
    'trace gain $\gamma^2$\hspace{0.5em}'}, ...
    'Interpreter','latex','Location','northwest','FontSize',12);
set(lg,'Box','on');
mvst_savepdf(gcf, fullfile(figdir,'perf_classes.pdf'));
end
