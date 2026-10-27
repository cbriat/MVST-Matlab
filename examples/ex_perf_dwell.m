function ex_perf_dwell()
%EX_PERF_DWELL  Dwell-time performance tubes and the deterministic peak-to-peak
%   (star-norm) sandwich (Chapter "Further Performance Measures"): the peak and
%   energy-to-peak certified bounds versus the minimum dwell-time on an impulsive
%   instance, and the star-norm sandwich on its continuous-time slice.
%
%   Impulsive instance (the stable-flow / expansive-jump instance of the Ch06
%   dwell illustration, with a disturbance and output added):
%     flow  dx=(A x+Bw w)dt, z=C x,  A=[-1 1;-1 -1] (Hurwitz, oscillatory),
%           Bw=[0;1], C=[1 0];
%     jump  x+=J x,  J=1.2 I (EXPANSIVE).
%   The reachability tube is bounded only for dwell-times above the critical
%   value T* (the tube monodromy radius rho(J e^{F T}) reaches 1); above T* the
%   certified peak and energy-to-peak bounds decrease as the guaranteed dwell
%   grows (more flow contraction between the expansive jumps).
%
%   Star-norm sandwich (continuous-time slice, no jumps; prop:perf:ptplower):
%     DC-gain lower bound ||C A^{-1} Bw||  <=  true star-norm (L1 of the impulse
%     response, = realized worst-case peak-to-peak)  <=  certified p2p gamma.
%
%   Outputs: figures/perf_dwell.pdf (peak and e2p versus T_min, T* marked),
%   macros in figures/data/perf_dwell.tex.
%
%   See also MVST_PERF, MVST_PERF_DWELL, MVST_LIFTSYM.

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

A = [-1 1;-1 -1]; Bw = [0;1]; C = [1 0]; J = 1.2*eye(2);
s = struct('A',A,'Bw',Bw,'C',C,'J',J);

% ---- critical dwell T* (tube monodromy radius = 1) ------------------------
LF = mvst_liftsym(A,{}); LJ = mvst_liftsym(zeros(2),{J});
rho = @(T) max(abs(eig(LJ*expm(LF*T))));
lo = 0.01; hi = 1;
for it = 1:40; m = 0.5*(lo+hi); if rho(m) < 1; hi = m; else; lo = m; end; end
Tstar = hi;

% ---- peak / energy-to-peak versus minimum dwell-time ----------------------
Tg = linspace(Tstar+0.05, 2, 24); pk = nan(size(Tg)); ep = nan(size(Tg));
for i = 1:numel(Tg)
    rp = mvst_perf(s,'peak',struct('dwell','min','Tmin',Tg(i),'deg',4));
    re = mvst_perf(s,'e2p', struct('dwell','min','Tmin',Tg(i),'deg',4));
    pk(i) = rp.gamma; ep(i) = re.gamma;
end
pkHalf = interp1(Tg,pk,0.5); pkOne = interp1(Tg,pk,1);
epHalf = interp1(Tg,ep,0.5); epOne = interp1(Tg,ep,1);

% ---- star-norm sandwich on the continuous-time slice (no jumps) ------------
ct = struct('A',A,'Bw',Bw,'C',C);
rp2p = mvst_perf(ct,'p2p',struct());
gLB = rp2p.gammaLB; gP2P = rp2p.gamma;
tg = linspace(0,40,40000); imp = arrayfun(@(t) C*expm(A*t)*Bw, tg);
starnorm = trapz(tg, abs(imp));                 % true star-norm = L1 of impulse response
% realized worst-case peak-to-peak by a sign-matched (bang-bang) disturbance
realized = realized_ptp(A,Bw,C);

fprintf('[ex_perf_dwell] T*=%.4f  peak(0.5)=%.4f peak(1)=%.4f  e2p(0.5)=%.4f e2p(1)=%.4f\n', ...
        Tstar, pkHalf, pkOne, epHalf, epOne);
fprintf('   star-norm sandwich (CT slice): gammaLB=%.4f <= realized=%.4f ~ star=%.4f <= p2p=%.4f\n', ...
        gLB, realized, starnorm, gP2P);

% ---- figure ---------------------------------------------------------------
figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
plot(Tg, pk, '-', 'Color', gc.order(1,:), 'LineWidth', 1.8); hold on
plot(Tg, ep, '-', 'Color', gc.order(2,:), 'LineWidth', 1.8);
xline(Tstar, '-.', 'Color', gc.order(3,:), 'LineWidth', 2.4);
grid on; box on
set(gca,'FontSize',12);
% Y-range and ticks pinned to [0, 1.4], step 0.2 (at 12pt the automatic ticks drop the top label).
ylim([0 1.4]); yticks(0:0.2:1.4);
xlabel('minimum dwell-time $T_{\min}$','Interpreter','latex');
ylabel('certified bound','Interpreter','latex');
% Legend entries are padded on the right (\hspace) so the frame clears the longest LaTeX entry.
legend({'peak (reachability) bound\hspace{0.8em}','energy-to-peak bound\hspace{0.8em}', ...
        'critical dwell-time $T^\star$\hspace{0.8em}'}, ...
       'Interpreter','latex','Location','northeast','FontSize',12);
mvst_savepdf(gcf, fullfile(figdir,'perf_dwell.pdf'));

% ---- data macros ----------------------------------------------------------
try
    w = mvst_texdata(fullfile(figdir,'data','perf_dwell.tex'));
    w.num('PerfDwTstar',Tstar);
    w.num('PerfDwPeakHalf',pkHalf); w.num('PerfDwPeakOne',pkOne);
    w.num('PerfDwEtpHalf',epHalf);  w.num('PerfDwEtpOne',epOne);
    w.num('PerfDwLB',gLB); w.num('PerfDwStar',starnorm);
    w.num('PerfDwSim',realized); w.num('PerfDwPtp',gP2P);
    w.close();
catch, end
end

function v = realized_ptp(A,Bw,C)
% Realized peak-to-peak under a unit-peak, sign-matched (bang-bang) disturbance.
% For the feedthrough-free SISO system the worst disturbance for the output peak
% at time T is w(t)=sign(h(T-t)), h(r)=C e^{A r}Bw the impulse response; the
% realized peak |z(T)| approaches the true star-norm (L1 kernel norm) from below.
    T = 40; N = 8000; h = T/N; tg = (0:N)*h;
    hker = arrayfun(@(r) C*expm(A*r)*Bw, tg);            % impulse response
    Ad = expm(A*h); Bd = (Ad - eye(size(A)))*(A\Bw);
    x = zeros(size(A,1),1); v = 0;
    for k = 1:N
        tau = T - (k-1)*h;
        ws  = sign(interp1(tg, hker, max(tau,0), 'linear', 0));  % w(t)=sign(h(T-t))
        x   = Ad*x + Bd*ws;
        v   = max(v, abs(C*x));
    end
end
