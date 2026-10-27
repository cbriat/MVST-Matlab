function ex_stab_levy()
%EX_STAB_LEVY  Linear system driven by LEVY noise: second-moment flow
%  operator, mean-square stability, and the "only the first two moments of
%  the Levy measure matter" phenomenon.  Stability chapter, Section
%  "Examples / Continuous-time systems" (Example ex:stab:levy).
%
%  SYSTEM.   dx = A0 x dt + sum_i A_i x dL_i,   L_i independent Levy
%  processes with FINITE SECOND MOMENT, mean rate mu_i = E[L_i(1)] and
%  variance rate c_i = E[[L_i](1)] = Var[L_i(1)].
%
%  SECOND MOMENT.  Taking expectations in the Ito formula for x x' (the
%  martingale parts average out; the quadratic variation contributes
%  c_i A_i X A_i'), the second moment X = E[x x'] obeys the LINEAR
%  matrix-valued flow
%      Xdot = F(X) = He[(A0 + sum_i mu_i A_i) X] + sum_i c_i A_i X A_i',
%  which is EXACTLY the canonical form He[A0bar X] + sum_l Abar_l X Abar_l'
%  with effective drift A0bar = A0 + sum_i mu_i A_i and channels
%  Abar_i = sqrt(c_i) A_i.  Hence the Levy triplet enters ONLY through
%  (mu_i, c_i): any two Levy noises with the same first two moments give
%  the SAME second-moment dynamics, hence the same mean-square stability.
%  Mean-square exponential stability <=> spectral abscissa of the lifted
%  operator (mvst_liftsym) is negative.
%
%  WHAT THIS SCRIPT SHOWS.
%   (a) Sample paths of three very different driving Levy noises, all with
%       mu = 0 and the SAME variance rate c: a Brownian motion, a
%       compensated Poisson process with frequent small jumps, and a
%       compensated Poisson process with rare large jumps.  This single
%       driver panel is the book figure stab_levy_combined.pdf.
%   (b) The mean-square stability boundary: the lifted spectral abscissa as
%       a function of the variance rate c, crossing zero at c*.
%
%  The boundary CURVE is a spectral quantity, drawn exactly from the reduced
%  matrix representation on the symmetric subspace; the CRITICAL variance rate
%  c* marked on it is found by BISECTION ON THE GENERIC MATRIX-SOS ANALYSIS
%  program MVST_STAB_SOS (the folded Levy generator is a constant-certificate CT
%  system, dwell type 'arbitrary': drift A0, channel sqrt(c) A1), and
%  cross-checked against the exact abscissa crossing (must agree to ~1e-2 rel.).
%  (The Monte-Carlo second-moment panel is not drawn: the
%  invariance is already the exact spectral statement, so only the driver
%  panel and the stability boundary are drawn.)

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();
datdir = fullfile(figdir,'data');

% ---------------------------------------------------------------- data
A0 = [0 1; -2 -0.5];            % Hurwitz drift (eig -0.25 +- 1.39i)
A1 = [0 1;  1  0];              % single Levy channel
n  = size(A0,1);
x0 = [1;0];  X0 = x0*x0';
fprintf('A0 eig %s   (Hurwitz)\n', mat2str(eig(A0)',3));

% -------------------------------------------- mean-square stability limit
% mu = 0, so A0bar = A0 and the channel is sqrt(c)*A1.  The boundary CURVE is the
% exact lifted spectral abscissa versus the variance rate c (a spectral reading).
absc = @(c) max(real(eig(mvst_liftsym(A0, {sqrt(c)*A1}))));
cgrid = linspace(0, 1.2, 241);
ag    = arrayfun(absc, cgrid);
% critical variance rate c* by BISECTION ON THE SOS FEASIBILITY (folded Levy is a
% constant-certificate CT system, dwell type 'arbitrary'): stable for small c.
sfeas = @(c) mvst_stab_sos(struct('A0',A0,'Ac',sqrt(c)*A1), ...
                           struct('type','arbitrary')).slack;
tolS = 1e-7; lo = 0; hi = 1.2;
assert(sfeas(lo) >  tolS, 'lower bracket infeasible');
assert(sfeas(hi) <= tolS, 'upper bracket feasible (raise hi)');
for k = 1:40
    mid = (lo+hi)/2;
    if sfeas(mid) > tolS, lo = mid; else, hi = mid; end
end
cstar = (lo+hi)/2;
% cross-check: exact abscissa crossing
lo2 = 0; hi2 = 1.2;
for k = 1:60
    mid = (lo2+hi2)/2;
    if absc(mid) < 0, lo2 = mid; else, hi2 = mid; end
end
cstarExact = (lo2+hi2)/2;
fprintf('critical variance rate c* = %.4f (SOS)  cross-check %.4f (exact abscissa), rel. dev. %.2e\n', ...
    cstar, cstarExact, abs(cstar-cstarExact)/cstarExact);
assert(abs(cstar-cstarExact)/cstarExact < 1e-2, ...
    'SOS critical variance rate disagrees with the exact abscissa crossing beyond 1e-2 relative');

% pick a STABLE operating point for the path comparison
c = 0.60*cstar;
M = mvst_liftsym(A0, {sqrt(c)*A1});
fprintf('operating point c = %.4f  -> lifted abscissa %.4f\n', c, max(real(eig(M))));

% ------------------------------------------------ three matched noises
% all with mu = 0 and variance rate c:
%   BM              : dL = sqrt(c) dW
%   Poisson (small) : dL = a(dN - lam dt),  a^2 lam = c,  lam large  a small
%   Poisson (large) : same, lam small, a large
T = 8; dt = 5e-4; nt = round(T/dt)+1; t = (0:nt-1)*dt;
lamS = 25;  aS = sqrt(c/lamS);
lamL = 0.4; aL = sqrt(c/lamL);
fprintf('matched noises: BM sigma=%.4f | Poisson small (lam=%.2f, a=%.4f) | Poisson large (lam=%.2f, a=%.4f)\n', ...
        sqrt(c), lamS, aS, lamL, aL);

% one representative driver sample path L(t) per matched noise (no Monte-Carlo
% second-moment estimate: the moment invariance is the exact spectral statement)
rng(7);
[~, ~, LBM] = mc_run('bm',   A0, A1, x0, c, [],   [],  T, dt, 1);
[~, ~, LPS] = mc_run('pois', A0, A1, x0, c, lamS, aS,  T, dt, 1);
[~, ~, LPL] = mc_run('pois', A0, A1, x0, c, lamL, aL,  T, dt, 1);

% --------------------------------------------------------------- figures
% One panel per figure (they are printed at full text width in the book).
% Decimate for plotting only; the curves are computed on the full grid.
d = 1:20:nt;

% (1) the three matched Levy drivers
fig1 = figure('Position',[100 100 560 330]); hold on; grid on; box on
plot(t(d), LBM(d), '-',  'Color', gc.sapphire, 'LineWidth', 1.0);
plot(t(d), LPS(d), '-',  'Color', gc.ruby,     'LineWidth', 1.0);
plot(t(d), LPL(d), '-',  'Color', gc.topaz,    'LineWidth', 1.2);
xlabel('$t$','Interpreter','latex');
ylabel('$L(t)$','Interpreter','latex');
legend({'Brownian','Poisson (frequent, small)','Poisson (rare, large)'}, ...
       'Location','best','Interpreter','latex','Box','off');
mvst_savepdf(fig1, fullfile(figdir,'stab_levy_paths.pdf'));

% (combined) single full-width panel: the three matched Levy drivers.  The
% Monte-Carlo second-moment panel is not drawn; the figure
% stab_levy_combined.pdf carries this driver panel alone.
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 8.5]);
axc = axes(figc); hold(axc,'on'); grid(axc,'on'); box(axc,'on');
plot(axc, t(d), LBM(d), '-', 'Color', gc.sapphire, 'LineWidth', 1.0);
plot(axc, t(d), LPS(d), '-', 'Color', gc.ruby,     'LineWidth', 1.0);
plot(axc, t(d), LPL(d), '-', 'Color', gc.topaz,    'LineWidth', 1.2);
xlabel(axc,'$t$','Interpreter','latex'); ylabel(axc,'$L(t)$','Interpreter','latex');
% Unframed legend placed above the axes (inside the axes the labels overlap the paths).
legend(axc,{'Brownian','Poisson (frequent, small)','Poisson (rare, large)'}, ...
       'Location','northoutside','NumColumns',2, ...
       'Interpreter','latex','Box','off');
set(axc,'FontSize',11);
mvst_savepdf(figc, fullfile(figdir,'stab_levy_combined.pdf'));

% (3) the mean-square stability boundary
fig3 = figure('Position',[100 100 560 330]); hold on; grid on; box on
plot(cgrid, ag, '-', 'Color', gc.sapphire, 'LineWidth', 1.4);
yline(0,'k-','LineWidth',0.8);
plot(cstar, 0, 'o', 'MarkerSize', 7, 'MarkerFaceColor', gc.ruby, ...
     'MarkerEdgeColor', gc.ruby);
xline(c,'--','Color',[.45 .45 .45],'LineWidth',1.0);
text(cstar, 0.12, sprintf('$c^\\star=%.3f$',cstar), 'Interpreter','latex', ...
     'HorizontalAlignment','center');
xlabel('variance rate $c$','Interpreter','latex');
ylabel('lifted spectral abscissa','Interpreter','latex');
mvst_savepdf(fig3, fullfile(figdir,'stab_levy_boundary.pdf'));
fprintf('  wrote stab_levy_{paths,combined,boundary}.pdf\n');

% ------------------------------------------------------------ data macros
if ~exist(datdir,'dir'); mkdir(datdir); end
fid = fopen(fullfile(datdir,'stab_levy.tex'),'w');
fprintf(fid,'%% Auto-generated by the mvst example scripts. Do not edit.\n');
fprintf(fid,'\\def\\StabLevyCstar{%.3f}\n', cstar);
fprintf(fid,'\\def\\StabLevyC{%.3f}\n', c);
fprintf(fid,'\\def\\StabLevyAbsc{%.3f}\n', max(real(eig(M))));
fprintf(fid,'\\def\\StabLevyLamSmall{%.2f}\n', lamS);
fprintf(fid,'\\def\\StabLevyASmall{%.3f}\n', aS);
fprintf(fid,'\\def\\StabLevyLamLarge{%.2f}\n', lamL);
fprintf(fid,'\\def\\StabLevyALarge{%.3f}\n', aL);
fclose(fid);
fprintf('  wrote data/stab_levy.tex\n');
end

% =====================================================================
function [trX, path1, Lpath] = mc_run(kind, A0, A1, x0, c, lam, a, T, dt, nreal)
%MC_RUN  Euler--Maruyama / jump simulation of dx = A0 x dt + A1 x dL,
%  with L a zero-mean Levy noise of variance rate c.  Returns the
%  Monte-Carlo E|x|^2, one sample path of x_1, and one path of L.
nt = round(T/dt)+1;
xs = repmat(x0(:), 1, nreal);
trX = zeros(1,nt); path1 = zeros(1,nt); Lpath = zeros(1,nt);
trX(1) = mean(sum(xs.^2,1)); path1(1) = xs(1,1);
sig = sqrt(c);
for k = 2:nt
    switch kind
        case 'bm'
            dW = sqrt(dt)*randn(1,nreal);
            xs = xs + (A0*xs)*dt + sig*(A1*xs).*dW;
            Lpath(k) = Lpath(k-1) + sig*dW(1);
        case 'pois'
            % compensated: dL = a(dN - lam dt), so drift -a*lam and jumps +a
            jump = double(rand(1,nreal) < lam*dt);
            xs = xs + (A0*xs)*dt - a*lam*(A1*xs)*dt + a*(A1*xs).*jump;
            Lpath(k) = Lpath(k-1) + a*(jump(1) - lam*dt);
    end
    trX(k) = mean(sum(xs.^2,1));
    path1(k) = xs(1,1);
end
end
