function ex_stab_rates_combined()
%EX_STAB_RATES_COMBINED  Native vertically-stacked consolidation of the three
%  stability-by-lifting plots of fig:stab:rates, replacing the former side-by-side
%  minipage tiling (ex_stab_rates.pdf, itself a two-row panel, next to
%  sd_stab_sos.pdf).  ALL curves are RECOMPUTED FROM SOURCE here (no stale PDFs
%  are loaded):
%    (a) Poissonian-jump diffusion (ex00_stability.m, section (b)): exact spectral
%        abscissa of the lifted flow operator versus the jump rate lambda_1;
%        the critical rate lambda* (off the plotted window) is returned by the
%        matrix-SOS feasibility bisection (MVST_STAB_SOS, arbitrary dwell-time).
%    (b) sampled-data system with the nominal gain (ex00_stability.m, section (c)):
%        exact monodromy spectral radius versus the sampling period, the maximum
%        stable period Tbar* returned by the matrix-SOS feasibility bisection
%        (periodic dwell-time), where the exact radius crosses 1.
%    (c) sampled-data worked instance (ex_stab_sd_sos.m): exact augmented
%        monodromy spectral radius versus the sampling period, with the maximum
%        stable period Tbar* from the timer-dependent matrix-SOS bisection marked.
%  Figure only: this script writes NO data macros.
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();
n = 2; nr = n*(n+1)/2;

fprintf('\n==== stability-by-lifting combined figure (all curves recomputed) ====\n');

% ===== panel (a): GBM + Poissonian jumps, abscissa versus rate =====
A  = [0 1;-2 -1];   C1 = 0.15*eye(2);   J1 = 1.15*eye(2);
L0 = mvst_liftsym(A,{C1});
LJ = mvst_liftsym(zeros(n),{J1}) - eye(nr);
lams = linspace(0, 2.5, 120); absc = zeros(size(lams));   % exact abscissa curve
for k = 1:numel(lams)
    absc(k) = max(real(eig(L0 + lams(k)*LJ)));
end
% critical Poisson rate lambda* by SOS feasibility bisection (crossing off window)
sfeasPoi = @(lam) mvst_stab_sos( ...
    struct('A0',A,'Ac',C1,'Jf',J1,'lamf',lam), struct('type','arbitrary')).slack;
lamCrit = sos_bisect(sfeasPoi, 0, 4.0, 'feaslo');
fprintf('(a) Poisson       : lambda* SOS=%.4f (off the plotted [0,2.5] window)\n', lamCrit);

% ===== panel (b): sampled-data with nominal gain, monodromy versus period =====
D = ex_data(); s = D.sd; TbarNom = s.Tbar;
out = mvst_igre_periodic_sd(s.A,s.B,s.Ac,s.Bc,s.Q,s.R,TbarNom);
Kn = out.K; m = size(s.B,2); na = n + m;
Ab0 = [s.A s.B; zeros(m,na)];
Abl = [s.Ac s.Bc; zeros(m,na)];
Jb  = [eye(n) zeros(n,m); Kn zeros(m,m)];
Lsd = mvst_liftsym(Ab0,{Abl});
Jlift = mvst_liftsym(zeros(na),{Jb});
Tbars = linspace(0.05, 2.0, 120); rsd = zeros(size(Tbars));   % exact radius curve
for k = 1:numel(Tbars)
    rsd(k) = max(abs(eig(Jlift*expm(Tbars(k)*Lsd))));
end
sfeasSd = @(Tb) mvst_stab_sos( ...
    struct('A0',Ab0,'Ac',Abl,'J',Jb), struct('type','periodic','Tbar',Tb,'deg',4)).slack;
TbarCrit = sos_bisect(sfeasSd, 0.05, 2.0, 'feaslo');
TbarCritExact = crossing(Tbars, rsd, 1);
fprintf('(b) sampled-data  : Tbar* SOS=%.4f exact=%.4f (rel. dev. %.2e), nominal Tbar=%.4f\n', ...
    TbarCrit, TbarCritExact, abs(TbarCrit-TbarCritExact)/TbarCritExact, TbarNom);

% ===== panel (c): sampled-data worked instance, augmented monodromy =====
Asd = [0 1; 0 0.1];  Bsd = [0;1];
Ksd = place(Asd,Bsd,[-1 -2]);  Ksd = -Ksd;
if max(real(eig(Asd+Bsd*Ksd))) >= -1e-9, Ksd = -place(Asd,Bsd,[-1 -2]); end
Atil = [Asd Bsd; zeros(1,3)];
Jtil = [eye(2) [0;0]; Ksd 0];
tol = 1e-6; lo = 0.05; hi = 3.0;
for it = 1:24
    mid = 0.5*(lo+hi);
    if slackval_sd(Atil,Jtil,mid) > tol, lo = mid; else, hi = mid; end
end
Tstar = lo;
Tgrid = linspace(0.4*Tstar, 1.6*Tstar, 50);
rad   = arrayfun(@(T) mono_radius_sd(Atil,Jtil,T), Tgrid);
radAt = mono_radius_sd(Atil,Jtil,Tstar);
fprintf('(c) sd instance   : Tbar* SOS=%.4f  radAt=%.4f\n', Tstar, radAt);

% ===== combined 3-row figure (native, readable) =====
figc = figure('Color','w','Units','centimeters','Position',[2 2 13 20]);
tlc  = tiledlayout(figc,3,1,'TileSpacing','compact','Padding','compact');

axa = nexttile(tlc); hold(axa,'on'); box(axa,'on'); grid(axa,'on');
plot(axa, lams, absc, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(axa, 0, ':','LineWidth',1.2,'Color',[0.4 0.4 0.4]);
xline(axa, lamCrit, '--','LineWidth',1.6,'Color',[0.5 0.5 0.5]);
xlabel(axa,'jump rate $\lambda_1$','Interpreter','latex');
ylabel(axa,'lifted-flow abscissa','Interpreter','latex');
lga = legend(axa, {'abscissa','$0$','$\lambda_1^\star$'}, ...
    'Interpreter','latex','Location','northwest'); set(lga,'Box','off');
xlim(axa, [lams(1) lams(end)]);   % keep the plotted lambda window (lambda* off-window)
title(axa,'(a) Poissonian-jump diffusion','Interpreter','latex','FontWeight','normal');
set(axa,'FontSize',11);

axb = nexttile(tlc); hold(axb,'on'); box(axb,'on'); grid(axb,'on');
plot(axb, Tbars, rsd, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(axb, 1, ':','LineWidth',1.2,'Color',[0.4 0.4 0.4]);
xline(axb, TbarCrit, '--','LineWidth',1.6,'Color',[0.5 0.5 0.5]);
xline(axb, TbarNom, ':','LineWidth',1.4,'Color',[0.5 0.5 0.5]);
xlabel(axb,'sampling period $\bar T$','Interpreter','latex');
ylabel(axb,'monodromy radius','Interpreter','latex');
lgb = legend(axb, {'$\rho(\bar T)$','$1$','$\bar T^\star$','nominal $\bar T$'}, ...
    'Interpreter','latex','Location','northwest'); set(lgb,'Box','off');
title(axb,'(b) Sampled-data system, nominal gain','Interpreter','latex','FontWeight','normal');
set(axb,'FontSize',11);

axc = nexttile(tlc); hold(axc,'on'); box(axc,'on'); grid(axc,'on');
plot(axc, Tgrid, rad, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(axc, 1, '--','Color',gc.order(2,:),'LineWidth',1.5,'HandleVisibility','off');
xline(axc, Tstar, ':','Color',gc.onyx,'LineWidth',1.5,'HandleVisibility','off');
plot(axc, Tstar, radAt, 'o','MarkerSize',7,'MarkerFaceColor',gc.order(3,:), ...
     'MarkerEdgeColor','k','HandleVisibility','off');
xlabel(axc,'sampling period $\bar T$','Interpreter','latex');
ylabel(axc,'augmented monodromy spectral radius','Interpreter','latex');
title(axc,sprintf('(c) Sampled-data worked instance, $\\bar T^\\star=%.3f$',Tstar), ...
    'Interpreter','latex','FontWeight','normal');
set(axc,'FontSize',11);

mvst_savepdf(figc, fullfile(figdir,'stab_rates_combined.pdf'));
fprintf('  wrote stab_rates_combined.pdf (figure only, no data macros)\n');
end

% =========================================================================
function s = slackval_sd(Atil,Jtil,Tbar)
%SLACKVAL_SD  Optimal slack of the sampled-data analysis (dwell type 'periodic',
%  degree-4 certificate), assembled and solved by MVST_STAB_SOS.  Feasible iff s>0.
sys = struct('A0', Atil, 'J', Jtil);
res = mvst_stab_sos(sys, struct('type','periodic','Tbar',Tbar,'deg',4));
s = res.slack;
end

% =========================================================================
function r = mono_radius_sd(Atil,Jtil,Tbar)
%MONO_RADIUS_SD  Exact augmented one-period monodromy spectral radius on the
%  lifted symmetric space: flow expm(Atil Tbar) then the sampling jump Jtil(.)Jtil'.
LF = mvst_liftsym(Atil, {});
MJ = mvst_liftsym(zeros(size(Atil)), {Jtil});
r = max(abs(eig(MJ*expm(LF*Tbar))));
end

% =========================================================================
function xc = sos_bisect(sfun, lo, hi, mode)
%SOS_BISECT  Bisection of a monotone SOS-feasibility boundary (see ex00_stability).
tol = 1e-7;
switch mode
    case 'feaslo'
        assert(sfun(lo) >  tol, 'sos_bisect: lo not feasible');
        assert(sfun(hi) <= tol, 'sos_bisect: hi not infeasible');
        for it = 1:30
            mid = 0.5*(lo+hi);
            if sfun(mid) > tol, lo = mid; else, hi = mid; end
        end
    case 'feashi'
        assert(sfun(lo) <= tol, 'sos_bisect: lo not infeasible');
        assert(sfun(hi) >  tol, 'sos_bisect: hi not feasible');
        for it = 1:30
            mid = 0.5*(lo+hi);
            if sfun(mid) > tol, hi = mid; else, lo = mid; end
        end
end
xc = 0.5*(lo+hi);
end

% =========================================================================
function xc = crossing(x, y, level)
%CROSSING  Last crossing of Y through LEVEL (linear interpolation).
xc = NaN;
for k = numel(x)-1:-1:1
    if (y(k) - level)*(y(k+1) - level) <= 0
        xc = x(k) + (level - y(k))*(x(k+1) - x(k))/(y(k+1) - y(k));
        return
    end
end
end
