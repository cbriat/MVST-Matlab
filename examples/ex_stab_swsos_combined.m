function ex_stab_swsos_combined()
%EX_STAB_SWSOS_COMBINED  Native vertically-stacked consolidation of the two
%  switched-stability panels of fig:stab:swsos, replacing the former side-by-side
%  minipage tiling (ex_stab_dwell.pdf + sw_stab_sos.pdf).  BOTH system instances
%  are RECOMPUTED FROM SOURCE here (no stale PDFs are loaded):
%    (a) first instance (Asw1,Asw2 of ex00_stability.m, section (d)): the exact
%        lifted two-mode monodromy spectral radius over one switching cycle versus
%        the dwell-time; the minimum dwell-time T* is returned by the mode-indexed
%        matrix-SOS feasibility bisection (MVST_STAB_SOS_MODES) and cross-checked
%        against the exact radius crossing of 1.
%    (b) second instance (A1,A2 of ex_stab_switched_sos.m): the exact two-cycle
%        monodromy spectral radius versus the dwell-time, with the minimum
%        dwell-time Tbar* returned by the mode-indexed matrix-SOS feasibility
%        bisection marked, where the exact radius crosses 1.
%  Figure only: this script writes NO data macros.
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

fprintf('\n==== switched-stability combined figure (both instances recomputed) ====\n');

% ===== panel (a): first instance, exact monodromy radius vs dwell =====
Asw1 = [-1 4; 0 -1];   Asw2 = [-1 0; 4 -1];        % Hurwitz, unstable average
Lsw1 = mvst_liftsym(Asw1,{}); Lsw2 = mvst_liftsym(Asw2,{});
Tds = linspace(0.02, 3.5, 160); rsw = zeros(size(Tds));   % exact radius curve
for k = 1:numel(Tds)
    rsw(k) = max(abs(eig(expm(Tds(k)*Lsw2)*expm(Tds(k)*Lsw1))));
end
% minimum dwell-time T* by mode-indexed SOS feasibility bisection, cross-checked
% by the exact radius crossing of 1.
sfeasSw = @(Td) mvst_stab_sos_modes( ...
    struct('A0',cat(3,Asw1,Asw2)), struct('type','switched','Tbar',Td,'deg',4)).slack;
TdCrit = sos_bisect(sfeasSw, 0.1, 3.5, 'feashi');
TdCritExact = crossing(Tds, rsw, 1);
fprintf('(a) first instance : T* SOS=%.4f  exact=%.4f  (rel. dev. %.2e)\n', ...
    TdCrit, TdCritExact, abs(TdCrit-TdCritExact)/TdCritExact);

% ===== panel (b): second instance, two-cycle monodromy + SOS min dwell =====
A1 = [-1 0; 2.5 -1];  A2 = [-1 2.5; 0 -1];         % double eig -1, shear split
tol = 1e-6; lo = 0.3; hi = 3.0;
for it = 1:22
    mid = 0.5*(lo+hi);
    if slackval_sw(A1,A2,mid) > tol, hi = mid; else, lo = mid; end
end
Tstar = hi;
Tg2  = linspace(0.5*Tstar, 1.8*Tstar, 50);
rad2 = arrayfun(@(T) mono_radius_sw(A1,A2,T), Tg2);
radAt    = mono_radius_sw(A1,A2,Tstar);
radBelow = mono_radius_sw(A1,A2,0.8*Tstar);
radAbove = mono_radius_sw(A1,A2,1.2*Tstar);
fprintf('(b) second instance: Tbar* SOS=%.4f  radAt=%.4f  radBelow=%.4f  radAbove=%.4f\n', ...
    Tstar, radAt, radBelow, radAbove);

% ===== combined 2-row figure (native, readable) =====
figc = figure('Color','w','Units','centimeters','Position',[2 2 13 15]);
tlc  = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');

axa = nexttile(tlc); hold(axa,'on'); box(axa,'on'); grid(axa,'on');
plot(axa, Tds, rsw, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(axa, 1, ':','LineWidth',1.2,'Color',[0.4 0.4 0.4]);
xline(axa, TdCrit, '--','LineWidth',1.6,'Color',[0.5 0.5 0.5]);
xlabel(axa,'dwell-time $T$','Interpreter','latex');
ylabel(axa,'monodromy radius over one cycle','Interpreter','latex');
title(axa,sprintf('(a) First instance, $T^\\star=%.3f$',TdCrit), ...
    'Interpreter','latex','FontWeight','normal');
lga = legend(axa, {'$\rho(T)$','$1$','$T^\star$'},'Interpreter','latex', ...
    'Location','northeast'); set(lga,'Box','off');
set(axa,'FontSize',11);

axb = nexttile(tlc); hold(axb,'on'); box(axb,'on'); grid(axb,'on');
plot(axb, Tg2, rad2, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(axb, 1, '--','Color',gc.order(2,:),'LineWidth',1.5,'HandleVisibility','off');
xline(axb, Tstar, ':','Color',gc.onyx,'LineWidth',1.5,'HandleVisibility','off');
plot(axb, Tstar, radAt, 'o','MarkerSize',7,'MarkerFaceColor',gc.order(3,:), ...
     'MarkerEdgeColor','k','HandleVisibility','off');
xlabel(axb,'dwell-time $\bar T$','Interpreter','latex');
ylabel(axb,'two-cycle monodromy spectral radius','Interpreter','latex');
title(axb,sprintf('(b) Second instance, $\\bar T^\\star=%.3f$',Tstar), ...
    'Interpreter','latex','FontWeight','normal');
set(axb,'FontSize',11);

mvst_savepdf(figc, fullfile(figdir,'stab_swsos_combined.pdf'));
fprintf('  wrote stab_swsos_combined.pdf (figure only, no data macros)\n');
end

% =========================================================================
function s = slackval_sw(A1,A2,Tbar)
%SLACKVAL_SW  Optimal common slack of the two-mode switched analysis (dwell type
%  'switched'), assembled and solved by MVST_STAB_SOS_MODES.  Feasible iff s > 0.
sys = struct('A0', cat(3, A1, A2));
res = mvst_stab_sos_modes(sys, struct('type','switched','Tbar',Tbar,'deg',4));
s = res.slack;
end

% =========================================================================
function r = mono_radius_sw(A1,A2,Tbar)
%MONO_RADIUS_SW  Exact two-mode monodromy spectral radius on the lifted symmetric
%  space over one switch cycle: flow mode 1 for Tbar, then flow mode 2 for Tbar.
L1 = mvst_liftsym(A1, {});  L2 = mvst_liftsym(A2, {});
r = max(abs(eig(expm(L2*Tbar)*expm(L1*Tbar))));
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
