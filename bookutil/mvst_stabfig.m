function mvst_stabfig(figdir, key, t, Xo, Xc, Uc, tc, trO, trC, disc, gc)
%MVST_STABFIG  Standard open-vs-closed stabilization figures (Stabilization
%  chapter Examples).  Writes THREE vector PDFs into FIGDIR:
%    stabz_<KEY>_paths.pdf    : 3 stacked panels, open-loop state sample paths
%       (top, red), closed-loop paths (middle, blue), closed-loop control
%       input (bottom, purple);
%    stabz_<KEY>_cov.pdf      : trace of the covariance, open loop RED, closed
%       loop BLUE, log-y;
%    stabz_<KEY>_combined.pdf : native single-figure merge of the two above,
%       the three path panels stacked in the left column and the semilog
%       covariance trace spanning the right column (2x1x3 tiledlayout).  This
%       is the figure the book uses; the two separate PDFs are kept so nothing
%       downstream breaks.
%  Arguments:
%    t      1 x nt time (DISC=false) or step (DISC=true) axis of the paths
%    Xo,Xc  n x nt x nreal open / closed state sample paths (row 1 is plotted)
%    Uc     m x nt x nreal closed-loop control input (row 1 is plotted)
%    tc,trO,trC  covariance axis and open / closed trace curves
%    disc   true for discrete-time (stairs, "step k"), false for CT (plot,
%           "time t"; the paths are decimated to keep the PDF light)
%    gc     mvst_gemcolors() struct (gc.order: 1=blue, 2=red, 4=purple)
%
%   See also MVST_STABTABLE, MVST_SAVEPDF, MVST_GEMCOLORS.
if disc; xu = '(k)'; xl = 'step $k$'; else; xu = '(t)'; xl = 'time $t$'; end

% decimate the plotted CT paths so the vector PDF stays light
if disc; d = 1:numel(t); else; d = 1:max(1,round(numel(t)/180)):numel(t); end
xr = [t(d(1)) t(d(end))];

% ----------------------------- paths figure (standalone) -------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 15]);
tl  = tiledlayout(fig,3,1,'Padding','compact','TileSpacing','compact');
pathpanels(tl, t, d, Xo, Xc, Uc, disc, xu, xl, xr, gc);
mvst_savepdf(fig, fullfile(figdir, ['stabz_' key '_paths.pdf']));

% ----------------------------- covariance figure (standalone) --------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
covpanel(axes(fig), tc, trO, trC, disc, xu, xl, gc);
mvst_savepdf(fig, fullfile(figdir, ['stabz_' key '_cov.pdf']));

% ----------------------------- combined figure (native merge) --------------
% Vertically STACKED (4x1) so every panel keeps the full text width and the
% fonts stay readable at \linewidth in the book: three path panels on top
% (open, closed, control) and the semilog covariance trace at the bottom.
fig = figure('Color','w','Units','centimeters','Position',[2 2 13.5 20]);
tl  = tiledlayout(fig,4,1,'Padding','compact','TileSpacing','compact');
a1 = nexttile(tl); a2 = nexttile(tl); a3 = nexttile(tl);
drawpaths(a1, a2, a3, t, d, Xo, Xc, Uc, disc, xu, xl, xr, gc);
a4 = covpanel(nexttile(tl), tc, trO, trC, disc, xu, xl, gc);
for a = [a1 a2 a3 a4]; set(a,'FontSize',11); end
mvst_savepdf(fig, fullfile(figdir, ['stabz_' key '_combined.pdf']));
end

% ===========================================================================
function pathpanels(tl, t, d, Xo, Xc, Uc, disc, xu, xl, xr, gc)
% The three stacked path panels inside an already-created 3x1 tiledlayout.
a1 = nexttile(tl); a2 = nexttile(tl); a3 = nexttile(tl);
drawpaths(a1, a2, a3, t, d, Xo, Xc, Uc, disc, xu, xl, xr, gc);
end

function drawpaths(a1, a2, a3, t, d, Xo, Xc, Uc, disc, xu, xl, xr, gc)
panel(a1, t(d), squeeze(Xo(1,d,:)), [gc.order(2,:) 0.35], disc, ...
      ['$x_1' xu '$'], 'open loop (unstable)');
panel(a2, t(d), squeeze(Xc(1,d,:)), [gc.order(1,:) 0.35], disc, ...
      ['$x_1' xu '$'], 'closed loop (stabilized)');
panel(a3, t(d), squeeze(Uc(1,d,:)), [gc.order(4,:) 0.35], disc, ...
      ['$u' xu '$'], 'control input');
xlabel(a3, xl, 'Interpreter','latex');
for a = [a1 a2 a3]; xlim(a, xr); end
end

function ax = panel(ax, t, Y, rgba, disc, ylab, ttl)
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
for r = 1:size(Y,2)
    if disc
        stairs(ax, t, Y(:,r), 'Color', rgba, 'LineWidth', 0.6);
    else
        plot(ax, t, Y(:,r), 'Color', rgba, 'LineWidth', 0.6);
    end
end
ylabel(ax, ylab, 'Interpreter','latex');
title(ax, ttl, 'Interpreter','latex');
end

function ax = covpanel(ax, tc, trO, trC, disc, xu, xl, gc)
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
if disc
    stairs(ax, tc, trC, '-','LineWidth',2.0,'Color',gc.order(1,:));
    stairs(ax, tc, trO, '-','LineWidth',2.0,'Color',gc.order(2,:));
else
    plot(ax, tc, trC, '-','LineWidth',2.0,'Color',gc.order(1,:));
    plot(ax, tc, trO, '-','LineWidth',2.0,'Color',gc.order(2,:));
end
set(ax,'YScale','log'); xlim(ax, [tc(1) tc(end)]);
xlabel(ax, xl, 'Interpreter','latex');
ylabel(ax, ['$\mathrm{tr}\,X' xu '$'], 'Interpreter','latex');
lg = legend(ax, {'closed loop','open loop'}, 'Interpreter','latex', ...
            'Location','northwest'); set(lg,'Box','off');
end
