%% EX_STAB_DWELLSCALAR  Example "Dwell-time acts in both directions".
%  Chapter "Stability", Section "Dwell-time stability".  Scalar impulsive
%  second moment X of  dx = a x dt,  x(t_k^+) = e^b x(t_k^-):
%    flow scales X by e^{2 a tau} over a lapse tau, each jump by e^{2 b},
%    so under constant dwell-time Tbar the pre-jump samples obey the
%    geometric law  X(t_{k+1}^-) = e^{2(a Tbar + b)} X(t_k^-).
%  Two regimes, with opposite stable sides of the threshold:
%    (a) a < 0 < b : stable iff Tbar > b/|a|   (flow must undo the jump);
%    (b) b < 0 < a : stable iff Tbar < |b|/a   (jumps must catch the flow).
%
%  Each regime is drawn as its own figure, and the per-cycle ratio curve
%  is a separate figure of its own (linear axes throughout):
%    ex_stab_dwellscalar_a.pdf        sample second moment, regime (a)
%    ex_stab_dwellscalar_b.pdf        sample second moment, regime (b)
%    ex_stab_dwellscalar_ratio_a.pdf  per-cycle ratio, regime (a)
%    ex_stab_dwellscalar_ratio_b.pdf  per-cycle ratio, regime (b)
%  plus ONE native full-width combined figure (2x2 tiles drawn inside
%  MATLAB at readable font sizes, panels (a)-(d) in the order above):
%    stab_dwellscalar_combined.pdf
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ---------- regime (a): contracting flow, expanding jump ----------
a1 = -1.0; b1 = 0.5;  th1 = b1/abs(a1);          % threshold b/|a| = 0.5
Ta_stab = 0.9;  Ta_unst = 0.3;                    % above / below threshold

% ---------- regime (b): expanding flow, contracting jump ----------
a2 =  1.0; b2 = -0.7; th2 = abs(b2)/a2;           % threshold |b|/a = 0.7
Tb_stab = 0.4;  Tb_unst = 1.0;                    % below / above threshold

Ncyc = 5;                                          % number of dwell-time cycles

% ===================== regime (a): two figures =====================
saw_figure(gc, a1, b1, Ta_unst, Ta_stab, Ncyc, ...
    '$a<0<b$ (flow contracts, jumps expand)', ...
    {sprintf('$\\bar T=%.1f<b/|a|$',Ta_unst), ...
     sprintf('$\\bar T=%.1f>b/|a|$',Ta_stab)}, ...
    fullfile(mvst_figdir(),'ex_stab_dwellscalar_a.pdf'));
ratio_figure(gc, a1, b1, th1, '$\bar T^\star=b/|a|$', ...
    fullfile(mvst_figdir(),'ex_stab_dwellscalar_ratio_a.pdf'));

% ===================== regime (b): two figures =====================
saw_figure(gc, a2, b2, Tb_stab, Tb_unst, Ncyc, ...
    '$b<0<a$ (jumps contract, flow expands)', ...
    {sprintf('$\\bar T=%.1f<|b|/a$',Tb_stab), ...
     sprintf('$\\bar T=%.1f>|b|/a$',Tb_unst)}, ...
    fullfile(mvst_figdir(),'ex_stab_dwellscalar_b.pdf'));
ratio_figure(gc, a2, b2, th2, '$\bar T^\star=|b|/a$', ...
    fullfile(mvst_figdir(),'ex_stab_dwellscalar_ratio_b.pdf'));

% ============ native combined 2x2 figure (full text width) ============
fig = figure('Color','w','Units','centimeters','Position',[2 2 15 13.5]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
saw_axes(ax, gc, a1, b1, Ta_unst, Ta_stab, Ncyc, ...
    {'(a) $a<0<b$','(flow contracts, jumps expand)'}, ...
    {sprintf('$\\bar T=%.1f<b/|a|$',Ta_unst), ...
     sprintf('$\\bar T=%.1f>b/|a|$',Ta_stab)}, 1);
ax = nexttile(tl);
saw_axes(ax, gc, a2, b2, Tb_stab, Tb_unst, Ncyc, ...
    {'(b) $b<0<a$','(jumps contract, flow expands)'}, ...
    {sprintf('$\\bar T=%.1f<|b|/a$',Tb_stab), ...
     sprintf('$\\bar T=%.1f>|b|/a$',Tb_unst)}, 1);
ax = nexttile(tl);
hl = ratio_axes(ax, gc, a1, b1, th1, '(c) $\bar T^\star=b/|a|$');
ax = nexttile(tl);
ratio_axes(ax, gc, a2, b2, th2, '(d) $\bar T^\star=|b|/a$');
% one shared legend for (c)-(d): identical entries in both panels
lg = legend(hl, {'ratio $e^{2(a\bar T+b)}$','unit level $1$','threshold'}, ...
    'Interpreter','latex','Orientation','horizontal','FontSize',10);
lg.Layout.Tile = 'south'; set(lg,'Box','off');
set(findall(fig,'Type','axes'),'FontSize',10);
mvst_savepdf(fig, fullfile(mvst_figdir(),'stab_dwellscalar_combined.pdf'));

% ========================================================================
function saw_figure(gc, a, b, Tfirst, Tsecond, Ncyc, ttl, leg, outpath)
%SAW_FIGURE  Sample second moment for two dwell-times, linear axes, legend
%  placed outside below the frame.  Colour is fixed by STABILITY, not call
%  order: the stable trajectory (per-cycle ratio e^{2(a T + b)} < 1) is blue
%  and solid, the unstable one red and dashed.  The two are started from two
%  different initial conditions (the stable one higher) so the convergence of
%  the stable dynamics is clearly visible alongside the growing unstable one.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
Ts = [Tfirst Tsecond]; h = gobjects(1,2);
for i = 1:2
    stable = exp(2*(a*Ts(i)+b)) < 1;
    if stable, col = gc.order(1,:); ls = '-';  X0 = 2.5;   % stable  = blue,  higher IC
    else,      col = gc.order(2,:); ls = '--'; X0 = 0.35;  % unstable = red, lower IC
    end
    h(i) = plot_saw(ax, a, b, Ts(i), Ncyc, col, ls, X0);
end
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'second moment $X(t)$','Interpreter','latex');
title(ax, ttl,'Interpreter','latex');
lg = legend(ax, h, leg, 'Interpreter','latex', ...
    'Orientation','horizontal','Location','southoutside');
set(lg,'Box','off');
mvst_savepdf(fig, outpath);
end

function ratio_figure(gc, a, b, thr, ttl, outpath)
%RATIO_FIGURE  Per-cycle ratio e^{2(a Tbar + b)} versus Tbar, own figure,
%  linear axes, legend outside below.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 7]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
Tg = linspace(0, 2*max(thr,0.5), 200);
hr = plot(ax, Tg, exp(2*(a*Tg+b)), '-','LineWidth',1.8,'Color',gc.order(1,:));
hy = yline(ax, 1, ':','LineWidth',1.2,'Color',[0.4 0.4 0.4]);
hx = xline(ax, thr, '--','LineWidth',1.4,'Color',[0.5 0.5 0.5]);
xlabel(ax,'dwell-time $\bar T$','Interpreter','latex');
ylabel(ax,'per-cycle ratio $e^{2(a\bar T+b)}$','Interpreter','latex');
title(ax, ttl,'Interpreter','latex');
lg = legend(ax, [hr hy hx], ...
    {'ratio $e^{2(a\bar T+b)}$','unit level $1$','threshold'}, ...
    'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off');
mvst_savepdf(fig, outpath);
end

function h = plot_saw(ax, a, b, Tbar, Ncyc, col, ls, X0)
%PLOT_SAW  Continuous flow between jumps, multiplicative reset at each jump,
%  started from X(0)=X0.  Returns a representative handle for the legend.
X = X0; t0 = 0; npt = 40; h = [];
for k = 1:Ncyc
    tt = linspace(0, Tbar, npt);
    hh = plot(ax, t0+tt, X*exp(2*a*tt), ls, 'LineWidth',1.8,'Color',col);
    if k == 1, h = hh; end
    Xm = X*exp(2*a*Tbar);                 % pre-jump value
    Xp = exp(2*b)*Xm;                      % post-jump value
    plot(ax, [t0+Tbar t0+Tbar], [Xm Xp], ':','LineWidth',1.0,'Color',col, ...
         'HandleVisibility','off');
    X = Xp; t0 = t0 + Tbar;
end
end

function saw_axes(ax, gc, a, b, Tfirst, Tsecond, Ncyc, ttl, leg, ncol)
%SAW_AXES  Same content as SAW_FIGURE, drawn into a given (tiled) axes; the
%  legend sits outside below the frame, NCOL columns, no box.
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
Ts = [Tfirst Tsecond]; h = gobjects(1,2);
for i = 1:2
    stable = exp(2*(a*Ts(i)+b)) < 1;
    if stable, col = gc.order(1,:); ls = '-';  X0 = 2.5;
    else,      col = gc.order(2,:); ls = '--'; X0 = 0.35;
    end
    h(i) = plot_saw(ax, a, b, Ts(i), Ncyc, col, ls, X0);
end
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'second moment $X(t)$','Interpreter','latex');
title(ax, ttl,'Interpreter','latex','FontSize',10);
lg = legend(ax, h, leg, 'Interpreter','latex','NumColumns',ncol, ...
    'Location','southoutside','FontSize',10);
set(lg,'Box','off');
end

function hl = ratio_axes(ax, gc, a, b, thr, ttl)
%RATIO_AXES  Same content as RATIO_FIGURE, drawn into a given (tiled) axes;
%  returns the three handles for a legend shared across tiles.
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
Tg = linspace(0, 2*max(thr,0.5), 200);
hr = plot(ax, Tg, exp(2*(a*Tg+b)), '-','LineWidth',1.8,'Color',gc.order(1,:));
hy = yline(ax, 1, ':','LineWidth',1.2,'Color',[0.4 0.4 0.4]);
hx = xline(ax, thr, '--','LineWidth',1.4,'Color',[0.5 0.5 0.5]);
% tight x-range, headroom above the curve so it never runs into the frame
xlim(ax,[0 Tg(end)]); ylim(ax,[0 1.12*max(exp(2*(a*Tg+b)))]);
xlabel(ax,'dwell-time $\bar T$','Interpreter','latex');
ylabel(ax,'per-cycle ratio $e^{2(a\bar T+b)}$','Interpreter','latex');
title(ax, ttl,'Interpreter','latex','FontSize',10);
hl = [hr hy hx];
end
