%% STAB_DWELLTYPES  Schematic jump patterns for the dwell-time families.
%  Chapter "Stability", Section "Dwell-time stability", subsection
%  "Dwell-time notions".  Each figure is a timeline on which the jump
%  instants t_k are drawn as vertical stems; dashed guides mark the
%  constraint(s) on the spacings Delta_k = t_{k+1} - t_k.  Purely
%  illustrative (no dynamics are simulated).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

% ---------- constant (periodic) dwell-time ----------
Tbar = 1.0; tk = Tbar*(1:8);
fig = timeline_fig(gc);
draw_timeline(gca, tk, gc.order(1,:));
annotate_span(gca, tk(3), tk(4), '$\bar T$', [0.5 0.5 0.5]);
title(gca,'Constant dwell-time: $\Delta_k\equiv\bar T$','Interpreter','latex');
mvst_savepdf(fig, fullfile(figdir,'stab_dwell_constant.pdf'));

% ---------- minimum dwell-time ----------
Tmin = 0.8;
d = [1.0 2.3 0.85 1.6 0.9 2.0 1.2];  tk = cumsum(d);   % all spacings >= Tmin
fig = timeline_fig(gc);
draw_timeline(gca, tk, gc.order(1,:));
yl = ylim(gca);
plot(gca, [0 tk(end)+1], (yl(1)+0.12*(yl(2)-yl(1)))*[1 1], 'k:'); % baseline
annotate_span(gca, tk(3), tk(3)+Tmin, '$\bar T$', [0.5 0.5 0.5]);
annotate_span(gca, tk(5), tk(5)+Tmin, '$\bar T$', [0.5 0.5 0.5]);
title(gca,'Minimum dwell-time: $\Delta_k\ge\bar T$','Interpreter','latex');
mvst_savepdf(fig, fullfile(figdir,'stab_dwell_minimum.pdf'));

% ---------- range dwell-time ----------
Tmn = 0.7; Tmx = 1.6;
d = [0.8 1.5 1.0 0.75 1.55 0.9 1.3];  tk = cumsum(d);  % spacings in [Tmn,Tmx]
fig = timeline_fig(gc);
draw_timeline(gca, tk, gc.order(1,:));
annotate_span(gca, tk(2), tk(2)+Tmn, '$T_{\min}$', [0.5 0.5 0.5]);
annotate_span(gca, tk(4), tk(4)+Tmx, '$T_{\max}$', [0.5 0.5 0.5]);
title(gca,'Range dwell-time: $\Delta_k\in[T_{\min},T_{\max}]$','Interpreter','latex');
mvst_savepdf(fig, fullfile(figdir,'stab_dwell_range.pdf'));

% ---------- other families (maximum / average / persistent) ----------
fig = figure('Color','w','Units','centimeters','Position',[2 2 15 11]);
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
Tmax = 1.2; d = [0.5 1.15 0.7 1.1 0.45 0.9 1.05]; tk = cumsum(d);  % all <= Tmax
draw_timeline(ax, tk, gc.order(1,:));
annotate_span(ax, tk(2), tk(2)+Tmax, '$\bar T$', [0.5 0.5 0.5]);
title(ax,'Maximum dwell-time: $\Delta_k\le\bar T$','Interpreter','latex');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
% average: bursts of close jumps compensated by long gaps (mean rate ~ 1/tau*)
d = [0.2 0.25 1.9 0.22 0.24 2.0 0.3]; tk = cumsum(d);
draw_timeline(ax, tk, gc.order(1,:));
title(ax,'Average dwell-time: $\kappa(t,s)\le N_0+(t-s)/\tau^\ast$', ...
      'Interpreter','latex');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
% persistent: a long jump-free interval (>= tau_D) recurring among busy runs
d = [0.3 0.35 0.3 2.4 0.32 0.28 0.34]; tk = cumsum(d);
draw_timeline(ax, tk, gc.order(1,:));
annotate_span(ax, tk(3), tk(4), '$\tau_D$', [0.5 0.5 0.5]);
title(ax,'Persistent dwell-time: recurring gaps $\ge\tau_D$', ...
      'Interpreter','latex');
xlabel(ax,'time $t$','Interpreter','latex');
mvst_savepdf(fig, fullfile(figdir,'stab_dwell_families.pdf'));

% ---------- native combined single column (six dwell-time families) ---------
% Vertically STACKED (6x1) so each timeline keeps the full text width and its
% labels stay readable at \linewidth in the book: constant, minimum, range,
% maximum, average, persistent (top to bottom).  Tall tiles (23 cm overall)
% leave vertical room inside each axes so that the span labels clear both
% the frame and the stems.
fig = figure('Color','w','Units','centimeters','Position',[2 2 15 23]);
tl = tiledlayout(fig,6,1,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
Tbar = 1.0; tk = Tbar*(1:8);
draw_timeline(ax, tk, gc.order(1,:));
annotate_span(ax, tk(3), tk(4), '$\bar T$', [0.5 0.5 0.5]);
title(ax,'(a) constant: $\Delta_k\equiv\bar T$','Interpreter','latex');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
Tmin = 0.8; d = [1.0 2.3 0.85 1.6 0.9 2.0 1.2]; tk = cumsum(d);
draw_timeline(ax, tk, gc.order(1,:));
annotate_span(ax, tk(3), tk(3)+Tmin, '$\bar T$', [0.5 0.5 0.5]);
annotate_span(ax, tk(5), tk(5)+Tmin, '$\bar T$', [0.5 0.5 0.5]);
title(ax,'(b) minimum: $\Delta_k\ge\bar T$','Interpreter','latex');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
Tmn = 0.7; Tmx = 1.6; d = [0.8 1.5 1.0 0.75 1.55 0.9 1.3]; tk = cumsum(d);
draw_timeline(ax, tk, gc.order(1,:));
annotate_span(ax, tk(2), tk(2)+Tmn, '$T_{\min}$', [0.5 0.5 0.5]);
annotate_span(ax, tk(4), tk(4)+Tmx, '$T_{\max}$', [0.5 0.5 0.5]);
title(ax,'(c) range: $\Delta_k\in[T_{\min},T_{\max}]$','Interpreter','latex');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
Tmax = 1.2; d = [0.5 1.15 0.7 1.1 0.45 0.9 1.05]; tk = cumsum(d);
draw_timeline(ax, tk, gc.order(1,:));
annotate_span(ax, tk(2), tk(2)+Tmax, '$\bar T$', [0.5 0.5 0.5]);
title(ax,'(d) maximum: $\Delta_k\le\bar T$','Interpreter','latex');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
d = [0.2 0.25 1.9 0.22 0.24 2.0 0.3]; tk = cumsum(d);
draw_timeline(ax, tk, gc.order(1,:));
title(ax,'(e) average: $\kappa(t,s)\le N_0+(t-s)/\tau^\ast$','Interpreter','latex');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on');
d = [0.3 0.35 0.3 2.4 0.32 0.28 0.34]; tk = cumsum(d);
draw_timeline(ax, tk, gc.order(1,:));
annotate_span(ax, tk(3), tk(4), '$\tau_D$', [0.5 0.5 0.5]);
title(ax,'(f) persistent: recurring gaps $\ge\tau_D$','Interpreter','latex');
xlabel(ax,'time $t$','Interpreter','latex');
mvst_savepdf(fig, fullfile(figdir,'stab_dwell_combined.pdf'));

% ========================================================================
function fig = timeline_fig(~)
fig = figure('Color','w','Units','centimeters','Position',[2 2 15 4.2]);
ax = axes(fig); hold(ax,'on'); box(ax,'on');
xlabel(ax,'time $t$','Interpreter','latex');
end

function draw_timeline(ax, tk, col)
%DRAW_TIMELINE  Baseline plus a vertical stem at each jump instant t_k.
%  The x-range starts slightly left of 0 so that the first stem and its
%  label stay clear of the left frame; the y-range leaves headroom above
%  the span guides (ANNOTATE_SPAN) and room below for the t_k labels.
tmax = tk(end) + 0.8; tmin = -0.3;
plot(ax, [tmin tmax], [0 0], '-','LineWidth',1.2,'Color',[0.3 0.3 0.3]);
for k = 1:numel(tk)
    plot(ax, [tk(k) tk(k)], [0 1], '-','LineWidth',2.0,'Color',col);
    plot(ax, tk(k), 1, '.','MarkerSize',14,'Color',col);
    text(ax, tk(k), -0.08, sprintf('$t_{%d}$',k), 'Interpreter','latex', ...
         'HorizontalAlignment','center','VerticalAlignment','top','FontSize',9);
end
set(ax,'YTick',[],'XTick',[]); ylim(ax,[-0.6 2.45]); xlim(ax,[tmin tmax]);
set(ax,'FontSize',10);
end

function annotate_span(ax, x1, x2, lab, col)
%ANNOTATE_SPAN  A double-headed guide between x1 and x2 with a label.
y = 1.45;
plot(ax, [x1 x2], [y y], '-','LineWidth',1.4,'Color',col);
plot(ax, [x1 x1], [y-0.12 y+0.12], '-','LineWidth',1.4,'Color',col);
plot(ax, [x2 x2], [y-0.12 y+0.12], '-','LineWidth',1.4,'Color',col);
text(ax, 0.5*(x1+x2), y+0.2, lab, 'Interpreter','latex', ...
     'HorizontalAlignment','center','VerticalAlignment','bottom', ...
     'FontSize',10,'Color',col);
end
