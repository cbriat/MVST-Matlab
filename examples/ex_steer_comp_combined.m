%% EX_STEER_COMP_COMBINED  Consolidated covariance-entry steering plots.
%  Native, vertically-stacked replacement for the former two-panel minipage
%  tiling (fig:steer:sd:comp): the covariance entries versus time for the two
%  structurally distinct steering cases with marked event times.
%    (a) sampled-data (sampling instants dotted),
%    (b) switched     (switching instants dotted).
%  Each case is recomputed by the generic steering wrapper MVST_STEER exactly
%  as in EX_STEER_SD / EX_STEER_SWITCHED.  Figure only: writes no data macros.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors(); D = ex_data();

Xf = [1.0 0;0 1.6];

% ---- (a) sampled-data --------------------------------------------------
s = D.sd; Ta = 3.0; Xa0 = [3 0;0 0.5];
ra = mvst_steer(s, Xf, Ta, struct('class','sd','X0',Xa0,'M0',-eye(2)));
tka = s.Tbar:s.Tbar:(Ta-1e-9);                       % sampling instants

% ---- (b) switched ------------------------------------------------------
Td = D.sw.Tdwell; Tb = 4.0; Xb0 = [3 0;0 0.5];
rb = mvst_steer(D.sw, Xf, Tb, struct('class','switched','X0',Xb0));
swb = Td:Td:(Tb-1e-9);                               % switching instants

% ---- combined figure ---------------------------------------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tl  = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(tl); comp_ax(ax1, ra.t, ra.Sigma, gc, tka, false);
title(ax1,'(a) sampled-data','Interpreter','latex','FontWeight','normal');
ax2 = nexttile(tl); comp_ax(ax2, rb.t, rb.Sigma, gc, swb, true);
title(ax2,'(b) switched','Interpreter','latex','FontWeight','normal');
% No legend: the color key (Sigma_11 gem blue, Sigma_12 gem red, Sigma_22 gem
% yellow, i.e. gc.order(1:3,:) in canonical order) is carried by line samples
% in the book caption.
mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_steer_comp_combined.pdf'));

% ========================================================================
function h = comp_ax(ax, tg, Sig, gc, events, showx)
%COMP_AX  Covariance entries versus time into AX (verbatim from
%  MVST_STEER_COMPONENTS, minus the figure/legend creation).
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
tg = tg(:).';
for e = events(:).'
    xline(ax, e, ':', 'Color',[0.6 0.6 0.6],'LineWidth',0.8,'HandleVisibility','off');
end
n = size(Sig,1);
cols = arrayfun(@(k) gc.order(k,:), 1:size(gc.order,1), 'UniformOutput', false);
h = gobjects(0); c = 0;
for i = 1:n
    for j = i:n
        c = c + 1; v = reshape(Sig(i,j,:),1,[]);
        h(end+1) = plot(ax, tg, v, '-','LineWidth',1.8,'Color',cols{mod(c-1,numel(cols))+1}); %#ok<AGROW>
    end
end
if showx, xlabel(ax,'$t$','Interpreter','latex'); end
ylabel(ax,'covariance entries','Interpreter','latex');
set(ax,'FontSize',11);
end
