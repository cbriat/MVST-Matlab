function mvst_steer_components(tg, Sig, ttl, fname, varargin) %#ok<INUSD>
%MVST_STEER_COMPONENTS  Covariance entries versus time.
%   MVST_STEER_COMPONENTS(TG,SIG,TTL,FNAME) plots the distinct entries of the
%   symmetric covariance SIG(:,:,k) against time TG, each entry in a different
%   colour, and saves a white-background vector PDF to FNAME.  This exposes
%   the time profile of the covariance, including the jumps or slope changes
%   that occur at sampling or switching instants.
%
%   MVST_STEER_COMPONENTS(...,EVENTS) additionally draws faint dotted vertical
%   guide lines at the times in EVENTS (e.g. sampling or switching instants)
%   to make the impulsive/switching structure visible.
%
%   See also MVST_STEER, MVST_STEER_ELLIPSE, MVST_SAVEPDF.
gc = mvst_gemcolors();
n = size(Sig,1); tg = tg(:).';
events = []; if ~isempty(varargin); events = varargin{1}(:).'; end
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
cols = arrayfun(@(k) gc.order(k,:), 1:size(gc.order,1), 'UniformOutput', false);
for e = events
    xline(ax, e, ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8, 'HandleVisibility','off');
end
h = gobjects(0); lbl = {}; c = 0;
for i = 1:n
    for j = i:n
        c = c + 1;
        v = reshape(Sig(i,j,:),1,[]);
        h(end+1) = plot(ax, tg, v, '-', 'LineWidth', 1.8, ...
                        'Color', cols{mod(c-1,numel(cols))+1}); %#ok<AGROW>
        lbl{end+1} = sprintf('$\\Sigma_{%d%d}$', i, j);          %#ok<AGROW>
    end
end
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'covariance entries','Interpreter','latex');
lg = legend(ax, h, lbl, 'Interpreter','latex', 'Orientation','horizontal', ...
            'Location','southoutside');
set(lg,'Box','off'); set(ax,'FontSize',11);
mvst_savepdf(fig, fname);
end
