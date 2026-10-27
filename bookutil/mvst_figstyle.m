function mvst_figstyle(fig)
%MVST_FIGSTYLE  Force a white background and black text on a figure.
%   MVST_FIGSTYLE(FIG) overrides any dark theme (e.g. MATLAB in dark mode)
%   so that exported figures have a white background, black axes/ticks,
%   black labels/titles, and legends with a white box and black text.
%   Call before exporting; MVST_SAVEPDF calls it automatically.
%
%   See also MVST_SAVEPDF, MVST_LATEXDEFAULTS.
if nargin < 1 || isempty(fig); fig = gcf; end

% R2025a+ : pin the figure to the light theme if the property exists.
try
    set(fig, 'Theme', 'light');
catch
end

set(fig, 'Color', 'w', 'InvertHardcopy', 'off');

ax = findall(fig, 'Type', 'axes');
for k = 1:numel(ax)
    a = ax(k);
    set(a, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
           'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.15, ...
           'MinorGridColor', [0.15 0.15 0.15]);
    set(get(a,'Title'),  'Color', 'k');
    set(get(a,'XLabel'), 'Color', 'k');
    set(get(a,'YLabel'), 'Color', 'k');
    set(get(a,'ZLabel'), 'Color', 'k');
end

lg = findall(fig, 'Type', 'Legend');
for k = 1:numel(lg)
    set(lg(k), 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.4 0.4 0.4]);
end

cb = findall(fig, 'Type', 'ColorBar');
for k = 1:numel(cb)
    set(cb(k), 'Color', 'k');
end
end
