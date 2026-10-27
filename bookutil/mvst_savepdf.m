function mvst_savepdf(fig, fname)
%MVST_SAVEPDF  Apply the white style and export a figure to a vector PDF.
%   MVST_SAVEPDF(FIG, FNAME) writes FIG to FNAME (a full path ending in
%   .pdf), forcing a white background and black text via MVST_FIGSTYLE so
%   that the result is correct even when MATLAB runs in dark mode.
%
%   See also MVST_FIGSTYLE, EXPORTGRAPHICS.
if ~(ischar(fname) || (isstring(fname) && isscalar(fname))) ...
        || strlength(string(fname)) == 0
    error('mvst:savepdf:fname', 'FNAME must be a non-empty file-path string.');
end
fname = char(fname);
d = fileparts(fname);
if ~isempty(d) && ~exist(d, 'dir')
    [ok, emsg] = mkdir(d);
    if ~ok
        error('mvst:savepdf:mkdir', ...
            'Cannot create output folder ''%s'': %s', d, emsg);
    end
end
mvst_figstyle(fig);
try
    exportgraphics(fig, fname, 'ContentType', 'vector', ...
                   'BackgroundColor', 'white');     % R2020a+
catch
    set(fig, 'Units', 'centimeters');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'centimeters', 'PaperSize', pos(3:4), ...
             'PaperPositionMode', 'manual', 'PaperPosition', [0 0 pos(3:4)]);
    % Legacy fallback for old MATLAB without EXPORTGRAPHICS; '-painters' is
    % deliberate: the '-vector' replacement postdates this branch's targets.
    print(fig, '-dpdf', '-painters', fname);   %#ok<PRTPT>
end
fprintf('  wrote %s\n', fname);
end
