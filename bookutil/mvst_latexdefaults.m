function mvst_latexdefaults()
%MVST_LATEXDEFAULTS  Use the LaTeX interpreter for all new text objects.
%   Call once at the start of a plotting script so that every label, title,
%   legend, and tick label is rendered with the LaTeX interpreter (the
%   MATLAB default is 'tex', which rejects \mathbb, \mathrm, ^\top, ...).
%
%   See also MVST_GEMCOLORS, MVST_FIGSTYLE, MVST_SAVEPDF.
set(groot, 'defaultTextInterpreter',            'latex');
set(groot, 'defaultAxesTickLabelInterpreter',   'latex');
set(groot, 'defaultLegendInterpreter',          'latex');
try
    set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');
catch
end
% Standard MATLAB "gem" color order (blue, red, yellow, purple, green, ...),
% so every new axes cycles colours in the canonical palette order and resets
% to blue on each new panel.  Explicit per-series colours below should index
% mvst_gemcolors().order(k,:) in plotting order (blue first, then red, ...).
try
    set(groot, 'defaultAxesColorOrder', mvst_gemcolors().order);
catch
end
end
