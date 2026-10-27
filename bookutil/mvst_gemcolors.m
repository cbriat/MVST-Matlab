function C = mvst_gemcolors()
%MVST_GEMCOLORS  MATLAB's built-in "gem" color order (with a fallback).
%   C = MVST_GEMCOLORS() returns the standard MATLAB gem palette obtained
%   from ORDEREDCOLORS("gem") (R2023b and later).  For backward
%   compatibility with the example scripts, the named fields
%   (.sapphire, .ruby, .topaz, .amethyst, .emerald, .teal, .garnet) are
%   simply positional handles into that order, and .order is the full Nx3
%   matrix (use it directly with COLORORDER).  On older MATLAB / Octave a
%   fixed gem-like palette is used instead.
%
%   Example:
%       gc = mvst_gemcolors();
%       plot(t, y, 'Color', gc.order(1,:));   % first (blue) series
%
%   See also MVST_LATEXDEFAULTS, ORDEREDCOLORS, COLORORDER.
try
    M = orderedcolors("gem");                 % standard built-in palette
catch
    M = [ 0.063 0.322 0.612;                  % fallback (gem-like)
          0.706 0.094 0.227;
          0.918 0.643 0.129;
          0.512 0.251 0.612;
          0.063 0.553 0.388;
          0.000 0.502 0.502;
          0.451 0.090 0.176 ];
end
C.order = M;
names = {'sapphire','ruby','topaz','amethyst','emerald','teal','garnet'};
for i = 1:min(numel(names), size(M,1))
    C.(names{i}) = M(i,:);
end
C.onyx = [0.149 0.149 0.176];
end
