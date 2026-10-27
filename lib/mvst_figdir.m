function d = mvst_figdir(varargin)
%MVST_FIGDIR  Output folder of the example scripts (figures and data macros).
%   D = MVST_FIGDIR() returns the folder to which the examples write their
%   vector PDFs and creates it (and its data/ subfolder) when missing.  The
%   folder is, in order of precedence,
%     1. the environment variable MVST_FIGDIR, when set and nonempty;
%     2. ./figures under the CURRENT folder otherwise.
%   D = MVST_FIGDIR(SUB1, SUB2, ...) returns fullfile(D, SUB1, SUB2, ...).
%
%   Example:
%       setenv('MVST_FIGDIR', fullfile(tempdir, 'mvst_out'));
%       ex_mjls_gains                      % writes under tempdir/mvst_out
%
%   See also MVST_SAVEPDF, MVST_TEXDATA.
d = getenv('MVST_FIGDIR');
if isempty(d), d = fullfile(pwd, 'figures'); end
if ~exist(d, 'dir'), mkdir(d); end
if ~exist(fullfile(d, 'data'), 'dir'), mkdir(fullfile(d, 'data')); end
if nargin > 0, d = fullfile(d, varargin{:}); end
end
