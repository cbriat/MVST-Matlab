function v = mvst_getf(s, f, d)
%MVST_GETF  Field of struct S named F, or default D if absent/empty.
%   V = MVST_GETF(S,F,D) returns S.(F) when S is a struct that carries a
%   non-empty field named F, and the default value D otherwise.  It is the
%   permissive field accessor used throughout the package to read optional
%   struct fields without an explicit ISFIELD guard at every call site.
%
%   INPUTS
%     S   struct (any other type yields the default);
%     F   field name (char);
%     D   default value returned when the field is absent or empty.
%   OUTPUT
%     V   S.(F) if present and non-empty, else D.
%
%   Example:
%       npts = mvst_getf(opts, 'npts', 60);   % opts.npts or 60
%
%   See also ISFIELD, GETFIELD.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
