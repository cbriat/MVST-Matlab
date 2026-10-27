function s = mvst_fmtmat(M, prec)
%MVST_FMTMAT  Format a numeric matrix as a LaTeX bmatrix string.
%   S = MVST_FMTMAT(M) returns '\begin{bmatrix} ... \end{bmatrix}' with the
%   entries of M (3 decimals by default).  S = MVST_FMTMAT(M,PREC) sets the
%   number of decimals.  Scalars are returned as a plain formatted number.
%
%   Example:
%       s = mvst_fmtmat([1 2; 3 4], 1);   % '\begin{bmatrix}...\end{bmatrix}'
%
%   See also MVST_TEXDATA.
if nargin < 2 || isempty(prec); prec = 3; end
if ~isnumeric(M) || ~ismatrix(M)
    error('mvst:fmtmat:input', ...
        'M must be a 2-D numeric matrix (got a %s with %d dimensions).', ...
        class(M), ndims(M));
end
f = ['%.' num2str(prec) 'f'];
if isscalar(M)
    s = sprintf(f, M); return;
end
[r,~] = size(M);
rows = cell(1,r);
for i = 1:r
    e = arrayfun(@(x) sprintf(f, x), M(i,:), 'UniformOutput', false);
    rows{i} = strjoin(e, ' & ');
end
% NOTE: strjoin translates escape sequences in its delimiter, so ' \\ '
% would collapse to a single backslash.  Use four backslashes so the LaTeX
% row separator '\\' survives the un-escaping.
s = ['\begin{bmatrix}' strjoin(rows, ' \\\\ ') '\end{bmatrix}'];
end
