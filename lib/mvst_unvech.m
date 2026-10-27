function M = mvst_unvech(v, n)
%MVST_UNVECH  Rebuild a symmetric n-by-n matrix from MVST_VECH output.
%   M = MVST_UNVECH(V,N) is the inverse of MVST_VECH: it places the length
%   n(n+1)/2 vector V into the upper triangle of an n-by-n matrix (row by row)
%   and mirrors it, returning a symmetric M.
%
%   INPUTS
%     V   n(n+1)/2 x 1 vector, as produced by MVST_VECH;
%     N   positive integer, the matrix dimension.
%   OUTPUT
%     M   symmetric n-by-n matrix.
%
%   Example:
%       M = mvst_unvech([1; 2; 3], 2);   % -> [1 2; 2 3]
%
%   See also MVST_VECH, MVST_SVECBASIS.
if ~(isnumeric(n) && isscalar(n) && isreal(n) && n >= 1 && n == floor(n))
    error('mvst:unvech:n', 'N must be a positive integer scalar.');
end
if numel(v) ~= n*(n+1)/2
    error('mvst:unvech:size', ...
        'V has %d entries but N=%d requires %d.', numel(v), n, n*(n+1)/2);
end
M = zeros(n); k = 0;
for i = 1:n
    for j = i:n; k = k+1; M(i,j) = v(k); M(j,i) = v(k); end
end
end
