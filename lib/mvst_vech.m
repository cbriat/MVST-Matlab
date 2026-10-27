function v = mvst_vech(M)
%MVST_VECH  Stack the independent (upper-triangular) entries of symmetric M.
%   V = MVST_VECH(M) returns the length n(n+1)/2 column vector of the
%   upper-triangular entries of the n-by-n (symmetric) matrix M, read row by
%   row: [M(1,1),M(1,2),...,M(1,n),M(2,2),...,M(n,n)]'.  It is the inverse of
%   MVST_UNVECH on symmetric matrices, i.e. M == mvst_unvech(mvst_vech(M),n)
%   whenever M is symmetric.
%
%   INPUT
%     M   square (symmetric) n-by-n matrix; only the upper triangle is read.
%   OUTPUT
%     V   n(n+1)/2 x 1 vector of the upper-triangular entries.
%
%   Example:
%       v = mvst_vech([1 2; 2 3]);   % -> [1; 2; 3]
%
%   See also MVST_UNVECH, MVST_SVECBASIS.
if ~ismatrix(M) || size(M,1) ~= size(M,2)
    error('mvst:vech:notSquare', 'M must be a square matrix.');
end
n = size(M,1); v = zeros(n*(n+1)/2,1); k = 0;
for i = 1:n
    for j = i:n; k = k+1; v(k) = M(i,j); end
end
end
