function U = mvst_svecbasis(n)
%MVST_SVECBASIS  Orthonormal basis matrix of the symmetric subspace.
%   U = MVST_SVECBASIS(N) returns U in R^(n^2 x n(n+1)/2) with orthonormal
%   columns (U'*U = I) such that vec(X) = U*svec(X) for symmetric X, where
%   svec is the scaled half-vectorization (off-diagonal entries scaled by
%   sqrt(2), so that svec is an isometry for the trace inner product).
%   Used to form the REDUCED matrix representation U'*M*U of a linear
%   operator on S^n from its full (Kronecker) representation M; the
%   spectrum of the operator on its domain S^n is that of the reduced
%   matrix, not of the full one.
%
%   Example:
%       U = mvst_svecbasis(3);   % 9 x 6, with U'*U = eye(6)
%
%   See also MVST_LIFTSYM, MVST_VECH, MVST_UNVECH.
if ~(isnumeric(n) && isscalar(n) && isreal(n) && n >= 1 && n == floor(n))
    error('mvst:svecbasis:n', 'N must be a positive integer scalar.');
end
m = n*(n+1)/2;
U = zeros(n*n, m);
col = 0;
for j = 1:n
    for i = j:n
        col = col + 1;
        if i == j
            U((j-1)*n + i, col) = 1;
        else
            U((j-1)*n + i, col) = 1/sqrt(2);
            U((i-1)*n + j, col) = 1/sqrt(2);
        end
    end
end
end
