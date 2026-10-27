function Mred = mvst_liftsym(A0, Cc)
%MVST_LIFTSYM  Reduced matrix representation of X -> He(A0 X) + sum C_l X C_l'.
%   MRED = MVST_LIFTSYM(A0, CC) returns the matrix representation, on the
%   symmetric subspace (scaled half-vectorization, see MVST_SVECBASIS), of
%   the flow operator
%       F(X) = A0*X + X*A0' + sum_l CC{l}*X*CC{l}',
%   i.e. MRED = U'*(kron(I,A0) + kron(A0,I) + sum_l kron(CC{l},CC{l}))*U.
%   The reduced matrix has dimension n(n+1)/2 and its spectrum is the
%   spectrum of the operator on its domain S^n; pass A0 = zeros(n) to lift
%   a purely completely positive (jump-type) operator.
%
%   Example:
%       M = mvst_liftsym([0 1; -1 0], {0.2*eye(2)});   % 3 x 3 reduced matrix
%
%   See also MVST_SVECBASIS, MVST_STABILITY.
if ~(ismatrix(A0) && size(A0,1) == size(A0,2))
    error('mvst:liftsym:notSquare', 'A0 must be a square matrix.');
end
if ~iscell(Cc)
    error('mvst:liftsym:cell', 'CC must be a cell array of n-by-n matrices.');
end
n = size(A0,1);
for l = 1:numel(Cc)
    if ~isequal(size(Cc{l}), [n n])
        error('mvst:liftsym:dim', ...
            'CC{%d} must be %d-by-%d to match A0.', l, n, n);
    end
end
M = kron(eye(n), A0) + kron(A0, eye(n));
for l = 1:numel(Cc)
    M = M + kron(Cc{l}, Cc{l});
end
U = mvst_svecbasis(n);
Mred = U' * M * U;
end
