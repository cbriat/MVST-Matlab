function [prog, Q] = mvst_sosposmatrvar(prog, n, d, vars)
%MVST_SOSPOSMATRVAR  SOS (positive-semidefinite) polynomial matrix variable.
%   [PROG, Q] = MVST_SOSPOSMATRVAR(PROG, N, D, VARS) declares an N-by-N matrix
%   Q(VARS) that is positive semidefinite BY CONSTRUCTION for every value of the
%   polynomial variables VARS, of degree D in VARS, and appends it to the
%   SOSTOOLS program PROG.  For N = 1 it returns a scalar sum-of-squares
%   polynomial of degree D (>= 0 everywhere) -- the workhorse multiplier used by
%   the positive-systems copositive-SOS gain/synthesis engines (gamma, alpha,
%   and every interval Positivstellensatz multiplier).
%
%   CLEAN PORT (attribution).  This reimplements the MatrixTools helper
%   `sosposmatrvar.m` (and the `sosposmatr.m` it calls) by
%       M. Peet, MatrixTools v.03 (matthew.peet@inria.fr),
%   distributed with the L-infinity positive-systems scripts of the NAHS
%   "Linf Positive Systems" package.  The original reaches into SOSTOOLS
%   internals (`prog.var.idx`) and reconstructs the Gram matrix from the
%   auto-generated `coeff_k` symbol names.  This port instead exposes the SOS
%   Gram matrix DIRECTLY through the modern SOSTOOLS API (`sossosvar` for the
%   scalar case, `sospolymatrixvar` + `sosmatrixineq` for the matrix case), so
%   there is no dependence on internal symbol bookkeeping.
%
%   See also MVST_SOSMATRVAR, MVST_POS_GAIN, MVST_POS_STABILIZE, SOSSOSVAR.

    Z = monomials(vars, 0:ceil(d/2));
    if n == 1
        % scalar SOS polynomial of degree d: Q = Z' * G * Z, G >= 0 (Gram).
        [prog, Q] = sossosvar(prog, Z);
    else
        % n x n symmetric polynomial matrix, constrained PSD over VARS.
        [prog, Q] = sospolymatrixvar(prog, Z, [n n], 'symmetric');
        prog = sosmatrixineq(prog, Q);
    end
end
