function [prog, P] = mvst_sosmatrvar(prog, Z, sizePQ)
%MVST_SOSMATRVAR  Free (sign-indefinite) polynomial matrix variable.
%   [PROG, P] = MVST_SOSMATRVAR(PROG, Z, [SP SQ]) declares an SP-by-SQ matrix P
%   whose entries are free (sign-indefinite) polynomials in the monomial basis Z
%   and appends it to the SOSTOOLS program PROG.  It is the free counterpart of
%   MVST_SOSPOSMATRVAR, used for the copositive certificate vectors
%   (lambda(tau), xi(tau)) and the free controller numerators (Uc(tau), Ud) of
%   the positive-systems copositive-SOS engines.
%
%   CLEAN PORT (attribution).  This reimplements the MatrixTools helper
%   `sosmatrvar.m` by M. Peet (MatrixTools; matthew.peet@inria.fr), distributed
%   with the L-infinity positive-systems scripts of the NAHS "Linf Positive
%   Systems" package.  The original fills the matrix entry-by-entry with
%   `sospolyvar(...,'wscoeff')`; this port uses the modern SOSTOOLS
%   `sospolymatrixvar`, which returns the whole polynomial matrix directly.
%
%   See also MVST_SOSPOSMATRVAR, MVST_POS_GAIN, MVST_POS_STABILIZE,
%   SOSPOLYMATRIXVAR.

    [prog, P] = sospolymatrixvar(prog, Z, sizePQ);
end
