function gam = mvst_l2gain_dt(A, Bw, C, Dw)
%MVST_L2GAIN_DT  Discrete-time L2-gain (H-infinity norm) by exact bisection.
%   GAM = MVST_L2GAIN_DT(A,BW,C,DW) returns the l2-gain (the H-infinity norm
%   over the unit circle) of the LTI system
%        x+ = A x + Bw w,   z = C x + Dw w,
%   with A Schur.  No Control System Toolbox is used.
%
%   METHOD (exact).  The bilinear map z = (1+s)/(1-s) sends the unit circle
%   onto the imaginary axis and preserves the H-infinity norm; the equivalent
%   continuous-time system is (A Schur implies -1 is not an eigenvalue of A)
%        Ac = (A+I)^{-1}(A-I),  Bc = (A+I)^{-1} Bw,
%        Cc = 2 C (A+I)^{-1},   Dc = Dw - C (A+I)^{-1} Bw,
%   with Ac Hurwitz, and its gain is obtained by the imaginary-axis Hamiltonian
%   bisection of MVST_L2GAIN_CT.  This is the Hamiltonian / symplectic-pencil
%   test of the discrete bounded-real lemma written in the s-plane; the result
%   is exact up to the CT bisection tolerance (1e-4 absolute).
%
%   Example:
%       gam = mvst_l2gain_dt(0.5, 1, 1, 0);     % = 1/(1-0.5) = 2
%
%   See also MVST_L2GAIN, MVST_L2GAIN_CT, MVST_L2GAIN_DT_MJLS, MVST_HINFNORM.
n = size(A,1);
if size(A,2) ~= n
    error('mvst:l2gain_dt:notSquare', 'A must be square; got %d x %d.', n, size(A,2));
end
if size(Bw,1) ~= n
    error('mvst:l2gain_dt:dim', 'BW must have %d rows to match A.', n);
end
if size(C,2) ~= n
    error('mvst:l2gain_dt:dim', 'C must have %d columns to match A.', n);
end
nz = size(C,1); nw = size(Bw,2);
if nargin < 4 || isempty(Dw) || (isscalar(Dw) && Dw == 0); Dw = zeros(nz, nw); end
if ~isequal(size(Dw), [nz nw])
    error('mvst:l2gain_dt:dim', 'DW must be %d x %d; got %d x %d.', nz, nw, size(Dw,1), size(Dw,2));
end
rho = max(abs(eig(A)));
if rho >= 1
    error('mvst:l2gain_dt:unstable', ...
        ['A must be Schur for a finite l2-gain; spectral radius is %.4g >= 1 ' ...
         '(the H-infinity norm is infinite).'], rho);
end
IA = eye(n) + A;
Ac = IA \ (A - eye(n));
Bc = IA \ Bw;
Cc = 2*(C / IA);
Dc = Dw - C*(IA \ Bw);
gam = mvst_l2gain_ct(Ac, Bc, Cc, Dc);
end
