function [passive, nu] = mvst_passivity_ct(A, Bw, C, Dw)
%MVST_PASSIVITY_CT  Passivity test of a square LTI system by a frequency sweep.
%   [PASSIVE,NU] = MVST_PASSIVITY_CT(A,BW,C,DW) tests passivity of the SQUARE
%   system (size(C,1)==size(Bw,2))
%        xdot = A x + Bw w,   z = C x + Dw w,
%   with A Hurwitz, by sweeping the transfer matrix
%        G(jw) = C ( jw I - A )^{-1} Bw + Dw
%   over a logspaced grid w in [1e-3,1e3] (plus w = 0 and a large w).
%
%   At each frequency the Hermitian part is formed and its smallest
%   (real) eigenvalue is m(w) = min(eig(G(jw)+G(jw)')).  The returned
%   passivity index is
%        NU = min_w m(w),
%   in the CONVENTION OF THE BOOK (Dissipativity chapter, Section
%   sec:diss:passivity: ISP supply He(w'z) - nu|w|^2, i.e. G+G' >= nu I).  A
%   POSITIVE nu certifies INPUT-STRICT passivity with index nu, a NEGATIVE nu
%   quantifies the SHORTAGE of passivity.  PASSIVE = (nu >= -1e-7) flags
%   (approximate) passivity.  CONVENTION CHANGE 2026-10-27: earlier versions
%   returned nu/2 (G+G' >= 2 nu I); the value now equals MVST_ISP_MARGIN_CT and
%   MVST_PASSIVITY.  The first output was formerly named ISP, a misnomer (it is
%   the passivity flag, not input-strict passivity).
%
%   This is the FREQUENCY-GRID reference; MVST_PASSIVITY gives the exact
%   (Hamiltonian) test and is the recommended entry point.
%
%   This is a frequency-domain certificate and uses no toolbox functions.
%
%   Example:
%       [passive,nu] = mvst_passivity_ct([-1 0;0 -2], [1;1], [1 1], 0.5);
%
%   See also MVST_PASSIVITY, MVST_ISP_MARGIN_CT, MVST_L2GAIN_CT, MVST_DISS_PT_CT.
if size(A,1) ~= size(A,2)
    error('mvst:passivity_ct:notSquare', 'A must be square.');
end
if size(Bw,1) ~= size(A,1)
    error('mvst:passivity_ct:dim', 'BW must have %d rows to match A.', size(A,1));
end
if size(C,2) ~= size(A,1)
    error('mvst:passivity_ct:dim', 'C must have %d columns to match A.', size(A,1));
end
n = size(A,1);
nu_in = size(Bw,2);
ny    = size(C,1);
assert(ny == nu_in, 'mvst:passivity_ct:square', ...
    'Passivity is defined for square systems (size(C,1)==size(Bw,2)).');

% Frequency grid: 0, a dense logspaced band, and a very large frequency.
w = [0, logspace(-3, 3, 4000), 1e6];

mmin = inf;
I = eye(n);
for k = 1:numel(w)
    G  = C*((1i*w(k)*I - A)\Bw) + Dw;       % ny x ny transfer matrix
    Gh = G + G';                            % Hermitian part (times 2)
    Gh = (Gh + Gh')/2;                      % enforce Hermitian numerically
    m  = min(real(eig(Gh)));
    if m < mmin; mmin = m; end
end

nu      = mmin;
passive = (nu >= -1e-7);
end
