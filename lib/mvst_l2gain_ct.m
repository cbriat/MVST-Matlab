function gam = mvst_l2gain_ct(A, Bw, C, Dw, Cc, Dwc)
%MVST_L2GAIN_CT  Continuous-time L2-gain (H-infinity norm) by bisection.
%   GAM = MVST_L2GAIN_CT(A,BW,C,DW) returns the L2-gain (i.e. the
%   H-infinity norm) of the LTI system
%        xdot = A x + Bw w,   z = C x + Dw w,
%   with A Hurwitz.  No Control System Toolbox is used: the test
%   "gain <= gamma" is the standard Hamiltonian imaginary-eigenvalue
%   test of the bounded-real lemma.
%
%   For a trial gamma > norm(Dw), with R = gamma^2 I - Dw'Dw (pos. def.),
%        H = [ A + Bw R^{-1} Dw' C ,            Bw R^{-1} Bw'        ;
%             -(C'C + C' Dw R^{-1} Dw' C) ,  -(A + Bw R^{-1} Dw' C)' ].
%   The gain is <= gamma iff H has NO eigenvalue on the imaginary axis,
%   i.e. no eigenvalue with |real part| < tol.  We bisect gamma on
%   [norm(Dw)+1e-6, hi] (growing hi until feasible) to tolerance 1e-4.
%
%   Works uniformly for SISO and MIMO via matrix operations.
%
%   GAM = MVST_L2GAIN_CT(A,BW,C,DW,CC,DWC) additionally accepts the
%   state-multiplicative and disturbance-multiplicative Brownian channels of
%   the diffusion system
%        dx = (A x + Bw w)dt + sum_l (Cc{l} x + Dwc{l} w) dW_l,   z = C x + Dw w,
%   as 1 x q cells CC (state noise) and DWC (disturbance noise; empty entries
%   or a shorter cell are padded with zeros).  When either cell is non-empty the
%   MEAN-SQUARE L2-gain is returned through the exact stochastic bounded-real
%   shooting Riccati (MVST_L2GAIN_CT_POISSON with no jumps): the deterministic
%   Hamiltonian test above no longer applies once multiplicative noise is
%   present.  With CC and DWC absent or empty the deterministic path is taken
%   unchanged (identical output to the four-argument call).
%
%   Example:
%       gam = mvst_l2gain_ct([-1 1; 0 -2], [1;0], [1 0], 0);
%       gam = mvst_l2gain_ct([0 1;-2 -1], [0;1], [1 0], 0.5, {0.15*eye(2)});
%
%   See also MVST_L2GAIN, MVST_L2GAIN_CT_POISSON, MVST_L2GAIN_CT_MJLS,
%   MVST_DISS_PT_CT, MVST_PASSIVITY_CT.
if nargin < 5; Cc  = {}; end
if nargin < 6; Dwc = {}; end
if ~iscell(Cc);  Cc  = tocell3(Cc);  end
if ~iscell(Dwc); Dwc = tocell3(Dwc); end
hasNoise = anynz(Cc) || anynz(Dwc);   % all-zero channels stay on the exact path
if hasNoise
    % Multiplicative noise present: route to the exact mean-square backend.
    % (Same Ci' Si Ci fold as the Poisson class, with no spontaneous jumps.)
    if size(A,1) ~= size(A,2)
        error('mvst:l2gain_ct:notSquare', 'A must be square; got %d x %d.', size(A,1), size(A,2));
    end
    n = size(A,1); m = size(Bw,2);
    if isempty(Dwc); Dwc = repmat({zeros(n,m)}, 1, numel(Cc)); end
    if numel(Dwc) < numel(Cc); Dwc(end+1:numel(Cc)) = {zeros(n,m)}; end
    % Pad the state channels to the disturbance channels as well (fix
    % 2026-10-04): with CC = {} and DWC nonempty the backend loops over
    % numel(CC) and silently dropped every disturbance-multiplicative channel.
    if numel(Cc) < numel(Dwc); Cc(end+1:numel(Dwc)) = {zeros(n)}; end
    if isempty(Dw) || (isscalar(Dw) && Dw == 0); Dw = zeros(size(C,1), m); end
    gam = mvst_l2gain_ct_poisson(A, Bw, C, Dw, Cc, Dwc, {}, [], {});
    return;
end
if size(A,1) ~= size(A,2)
    error('mvst:l2gain_ct:notSquare', 'A must be square; got %d x %d.', size(A,1), size(A,2));
end
if size(Bw,1) ~= size(A,1)
    error('mvst:l2gain_ct:dim', 'BW must have %d rows to match A.', size(A,1));
end
if size(C,2) ~= size(A,1)
    error('mvst:l2gain_ct:dim', 'C must have %d columns to match A.', size(A,1));
end
% The Hamiltonian bounded-real test below is only valid for a Hurwitz A: the
% H-infinity norm of a non-Hurwitz system is infinite, and the imaginary-axis
% eigenvalue test would silently return a finite garbage value.  Reject such
% input explicitly (behaviour change on ill-posed input only; every legitimate
% caller feeds a stable A, as the synthesis validation gate enforces upstream).
absc = max(real(eig(A)));
if absc >= 0
    error('mvst:l2gain_ct:unstable', ...
        ['A must be Hurwitz for a finite L2-gain; spectral abscissa is ' ...
         '%.4g >= 0 (the H-infinity norm is infinite).'], absc);
end
nz  = size(C,1);    % number of performance outputs z (rows of C, Dw)
nw  = size(Bw,2);   % number of disturbance inputs w (cols of Bw, Dw)

% Normalise the feedthrough Dw to an nz x nw matrix.  Callers routinely pass
% the scalar 0 as shorthand for "no feedthrough"; expand it here so the
% bounded-real R block below is sized by the disturbance count nw (NOT nz).
if isempty(Dw) || (isscalar(Dw) && Dw == 0)
    Dw = zeros(nz, nw);
end
if ~isequal(size(Dw), [nz nw])
    error('mvst:l2gain_ct:dim', ...
        'DW must be %d x %d (rows of C by cols of BW); got %d x %d.', ...
        nz, nw, size(Dw,1), size(Dw,2));
end
tolImag = 1e-7;     % threshold on |Re(eig(H))| for "on the imaginary axis"
tolBis  = 1e-4;     % bisection tolerance on gamma

dnorm = norm(Dw);   % largest singular value of the feedthrough

% Feasibility test: is the gain <= gamma ?
    function ok = feasible(g)
        if g <= dnorm
            ok = false; return;
        end
        R = g^2*eye(nw) - Dw.'*Dw;        % nw x nw, must be pos. def.
        [~,pchol] = chol((R+R.')/2);
        if pchol ~= 0
            ok = false; return;
        end
        RiDtC = R\(Dw.'*C);               % R^{-1} Dw' C   (nz x n)
        Acl   = A + Bw*RiDtC;             % A + Bw R^{-1} Dw' C
        H = [ Acl,                         Bw*(R\Bw.') ; ...
             -(C.'*C + C.'*Dw*RiDtC),     -Acl.' ];
        ev = eig(H);
        ok = all(abs(real(ev)) >= tolImag);
    end

% Lower bound just above the feedthrough norm.
lo = dnorm + 1e-6;

% Grow an upper bound until it is feasible.
hi = max(lo*2, 1.0);
maxGrow = 60;
for kk = 1:maxGrow
    if feasible(hi); break; end
    hi = hi*2;
    if kk == maxGrow
        warning('mvst:l2gain_ct:nohi', ...
            'Could not bracket a feasible gamma; returning upper bound.');
        gam = hi; return;
    end
end

% Bisection: keep hi feasible, lo infeasible.
while (hi - lo) > tolBis
    mid = 0.5*(lo + hi);
    if feasible(mid)
        hi = mid;
    else
        lo = mid;
    end
end
gam = hi;
end

% ------------------------------------------------------------------------
function c = tocell3(X)
%TOCELL3  Slice a 3-D array along its 3rd dimension into a row cell.
if iscell(X)
    c = X;
elseif isempty(X)
    c = {};
else
    p = size(X, 3);
    c = cell(1, p);
    for k = 1:p; c{k} = X(:,:,k); end
end
end

% ------------------------------------------------------------------------
function tf = anynz(C)
%ANYNZ  True if a channel cell carries any structurally nonzero matrix.
%   An explicit all-zero channel is treated as absent, so it stays on the
%   deterministic exact path rather than switching to the shooting Riccati.
if ~iscell(C); tf = ~isempty(C) && nnz(C) > 0; return; end
tf = any(cellfun(@(x) ~isempty(x) && nnz(x) > 0, C));
end
