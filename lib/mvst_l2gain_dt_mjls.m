function [gam, S] = mvst_l2gain_dt_mjls(Ac, Bw, Cc, Dw, Ptr, N, rtol)
%MVST_L2GAIN_DT_MJLS  Coupled mean-square l2-gain of a discrete-time MJLS.
%   GAM = MVST_L2GAIN_DT_MJLS(AC,BW,CC,DW,PTR) returns the mean-square
%   l2-gain of the Markov jump linear system
%        x+ = A_i x + Bw w,   z = C_i x + Dw w,   mode i -> j w.p. PTR(i,j),
%   with 1 x M cells AC, CC of mode matrices, common input Bw and feedthrough
%   Dw, and a row-stochastic M x M transition matrix PTR.  The gain is the
%   smallest gamma for which the COUPLED bounded-real Riccati recursion
%        E_i = sum_j PTR(i,j) S_j,   Th_i = gamma^2 I - Dw'Dw - Bw'E_i Bw,
%        S_i <- A_i'E_i A_i + C_i'C_i + L_i Th_i^{-1} L_i',
%        L_i = A_i'E_i Bw + C_i'Dw,
%   started at S_i = 0, keeps every pivot Th_i positive definite and stays
%   bounded (it then converges monotonically to a fixed point).  The test is
%   decided by iterating at most N steps (default 5000) with early exit on
%   convergence; gamma is bisected to the relative tolerance RTOL (default
%   1e-4).  Like the CT shooting backend it can under-report the gain by an
%   amount that vanishes as N grows.
%
%   [GAM,S] = MVST_L2GAIN_DT_MJLS(...) also returns the 1 x M cell S of
%   fixed-point storage matrices computed at 1.001*GAM (the certificate:
%   V(x,i) = x'S_i x, S_i >= 0), or {} when that level is not certified.
%
%   Example:
%       gam = mvst_l2gain_dt_mjls({0.5,0.4}, 1, {1,1}, 0, [0.9 0.1;0.2 0.8]);
%
%   See also MVST_L2GAIN, MVST_L2GAIN_DT, MVST_L2GAIN_CT_MJLS.
if ~iscell(Ac) || ~iscell(Cc)
    error('mvst:l2gain_dt_mjls:modes', 'AC and CC must be 1 x M cells.');
end
M = numel(Ac);
if ~isequal(size(Ptr), [M M]) || any(Ptr(:) < -1e-12) || any(abs(sum(Ptr,2)-1) > 1e-8)
    error('mvst:l2gain_dt_mjls:transition', ...
        'PTR must be a row-stochastic %d x %d matrix.', M, M);
end
if nargin < 6 || isempty(N);    N = 5000;    end
if nargin < 7 || isempty(rtol); rtol = 1e-4; end
n = size(Ac{1},1); m = size(Bw,2); nz = size(Cc{1},1);
if isempty(Dw) || (isscalar(Dw) && Dw == 0); Dw = zeros(nz, m); end
dnorm = norm(Dw);
lo = dnorm + 1e-6; hi = max(2*lo, 1);
for kk = 1:60
    if feasible(hi); break; end
    hi = 2*hi;
    if kk == 60
        warning('mvst:l2gain_dt_mjls:nohi', 'No feasible gamma bracketed.');
        gam = hi; S = {}; return;
    end
end
while (hi - lo) > rtol*max(1, hi)
    mid = 0.5*(lo + hi);
    if feasible(mid); hi = mid; else; lo = mid; end
end
gam = hi;
if nargout > 1
    [ok, S] = feasible(1.001*gam);
    if ~ok; S = {}; end
end

    function [ok, Sf] = feasible(g)
        Sk = repmat({zeros(n)}, 1, M); ok = true;
        for k = 1:N
            Sn = cell(1, M); dmax = 0; smax = 0;
            for i = 1:M
                E = zeros(n);
                for j = 1:M; E = E + Ptr(i,j)*Sk{j}; end
                Th = g^2*eye(m) - Dw.'*Dw - Bw.'*E*Bw;
                Th = (Th + Th.')/2;
                [~, pc] = chol(Th);
                if pc ~= 0; ok = false; Sf = {}; return; end
                L = Ac{i}.'*E*Bw + Cc{i}.'*Dw;
                Si = Ac{i}.'*E*Ac{i} + Cc{i}.'*Cc{i} + L*(Th\L.');
                Sn{i} = (Si + Si.')/2;
                dmax = max(dmax, norm(Sn{i} - Sk{i}));
                smax = max(smax, norm(Sn{i}));
            end
            Sk = Sn;
            if smax > 1e8; ok = false; Sf = {}; return; end
            if dmax <= 1e-11*(1 + smax); break; end
        end
        Sf = Sk;
    end
end
