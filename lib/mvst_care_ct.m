function [P, K, info] = mvst_care_ct(md, opts)
%MVST_CARE_CT  Infinite-horizon CT I-GRE (stochastic ARE) for one mode.
%   [P,K,info] = MVST_CARE_CT(MD) computes the stabilizing solution of the
%   algebraic I-GRE G(P)=0 (see MVST_RIC_CT) and the optimal gain K, by the
%   monotone-limit construction of the monograph: integrate the differential
%   I-GRE backward from P(T)=0 over a growing horizon until the residual is
%   negligible.  This needs only ODE integration (no toolboxes) and mirrors
%   the theory (P_T increases monotonically to P_inf).
%
%   OPTS (optional struct) fields: tol (default 1e-9), Tstart (5),
%   Tmax (400), reltol (1e-9).  INFO returns .residual, .horizon, .iters.
%
%   Example:
%       md = mvst_mode([-1 0.5; 0 -2], [1;1], eye(2), 1);
%       [P,K] = mvst_care_ct(md);
%
%   See also MVST_DRE_CT, MVST_RIC_CT, MVST_MODE, MVST_SOLVE_INFINITE.
if nargin < 2; opts = struct(); end
if ~isstruct(md) || ~isfield(md, 'A')
    error('mvst:care_ct:mode', ...
        'MD must be a mode struct from MVST_MODE (missing field .A).');
end
if size(md.A,1) ~= size(md.A,2)
    error('mvst:care_ct:notSquare', ...
        'MD.A must be square; got %d x %d.', size(md.A,1), size(md.A,2));
end
tol    = getf(opts,'tol',1e-9);
T      = getf(opts,'Tstart',5);
Tmax   = getf(opts,'Tmax',400);
reltol = getf(opts,'reltol',1e-9);

P = zeros(md.n);  res = inf;  it = 0;
while res > tol && T <= Tmax
    it = it + 1;
    [~,~,P] = mvst_dre_ct(md, T, zeros(md.n), reltol);
    res = norm(mvst_ric_ct(P, md), 'fro');
    T = T * 1.6;
end
[~, K] = mvst_ric_ct(P, md);
info = struct('residual', res, 'horizon', T/1.6, 'iters', it);
if res > tol
    warning('mvst:care_ct:converge', ...
        'I-GRE residual %.2e exceeds tol %.1e (horizon %.1f). Check stabilizability/detectability.', ...
        res, tol, T/1.6);
end
end

function v = getf(s,f,d)
if isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
