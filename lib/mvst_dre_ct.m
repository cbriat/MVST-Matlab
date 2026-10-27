function [t, Pt, P0, Kt] = mvst_dre_ct(md, T, Pterm, reltol)
%MVST_DRE_CT  Backward solution of the continuous-time differential I-GRE.
%   [t,Pt,P0,Kt] = MVST_DRE_CT(MD,T,PTERM) integrates Pdot = -G(P) backward
%   on [0,T] from P(T)=PTERM, where G is the I-GRE residual (MVST_RIC_CT).
%   Returns increasing times t, the trajectory Pt(:,:,k)=P(t(k)), the
%   horizon-T dual variable P0=P(0), and the optimal gains Kt(:,:,k).
%   MVST_DRE_CT(...,RELTOL) sets the ODE relative tolerance (default 1e-8).
%
%   Example:
%       md = mvst_mode([-1 0.5; 0 -2], [1;1], eye(2), 1);
%       [t,Pt,P0] = mvst_dre_ct(md, 5);
%
%   See also MVST_RIC_CT, MVST_CARE_CT, MVST_DRE_CT_TV, MVST_SOLVE_FINITE.
if ~isstruct(md) || ~isfield(md, 'A')
    error('mvst:dre_ct:mode', ...
        'MD must be a mode struct from MVST_MODE (missing field .A).');
end
if nargin < 4 || isempty(reltol); reltol = 1e-8; end
if nargin < 3 || isempty(Pterm);  Pterm  = zeros(md.n); end
n = md.n;
if ~isequal(size(Pterm), [n n])
    error('mvst:dre_ct:dim', 'PTERM must be %d x %d to match the mode.', n, n);
end
odefun = @(tt,pv) reshape(-mvst_ric_ct(reshape(pv,n,n), md), [], 1);
% The norm-cap event stops cleanly on finite escape (e.g. an indefinite
% terminal condition tried by a steering Newton step) instead of grinding
% the step size down and spamming integration-tolerance warnings; callers
% see a truncated, clearly wrong trajectory and reject the step.
opts = odeset('RelTol', reltol, 'AbsTol', reltol*1e-2, ...
              'Events', @(tt,pv) mvst_odecap(pv, 1e8));
[ts, Pv] = ode45(odefun, [T 0], Pterm(:), opts);
[t, idx] = sort(ts);  Pv = Pv(idx,:);
Pt = zeros(n,n,numel(t));
Kt = zeros(md.m,n,numel(t));
for k = 1:numel(t)
    Pk = reshape(Pv(k,:), n, n);  Pk = (Pk+Pk.')/2;
    Pt(:,:,k) = Pk;
    [~, Kt(:,:,k)] = mvst_ric_ct(Pk, md);
end
P0 = Pt(:,:,1);
end
