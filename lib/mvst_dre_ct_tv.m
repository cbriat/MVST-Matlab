function [t, Pt, Kt] = mvst_dre_ct_tv(md, tspan, Pterm, reltol)
%MVST_DRE_CT_TV  Backward time-varying differential I-GRE.
%   [t,Pt,Kt] = MVST_DRE_CT_TV(MD,TSPAN,PTERM) integrates the time-varying
%   differential I-GRE  Pdot = -G(t,P)  backward on TSPAN = [t0 T] from the
%   terminal condition P(T) = PTERM, where G is the time-varying I-GRE
%   residual (MVST_RIC_CT_TV).  Returns increasing times t, the trajectory
%   Pt(:,:,k) = P(t(k)), and the optimal gains Kt(:,:,k) = K(t(k)).
%   MVST_DRE_CT_TV(...,RELTOL) sets the ODE relative tolerance (default 1e-8).
%
%   See also MVST_RIC_CT_TV, MVST_DRE_CT, MVST_SOLVE_FINITE_TV, MVST_EVALC.
if ~isstruct(md) || ~isfield(md, 'A')
    error('mvst:dre_ct_tv:mode', ...
        'MD must be a mode struct from MVST_MODE (missing field .A).');
end
if numel(tspan) < 2
    error('mvst:dre_ct_tv:tspan', 'TSPAN must be [t0 T].');
end
if nargin < 4 || isempty(reltol); reltol = 1e-8; end
n = md.n;
if nargin < 3 || isempty(Pterm); Pterm = zeros(n); end
odef = @(tt,pv) reshape(-mvst_ric_ct_tv(tt, sym(pv,n), md), [], 1);
opts = odeset('RelTol', reltol, 'AbsTol', reltol*1e-2);
[ts, Pv] = ode45(odef, [tspan(2) tspan(1)], Pterm(:), opts);
[t, idx] = sort(ts);  Pv = Pv(idx,:);
Pt = zeros(n,n,numel(t));
Kt = zeros(md.m,n,numel(t));
for k = 1:numel(t)
    Pk = sym(Pv(k,:).', n);  Pt(:,:,k) = Pk;
    [~, Kt(:,:,k)] = mvst_ric_ct_tv(t(k), Pk, md);
end
end

function M = sym(v, n)
M = reshape(v, n, n);  M = (M + M.')/2;
end
