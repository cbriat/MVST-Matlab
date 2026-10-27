function res = mvst_solve_finite_tv(sys, T, opts)
%MVST_SOLVE_FINITE_TV  Finite-horizon linear-cost optimal control, TIME-VARYING.
%   RES = MVST_SOLVE_FINITE_TV(SYS,T) is the time-varying counterpart of
%   MVST_SOLVE_FINITE: each mode's coefficients may be function handles of
%   time (evaluated through MVST_RIC_CT_TV).  It integrates the (coupled)
%   differential I-GRE backward from P(T)=ZT on [0,T].
%
%   Single mode (M=1) -> RES has fields .t, .P (n x n x nt), .K (m x n x nt),
%       .P0=P(0), .K0=K(0), .cost=trace(P0*X0(0)).
%   Multi-mode (M>1) -> RES.P,.K are 1xM cells of 3-D arrays, RES.P0/.K0 are
%       1xM cells, RES.cost = sum_i pi0(i) tr(P0_i X0_i(0)).
%
%   See also MVST_SOLVE_FINITE, MVST_SOLVE_INFINITE_TV, MVST_DRE_CT_TV, MVST_RIC_CT_TV.
if nargin < 3; opts = struct(); end
if ~isstruct(sys) || ~isfield(sys, 'modes')
    error('mvst:solve_finite_tv:sys', ...
        'SYS must be a thin system struct with a .modes cell (see MVST_SYSTEM).');
end
if ~(isscalar(T) && isreal(T) && T > 0)
    error('mvst:solve_finite_tv:horizon', 'T must be a positive scalar.');
end
M  = numel(sys.modes);
n  = sys.modes{1}.n;
rt = mvst_getf(opts,'reltol',1e-8);
ZT = mvst_getf(sys,'ZT',zeros(n));

if M == 1
    md = sys.modes{1};
    [t, Pt, Kt] = mvst_dre_ct_tv(md, [0 T], ZT, rt);
    res.t = t; res.P = Pt; res.K = Kt; res.P0 = Pt(:,:,1); res.K0 = Kt(:,:,1);
    res.cost = trace(res.P0 * mvst_evalc(md.X0,0));
    return;
end

Pi  = mvst_getf(sys,'Pi',zeros(M));
pi0 = mvst_getf(sys,'pi0',[1 zeros(1,M-1)]);
oo  = odeset('RelTol',rt,'AbsTol',rt*1e-2);
odef = @(tt,pv) coupled_rhs_tv(tt, pv, sys, Pi, n, M);
[ts, Y] = ode45(odef, [T 0], repmat(ZT(:), M, 1), oo);
[t, idx] = sort(ts); Y = Y(idx,:); nt = numel(t);
P = cell(1,M); K = cell(1,M); P0 = cell(1,M); K0 = cell(1,M);
for i = 1:M; P{i} = zeros(n,n,nt); K{i} = zeros(sys.m,n,nt); end
for k = 1:nt
    for i = 1:M
        Pk = reshape(Y(k,(i-1)*n*n+1:i*n*n), n, n); Pk = (Pk+Pk.')/2;
        P{i}(:,:,k) = Pk;
        [~, K{i}(:,:,k)] = mvst_ric_ct_tv(t(k), Pk, sys.modes{i});
    end
end
cost = 0;
for i = 1:M
    P0{i} = P{i}(:,:,1); K0{i} = K{i}(:,:,1);
    cost = cost + pi0(i)*trace(P0{i}*mvst_evalc(sys.modes{i}.X0,0));
end
res.t = t; res.P = P; res.K = K; res.P0 = P0; res.K0 = K0; res.cost = cost;
end

function dpv = coupled_rhs_tv(tt, pv, sys, Pi, n, M)
P = cell(1,M);
for i = 1:M; Pi_ = reshape(pv((i-1)*n*n+1:i*n*n), n, n); P{i} = (Pi_+Pi_.')/2; end
dpv = zeros(M*n*n,1);
for i = 1:M
    C = zeros(n);
    for j = 1:M; C = C + Pi(i,j)*P{j}; end
    G = mvst_ric_ct_tv(tt, P{i}, sys.modes{i}, C);
    dpv((i-1)*n*n+1:i*n*n) = reshape(-G, [], 1);
end
end
