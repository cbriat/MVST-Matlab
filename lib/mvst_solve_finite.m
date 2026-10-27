function res = mvst_solve_finite(sys, T, opts)
%MVST_SOLVE_FINITE  Finite-horizon linear-cost optimal control (CT).
%   RES = MVST_SOLVE_FINITE(SYS,T) solves the finite-horizon problem on
%   [0,T] by integrating the (coupled) differential I-GRE backward from
%   P(T)=ZT.  Works for a single mode and for a Markov-jump system
%   (SYS.Pi the generator coupling the modes).
%
%   Single mode (M=1) -> RES has fields:
%       .t  (times),  .P (n x n x nt),  .K (m x n x nt),
%       .P0 = P(0),   .K0 = K(0),       .cost = trace(P0*X0).
%   Multi-mode (M>1) -> RES.P and RES.K are 1xM cell arrays of such 3-D
%   arrays, RES.P0/.K0 are 1xM cells, and RES.cost = sum_i pi0(i) tr(P0_i X0_i).
%
%   Example:
%       S = mvst.System.ct([-1 0;0 -2], [1;1], 'Q', eye(2), 'R', 1);
%       r = mvst_solve_finite(S.thin(), 3);   r.cost
%
%   See also MVST_SOLVE_INFINITE, MVST_SOLVE_FINITE_TV, MVST_DRE_CT, MVST_RIC_CT.
if nargin < 3; opts = struct(); end
if ~isstruct(sys) || ~isfield(sys, 'modes')
    error('mvst:solve_finite:sys', ...
        'SYS must be a thin system struct with a .modes cell (see MVST_SYSTEM).');
end
if ~(isscalar(T) && isreal(T) && T > 0)
    error('mvst:solve_finite:horizon', 'T must be a positive scalar.');
end
M  = numel(sys.modes);
n  = sys.modes{1}.n;
rt = mvst_getf(opts,'reltol',1e-8);
ZT = mvst_getf(sys,'ZT',zeros(n));

if M == 1
    md = sys.modes{1};
    [t, Pt, P0, Kt] = mvst_dre_ct(md, T, ZT, rt);
    res.t = t; res.P = Pt; res.K = Kt; res.P0 = P0; res.K0 = Kt(:,:,1);
    res.cost = trace(P0 * md.X0);
    return;
end

% --- coupled Markov-jump finite horizon ---
Pi  = mvst_getf(sys,'Pi',zeros(M));
pi0 = mvst_getf(sys,'pi0',[1 zeros(1,M-1)]);
odefun = @(tt,pv) coupled_rhs(pv, sys, Pi, n, M);
% Norm-cap event: stop cleanly on finite escape (bad terminal multiplier
% during a steering Newton line search) instead of raising step-size
% warnings; the truncated trajectory yields a rejected Newton step.
oo = odeset('RelTol',rt,'AbsTol',rt*1e-2, ...
            'Events',@(tt,pv) mvst_odecap(pv, 1e8));
y0 = repmat(ZT(:), M, 1);
[ts, Y] = ode45(odefun, [T 0], y0, oo);
[t, idx] = sort(ts); Y = Y(idx,:);
nt = numel(t);
P = cell(1,M); K = cell(1,M); P0 = cell(1,M); K0 = cell(1,M);
for i = 1:M; P{i} = zeros(n,n,nt); K{i} = zeros(sys.m,n,nt); end
for k = 1:nt
    for i = 1:M
        Pk = reshape(Y(k,(i-1)*n*n+1:i*n*n), n, n); Pk = (Pk+Pk.')/2;
        P{i}(:,:,k) = Pk;
        [~, K{i}(:,:,k)] = mvst_ric_ct(Pk, sys.modes{i});
    end
end
cost = 0;
for i = 1:M
    P0{i} = P{i}(:,:,1); K0{i} = K{i}(:,:,1);
    cost = cost + pi0(i)*trace(P0{i}*sys.modes{i}.X0);
end
res.t = t; res.P = P; res.K = K; res.P0 = P0; res.K0 = K0; res.cost = cost;
end

function dpv = coupled_rhs(pv, sys, Pi, n, M)
P = cell(1,M);
for i = 1:M; Pi_ = reshape(pv((i-1)*n*n+1:i*n*n), n, n); P{i} = (Pi_+Pi_.')/2; end
dpv = zeros(M*n*n,1);
for i = 1:M
    C = zeros(n);
    for j = 1:M; C = C + Pi(i,j)*P{j}; end
    G = mvst_ric_ct(P{i}, sys.modes{i}, C);
    dpv((i-1)*n*n+1:i*n*n) = reshape(-G, [], 1);
end
end
