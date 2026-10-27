function res = mvst_solve_infinite_tv(sys, Tinf, opts)
%MVST_SOLVE_INFINITE_TV  Infinite-horizon dual variable, TIME-VARYING.
%   RES = MVST_SOLVE_INFINITE_TV(SYS,TINF) returns the infinite-horizon dual
%   P_inf(t) of a time-varying system whose coefficients settle as t grows
%   (e.g. A(t)=A0+A1*exp(-t) -> A0).  The construction:
%     1. freeze the coefficients at the large time TINF (where the system is
%        essentially LTI) and solve the LTI infinite-horizon I-GRE -> P_LTI;
%     2. integrate the time-varying differential I-GRE backward from TINF
%        seeded with the terminal value P_LTI.
%   The result P_inf(t) on [0,TINF] is the genuine infinite-horizon dual of
%   the time-varying system; it settles to the asymptotic value P_LTI.
%   Seeding the terminal with P_LTI (rather than 0) removes the
%   finite-horizon terminal transient, so P_inf(t) is meaningful up to TINF.
%
%   RES has fields .t, .P (single mode: n x n x nt; multi-mode: 1xM cell),
%   .K (same shape), .Plti (1xM cell, the asymptotic LTI dual), .Klti.
%
%   See also MVST_SOLVE_INFINITE, MVST_SOLVE_FINITE_TV, MVST_DRE_CT_TV.
if nargin < 3; opts = struct(); end
if ~isstruct(sys) || ~isfield(sys, 'modes')
    error('mvst:solve_infinite_tv:sys', ...
        'SYS must be a thin system struct with a .modes cell (see MVST_SYSTEM).');
end
if ~(isscalar(Tinf) && isreal(Tinf) && Tinf > 0)
    error('mvst:solve_infinite_tv:horizon', 'TINF must be a positive scalar.');
end
M  = numel(sys.modes);
n  = sys.modes{1}.n;
rt = mvst_getf(opts,'reltol',1e-8);

% --- asymptotic LTI system: freeze coefficients at t = Tinf ---
sys0 = sys;
for i = 1:M; sys0.modes{i} = snapshot(sys.modes{i}, Tinf); end
ri   = mvst_solve_infinite(sys0, opts);
Plti = ri.P;  Klti = ri.K;

oo = odeset('RelTol',rt,'AbsTol',rt*1e-2);
if M == 1
    [t, Pt, Kt] = mvst_dre_ct_tv(sys.modes{1}, [0 Tinf], Plti{1}, rt);
    res.t = t; res.P = Pt; res.K = Kt; res.Plti = Plti; res.Klti = Klti;
    return;
end

Pi = mvst_getf(sys,'Pi',zeros(M));
odef = @(tt,pv) coupled_rhs_tv(tt, pv, sys, Pi, n, M);
y0 = zeros(M*n*n,1);
for i = 1:M; y0((i-1)*n*n+1:i*n*n) = Plti{i}(:); end
[ts, Y] = ode45(odef, [Tinf 0], y0, oo);
[t, idx] = sort(ts); Y = Y(idx,:); nt = numel(t);
P = cell(1,M); K = cell(1,M);
for i = 1:M; P{i} = zeros(n,n,nt); K{i} = zeros(sys.m,n,nt); end
for k = 1:nt
    for i = 1:M
        Pk = reshape(Y(k,(i-1)*n*n+1:i*n*n), n, n); Pk = (Pk+Pk.')/2;
        P{i}(:,:,k) = Pk;
        [~, K{i}(:,:,k)] = mvst_ric_ct_tv(t(k), Pk, sys.modes{i});
    end
end
res.t = t; res.P = P; res.K = K; res.Plti = Plti; res.Klti = Klti;
end

function md0 = snapshot(md, t)
%SNAPSHOT  Freeze a (possibly time-varying) mode at time t into an LTI mode.
md0 = md;
md0.A = mvst_evalc(md.A,t);  md0.B = mvst_evalc(md.B,t);
md0.Q = (mvst_evalc(md.Q,t)+mvst_evalc(md.Q,t).')/2;
md0.R = (mvst_evalc(md.R,t)+mvst_evalc(md.R,t).')/2;
md0.S = mvst_evalc(md.S,t);
md0.Ac = mvst_evalc(md.Ac,t);  md0.Bc = mvst_evalc(md.Bc,t);
md0.Nc = size(md0.Ac,3);
md0.X0 = mvst_evalc(md.X0,t);
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
