function res = mvst_simulate_mjls(sys, K, tspan, reltol)
%MVST_SIMULATE_MJLS  Closed-loop second-moment simulation for an MJLS.
%   RES = MVST_SIMULATE_MJLS(SYS,K,TSPAN) propagates the mode-indexed
%   second moments X_i(t) = E[x x' 1_{sigma=i}] of the closed-loop
%   Markov-jump diffusion under mode-dependent feedback u = K{i} x:
%
%     dX_i/dt = Acl_i X_i + X_i Acl_i' + sum_l Lc_il X_i Lc_il'
%               + sum_j pi_ji X_j + p_i(t) W_i,
%
%   with Acl_i = A_i + B_i K{i}, Lc_il = Ac_il + Bc_il K{i}.  The additive
%   process noise of mode i contributes only while the chain occupies that
%   mode, so it enters weighted by the mode-occupation probability p_i(t),
%   which satisfies dp/dt = Pi' p with p(0) = pi0.  At steady state p -> the
%   stationary distribution, so the running-average cost rate converges to
%   the analytic optimum rho* = sum_i pi_stat(i) trace(P_i W_i) returned by
%   MVST_SOLVE_AVERAGE.  Initial condition X_i(0) = pi0(i) * MD_i.X0.  The
%   total state covariance is E[x x'] = sum_i X_i.  RES has fields .t, .X
%   (1xM cell of n x n x nt), .Xtot (n x n x nt), .p (nt x M), .inst,
%   .runavg, .cost.
%
%   See also MVST_SIMULATE_CT, MVST_SAMPLE_MJLS, MVST_SOLVE_AVERAGE.
if nargin < 4 || isempty(reltol); reltol = 1e-8; end
M = numel(sys.modes); n = sys.n; Pi = sys.Pi; pi0 = sys.pi0;

Acl = cell(1,M); Lc = cell(1,M); L = cell(1,M);
for i = 1:M
    md = sys.modes{i};
    Acl{i} = md.A + md.B*K{i};
    Lc{i}  = zeros(n,n,md.Nc);
    for l = 1:md.Nc; Lc{i}(:,:,l) = md.Ac(:,:,l) + md.Bc(:,:,l)*K{i}; end
    Li = md.Q + md.S*K{i} + (md.S*K{i}).' + K{i}.'*md.R*K{i};
    L{i} = (Li+Li.')/2;
end

np = M*n*n;                                  % length of the stacked moments
y0 = zeros(np+M+1,1);
for i = 1:M
    Xi0 = pi0(i)*sys.modes{i}.X0;
    y0((i-1)*n*n+1:i*n*n) = Xi0(:);
end
y0(np+1:np+M) = pi0(:);                      % mode-occupation probabilities p_i(t)
odef = @(t,y) derivM(y, sys, Acl, Lc, L, Pi, n, M, np);
opts = odeset('RelTol',reltol,'AbsTol',reltol*1e-2);
[t, Y] = ode45(odef, tspan, y0, opts);

nt = numel(t);
X = repmat({zeros(n,n,nt)}, 1, M); Xtot = zeros(n,n,nt); inst = zeros(nt,1);
for k = 1:nt
    c = 0;
    for i = 1:M
        Xi = reshape(Y(k,(i-1)*n*n+1:i*n*n), n, n); Xi = (Xi+Xi.')/2;
        X{i}(:,:,k) = Xi; Xtot(:,:,k) = Xtot(:,:,k) + Xi;
        c = c + trace(L{i}*Xi);
    end
    inst(k) = c;
end
g = Y(:,end);
runavg = g ./ max(t - t(1), eps); runavg(1) = inst(1);
% Blank the burn-in: dividing by a near-zero elapsed time makes the first
% samples spike; NaNs are skipped by plot so the curve starts cleanly.
runavg(t - t(1) < 0.02*(t(end) - t(1))) = NaN;
res.t = t; res.X = X; res.Xtot = Xtot; res.p = Y(:,np+1:np+M);
res.inst = inst; res.runavg = runavg; res.cost = g;
end

function dy = derivM(y, sys, Acl, Lc, L, Pi, n, M, np)
X = cell(1,M);
for i = 1:M; Xi = reshape(y((i-1)*n*n+1:i*n*n), n, n); X{i} = (Xi+Xi.')/2; end
p = y(np+1:np+M);
dy = zeros(np+M+1,1); c = 0;
for i = 1:M
    md = sys.modes{i};
    dX = Acl{i}*X{i} + X{i}*Acl{i}.' + p(i)*md.W;   % noise weighted by occupation prob.
    for l = 1:md.Nc; dX = dX + Lc{i}(:,:,l)*X{i}*Lc{i}(:,:,l).'; end
    for j = 1:M; dX = dX + Pi(j,i)*X{j}; end          % Markov coupling (adjoint)
    dy((i-1)*n*n+1:i*n*n) = dX(:);
    c = c + trace(L{i}*X{i});
end
dy(np+1:np+M) = Pi.'*p;                                % mode probabilities: dp/dt = Pi' p
dy(end) = c;
end
