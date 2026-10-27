function res = mvst_simulate_ct(md, K, tspan, reltol)
%MVST_SIMULATE_CT  Closed-loop second-moment (covariance) simulation (CT).
%   RES = MVST_SIMULATE_CT(MD,K,TSPAN) propagates the closed-loop state
%   covariance Sigma(t)=E[x x'] under the feedback u = K x for the diffusion
%   of mode MD with additive process noise MD.W:
%
%     dSigma/dt = Acl Sigma + Sigma Acl' + sum_l Lc_l Sigma Lc_l' + W,
%     Acl = A + B K,   Lc_l = Ac_l + Bc_l K,   Sigma(0) = MD.X0,
%
%   and accumulates the cost rate tr(L Sigma) with
%   L = Q + S K + (S K)' + K' R K.  TSPAN = [t0 tf] (or a time vector).
%   RES has fields .t, .Sigma (n x n x numel(t)), .inst (instantaneous cost
%   rate), .runavg (running average cost rate), .cost (accumulated cost).
%
%   See also MVST_SIMULATE_MJLS, MVST_SIMULATE_IMPULSIVE, MVST_SAMPLE_CT.
if nargin < 4 || isempty(reltol); reltol = 1e-8; end
n = md.n;
Acl = md.A + md.B*K;
Lc  = zeros(n,n,md.Nc);
for l = 1:md.Nc; Lc(:,:,l) = md.Ac(:,:,l) + md.Bc(:,:,l)*K; end
W = md.W;
L = md.Q + md.S*K + (md.S*K).' + K.'*md.R*K;  L = (L+L.')/2;

odef = @(t,y) deriv(y, Acl, Lc, W, L, n, md.Nc);
y0 = [md.X0(:); 0];
opts = odeset('RelTol',reltol,'AbsTol',reltol*1e-2);
[t, Y] = ode45(odef, tspan, y0, opts);

K1 = numel(t);
Sigma = zeros(n,n,K1); inst = zeros(K1,1);
for k = 1:K1
    Sg = reshape(Y(k,1:n*n), n, n); Sg = (Sg+Sg.')/2;
    Sigma(:,:,k) = Sg; inst(k) = trace(L*Sg);
end
g = Y(:,end);
runavg = g ./ max(t - t(1), eps);  runavg(1) = inst(1);
% Blank the burn-in: dividing by a near-zero elapsed time makes the first
% samples spike; NaNs are skipped by plot so the curve starts cleanly.
runavg(t - t(1) < 0.02*(t(end) - t(1))) = NaN;
res.t = t; res.Sigma = Sigma; res.inst = inst; res.runavg = runavg; res.cost = g;
end

function dy = deriv(y, Acl, Lc, W, L, n, Nc)
Sg = reshape(y(1:n*n), n, n); Sg = (Sg+Sg.')/2;
dS = Acl*Sg + Sg*Acl.' + W;
for l = 1:Nc; dS = dS + Lc(:,:,l)*Sg*Lc(:,:,l).'; end
dy = [dS(:); trace(L*Sg)];
end
