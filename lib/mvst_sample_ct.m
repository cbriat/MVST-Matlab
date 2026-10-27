function out = mvst_sample_ct(md, K, x0, T, dt, nreal, seed)
%MVST_SAMPLE_CT  Euler--Maruyama sample paths of a controlled diffusion.
%   OUT = MVST_SAMPLE_CT(MD,K,X0,T,DT,NREAL) simulates NREAL realizations of
%       dx = (A+BK) x dt + sum_l (Ac_l+Bc_l K) x dW_l + Gw dW0,   u = K x,
%   on [0,T] with step DT, where Gw Gw' = MD.W is the additive process noise.
%   OUT has fields .t (1 x nt), .x (n x nt x NREAL), .u (m x nt x NREAL),
%   .xmean, .umean.  Pass NREAL=1 (and MD.W=0, no noise channels) for the
%   deterministic case.  Optional SEED makes the run reproducible.
%
%   See also MVST_SAMPLE_MJLS, MVST_SIMULATE_CT.
if nargin < 7 || isempty(seed); seed = 1; end
if nargin < 6 || isempty(nreal); nreal = 8; end
rng(seed);
n = md.n; m = md.m; Nc = md.Nc;
Acl = md.A + md.B*K;
Lc = zeros(n,n,Nc);
for l = 1:Nc; Lc(:,:,l) = md.Ac(:,:,l) + md.Bc(:,:,l)*K; end
Gw = chol_psd(md.W);                          % n x n, Gw*Gw' = W

nt = round(T/dt) + 1; t = (0:nt-1)*dt;
x = zeros(n,nt,nreal); u = zeros(m,nt,nreal);
for r = 1:nreal
    xr = x0(:);
    for k = 1:nt
        x(:,k,r) = xr; u(:,k,r) = K*xr;
        if k < nt
            inc = Acl*xr*dt + Gw*(sqrt(dt)*randn(n,1));
            if Nc > 0
                dW = sqrt(dt)*randn(Nc,1);
                for l = 1:Nc; inc = inc + (Lc(:,:,l)*xr)*dW(l); end
            end
            xr = xr + inc;
        end
    end
end
out.t = t; out.x = x; out.u = u; out.xmean = mean(x,3); out.umean = mean(u,3);
end

function G = chol_psd(W)
W = (W+W.')/2;
[V,Dd] = eig(W);
G = V*diag(sqrt(max(real(diag(Dd)),0)));
end
