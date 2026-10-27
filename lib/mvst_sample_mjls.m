function out = mvst_sample_mjls(sys, K, x0, T, dt, nreal, seed)
%MVST_SAMPLE_MJLS  Sample paths of a Markov-jump controlled diffusion.
%   OUT = MVST_SAMPLE_MJLS(SYS,K,X0,T,DT,NREAL) simulates NREAL realizations
%   of the closed loop dx = (A_s+B_s K_s) x dt + sum_l (Ac_{s,l}+Bc_{s,l}K_s)
%   x dW_l, with the mode s(t) a Markov chain with generator SYS.Pi and the
%   feedback u = K{s} x.  OUT has fields .t, .x (n x nt x NREAL),
%   .u (m x nt x NREAL), .mode (nt x NREAL).
%
%   See also MVST_SAMPLE_CT, MVST_SIMULATE_MJLS.
if nargin < 7 || isempty(seed); seed = 1; end
if nargin < 6 || isempty(nreal); nreal = 6; end
rng(seed);
n = sys.n; m = sys.m; Pi = sys.Pi; pi0 = sys.pi0;
nt = round(T/dt)+1; t = (0:nt-1)*dt;
x = zeros(n,nt,nreal); u = zeros(m,nt,nreal); mode = zeros(nt,nreal);

for r = 1:nreal
    xr = x0(:); s = pick(pi0);
    for k = 1:nt
        x(:,k,r) = xr; u(:,k,r) = K{s}*xr; mode(k,r) = s;
        if k < nt
            md = sys.modes{s};
            Acl = md.A + md.B*K{s};
            inc = Acl*xr*dt;
            for l = 1:md.Nc
                Lcl = md.Ac(:,:,l) + md.Bc(:,:,l)*K{s};
                inc = inc + (Lcl*xr)*(sqrt(dt)*randn);
            end
            xr = xr + inc;
            % Markov mode update (Euler on the generator).
            pr = Pi(s,:)*dt; pr(s) = 0; pr = max(pr,0);
            if rand < sum(pr); s = pick(pr/sum(pr)); end
        end
    end
end
out.t = t; out.x = x; out.u = u; out.mode = mode;
end

function i = pick(p)
p = p(:).'/sum(p); i = find(rand <= cumsum(p), 1, 'first');
if isempty(i); i = numel(p); end
end
