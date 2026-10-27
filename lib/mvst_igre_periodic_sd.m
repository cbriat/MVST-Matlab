function out = mvst_igre_periodic_sd(A, B, Ac, Bc, Q, R, Tbar, opts)
%MVST_IGRE_PERIODIC_SD  Periodic-sampling sampled-data optimal control.
%   OUT = MVST_IGRE_PERIODIC_SD(A,B,Ac,Bc,Q,R,Tbar) solves the impulsive
%   I-GRE for the periodic sampled-data control of the diffusion
%       dx = (A x + B u)dt + sum_l (Ac(:,:,l) x + Bc(:,:,l) u) dW_l,
%   with zero-order-hold input updated every Tbar by u(t_k^+) = K x(t_k),
%   and running cost x'Qx + u'Ru.  The augmented state is z = [x; u]; the
%   input is held during the flow (no flow control) and the gain K is the
%   jump decision.  The fixed point of
%       P0 = FlowBack( JumpRiccati(P0), Tbar )
%   is computed, where the flow step integrates the control-free augmented
%   Lyapunov--Riccati ODE backward over one period, and the jump step
%   eliminates the input block by its Schur complement (the discrete-time
%   jump Riccati), giving K = -P0(u,u)^{-1} P0(u,x).
%
%   OUT has fields .P0 (post-sample (n+m)x(n+m) value), .K (m x n gain),
%   .iters, .residual, and .Pfun(tau) (within-period P(tau), tau in [0,Tbar]).
%
%   See also MVST_SOLVE_FINITE, MVST_L2GAIN_SD_HYBRID, MVST_STEER.
if nargin < 8; opts = struct(); end
if size(A,1) ~= size(A,2)
    error('mvst:igre_periodic_sd:notSquare', ...
        'A must be square; got %d x %d.', size(A,1), size(A,2));
end
if size(B,1) ~= size(A,1)
    error('mvst:igre_periodic_sd:dim', 'B must have %d rows to match A.', size(A,1));
end
if ~(isscalar(Tbar) && isreal(Tbar) && Tbar > 0)
    error('mvst:igre_periodic_sd:period', 'TBAR must be a positive scalar.');
end
n = size(A,1); m = size(B,2); Nc = size(Ac,3);

Ab0 = [A B; zeros(m,n+m)];
Abl = zeros(n+m,n+m,Nc);
for l = 1:Nc; Abl(:,:,l) = [Ac(:,:,l) Bc(:,:,l); zeros(m,n+m)]; end
Zf  = blkdiag((Q+Q.')/2, (R+R.')/2);

tol   = mvst_getf(opts,'tol',1e-9);
maxit = mvst_getf(opts,'maxit',200);
rt    = mvst_getf(opts,'reltol',1e-9);
oo    = odeset('RelTol',rt,'AbsTol',rt*1e-2);
flowfun = @(t,pv) reshape(-flowrhs(reshape(pv,n+m,n+m), Ab0, Abl, Zf, Nc), [], 1);

P0 = blkdiag(eye(n), eye(m)); delta = inf; it = 0;  % PD start (invertible u-block)
while it < maxit && delta > tol
    it = it + 1;
    Ppre = jumpriccati(P0, n, m);
    [~, Pv] = ode45(flowfun, [Tbar 0], Ppre(:), oo);   % integrate Tbar -> 0
    P0new = reshape(Pv(end,:), n+m, n+m); P0new = (P0new+P0new.')/2;
    delta = norm(P0new - P0, 'fro');
    P0 = P0new;
end
[K, ~] = jumpriccati_gain(P0, n, m);

out.P0 = P0; out.K = K; out.iters = it; out.residual = delta;
out.Pfun = @(tau) flowback(jumpriccati(P0,n,m), Tbar, tau, flowfun, n, m);
if delta > tol
    warning('mvst:igre_periodic_sd:converge', ...
        'Periodic I-GRE fixed-point residual %.2e exceeds tol %.1e.', delta, tol);
end
end

function dP = flowrhs(P, Ab0, Abl, Zf, Nc)
dP = Ab0.'*P + P*Ab0 + Zf;
for l = 1:Nc; dP = dP + Abl(:,:,l).'*P*Abl(:,:,l); end
dP = (dP+dP.')/2;
end

function [K, sc] = jumpriccati_gain(P0, n, m)
ix = 1:n; iu = n+1:n+m;
Puu = (P0(iu,iu)+P0(iu,iu).')/2;
Pux = P0(iu,ix);
Pxx = (P0(ix,ix)+P0(ix,ix).')/2;
Pxu = P0(ix,iu);
K  = -(Puu \ Pux);
sc = Pxx - Pxu*(Puu \ Pux);  sc = (sc+sc.')/2;     % Schur complement
end

function Ppre = jumpriccati(P0, n, m)
[~, sc] = jumpriccati_gain(P0, n, m);
Ppre = blkdiag(sc, zeros(m));
end

function P = flowback(Pterm, Tbar, tau, flowfun, n, m)
% P(tau) within a period: integrate the flow ODE backward from Tbar to tau.
if abs(tau - Tbar) < eps; P = (Pterm+Pterm.')/2; return; end
[~, Pv] = ode45(flowfun, [Tbar tau], Pterm(:));
P = reshape(Pv(end,:), n+m, n+m); P = (P+P.')/2;
end
