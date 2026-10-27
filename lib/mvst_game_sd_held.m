function out = mvst_game_sd_held(A, Bu, Bw, Q, Ru, gam, Tp, Dwc, Nmax, tol)
%MVST_GAME_SD_HELD  Periodic HJI-GRE of the held-input sampled-data game.
%   OUT = MVST_GAME_SD_HELD(A,BU,BW,Q,RU,GAM,TP,DWC,NMAX,TOL) solves, by
%   backward iteration over sampling periods from a zero terminal value, the
%   periodic I-HJI-GRE (Definition def:game:GRE, Theorem th:game:periodic of
%   the Game chapter) of the zero-sum game on the held-input sampled-data
%   class (Section sec:pre:cov:sd):
%
%     dx = (A x + Bu u + Bw w)dt + Dwc w dW,  du = 0   on (t_k, t_k + TP),
%     x(t_k^+) = x(t_k),  u(t_k^+) = u_d(k)              at t_k = t_0 + k TP,
%
%   with running cost x'Qx + u'Ru u - GAM^2 |w|^2 and no cost at the samples.
%   On the flow the state is chi = (x,u); the Isaacs block is the maximizer
%   block alone (u is a state on the flow), so the flow Riccati eliminates w
%   only (backward in tau, with P(TP^-) given):
%
%     dP/dtau + Abar'P + P Abar + Qbar + P Bwbar Theta^{-1} Bwbar' P = 0,
%     Theta = GAM^2 I - Dwbar' P Dwbar,
%     Abar = [A Bu; 0 0], Bwbar = [Bw; 0], Dwbar = [Dwc; 0], Qbar = blkdiag(Q,Ru).
%
%   At each sample the Isaacs jump block is the minimizer block alone, and the
%   jump condition (Schur step over u_d) reads
%
%     P(TP^-) = blkdiag( Pxx - Pxu Puu^{-1} Pux , 0 ),   P = P(0^+).
%
%   Saddle gains: u_d = Ku x(t_k) with Ku = -Puu^{-1} Pux (at P(0^+)), and
%   w = Kw(tau) chi with Kw(tau) = Theta(tau)^{-1} Bwbar' P(tau).
%
%   DWC may be empty (no diffusion).  The iteration stops when the relative
%   change of the pre-sample value S_N (x-block value at t_0^-, horizon N TP)
%   falls below TOL (default 1e-10), or after NMAX periods (default 5000).
%   OUT fields:
%     ok      true if converged with Theta > 0 and Puu > 0 throughout
%     reason  'converged' | 'escape' | 'isaacs_flow' | 'isaacs_jump' | 'blowup' | 'maxiter'
%     N       number of periods iterated
%     S       converged pre-sample value (n x n), Shist(:,:,N) all iterates
%     Pplus   P(0^+) ((n+mu) x (n+mu)),  Pminus = P(TP^-)
%     Ku      sample gain (mu x n),  Kw0 = Kw(0^+) (mw x (n+mu))
%     thetaMin minimum of the least eigenvalue of Theta over the last period
%     jumpRes residual of the jump condition at the fixed point
%
%   See also MVST_GAME_GAMMA_SD_HELD, MVST_GAME_DRE_CT.
if nargin < 8; Dwc = []; end
if nargin < 9 || isempty(Nmax); Nmax = 5000; end
if nargin < 10 || isempty(tol); tol = 1e-10; end
n = size(A,1); mu = size(Bu,2); mw = size(Bw,2); N3 = n + mu;
if isempty(Dwc); Dwc = zeros(n,mw); end
Ab  = [A Bu; zeros(mu, N3)];
Bwb = [Bw; zeros(mu, mw)];
Dwb = [Dwc; zeros(mu, mw)];
Qb  = blkdiag(Q, Ru);
cap = 1e8;
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
Pm = zeros(N3);                       % P(TP^-) of the last period: zero terminal
Shist = zeros(n,n,0); out.ok = false; out.reason = 'maxiter';
Sprev = [];
for k = 1:Nmax
    [s, Y] = ode45(@(s,p) rhs(p, gam, Ab, Bwb, Dwb, Qb, N3, mw), [0 Tp], Pm(:), opts);
    if s(end) < Tp*(1-1e-12) || any(~isfinite(Y(end,:)))
        out.reason = 'escape'; if any(isnan(Y(end,:))); out.reason = 'isaacs_flow'; end
        break
    end
    Pp = reshape(Y(end,:), N3, N3); Pp = (Pp+Pp.')/2;          % P(0^+)
    Puu = Pp(n+1:end, n+1:end);
    [~, pc] = chol((Puu+Puu.')/2);
    if pc ~= 0; out.reason = 'isaacs_jump'; break; end
    S = Pp(1:n,1:n) - Pp(1:n,n+1:end)*(Puu\Pp(n+1:end,1:n)); S = (S+S.')/2;
    if norm(S) > cap; out.reason = 'blowup'; break; end
    Shist(:,:,k) = S; %#ok<AGROW>
    Pm = blkdiag(S, zeros(mu));
    if ~isempty(Sprev) && norm(S - Sprev) <= tol*max(1, norm(S))
        out.ok = true; out.reason = 'converged'; break
    end
    Sprev = S;
end
out.N = size(Shist,3); out.Shist = Shist;
if out.ok
    % one more period from the fixed point, for the gains and the residuals
    [s, Y] = ode45(@(s,p) rhs(p, gam, Ab, Bwb, Dwb, Qb, N3, mw), [0 Tp], Pm(:), opts);
    Pp = reshape(Y(end,:), N3, N3); Pp = (Pp+Pp.')/2;
    Puu = Pp(n+1:end, n+1:end);
    S = Pp(1:n,1:n) - Pp(1:n,n+1:end)*(Puu\Pp(n+1:end,1:n)); S = (S+S.')/2;
    out.S = S; out.Pplus = Pp; out.Pminus = Pm;
    out.jumpRes = norm(blkdiag(S, zeros(mu)) - Pm);
    out.Ku  = -(Puu\Pp(n+1:end,1:n));
    Th0 = gam^2*eye(mw) - Dwb.'*Pp*Dwb;
    out.Kw0 = Th0\(Bwb.'*Pp);
    th = zeros(numel(s),1);
    for i = 1:numel(s)
        Pi = reshape(Y(i,:), N3, N3); Pi = (Pi+Pi.')/2;
        th(i) = min(eig(gam^2*eye(mw) - Dwb.'*Pi*Dwb));
    end
    out.thetaMin = min(th);
    out.tau = Tp - s;                  % timer value tau for each ODE node
    out.Ptau = Y;                      % P(tau) rows (vectorized)
end
end

function dp = rhs(p, gam, Ab, Bwb, Dwb, Qb, N3, mw)
%RHS  dP/ds with s = TP - tau (backward timer): Abar'P + P Abar + Qbar + P Bw Theta^{-1} Bw' P.
P = reshape(p, N3, N3); P = (P+P.')/2;
Th = gam^2*eye(mw) - Dwb.'*P*Dwb; Th = (Th+Th.')/2;
[~, pc] = chol(Th);
if pc ~= 0; dp = nan(N3*N3,1); return; end         % Isaacs fails on the flow
v = P*Bwb;
G = Ab.'*P + P*Ab + Qb + v*(Th\v.'); G = (G+G.')/2;
dp = G(:);
end
