function ex_stabz_poissample()
%EX_STABZ_POISSAMPLE  Sampled-data system with Poissonian sampling: mean-square
%  stabilization by design and certify.  Synthesis mirror of ex:stab:poissample.
%  The plant with multiplicative noise
%     dx = (A x + B u) dt + Ac x dW1,   u = K x(t_k) held on (t_k,t_{k+1}],
%  is sampled at exponentially distributed instants of rate lam (Poissonian
%  sampling).  The state is extended to z = col(x,u): the flow holds the input
%  (Atil = [A B; 0 0], Actil = [Ac 0; 0 0]) and the sample refreshes it through
%  the non-invertible map Jtil = [I 0; K 0].  Averaging over the Poisson
%  sampling sends the refresh into the flow, giving the rate-weighted lifted
%  operator F(Z) = He(Atil Z) + Actil Z Actil' + lam (Jtil Z Jtil' - Z).
%  The gain sits inside Jtil, so the synthesis is NOT convex: the gain is
%  DESIGNED from the continuous flow synthesis condition on the plant (solved as
%  a matrix sum-of-squares program with arbitrary dwell-time through
%  mvst_stabz_sos) and the closed loop is then CERTIFIED by the rate-weighted
%  analysis condition (augmented rate-weighted lifted abscissa < 0).
%  Open/closed measure: that augmented rate-weighted spectral abscissa.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

n = 2; m = 1;
A = [0 1; 0 0]; B = [0; 1]; Ac = 0.2*eye(n); lam = 8.0;   % master double integrator (eq poissample data)

rc = mvst_stabz_sos(struct('A0',A,'B0',B,'Ac',Ac), struct('type','arbitrary'));
Kc = rc.Kfun(0);                               % continuous design (flow SOS)
[absOpen,   Mo] = poisop(A, B, Ac, zeros(m,n), lam, n, m);   % open loop: Jtil=[I 0;0 0]
[absClosed, Mc] = poisop(A, B, Ac, Kc,          lam, n, m);

% ------------------- sample paths (random Poisson sampling) -----------------
T = 4; dt = 2e-3; nreal = 16; x0 = [1; -1];
[t, Xo, ~ ] = sim_poisson(A, B, Ac, zeros(m,n), lam, x0, T, dt, nreal, 31);
[~, Xc, Uc] = sim_poisson(A, B, Ac, Kc,          lam, x0, T, dt, nreal, 31);

% ------------------- covariance trace (exact rate-weighted operator) --------
n2 = n + m;
[tc, trO] = poistr(Mo, T, n2, n);
[~,  trC] = poistr(Mc, T, n2, n);

mvst_stabfig(figdir, 'poissample', t, Xo, Xc, Uc, tc, trO, trC, false, gc);
mvst_stabtable(figdir, 'poissample');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_poissample.tex'));
w.num('StabzPoisOpen', absOpen); w.num('StabzPoisClosed', absClosed);
w.mat('StabzPoisK', Kc); w.num('StabzPoisLam', lam); w.close();
fprintf('Poisson-sampling stabilization: open abscissa %.4f, closed %.4f\n', absOpen, absClosed);
end

% ===========================================================================
function [ab, M] = poisop(A, B, Ac, K, lam, n, m)
% Augmented rate-weighted lifted flow operator on z = col(x,u):
%   F(Z) = He(Atil Z) + Actil Z Actil' + lam (Jtil Z Jtil' - Z).
Atil  = [A B; zeros(m,n) zeros(m)];
Actil = [Ac zeros(n,m); zeros(m,n) zeros(m)];
Jtil  = [eye(n) zeros(n,m); K zeros(m)];
M  = mvst_liftsym(Atil, {Actil, sqrt(lam)*Jtil}) - lam*eye(size(mvst_svecbasis(n+m),2));
ab = max(real(eig(M)));
end

function [t, tr] = poistr(M, T, n2, n)
% Trace of the x-block of the augmented second moment propagated exactly by
% the linear rate-weighted operator M: Z(t) = expm(t M) Z0, Z0 = diag(I_n,0).
Ud = mvst_svecbasis(n2); Z0 = zeros(n2); Z0(1:n,1:n) = eye(n);
s0 = Ud'*reshape(Z0,[],1); t = linspace(0,T,300); tr = zeros(size(t));
for k = 1:numel(t)
    Z = reshape(Ud*(expm(t(k)*M)*s0), n2, n2);
    tr(k) = trace(Z(1:n,1:n));
end
end

function [t, X, U] = sim_poisson(A, B, Ac, K, lam, x0, T, dt, nreal, seed)
% Sample paths of the actual sampled-data jump diffusion with Poisson
% sampling: Euler flow holding u; at each Poisson instant u <- K x.
rng(seed); n = numel(x0); m = size(K,1); nt = round(T/dt)+1; t = (0:nt-1)*dt;
X = zeros(n,nt,nreal); U = zeros(m,nt,nreal);
for r = 1:nreal
    x = x0; u = zeros(m,1);
    for k = 1:nt
        X(:,k,r) = x; U(:,k,r) = u;
        if k < nt
            x = x + (A*x + B*u)*dt + (Ac*x)*(sqrt(dt)*randn);
            if rand < lam*dt; u = K*x; end     % Poisson sample: refresh held input
        end
    end
end
end
