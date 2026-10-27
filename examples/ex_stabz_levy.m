function ex_stabz_levy()
%EX_STABZ_LEVY  Diffusion driven by a Levy noise: mean-square state-feedback
%  stabilization.  Synthesis mirror of the analysis example ex:stab:levy.
%     dx = (A0 x + B0 u) dt + A1 x dL,   u = K x,
%  with a Levy noise L of mean rate mu and variance rate c.  The second moment
%  sees the noise only through (mu,c): effective drift A0 + mu A1 and channel
%  sqrt(c) A1.  K is obtained by solving the chapter synthesis condition on
%  that effective diffusion as a matrix sum-of-squares program (arbitrary
%  dwell-time) through mvst_stabz_sos, with flow channel Ac = sqrt(c) A1.
%  The sample paths are driven by a compensated-Poisson Levy noise (mu=0,
%  c = lam*a^2), so the paths are genuine jump paths; the covariance trace is
%  exact through the lifted effective operator.  Open/closed measure: lifted
%  effective-flow spectral abscissa.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

% ------------------- system + Levy noise (mu, c) ----------------------------
n = 2; m = 1;
A0 = [0 1; 2 0.1]; B0 = [0; 1]; A1 = 0.3*eye(n);
B1 = 0.1*[0; 1];           % input selector actuating the Levy channel (eq:stabz:levy data)
mu = 0;                    % centered Levy driver
lam = 6; a = sqrt(0.09/lam);   % compensated Poisson: c = lam*a^2 = 0.09
c = lam*a^2;

Aeff = A0 + mu*A1; Ach = sqrt(c)*A1; Bch = sqrt(c)*B1;
Mopen   = mvst_liftsym(Aeff, {Ach});
absOpen = max(real(eig(Mopen)));

sys = struct('A0', Aeff, 'B0', B0, 'Ac', Ach, 'Bc', Bch);
res = mvst_stabz_sos(sys, struct('type','arbitrary'));
K   = res.Kfun(0);
Mcl       = mvst_liftsym(Aeff + B0*K, {Ach + Bch*K});   % control routed into the channel
absClosed = max(real(eig(Mcl)));

% ------------------- sample paths (compensated-Poisson Levy) ----------------
T = 4; dt = 2e-3; nreal = 16; x0 = [1; -1];
[t, Xo, ~ ] = sim_levy(A0, B0, A1, B1, zeros(m,n), mu, lam, a, x0, T, dt, nreal, 7);
[~, Xc, Uc] = sim_levy(A0, B0, A1, B1, K,          mu, lam, a, x0, T, dt, nreal, 7);

% ------------------- covariance trace (exact, effective operator) -----------
Us = mvst_svecbasis(n); s0 = Us'*reshape(eye(n),[],1);
tc = linspace(0,T,300); trO = zeros(size(tc)); trC = zeros(size(tc));
for k = 1:numel(tc)
    trO(k) = trace(reshape(Us*(expm(tc(k)*Mopen)*s0), n, n));
    trC(k) = trace(reshape(Us*(expm(tc(k)*Mcl )*s0), n, n));
end

mvst_stabfig(figdir, 'levy', t, Xo, Xc, Uc, tc, trO, trC, false, gc);
mvst_stabtable(figdir, 'levy');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_levy.tex'));
w.num('StabzLevyOpen', absOpen); w.num('StabzLevyClosed', absClosed);
w.mat('StabzLevyK', K); w.close();
fprintf('Levy stabilization: open abscissa %.4f, closed %.4f\n', absOpen, absClosed);
end

% ===========================================================================
function [t, X, U] = sim_levy(A0, B0, A1, B1, K, mu, lam, a, x0, T, dt, nreal, seed)
rng(seed); n = numel(x0); m = size(K,1);
nt = round(T/dt)+1; t = (0:nt-1)*dt;
X = zeros(n,nt,nreal); U = zeros(m,nt,nreal); Acl = A0 + B0*K; Achl = A1 + B1*K;
for r = 1:nreal
    x = x0;
    for k = 1:nt
        X(:,k,r) = x; U(:,k,r) = K*x;
        if k < nt
            dN = a*poissrnd(lam*dt);              % Poisson jump part
            dL = mu*dt + dN - lam*a*dt;           % compensated (mean mu*dt)
            x = x + Acl*x*dt + (Achl*x)*dL;       % channel actuated by B1
        end
    end
end
end
