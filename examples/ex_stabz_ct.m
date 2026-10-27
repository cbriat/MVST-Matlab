function ex_stabz_ct()
%EX_STABZ_CT  Continuous-time (Brownian) diffusion: mean-square state-feedback
%  stabilization.  Synthesis mirror of the analysis example ex:stab:poisson
%  (Stochastic systems with Brownian motions and Poissonian jumps).  The
%  controlled diffusion
%     dx = (A0 x + B0 u) dt + A1 x dW1,   u = K x,
%  is open-loop mean-square unstable; the gain K is obtained by solving the
%  chapter synthesis condition (constant-data Brownian slice of
%  Theorem th:stab:stabzpoisson / eq:stab:stabzpoisson) as a matrix
%  sum-of-squares program (arbitrary dwell-time) through mvst_stabz_sos.
%
%  Produces two figures and one data file:
%    figures/stabz_ct_paths.pdf : open-loop state sample paths (top),
%                                 closed-loop paths (middle), control u (bottom);
%    figures/stabz_ct_cov.pdf   : tr X(t), open loop RED vs closed loop BLUE;
%    figures/data/ex_stabz_ct.tex : \StabzCtOpen, \StabzCtClosed (lifted-flow
%                                   spectral abscissa, open / closed).
%  Sample paths are Euler--Maruyama (mvst_sample_ct); the covariance trace is
%  propagated EXACTLY through the reduced flow operator (mvst_liftsym).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

% ------------------- system: unstable controlled diffusion -------------------
n = 2; m = 1;
A0 = [0 1; 2 0.1];        % unstable drift
B0 = [0; 1];              % single input
A1 = 0.3*eye(n);          % state-multiplicative Brownian channel
B1 = 0.1*[0; 1];          % input selector actuating the Brownian channel (eq:stabz:ct:data)

Mopen   = mvst_liftsym(A0, {A1});          % lifted open-loop flow operator
absOpen = max(real(eig(Mopen)));           % > 0 : mean-square unstable

% ------------------- chapter synthesis (matrix SOS, SOSTOOLS + SeDuMi) -------
sys = struct('A0', A0, 'B0', B0, 'Ac', A1, 'Bc', B1);
res = mvst_stabz_sos(sys, struct('type','arbitrary'));
K   = res.Kfun(0);
Mcl       = mvst_liftsym(A0 + B0*K, {A1 + B1*K});   % control routed into the channel
absClosed = max(real(eig(Mcl)));           % < 0 : mean-square stable

% ------------------- sample paths: open (K=0) and closed --------------------
md = mvst_mode(A0, B0, eye(n), eye(m), 'Ac', A1, 'Bc', B1, 'X0', eye(n));
T = 4; dt = 2e-3; nreal = 16; x0 = [1; -1];
so = mvst_sample_ct(md, zeros(m,n), x0, T, dt, nreal, 7);   % open loop
sc = mvst_sample_ct(md, K,          x0, T, dt, nreal, 7);   % closed loop

% ------------------- covariance trace (exact) -------------------------------
Us = mvst_svecbasis(n); s0 = Us' * reshape(eye(n), [], 1);
tc = linspace(0, T, 300); trO = zeros(size(tc)); trC = zeros(size(tc));
for k = 1:numel(tc)
    trO(k) = trace(reshape(Us*(expm(tc(k)*Mopen)*s0), n, n));
    trC(k) = trace(reshape(Us*(expm(tc(k)*Mcl )*s0), n, n));
end

% ------------------- the two standard figures -------------------------------
mvst_stabfig(figdir, 'ct', so.t, so.x, sc.x, sc.u, tc, trO, trC, false, gc);
mvst_stabtable(figdir, 'ct');   %% generated summary-table row

% ------------------- data macros --------------------------------------------
w = mvst_texdata(fullfile(figdir,'data','ex_stabz_ct.tex'));
w.num('StabzCtOpen',   absOpen);
w.num('StabzCtClosed', absClosed);
w.mat('StabzCtK', K);
w.close();
fprintf('CT stabilization: open abscissa %.4f, closed %.4f\n', absOpen, absClosed);
end
