function ex_stabz_mjls()
%EX_STABZ_MJLS  Markov jump linear system: mean-square stabilization by
%  mode-dependent state feedback.  Synthesis mirror of ex:stab:mjls.
%     dx = A_{s(t)} x dt + Bc_{s(t)} x dW1 + B_{s(t)} u dt,   u = K_{s(t)} x,
%  with s(t) a continuous-time Markov chain of generator Pi.  Each mode gain
%  K_i is obtained by solving the per-mode continuous-time synthesis condition
%  (flow channel Ac = Bc_i) as a matrix sum-of-squares program (arbitrary
%  dwell-time) through mvst_stabz_sos; the closed MJLS is verified mean-square stable through
%  the coupled lifted generator.  Sample paths simulate the chain and the
%  mode-dependent diffusion; the covariance trace is the Monte-Carlo mean of
%  ||x||^2 over many realizations (open-loop unstable, closed-loop stable).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

% ------------------- 2-mode MJLS --------------------------------------------
n = 2; m = 1; M = 2;
A = cat(3, [0 1; 2 0.1], [0.5 1; -3 0.2]);        % both unstable
B = cat(3, [0;1], [0;1]);
Bc = cat(3, 0.2*eye(n), 0.2*eye(n));              % state-multiplicative noise
Pi = [-1 1; 1 -1];                                % generator (rate 1)

radOpen = mjls_rad(A, zeros(n,m,M), Bc, Pi, n, M);   % lifted coupled abscissa
K = zeros(m,n,M);
for i = 1:M
    sys = struct('A0', A(:,:,i), 'B0', B(:,:,i), 'Ac', Bc(:,:,i));
    res = mvst_stabz_sos(sys, struct('type','arbitrary'));
    K(:,:,i) = res.Kfun(0);
end
Acl = A; for i = 1:M; Acl(:,:,i) = A(:,:,i) + B(:,:,i)*K(:,:,i); end
radClosed = mjls_rad(Acl, zeros(n,m,M), Bc, Pi, n, M);

% ------------------- sample paths + Monte-Carlo covariance ------------------
T = 5; dt = 2e-3; nshow = 16; nmc = 400; x0 = [1; -1];
[t, Xo, ~,  trO] = sim_mjls(A,  B, Bc, zeros(m,n,M), Pi, x0, T, dt, nshow, nmc, 11);
[~, Xc, Uc, trC] = sim_mjls(A,  B, Bc, K,            Pi, x0, T, dt, nshow, nmc, 11);

mvst_stabfig(figdir, 'mjls', t, Xo, Xc, Uc, t, trO, trC, false, gc);
mvst_stabtable(figdir, 'mjls');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_mjls.tex'));
w.num('StabzMjlsOpen', radOpen); w.num('StabzMjlsClosed', radClosed);
w.mat('StabzMjlsKone', K(:,:,1)); w.mat('StabzMjlsKtwo', K(:,:,2)); w.close();
fprintf('MJLS stabilization: open abscissa %.4f, closed %.4f\n', radOpen, radClosed);
end

% ===========================================================================
function [t, X, U, tr] = sim_mjls(A, B, Bc, K, Pi, x0, T, dt, nshow, nmc, seed)
rng(seed); n = numel(x0); m = size(K,1); M = size(A,3);
nt = round(T/dt)+1; t = (0:nt-1)*dt;
X = zeros(n,nt,nshow); U = zeros(m,nt,nshow); sq = zeros(1,nt);
rates = -diag(Pi);
for r = 1:nmc
    x = x0; s = 1 + (rand > 0.5);                 % random initial mode
    for k = 1:nt
        if r <= nshow; X(:,k,r) = x; U(:,k,r) = K(:,:,s)*x; end
        sq(k) = sq(k) + (x.'*x);
        if k < nt
            Acl = A(:,:,s) + B(:,:,s)*K(:,:,s);
            x = x + Acl*x*dt + (Bc(:,:,s)*x)*(sqrt(dt)*randn);
            if rand < rates(s)*dt                  % Markov transition
                pj = Pi(s,:); pj(s) = 0; pj = pj/sum(pj);
                s = find(rand <= cumsum(pj), 1);
            end
        end
    end
end
tr = sq / nmc;
end

function rad = mjls_rad(A, ~, Bc, Pi, n, M)
% Spectral abscissa of the coupled lifted second-moment generator on (S^n)^M.
Us = mvst_svecbasis(n); d = size(Us,2); G = zeros(M*d);
for i = 1:M
    Mi = Us' * (kron(eye(n),A(:,:,i)) + kron(A(:,:,i),eye(n)) ...
              + kron(Bc(:,:,i),Bc(:,:,i))) * Us;
    G((i-1)*d+(1:d), (i-1)*d+(1:d)) = Mi + Pi(i,i)*eye(d);
    for j = 1:M
        if j ~= i
            G((i-1)*d+(1:d), (j-1)*d+(1:d)) = Pi(j,i)*eye(d);
        end
    end
end
rad = max(real(eig(G)));
end
