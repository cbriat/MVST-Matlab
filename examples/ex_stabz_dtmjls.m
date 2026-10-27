function ex_stabz_dtmjls()
%EX_STABZ_DTMJLS  Discrete-time Markov jump linear system: mean-square
%  stabilization by mode-dependent state feedback.  Synthesis mirror of
%  ex:stab:dtmjls.
%     x(k+1) = A_s x(k) + B_s u(k) + nu(k) Bc_s x(k),   u = K_s x,
%  s a discrete-time Markov chain of transition matrix P.  Each K_i solves the
%  DT synthesis condition (per-mode, contraction margin rho) as a matrix
%  sum-of-squares program (discrete dwell-time) through mvst_stabz_sos,
%  verified through the coupled lifted transition radius.  Sample
%  paths simulate the chain; the covariance trace is a Monte-Carlo mean.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

n = 2; m = 2; M = 2;
A  = cat(3, [1.1 0.3; 0 1.05], [1.05 -0.4; 0.2 1.1]);
B  = cat(3, eye(n), eye(n));              % full actuation (mode-dependent gains)
Bc = cat(3, 0.15*eye(n), 0.15*eye(n));
P  = [0.7 0.3; 0.4 0.6];                 % Markov transition matrix

radOpen = dtmjls_rad(A, zeros(n,m,M), Bc, P, n, M);
% Per-mode DT synthesis with a contraction margin rho, tightened until the
% COUPLED MJLS radius is certified below 1 (the design condition).
rho = 0.6; radClosed = inf;
while radClosed >= 1 && rho > 0.15
    K = zeros(m,n,M);
    for i = 1:M
        sys = struct('A0', A(:,:,i), 'B0', B(:,:,i), 'Ac', Bc(:,:,i));
        res = mvst_stabz_sos(sys, struct('type','discrete','rho',rho));
        K(:,:,i) = res.Kfun(0);
    end
    Acl = A; for i = 1:M; Acl(:,:,i) = A(:,:,i) + B(:,:,i)*K(:,:,i); end
    radClosed = dtmjls_rad(Acl, zeros(n,m,M), Bc, P, n, M);
    rho = rho - 0.1;
end

Nstep = 18; nshow = 16; nmc = 400; x0 = [1; -1]; kk = 0:Nstep;
[Xo, ~,  trO] = sim_dtmjls(A, B, Bc, zeros(m,n,M), P, x0, Nstep, nshow, nmc, 51);
[Xc, Uc, trC] = sim_dtmjls(A, B, Bc, K,            P, x0, Nstep, nshow, nmc, 51);

mvst_stabfig(figdir, 'dtmjls', kk, Xo, Xc, Uc, kk, trO, trC, true, gc);
mvst_stabtable(figdir, 'dtmjls');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_dtmjls.tex'));
w.num('StabzDtmOpen', radOpen); w.num('StabzDtmClosed', radClosed);
w.mat('StabzDtmKone', K(:,:,1)); w.mat('StabzDtmKtwo', K(:,:,2)); w.close();
fprintf('DT-MJLS stabilization: open radius %.4f, closed %.4f\n', radOpen, radClosed);
end

% ===========================================================================
function [X, U, tr] = sim_dtmjls(A, B, Bc, K, P, x0, Nstep, nshow, nmc, seed)
rng(seed); n = numel(x0); m = size(K,1); M = size(A,3);
X = zeros(n,Nstep+1,nshow); U = zeros(m,Nstep+1,nshow); sq = zeros(1,Nstep+1);
for r = 1:nmc
    x = x0; s = 1 + (rand > 0.5);
    for j = 1:Nstep+1
        if r <= nshow; X(:,j,r) = x; U(:,j,r) = K(:,:,s)*x; end
        sq(j) = sq(j) + x.'*x;
        if j <= Nstep
            Acl = A(:,:,s) + B(:,:,s)*K(:,:,s);
            x = Acl*x + randn*(Bc(:,:,s)*x);
            s = find(rand <= cumsum(P(s,:)), 1);
        end
    end
end
tr = sq / nmc;
end

function rad = dtmjls_rad(A, ~, Bc, P, n, M)
% Coupled DT MJLS second-moment map on (S^n)^M: X_i^+ = sum_j p_{ji} T_j(X_j),
% so the (i,j) block is p_{ji} * T_j (incoming from source mode j).
Us = mvst_svecbasis(n); d = size(Us,2); T = cell(1,M); G = zeros(M*d);
for j = 1:M
    T{j} = Us' * (kron(A(:,:,j),A(:,:,j)) + kron(Bc(:,:,j),Bc(:,:,j))) * Us;
end
for i = 1:M
    for j = 1:M
        G((i-1)*d+(1:d), (j-1)*d+(1:d)) = P(j,i)*T{j};
    end
end
rad = max(abs(eig(G)));
end
