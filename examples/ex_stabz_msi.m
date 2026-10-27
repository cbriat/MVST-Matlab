function ex_stabz_msi()
%EX_STABZ_MSI  Stochastic impulses and switching: mean-square stabilization.
%  Synthesis mirror of ex:stab:msi.  Markovian switching between two unstable
%  continuous modes, each carrying a mode gain, with an impulse reset J applied
%  at every switch:
%     dx = (A_s x + B u) dt + Ac x dW1,  u = K_s x;   x^+ = J x at switches.
%  Each K_i is obtained by solving the continuous-time flow synthesis condition
%  as a matrix sum-of-squares program (arbitrary dwell-time) through
%  mvst_stabz_sos (as in the switched
%  example); the closed loop is verified mean-square stable through the coupled
%  lifted generator with the reset.  Sample paths simulate the chain and the
%  resets; the covariance trace is a Monte-Carlo mean.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

n = 2; m = 1; M = 2;
A  = cat(3, [0 1; 2 0.1], [0.3 1; -3 0.2]);
B  = [0;1]; Ac = 0.2*eye(n); J = 0.85*eye(n); Pi = [-1.5 1.5; 1.5 -1.5];
B1  = 0.1*[0; 1];         % input selector actuating the diffusion channel B_{l,i}
BJ0 = [0; 1];             % jump input selector B_{J,r,ell}

% ---- reset feedback u_J = K_J x: DT synthesis on the (contractive) reset ----
resJ = mvst_stabz_sos(struct('A0', J, 'B0', BJ0), struct('type','discrete'));
KJ   = resJ.Kfun(0);
Jcl  = J + BJ0*KJ;                                  % closed-loop reset map

% ---- open-loop coupled abscissa (K=0, uncontrolled reset J) -----------------
absOpen = msi_absc(A, zeros(m,n,M), B, Ac, B1, J, Pi, n, M);

% ---- per-mode flow feedback with the actuated diffusion channel -------------
K = zeros(m,n,M);
for i = 1:M
    res = mvst_stabz_sos(struct('A0',A(:,:,i),'B0',B,'Ac',Ac,'Bc',B1), struct('type','arbitrary'));
    K(:,:,i) = res.Kfun(0);
end
absClosed = msi_absc(A, K, B, Ac, B1, Jcl, Pi, n, M);

T = 6; dt = 2e-3; nshow = 16; nmc = 400; x0 = [1;-1];
[t, Xo, ~,  trO] = sim_msi(A, B, Ac, B1, J,   zeros(m,n,M), Pi, x0, T, dt, nshow, nmc, 81);
[~, Xc, Uc, trC] = sim_msi(A, B, Ac, B1, Jcl, K,            Pi, x0, T, dt, nshow, nmc, 81);

mvst_stabfig(figdir, 'msi', t, Xo, Xc, Uc, t, trO, trC, false, gc);
mvst_stabtable(figdir, 'msi');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_msi.tex'));
w.num('StabzMsiOpen', absOpen); w.num('StabzMsiClosed', absClosed);
w.mat('StabzMsiKone', K(:,:,1)); w.mat('StabzMsiKtwo', K(:,:,2));
w.mat('StabzMsiKJ', KJ); w.close();
fprintf('MSI stabilization: open abscissa %.4f, closed %.4f, KJ=[%.3f %.3f]\n', ...
        absOpen, absClosed, KJ(1), KJ(2));
end

% ===========================================================================
function a = msi_absc(A, K, B, Ac, B1, Jre, Pi, n, M)
% Coupled lifted generator: per-mode flow (channel actuated by B1), mode
% coupling with the reset Jre weighting the incoming rates (the switch carries
% an impulse reset; Jre is the closed reset J + B_{J,0} K_J in closed loop).
Us = mvst_svecbasis(n); d = size(Us,2); MJ = Us'*kron(Jre,Jre)*Us; G = zeros(M*d);
for i = 1:M
    Mi = mvst_liftsym(A(:,:,i)+B*K(:,:,i), {Ac + B1*K(:,:,i)});
    G((i-1)*d+(1:d),(i-1)*d+(1:d)) = Mi + Pi(i,i)*eye(d);
    for j = 1:M
        if j ~= i; G((i-1)*d+(1:d),(j-1)*d+(1:d)) = Pi(j,i)*MJ; end
    end
end
a = max(real(eig(G)));
end

function [t, X, U, tr] = sim_msi(A, B, Ac, B1, Jre, K, Pi, x0, T, dt, nshow, nmc, seed)
rng(seed); n = numel(x0); m = size(K,1); M = size(A,3);
nt = round(T/dt)+1; t = (0:nt-1)*dt; rates = -diag(Pi);
X = zeros(n,nt,nshow); U = zeros(m,nt,nshow); sq = zeros(1,nt);
for r = 1:nmc
    x = x0; s = 1 + (rand > 0.5);
    for k = 1:nt
        if r <= nshow; X(:,k,r) = x; U(:,k,r) = K(:,:,s)*x; end
        sq(k) = sq(k) + x.'*x;
        if k < nt
            Acl = A(:,:,s) + B*K(:,:,s);
            x = x + Acl*x*dt + ((Ac + B1*K(:,:,s))*x)*(sqrt(dt)*randn);
            if rand < rates(s)*dt
                x = Jre*x;                        % impulse reset (closed) at the switch
                pj = Pi(s,:); pj(s) = 0; pj = pj/sum(pj);
                s = find(rand <= cumsum(pj), 1);
            end
        end
    end
end
tr = sq / nmc;
end
