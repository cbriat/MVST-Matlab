function ex_stabz_dtdelay()
%EX_STABZ_DTDELAY  Discrete-time stochastic system with delay: mean-square
%  state-feedback stabilization.  Synthesis mirror of ex:stab:dtdelay.
%     x(k+1) = A0 x(k) + B0 u(k) + nu(k) A1 x(k-tau),   u = K x.
%  Stabilization is delay-independent: a gain feasible for the DT synthesis
%  condition (solved as a matrix sum-of-squares program with discrete
%  dwell-time through mvst_stabz_sos, flow channel Ac = A1) stabilizes for
%  every integer delay tau.  Sample paths use a state-history buffer; the
%  covariance trace is the exact delayed second-moment recursion.  Open/closed
%  measure: lifted zero-delay transition spectral radius.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

n = 2; m = 1;
A0 = [1.1 0.3; 0 1.05]; B0 = [0; 1]; A1 = 0.25*eye(n); tau = 3;

radOpen  = max(abs(eig(mvst_liftsym(zeros(n), {A0, A1}))));
sys = struct('A0', A0, 'B0', B0, 'Ac', A1);
res = mvst_stabz_sos(sys, struct('type','discrete'));
K   = res.Kfun(0);
radClosed = max(abs(eig(mvst_liftsym(zeros(n), {A0+B0*K, A1}))));

Nstep = 18; nreal = 16; x0 = [1; -1]; kk = 0:Nstep;
[Xo, ~ ] = sim_dtdelay(A0, B0, A1, zeros(m,n), tau, x0, Nstep, nreal, 41);
[Xc, Uc] = sim_dtdelay(A0, B0, A1, K,          tau, x0, Nstep, nreal, 41);

trO = covtr_dtdelay(A0,        A1, tau, Nstep, eye(n));
trC = covtr_dtdelay(A0+B0*K,   A1, tau, Nstep, eye(n));

mvst_stabfig(figdir, 'dtdelay', kk, Xo, Xc, Uc, kk, trO, trC, true, gc);
mvst_stabtable(figdir, 'dtdelay');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_dtdelay.tex'));
w.num('StabzDtdOpen', radOpen); w.num('StabzDtdClosed', radClosed);
w.mat('StabzDtdK', K); w.close();
fprintf('DT-delay stabilization: open radius %.4f, closed %.4f\n', radOpen, radClosed);
end

% ===========================================================================
function [X, U] = sim_dtdelay(A0, B0, A1, K, tau, x0, Nstep, nreal, seed)
rng(seed); n = numel(x0); m = size(K,1); Acl = A0 + B0*K;
X = zeros(n,Nstep+1,nreal); U = zeros(m,Nstep+1,nreal);
for r = 1:nreal
    buf = repmat(x0,1,tau+1); x = x0;
    for j = 1:Nstep+1
        X(:,j,r) = x; U(:,j,r) = K*x;
        if j <= Nstep
            xd = buf(:,1); nu = randn;
            x = Acl*x + nu*(A1*xd);
            buf = [buf(:,2:end), x];      %#ok<AGROW>
        end
    end
end
end

function tr = covtr_dtdelay(Acl, A1, tau, Nstep, X0)
n = size(X0,1); Xb = repmat(X0(:),1,tau+1); X = X0; tr = zeros(1,Nstep+1);
for j = 1:Nstep+1
    tr(j) = trace(X);
    if j <= Nstep
        Xd = reshape(Xb(:,1), n, n);
        X = Acl*X*Acl.' + A1*Xd*A1.'; X = (X+X.')/2;
        Xb = [Xb(:,2:end), X(:)];         %#ok<AGROW>
    end
end
end
