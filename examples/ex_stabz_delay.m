function ex_stabz_delay()
%EX_STABZ_DELAY  Delayed diffusion: mean-square state-feedback stabilization.
%  Synthesis mirror of the analysis example ex:stab:delay.
%     dx = (A0 x + B0 u) dt + A1 x(t-h) dW1,   u = K x.
%  Stabilization is DELAY-INDEPENDENT: a gain that stabilizes the zero-delay
%  synthesis condition (solved as a matrix sum-of-squares program with
%  arbitrary dwell-time through mvst_stabz_sos, flow channel Ac = A1)
%  stabilizes the closed loop for every delay h >= 0.  Sample paths are
%  Euler--Maruyama with a state-history buffer; the covariance trace is
%  propagated from the delayed second-moment equation
%     dX/dt = He((A0+B0K) X) + A1 X(t-h) A1'
%  with a matrix-history buffer.  Open/closed measure: lifted zero-delay flow
%  abscissa.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

n = 2; m = 1;
A0 = [0 1; 2 0.1]; B0 = [0; 1]; A1 = 0.3*eye(n); h = 0.4;

Mopen   = mvst_liftsym(A0, {A1});          % zero-delay lifted operator
absOpen = max(real(eig(Mopen)));
sys = struct('A0', A0, 'B0', B0, 'Ac', A1);
res = mvst_stabz_sos(sys, struct('type','arbitrary'));
K   = res.Kfun(0);
Mcl       = mvst_liftsym(A0 + B0*K, {A1});
absClosed = max(real(eig(Mcl)));

% ------------------- sample paths (delayed Euler + history buffer) ----------
T = 4; dt = 2e-3; nreal = 16; x0 = [1; -1];
[t, Xo, ~ ] = sim_delay(A0, B0, A1, zeros(m,n), h, x0, T, dt, nreal, 7);
[~, Xc, Uc] = sim_delay(A0, B0, A1, K,          h, x0, T, dt, nreal, 7);

% ------------------- covariance trace (delayed 2nd-moment ODE) --------------
tc = linspace(0,T,300);
trO = covtr_delay(A0,          A1, h, tc, eye(n));
trC = covtr_delay(A0 + B0*K,   A1, h, tc, eye(n));

mvst_stabfig(figdir, 'delay', t, Xo, Xc, Uc, tc, trO, trC, false, gc);
mvst_stabtable(figdir, 'delay');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_delay.tex'));
w.num('StabzDelayOpen', absOpen); w.num('StabzDelayClosed', absClosed);
w.mat('StabzDelayK', K); w.close();
fprintf('Delay stabilization: open abscissa %.4f, closed %.4f\n', absOpen, absClosed);
end

% ===========================================================================
function [t, X, U] = sim_delay(A0, B0, A1, K, h, x0, T, dt, nreal, seed)
rng(seed); n = numel(x0); m = size(K,1);
nt = round(T/dt)+1; t = (0:nt-1)*dt; hd = round(h/dt); Acl = A0 + B0*K;
X = zeros(n,nt,nreal); U = zeros(m,nt,nreal);
for r = 1:nreal
    x = x0; hist = repmat(x0,1,hd+1);       % history for x(t-h)
    for k = 1:nt
        X(:,k,r) = x; U(:,k,r) = K*x;
        if k < nt
            xd = hist(:,1);                 % delayed state
            x = x + Acl*x*dt + (A1*xd)*(sqrt(dt)*randn);
            hist = [hist(:,2:end), x];       %#ok<AGROW>
        end
    end
end
end

function tr = covtr_delay(Acl, A1, h, tc, X0)
% Euler on dX/dt = He(Acl X) + A1 X(t-h) A1', X(s)=X0 for s<=0.
n = size(X0,1); dt = 5e-4; T = tc(end); nt = round(T/dt)+1; hd = round(h/dt);
Xh = repmat(X0(:),1,hd+1); X = X0; tt = (0:nt-1)*dt; trr = zeros(1,nt);
for k = 1:nt
    trr(k) = trace(X);
    if k < nt
        Xd = reshape(Xh(:,1), n, n);
        dX = Acl*X + X*Acl.' + A1*Xd*A1.';
        X = X + dX*dt; X = (X+X.')/2;
        Xh = [Xh(:,2:end), X(:)];            %#ok<AGROW>
    end
end
tr = interp1(tt, trr, tc, 'linear', 'extrap');
end
