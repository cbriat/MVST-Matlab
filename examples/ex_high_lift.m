%% EX_HIGH_LIFT  Lifted certificates for a matrix-valued impulsive system.
%  Chapter "Higher-order systems", numerical illustration of the
%  tensor-power lift (Lemmas lem:high:lift and lem:high:cp) and of the
%  Stein certificate of the lifted monodromy.  It does NOT compute the
%  linear certificate <P_m, X^{otimes m}> of Theorem th:high:lifted (degree
%  m in X).  The underlying system is a single-mode impulsive
%  covariance system on the Hermitian 2x2 cone,
%      flow   Fdot(X) = A0 X + X A0' + C X C'      (completely positive),
%      jump   J(X)    = J0 X J0'                     at each period Tbar,
%  and for m = 1,2,3 the reduced tensor-power lift (mvst_liftsym /
%  mvst_svecbasis, compressed to the symmetric power) builds the pair
%  (Fdot^[m], J^[m]).  The script
%    * verifies the semigroup identity e^{t Fdot^[m]} = (e^{t Fdot})^{otimes m}
%      of Lemma lem:high:lift to machine precision,
%    * verifies that the lift stays completely positive (Lemma lem:high:cp)
%      by a Choi-matrix test on the flow map and on the tensor-power jump,
%    * solves, at each order and at the periodic operating dwell-time Top,
%      the Stein equation Mm' P_m Mm - P_m = -I of the lifted monodromy
%      Mm = J^[m] e^{Top Fdot^[m]} (flow for Top, then jump), written in an
%      orthonormal basis of the reduced symmetric power, and checks P_m > 0.
%      With xi the coordinate vector of X in the reduced basis, the form
%      (xi^{otimes m})' P_m xi^{otimes m} is quadratic in the lifted
%      coordinate vector, of degree 2m in X, and decreases from one period
%      to the next: the certificate eq:high:momentsdp:cert of the covariance
%      lift of the vectorized system, with k = 0 and the coupling row only
%      (no flow row, no multiplier),
%    * reports the spectral radius rho_m of the lifted monodromy and checks
%      rho_m = rho_1^m (a spectral radius, not a certified contraction).
%  Data macros are written to figures/data/high_lift.tex.  No figure.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)

% ---------------- base impulsive instance -----------------------------
A0 = [-1 1; -1 -1];          % Hurwitz flow drift (eig -1 +- i)
C  = 0.2*eye(2);             % one completely positive diffusion channel
J0 = [1.2 0.4; 0 1.1];       % expanding, non-normal reset
n  = 2;
Lred = mvst_liftsym(A0, {C});          % 3x3 reduced flow generator
Jred = mvst_liftsym(zeros(n), {J0});   % 3x3 reduced jump
d    = size(Lred,1);                    % = 3

% exact periodic (minimum) dwell-time from the base monodromy
rhoM  = @(T) max(abs(eig(Jred*expm(Lred*T))));
Tstar = bisect(@(T) rhoM(T)-1, 1e-4, 5, 1e-12);
Top   = 0.30;                            % operating dwell-time > Tstar

% flow map completely positive?  (Choi of e^{Top Fdot} on the base)
Mfull   = kron(eye(n),A0) + kron(A0,eye(n)) + kron(C,C);
choiF   = choi_super(expm(Mfull*Top), n);
choiFmin = min(real(eig(choiF)));

dims  = zeros(1,3); rho = zeros(1,3);
semiR = zeros(1,3); choiJ = zeros(1,3); mineigP = zeros(1,3);
Tchk  = 0.37;
for m = 1:3
    W  = symb(d,m);  dm = size(W,2);  dims(m) = dm;
    Gm = W'*ksum(Lred,m)*W;                       % lifted flow generator
    Jm = W'*kpow(Jred,m)*W;                        % lifted tensor-power jump
    % Lemma lem:high:lift: e^{t Gm} = (e^{t Lred})^{otimes m} on the sym. power
    semiR(m) = norm(expm(Gm*Tchk) - W'*kpow(expm(Lred*Tchk),m)*W, 'fro');
    % Lemma lem:high:cp: tensor-power jump completely positive (Choi >= 0)
    choiJ(m) = min(real(eig(choi_conj(kpow(J0,m)))));
    % periodic monodromy and Stein certificate (Stein equation)
    Mm = Jm*expm(Gm*Top);  rho(m) = max(abs(eig(Mm)));
    P  = dlyap(Mm', eye(dm));  P = (P+P')/2;  mineigP(m) = min(eig(P));
end

% norm-power identity ||X^{otimes m}|| = ||X||^m on the cone (trace)
Xt = [2 0.5; 0.5 1];  Xt = Xt*Xt';   % a PSD test state
normR = abs(trace(kpow(Xt,3)) - trace(Xt)^3);

% ---------------- soundness verification ------------------------------
assert(all(semiR < 1e-9),  'semigroup identity failed');
assert(all(choiJ > -1e-9), 'lifted jump not completely positive');
assert(choiFmin > -1e-9,   'flow map not completely positive');
assert(all(mineigP > 0),   'Stein certificate P_m not positive definite');
assert(normR < 1e-9,       'norm-power identity failed');
assert(abs(rho(2)-rho(1)^2) < 1e-9 && abs(rho(3)-rho(1)^3) < 1e-9, ...
    'lifted monodromy is not the tensor power of the base');

fprintf('EX_HIGH_LIFT: Tstar=%.6f, dims=[%d %d %d]\n', Tstar, dims);
fprintf('  rho_m(Top=%.2f) = [%.6f %.6f %.6f]\n', Top, rho);
fprintf('  max semigroup resid=%.2e, min Choi(J^m)=%.2e, Choi(flow)=%.2e\n', ...
    max(semiR), min(choiJ), choiFmin);

% ---------------- data macros -----------------------------------------
w = mvst_texdata(fullfile(mvst_figdir(),'data','high_lift.tex'));
w.num('HighLiftN', n, 0);
w.num('HighLiftDimOne',   dims(1), 0);
w.num('HighLiftDimTwo',   dims(2), 0);
w.num('HighLiftDimThree', dims(3), 0);
w.num('HighLiftTstar', Tstar, 3);
w.num('HighLiftTop',   Top,   2);
w.num('HighLiftRhoOne',   rho(1), 3);
w.num('HighLiftRhoTwo',   rho(2), 3);
w.num('HighLiftRhoThree', rho(3), 3);
w.raw('HighLiftSemiResid', sprintf('%.0e', max(semiR)));
w.close();

% ======================= local functions ==============================
function W = symb(d,m)
%SYMB Orthonormal basis of the symmetric subspace of (R^d)^{otimes m}.
if m==1, W=eye(d); return; end
P = zeros(d^m); pm = perms(1:m);
for p = 1:size(pm,1), P = P + permop(d,m,pm(p,:)); end
P = P/size(pm,1);
[U,S] = svd(P); W = U(:, diag(S) > 0.5);
end
function Q = permop(d,m,sigma)
%PERMOP Operator on (R^d)^{otimes m} permuting the m tensor factors by SIGMA.
N=d^m; Q=zeros(N); idx=cell(1,m);
for lin=1:N
    [idx{:}] = ind2sub(repmat(d,1,m), lin);  v = cell2mat(idx);
    w = zeros(1,m); w(sigma) = v;  wc = num2cell(w);
    Q(sub2ind(repmat(d,1,m), wc{:}), lin) = 1;
end
end
function G = ksum(L,m)
%KSUM m-fold Kronecker sum L (+) ... (+) L on (R^d)^{otimes m}.
d=size(L,1); G=zeros(d^m);
for i=1:m
    t=1; for k=1:m, if k==i, t=kron(t,L); else, t=kron(t,eye(d)); end; end
    G=G+t;
end
end
function K = kpow(M,m)
%KPOW m-fold Kronecker power of M.
K=1; for k=1:m, K=kron(K,M); end
end
function C = choi_conj(N)
%CHOI_CONJ Choi matrix of the completely positive map X -> N X N'.
d=size(N,1); C=zeros(d^2);
for i=1:d, for j=1:d
    E=zeros(d); E(i,j)=1;  C = C + kron(E, N*E*N');
end, end
C=(C+C')/2;
end
function C = choi_super(Msup, d)
%CHOI_SUPER Choi matrix of a superoperator given by its vec-matrix MSUP.
C=zeros(d^2);
for i=1:d, for j=1:d
    E=zeros(d); E(i,j)=1;  PhiE = reshape(Msup*E(:), d, d);
    C = C + kron(E, PhiE);
end, end
C=(C+C')/2;
end
function x = bisect(f,a,b,tol)
%BISECT Sign-change bisection of F on [A,B].
fa=f(a); fb=f(b);
if fa*fb>0, x=NaN; return; end
for it=1:200
    x=(a+b)/2; fx=f(x);
    if abs(fx)<tol || (b-a)/2<tol, return; end
    if fa*fx<0, b=x; else, a=x; fa=fx; end
end
end
