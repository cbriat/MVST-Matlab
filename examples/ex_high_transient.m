%% EX_HIGH_TRANSIENT  Lifted Gramian sets of the transient-growth instance.
%  Chapter "Higher-order systems", Section sec:high:emergence, paragraph
%  after Corollary cor:high:tracemono (and the closing sentence of the
%  remark before Example ex:high:loewner).  CHAPTER CLAIM: for the instance
%      Fdot(X) = He(A X) = A X + X A',   J(X) = J X J',
%      A = [-0.5 12; 0 -0.7],   J = [0.8 0.2; 0 0.5],   Tbar = 0.3,
%  over the minimum dwell-time family S_min[Tbar], "no greatest element of
%  the computed lifted Gramians was found up to order six", and the
%  product-state ratios decrease across m = 1..6.  STATUS: numerical
%  evidence over a FINITE grid of admissible sequences, NOT a proof; the
%  chapter states it as such ("These observations are not proved here").
%
%  LIFTED GRAMIAN.  By Lemma lem:high:moment (eq:high:momentop),
%      Pbar^[m]_sigma(0) = int_0^inf G_s^{otimes m} ds + sum_k G_{t_k^-}^{otimes m},
%      G_s = S_sigma^*(s,0) I_n.
%  For a single-Kraus jump every transition is a congruence: on the flow
%  interval [t_j, t_{j+1}) one has S_sigma(s,0)Y = M Y M' with
%  M = e^{A(s-t_j)} N_j and N_j = J e^{A Delta_{j-1}} ... J e^{A t_1}
%  (N_0 = I_n), so G_s = M'M and G_s^{otimes m} = (M^{otimes m})' M^{otimes m}.
%  Since (e^{Au})^{otimes m} = e^{A^[m] u} (A^[m] the m-fold Kronecker sum,
%  Lemma lem:high:lift), the interval of length D starting from N
%  contributes EXACTLY
%      (N^{otimes m})' [ W_m(D) + E' E ] N^{otimes m},
%      E = e^{A^[m] D},   W_m(D) = W_inf - E' W_inf E,
%      A^[m]' W_inf + W_inf A^[m] = -I_{n^m}
%  (integral term plus the jump-sum term at the end of the interval), and
%  the periodic tail is summed in closed form by a Stein equation.  No
%  time discretization enters the Gramians.
%
%  ADMISSIBLE SEQUENCES (Definition def:high:families, time measured from
%  the current instant, eq:high:futures: sigma = {t_1 < t_2 < ...} in
%  [0,inf), only Delta_k = t_{k+1} - t_k >= Tbar is imposed, t_1 = 0 means
%  a jump at the current instant).  GRID: sigma = (t_1; Delta_1, Delta_2;
%  periodic tail Delta_k = Dtail for k >= 3) with
%      t_1   in T1grid = {0, 0.1, 0.3, 0.6, 1, 1.5, 2.5, 4, 8}        (9)
%      Delta_1, Delta_2, Dtail in Dgrid = {0.3, 0.4, 0.6, 1, 1.5, 2.5, 4, 8} (8)
%  for 9*8*8*8 = 4608 sequences, all in S_min[Tbar].  The grid covers the
%  transient peak of ||e^{At}||^2 (about 55 at t = 1.67) on both sides.
%
%  GREATEST-ELEMENT TEST at each order m = 1..6 (Loewner order of the full
%  space H^{n^m}, n^m = 2..64).  Pairwise domination P_i >= P_j is tested
%  as lambda_min(P_i - P_j) >= -tol*lambda_max(P_j), tol = 1e-9 (relative
%  tolerance: the Gramians grow like 55^m).  A grid Gramian dominating all
%  others has the largest trace, so the test is run, against all 4608 grid
%  Gramians, for every trace maximizer (trace within the relative tolerance
%  of the maximum): this is the complete pairwise test for the property
%  "some grid Gramian dominates all others", at linear cost.  The worst
%  relative violation of the best candidate is reported (negative means no
%  greatest element), with the sequence that breaks it.  Independently,
%  the maximizers of <P, Z> over the grid are computed for 200 positive
%  definite test directions Z (I_{n^m}, the product states Y^{otimes m} of
%  the Y below, and random Wishart matrices): each maximizer for Z > 0 is a
%  maximal element of the grid set, and two distinct maximizers already
%  exclude a greatest element.
%
%  PRODUCT-STATE RATIO (informative, eq:high:betaeffdom on the grid): with
%  Y on a polar grid of the real base slice {Y >= 0, tr Y = 1} (97 points;
%  real Y suffice since the data are real), the smallest 1+eps for which
%  eq:high:betaeffdom holds on the grid at order m is
%      r_m = min_{sigma*} max_{sigma, Y}
%            ( <P^[m]_sigma, Y^{otimes m}> / <P^[m]_{sigma*}, Y^{otimes m}> )^{1/m}.
%  The ratio Pi_*(Y)/Pi_{sigma*}(Y) of eq:high:betaeff does not depend on m;
%  its grid value beta_eff (peaks eq:high:peak sampled with step 0.005,
%  rigorous stopping rule below) is reported as well.  Neither quantity is
%  gated.  Observed: beta_eff = 1 on the grid, while r_m is
%  NOT monotone in m (about 1.74, 1.67, 1.81, 1.83, 1.77, 1.69).
%
%  GATES (stop on violation): (G1) rho(J e^{A D}) < 1 for every D in Dgrid
%  (the tail series converge); (G2) the instance is uniformly globally
%  exponentially stable over S_min[Tbar], certified by a constant linear
%  certificate P > 0 with A'P + PA < 0 and J'e^{A'Tbar} P e^{A Tbar}J - P < 0
%  (MathWorks LMI Lab feasp, eigenvalues re-checked); (G3) the closed-form
%  Gramians satisfy the power-moment identity eq:high:moment, checked against
%  an independent forward propagation in svec coordinates
%  (mvst_liftsym/mvst_svecbasis) with adaptive quadrature, m = 1..3;
%  (G4) every grid Gramian is positive definite; (G5, observation /
%  regression) no grid Gramian dominates all others at any m = 1..6: a
%  failure here CONTRADICTS the chapter's claim.  Runtime: about 30 s.
%  No figure, no data macro.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)

% ---------------- instance, exactly as displayed in the chapter -----------
A    = [-0.5 12; 0 -0.7];
J    = [0.8 0.2; 0 0.5];
Tbar = 0.3;
n    = 2;
mmax = 6;
tol  = 1e-9;                      % relative Loewner tolerance

% ---------------- grid of admissible sequences ----------------------------
T1grid = [0 0.1 0.3 0.6 1.0 1.5 2.5 4.0 8.0];
Dgrid  = [0.3 0.4 0.6 1.0 1.5 2.5 4.0 8.0];
[i1,i2,i3,i4] = ndgrid(1:numel(T1grid), 1:numel(Dgrid), ...
                       1:numel(Dgrid), 1:numel(Dgrid));
SEQ = [T1grid(i1(:))', Dgrid(i2(:))', Dgrid(i3(:))', Dgrid(i4(:))'];
Ns  = size(SEQ,1);                % 4608 sequences (t_1, D_1, D_2, Dtail)
assert(all(all(SEQ(:,2:4) >= Tbar)), 'grid sequence outside S_min[Tbar]');

% ---------------- (G1) convergence of the periodic tails ------------------
rhoTail = arrayfun(@(D) max(abs(eig(J*expm(A*D)))), Dgrid);
assert(all(rhoTail < 1), 'G1: a periodic tail does not converge');

% ---------------- (G2) stability over S_min[Tbar] -------------------------
E0 = expm(A*Tbar);
setlmis([]);  Pv = lmivar(1,[n 1]);
lmiterm([1 1 1 Pv], A', 1, 's');                       % A'P + PA < 0
lmiterm([2 1 1 Pv], (E0*J)', E0*J); lmiterm([2 1 1 Pv], -1, 1);
lmiterm([-3 1 1 Pv], 1, 1); lmiterm([3 1 1 0], 1);     % P > I
lmis = getlmis;
[tminS, xS] = feasp(lmis, [1e-6 200 0 0 1]);
Pc = dec2mat(lmis, xS, Pv);  Pc = (Pc+Pc')/2;
stabMargin = max([eig(A'*Pc + Pc*A); eig(J'*E0'*Pc*E0*J - Pc)]);
assert(tminS < 0 && stabMargin < 0 && min(eig(Pc)) > 0, ...
    'G2: constant minimum dwell-time certificate not found at Tbar');

% ---------------- transfer matrices of the grid sequences -----------------
% N1 = J e^{A t_1}, N2 = J e^{A D_1} N1, N3 = J e^{A D_2} N2
eA  = @(t) expm(A*t);
N1 = zeros(n,n,Ns); N2 = N1; N3 = N1;
for s = 1:Ns
    N1(:,:,s) = J*eA(SEQ(s,1));
    N2(:,:,s) = J*eA(SEQ(s,2))*N1(:,:,s);
    N3(:,:,s) = J*eA(SEQ(s,3))*N2(:,:,s);
end

% ---------------- product-state directions (real base slice) --------------
[rr,th] = ndgrid([0.25 0.5 0.75 1], (0:23)*pi/12);
bx = [0; rr(:).*cos(th(:))];  bz = [0; rr(:).*sin(th(:))];   % 97 points
NY = numel(bx);
Ys = cell(1,NY);
for y = 1:NY
    Ys{y} = [(1+bz(y))/2, bx(y)/2; bx(y)/2, (1-bz(y))/2];
end
rng(20261025);                     % reproducible random test directions

% ---------------- main loop over the lifting order ------------------------
res = struct('m',{},'dim',{},'nCand',{},'bestViol',{},'bestSeq',{}, ...
             'witSeq',{},'nMax',{},'r',{},'minEig',{});
Pstore = cell(1,3);                % orders 1..3 kept for gate G3
tStart = tic;
for m = 1:mmax
    d   = n^m;
    Am  = ksum(A,m);
    Winf = lyap(Am', eye(d));  Winf = (Winf+Winf')/2;
    Qint = @(D) qint(kpow(eA(D),m), Winf);           % W_m(D) + E'E
    % tail sums S(Dtail) = sum_i L'^i Q L^i, L = (J e^{A Dtail})^{otimes m}
    Stail = cell(1,numel(Dgrid));
    for k = 1:numel(Dgrid)
        L = kpow(J*eA(Dgrid(k)), m);
        S = dlyap(L', Qint(Dgrid(k)));  Stail{k} = (S+S')/2;
    end
    Qt1 = cell(1,numel(T1grid));  for k=1:numel(T1grid), Qt1{k}=Qint(T1grid(k)); end
    Qd  = cell(1,numel(Dgrid));   for k=1:numel(Dgrid),  Qd{k} =Qint(Dgrid(k));  end
    % lifted Gramians of all grid sequences, stored as columns vec(P)
    PV = zeros(d*d, Ns);
    for s = 1:Ns
        K1 = kpow(N1(:,:,s),m); K2 = kpow(N2(:,:,s),m); K3 = kpow(N3(:,:,s),m);
        P = Qt1{i1(s)} + K1'*Qd{i2(s)}*K1 + K2'*Qd{i3(s)}*K2 ...
            + K3'*Stail{i4(s)}*K3;
        P = (P+P')/2;
        PV(:,s) = P(:);
    end
    if m <= 3, Pstore{m} = PV; end
    % (G4) positivity of every grid Gramian
    minEig = inf;
    for s = 1:Ns, minEig = min(minEig, min(eig(reshape(PV(:,s),d,d)))); end
    % greatest-element test: every trace maximizer against all grid Gramians
    tr  = sum(PV(1:d+1:end,:),1);                     % traces
    cand = find(tr >= max(tr)*(1-tol));
    bestViol = -inf; bestSeq = NaN; witSeq = NaN;
    for c = cand
        Pc_ = reshape(PV(:,c),d,d);  worst = inf;  wj = NaN;
        for j = 1:Ns
            Pj = reshape(PV(:,j),d,d);
            v  = min(eig(Pc_-Pj)) / max(eig(Pj));
            if v < worst, worst = v; wj = j; end
        end
        if worst > bestViol, bestViol = worst; bestSeq = c; witSeq = wj; end
    end
    % maximizers of <P,Z> over positive definite test directions Z
    Zs = zeros(d*d, 200);
    Zs(:,1) = reshape(eye(d), [], 1);
    for y = 1:NY, Zs(:,1+y) = reshape(kpow(Ys{y},m), [], 1); end
    for z = NY+2:200
        G = randn(d,d+2); Z = G*G'/(d+2) + 1e-3*eye(d); Zs(:,z) = Z(:);
    end
    [~, amax] = max(PV'*Zs, [], 1);
    nMax = numel(unique(amax));
    % product-state ratio r_m (informative)
    V  = PV' * Zs(:,2:NY+1);                         % Ns x NY, <P,Y^{otimes m}>
    Vmax = max(V, [], 1);
    r  = min(max(Vmax ./ V, [], 2))^(1/m);
    res(m) = struct('m',m,'dim',d,'nCand',numel(cand),'bestViol',bestViol, ...
        'bestSeq',bestSeq,'witSeq',witSeq,'nMax',nMax,'r',r,'minEig',minEig);
end
tRun = toc(tStart);

% ---------------- (G3) power-moment identity, independent check -----------
% forward propagation of svec(X) with the reduced flow and jump operators,
% phi(s) = tr X(s) = <I_n, X(s)>, quadrature of phi^m over each interval
U    = mvst_svecbasis(n);
Lf   = mvst_liftsym(A, {});                  % svec matrix of X -> AX + XA'
Lj   = mvst_liftsym(zeros(n), {J});          % svec matrix of X -> J X J'
uI   = U'*reshape(eye(n),[],1);              % svec(I_n)
[Vf, Df] = eig(Lf);  lamf = diag(Df);        % distinct eigenvalues -1,-1.2,-1.4
assert(cond(Vf) < 1e8, 'reduced flow generator not diagonalizable');
cI   = uI'*Vf;                               % row vector
chkSeq = [1, Ns, 1+floor(Ns/3), 1+floor(2*Ns/3)];
Ychk = {[0.7 0.2; 0.2 0.3], [0.5 -0.45; -0.45 0.5]};
momRes = 0;
for m = 1:3
    for s = chkSeq
        for y = 1:numel(Ychk)
            Y = Ychk{y};
            lhs = Pstore{m}(:,s)' * reshape(kpow(Y,m),[],1);
            dts = [SEQ(s,1), SEQ(s,2), SEQ(s,3), repmat(SEQ(s,4),1,400)];
            x = U'*Y(:);  rhs = 0;
            for k = 1:numel(dts)
                D = dts(k);
                if D > 0                              % phi(t) = c' e^{Lam t} w
                    w = Vf\x;
                    f = @(t) (real((cI.*w.')*exp(lamf*t))).^m;
                    rhs = rhs + integral(f, 0, D, 'RelTol',1e-12,'AbsTol',0);
                end
                x = expm(Lf*D)*x;                    % pre-jump state
                rhs = rhs + (uI'*x)^m;                % jump-sum term
                x = Lj*x;                             % jump
                if k > 3 && (uI'*x)^m < 1e-18*rhs, break; end
            end
            momRes = max(momRes, abs(lhs-rhs)/abs(rhs));
        end
    end
end

% ---------------- peak functionals and beta_eff on the grid (informative) ---
% Pi_sigma(Y) = sup_s phi_sigma(s;Y) (eq:high:peak), sampled with step
% hs on each flow interval; the propagation of a sequence stops once
% gmax*phi(t_k^+;Y) < Pi_sigma(Y) for every Y, gmax = sup_u ||e^{Au}||^2,
% since phi(t_k^+ + u;Y) <= ||e^{Au}||^2 phi(t_k^+;Y) until the next jump
% and the jumps do not increase the trace (J'J <= I_n).
assert(max(eig(J'*J)) <= 1, 'J''J <= I_n fails');
gmax = max(arrayfun(@(u) norm(eA(u))^2, 0:0.005:20));
hs   = 0.005;
Xb   = zeros(3,NY);  for y=1:NY, Xb(:,y) = U'*Ys{y}(:); end
trY  = uI'*Xb;                                     % = 1 on the base slice
PI   = zeros(Ns,NY);
for s = 1:Ns
    x = Xb;  pk = trY;  k = 0;
    while true
        k = k + 1;
        if k <= 3, D = SEQ(s,k); else, D = SEQ(s,4); end
        if D > 0
            tt  = 0:hs:D;  if tt(end) < D, tt = [tt D]; end %#ok<AGROW>
            ph  = (exp(lamf*tt).'.*cI) * (Vf\x);     % numel(tt) x NY
            pk  = max(pk, max(real(ph),[],1));
        end
        x = Lj*(expm(Lf*D)*x);                       % flow D, then jump
        if k >= 3 && all(gmax*(uI'*x) < pk), break; end
    end
    PI(s,:) = pk;
end
PIstar  = max(PI,[],1);                            % family peak on the grid
[betaEffGrid, sStar] = min(max(PIstar./PI,[],2));

% ---------------- gates -----------------------------------------------------
assert(momRes < 1e-8, 'G3: power-moment identity residual %.2e', momRes);
assert(all([res.minEig] > 0), 'G4: a grid Gramian is not positive definite');

% ---------------- report ----------------------------------------------------
fprintf('EX_HIGH_TRANSIENT: A=[-0.5 12;0 -0.7], J=[0.8 0.2;0 0.5], Tbar=%.1f\n', Tbar);
fprintf('  grid: %d sequences (t_1 in %d values, Delta_1, Delta_2, Dtail in %d values)\n', ...
    Ns, numel(T1grid), numel(Dgrid));
fprintf('  G1 tail spectral radii in [%.4f, %.4f]; G2 constant certificate at Tbar: margin %.3e\n', ...
    min(rhoTail), max(rhoTail), stabMargin);
fprintf('  G3 power-moment identity, max relative residual %.2e; run time %.1f s\n', ...
    momRes, tRun);
fprintf('   m  n^m  #trace-max  best viol. (rel.)  breaking seq. (t1,D1,D2,Dtail)   #maximizers/200   r_m\n');
for m = 1:mmax
    q = res(m);
    fprintf('  %2d  %3d  %6d      %+.3e        (%g, %g, %g, %g) vs (%g, %g, %g, %g)   %4d   %.4f\n', ...
        q.m, q.dim, q.nCand, q.bestViol, SEQ(q.bestSeq,:), SEQ(q.witSeq,:), q.nMax, q.r);
end
noGreatest = all([res.bestViol] < -tol);
if noGreatest
    fprintf(['  OBSERVATION (consistent with the chapter): no grid Gramian dominates ' ...
             'all others at any order m = 1..%d.\n'], mmax);
else
    bad = find([res.bestViol] >= -tol);
    fprintf(2, ['  *** CONTRADICTS THE CHAPTER: a grid Gramian dominating all others ' ...
                'exists at order(s) %s ***\n'], mat2str(bad));
end
rDec = all(diff([res.r]) < 0);
fprintf('  product-state ratio r_m decreasing in m: %d  (r_m = %s)\n', ...
    rDec, mat2str([res.r],5));
fprintf(['  grid beta_eff = min_sigma* max_Y Pi_*(Y)/Pi_sigma*(Y) = %.6f, ' ...
         'attained at sigma* = (%g, %g, %g, %g)\n'], betaEffGrid, SEQ(sStar,:));
% (G5) observation / regression gate
if ~noGreatest
    error('mvst:high:gate', ...
        'G5: a greatest element was found on the grid (contradicts Ch13 claim).');
end
fprintf('  PASS: gates G1-G5.\n');

% ======================= local functions ==============================
function Q = qint(E, Winf)
%QINT Integral term W_m(D) = Winf - E'Winf E plus the jump-sum term E'E.
Q = Winf - E'*Winf*E + E'*E;  Q = (Q+Q')/2;
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
