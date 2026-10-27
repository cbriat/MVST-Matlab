%% EX_HIGH_HIERARCHY  Minimum dwell-time hierarchy from the tensor-power lift.
%  Chapter "Higher-order systems", numerical illustration of
%  Proposition prop:high:monotone for the order-m programs WITHOUT
%  multiplier on the coupling inequality (the case k=0, in the
%  mode-dependent form below).  The convergence of these estimates as m
%  grows is NOT established: Theorem th:high:approach proves convergence
%  for the multiplier-augmented estimates T_M of eq:high:TM, which this
%  script does not compute; the decrease seen here is observed
%  numerically.  The base is
%  the switched vector system xdot = A_sigma x (worst-case, rank-one reading
%  of the matrix system, prob:high:switched) with the two Hurwitz modes
%      A1 = [-1 -1; 1 -1],   A2 = [-1 -10; 1/10 -1],
%  both with spectrum {-1 +- i}.  Their product A1*A2 has negative real
%  eigenvalues, so by the Shorten-Narendra pencil criterion no COMMON
%  quadratic Lyapunov function exists: the order-one certificate is
%  conservative.
%
%  For m = 1,2,3 the order-m lift replaces the quadratic x'Px by the
%  homogeneous degree-2m form (x^{otimes m})' P (x^{otimes m}), whose
%  generator on the symmetric power is the m-fold Kronecker sum
%  A_i (+) ... (+) A_i (Corollary cor:high:hplf).  The certified minimum
%  dwell-time estimate
%      T_m := inf{ T : order-m program feasible at dwell-time T }
%  (mode-dependent Geromel-Colaneri LMI without multiplier, solved with the
%  MathWorks LMI Lab feasp, feasible iff feasp returns tmin < 0) is obtained
%  by bisection at each order on the bracket [0, 5]: T = 0 is checked
%  infeasible (at T = 1e-6; if it is feasible, T_m = 0), T = 5 is checked
%  feasible, and the bracket is halved until its width is below 1e-6; the
%  reported T_m is the upper end of the final bracket, a dwell-time at which
%  the program is feasible.  Bisection is valid because feasibility is
%  monotone in T: if E_i(T)'Q_j E_i(T) < Q_i and S_i'Q_i + Q_i S_i < 0, then
%  E_i(T+s)'Q_j E_i(T+s) < e^{S_i's} Q_i e^{S_i s} <= Q_i for s >= 0.
%  The exact minimum dwell-time is T* = 0: the system is stable under
%  arbitrary switching, certified by a COMMON homogeneous polynomial
%  Lyapunov function of finite degree, found here at order m* (degree 2 m*).
%
%  GATES (assert, stop on violation): soundness, every certified T_m is an
%  upper bound on the exact value (T_m >= T* = 0); T_2, T_3 <= T_1 (orders
%  2 and 3 are multiples of 1, Proposition prop:high:monotone); T_3 <= T_2
%  is not covered by that proposition and is checked as an observed
%  property (regression gate).
%
%  Figure figures/high_hierarchy.pdf (the estimates T_1, T_2, T_3 against
%  the order m, with the exact value T* = 0 as a reference line), data
%  figures/data/high_hierarchy.tex.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

A1 = [-1 -1; 1 -1];
A2 = [-1 -10; 1/10 -1];
n  = 2;

% ---- Shorten-Narendra: no common quadratic (eig(A1 A2) negative real) ----
eProd = sort(real(eig(A1*A2)));            % both negative here

% ---- exact minimum dwell-time --------------------------------------------
% The exact minimum dwell-time is T* = 0: the base is stable under ARBITRARY
% switching, certified rigorously below by the existence of a COMMON
% homogeneous polynomial Lyapunov function of finite degree (found at order
% m*), so every dwell-time, including zero, is admissible.
Tstar = 0.0;

% ---- order-m minimum dwell-time by bisection (mode-dependent LMI) --------
% Bisection on the bracket [0,5] with tolerance 1e-6; the reported estimate is the feasible upper end.
Tlo = 0;  Thi = 5;  Ttol = 1e-6;      % bisection bracket and tolerance
Tm   = zeros(1,3);  dims = zeros(1,3);
for m = 1:3
    W = symb(n,m);  dims(m) = size(W,2);
    S = {W'*ksum(A1,m)*W, W'*ksum(A2,m)*W};
    if mdt_feas(S,1e-6) < 0            % common certificate already suffices
        Tm(m) = 0;
    else
        lo=Tlo; hi=Thi;
        assert(mdt_feas(S,hi) < 0, 'order-%d infeasible at T=%.1f', m, hi);
        while hi-lo > Ttol
            mid=(lo+hi)/2;
            if mdt_feas(S,mid) < 0, hi=mid; else, lo=mid; end
        end
        Tm(m) = hi;                    % feasible end of the final bracket
    end
end

% ---- first order with a COMMON homogeneous polynomial certificate --------
mstar = NaN;
for m = 1:6
    W = symb(n,m);  S = {W'*ksum(A1,m)*W, W'*ksum(A2,m)*W};
    if comm_lyap(S) < 0, mstar = m; break; end
end
degstar = 2*mstar;

% ---- soundness gates -----------------------------------------------------
if ~all(Tm >= Tstar - 1e-9)
    error('mvst:high:gate', 'SOUNDNESS: some T_m < exact value T*=%.3f.', Tstar);
end
if ~(Tm(1) >= Tm(2)-Ttol && Tm(1) >= Tm(3)-Ttol)
    error('mvst:high:gate', 'MONOTONICITY: T_2 or T_3 > T_1: %s.', mat2str(Tm,4));
end
if ~(Tm(2) >= Tm(3)-Ttol)
    error('mvst:high:gate', 'REGRESSION: observed T_3 <= T_2 lost: %s.', mat2str(Tm,4));
end
assert(all(eProd < 0), 'Shorten-Narendra: eig(A1 A2) not both negative.');
assert(~isnan(mstar), 'no common polynomial certificate found up to order 6.');

fprintf('EX_HIGH_HIERARCHY: eig(A1 A2)=[%.4f %.4f]\n', eProd);
fprintf('  T_1=%.6f T_2=%.6f T_3=%.6f  (exact T*=0), common poly at m*=%d (deg %d)\n', ...
    Tm, mstar, degstar);

% ---- figure: certified estimates T_m against the order m -----------------
% Certified estimates T_m against the order m, with the exact value T* = 0 as reference line.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'FontSize',12);
hr = yline(ax, Tstar, '--', 'LineWidth',1.6, 'Color', gc.order(2,:));
hm = plot(ax, 1:3, Tm, 'o', 'MarkerSize',9, 'LineWidth',1.8, ...
    'MarkerEdgeColor', gc.order(1,:), 'MarkerFaceColor', gc.order(1,:));
for m = 1:3                              % value labels above the markers
    text(ax, m, Tm(m)+0.045, sprintf('$%.3f$', Tm(m)), ...
        'Interpreter','latex','FontSize',12, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
end
xlim(ax, [0.5 3.5]);  ylim(ax, [-0.05 0.6]);
set(ax,'XTick',1:3);
xlabel(ax, 'lifting order $m$', 'Interpreter','latex');
ylabel(ax, 'minimum dwell-time', 'Interpreter','latex');
lg = legend(ax, [hm hr], ...
    {'certified estimate $\bar T_m$', ...
     sprintf('exact value $\\bar T^\\star=%g$', Tstar)}, ...
    'Interpreter','latex','Location','northeast','FontSize',12);
set(lg,'Box','off');
mvst_savepdf(fig, fullfile(mvst_figdir(),'high_hierarchy.pdf'));

% ---- data macros ---------------------------------------------------------
w = mvst_texdata(fullfile(mvst_figdir(),'data','high_hierarchy.tex'));
w.num('HighHierTone',   Tm(1), 3);
w.num('HighHierTtwo',   Tm(2), 3);
w.num('HighHierTthree', Tm(3), 3);
w.num('HighHierTstar',  Tstar, 3);
w.num('HighHierMstar',  mstar, 0);
w.num('HighHierDegStar',degstar, 0);
w.num('HighHierProdOne',eProd(1), 3);
w.num('HighHierProdTwo',eProd(2), 3);
w.num('HighHierDimOne',   dims(1), 0);
w.num('HighHierDimTwo',   dims(2), 0);
w.num('HighHierDimThree', dims(3), 0);
w.close();

% ======================= local functions ==============================
function tmin = comm_lyap(S)
%COMM_LYAP Feasibility (feasp tmin) of a COMMON quadratic S_i'Q+Q S_i<0, Q>I.
d=size(S{1},1);
setlmis([]);  Q=lmivar(1,[d 1]);
for i=1:numel(S), lmiterm([i 1 1 Q],S{i}',1,'s'); end
k=numel(S)+1;  lmiterm([-k 1 1 Q],1,1);  lmiterm([k 1 1 0],1);
tmin = feasp(getlmis, [1e-6 200 0 0 1]);
end
function tmin = mdt_feas(S,T)
%MDT_FEAS Feasibility (feasp tmin) of the mode-dependent minimum dwell-time
%  LMI at dwell-time T: S_i'Q_i+Q_i S_i<0, E_i'Q_j E_i - Q_i<0 (i~=j), Q_i>I.
d=size(S{1},1); N=numel(S);  E=cell(1,N);
for i=1:N, E{i}=expm(S{i}*T); end
setlmis([]);  Q=cell(1,N);
for i=1:N, Q{i}=lmivar(1,[d 1]); end
c=0;
for i=1:N, c=c+1; lmiterm([c 1 1 Q{i}],S{i}',1,'s'); end
for i=1:N, for j=1:N
    if i==j, continue; end
    c=c+1; lmiterm([c 1 1 Q{j}],E{i}',E{i}); lmiterm([c 1 1 Q{i}],-1,1);
end, end
for i=1:N, c=c+1; lmiterm([-c 1 1 Q{i}],1,1); lmiterm([c 1 1 0],1); end
tmin = feasp(getlmis, [1e-6 200 0 0 1]);
end
function W = symb(d,m)
if m==1, W=eye(d); return; end
P=zeros(d^m); pm=perms(1:m);
for p=1:size(pm,1), P=P+permop(d,m,pm(p,:)); end
P=P/size(pm,1);
[U,S]=svd(P); W=U(:, diag(S)>0.5);
end
function Q=permop(d,m,sigma)
N=d^m; Q=zeros(N); idx=cell(1,m);
for lin=1:N
    [idx{:}]=ind2sub(repmat(d,1,m),lin); v=cell2mat(idx);
    w=zeros(1,m); w(sigma)=v; wc=num2cell(w);
    Q(sub2ind(repmat(d,1,m),wc{:}),lin)=1;
end
end
function G=ksum(L,m)
d=size(L,1); G=zeros(d^m);
for i=1:m
    t=1; for k=1:m, if k==i, t=kron(t,L); else, t=kron(t,eye(d)); end; end
    G=G+t;
end
end
