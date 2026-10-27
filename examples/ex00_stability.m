%% EX00_STABILITY  Stability tests for the four running system classes.
%  Chapter "Stability Analysis".  The stability of each system is now CERTIFIED
%  through the generic matrix sum-of-squares analysis programs (MVST_STAB_SOS for
%  the single-certificate impulsive systems, MVST_STAB_SOS_MODES for the
%  product-cone / mode-indexed ones), and every CRITICAL PARAMETER is found by
%  BISECTION ON THE SOS FEASIBILITY.  The plotted curves remain exact spectral
%  readings of the impulsive Lyapunov operator, computed from the REDUCED matrix
%  representations on the symmetric subspace (mvst_liftsym/mvst_svecbasis); the
%  exact spectral crossing is the CROSS-CHECK of the SOS critical value (they
%  must agree to ~1e-2 relative):
%    (a) MJLS with diffusion: coupled lifted generator abscissa, cross-checked
%        by the product-cone SOS certificate (feasible iff abscissa < 0);
%    (b) GBM + Poissonian jumps: exact abscissa versus the jump rate lambda;
%        critical rate lambda* by SOS bisection (arbitrary dwell-time);
%    (c) sampled-data: exact monodromy radius versus the sampling period;
%        maximum stable period Tbar* by SOS bisection (periodic dwell-time),
%        plus a ranged-dwell certificate of an admissible interval;
%    (d) switched system: exact monodromy radius over one switching cycle;
%        minimum dwell time T* by mode-indexed SOS bisection.
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ---------- (a) MJLS with diffusion: coupled lifted generator ----------
A1 = [0 1;-2 -1];   A2 = [0 1;-3 -1.5];
Cd = 0.15*eye(2);   Pig = [-1 1; 2 -2];
n = 2; nr = n*(n+1)/2;                          % reduced dimension
Lc = {mvst_liftsym(A1,{Cd}), mvst_liftsym(A2,{Cd})};
Lmj = zeros(2*nr);
for i = 1:2
    for j = 1:2
        blk = Pig(j,i)*eye(nr);                 % coefficient of X_j in Xdot_i
        if i == j; blk = blk + Lc{i}; end
        Lmj((i-1)*nr+1:i*nr, (j-1)*nr+1:j*nr) = blk;
    end
end
absMjls = max(real(eig(Lmj)));
% product-cone SOS certificate (coupled constant certificates {P_i}): feasible
% iff the coupled abscissa is negative (mean-square stable MJLS).
resMjls = mvst_stab_sos_modes( ...
    struct('A0',cat(3,A1,A2),'Ac',cat(4,Cd,Cd),'Pi',Pig), struct('type','ctmjls'));
fprintf('Stability (MJLS): coupled lifted abscissa %.4f  (SOS feasible=%d, slack %.2e)\n', ...
    absMjls, resMjls.feas, resMjls.slack);
assert(resMjls.feas == (absMjls < 0), ...
    'MJLS product-cone SOS feasibility inconsistent with the coupled abscissa sign');

% ---------- (a') DT-MJLS analysis: both modes UNSTABLE, chain stabilizes ----
% Two FULL (non-triangular) modes, each with spectral radius 1.5 (an eigenvalue
% outside the unit disk, so each mode is discrete-time UNSTABLE), a common
% state-multiplicative channel 0.2 I2, and a fast-alternating Markov chain.  The
% chain, not the modes, carries mean-square stability: the coupled second-moment
% lift is a contraction (radius < 1) even though every mode expands.  Under the
% i.i.d. chain the same modes are mean-square UNSTABLE, the converse warning.
Ad1 = [1 0.5; 0.5 1];   Ad2 = [1 -0.5; -0.5 1];   % full, non-triangular, unstable
Cdmj = 0.2*eye(n);                                 % common multiplicative channel
Pdt  = [0.1 0.9; 0.9 0.1];                         % fast-alternating transition
rhoDt = [max(abs(eig(Ad1))), max(abs(eig(Ad2)))];  % per-mode spectral radii (1.5)
radDtMjls    = dtmjls_msrad(cat(3,Ad1,Ad2), Cdmj, Pdt);
radDtMjlsIid = dtmjls_msrad(cat(3,Ad1,Ad2), Cdmj, [0.5 0.5; 0.5 0.5]);
% coupled DT-MJLS SOS certificate: feasible iff the coupled radius is below 1.
resDt = mvst_stab_sos_modes( ...
    struct('A0',cat(3,Ad1,Ad2),'Ac',cat(4,Cdmj,Cdmj),'P',Pdt), struct('type','dtmjls'));
fprintf(['Stability (DT-MJLS): rho(A1)=%.3f rho(A2)=%.3f (both unstable), coupled ', ...
    'MS radius %.4f (SOS feas=%d, slack %.2e); i.i.d.-chain radius %.4f\n'], ...
    rhoDt(1), rhoDt(2), radDtMjls, resDt.feas, resDt.slack, radDtMjlsIid);
assert(all(rhoDt > 1), 'DT-MJLS modes must both be unstable (spectral radius > 1)');
assert(resDt.feas == (radDtMjls < 1), ...
    'DT-MJLS SOS feasibility inconsistent with the coupled second-moment radius');
assert(radDtMjls < 1 && radDtMjlsIid > 1, ...
    'DT-MJLS story broken: alternating chain must stabilize while i.i.d. chain does not');
% Converse warning: the nilpotent (Schur) pair S1, S2 under the same chain Pdt,
% modes alone (no multiplicative channel), is mean-square UNSTABLE.
S1 = [0 2; 0 0];   S2 = [0 0; 2 0];
radDtMjlsNil = dtmjls_msrad(cat(3,S1,S2), zeros(n), Pdt);
fprintf('Stability (DT-MJLS, nilpotent pair): coupled MS radius %.4f\n', radDtMjlsNil);
assert(radDtMjlsNil > 1, 'nilpotent-pair DT-MJLS must be mean-square unstable');

% ---------- (b) GBM + Poissonian jumps: abscissa versus rate ----------
A  = [0 1;-2 -1];   C1 = 0.15*eye(2);   J1 = 1.15*eye(2);
lamNom = 0.5;
L0 = mvst_liftsym(A,{C1});
LJ = mvst_liftsym(zeros(n),{J1}) - eye(nr);
lams = linspace(0, 2.5, 120); absc = zeros(size(lams));   % exact abscissa curve
for k = 1:numel(lams)
    absc(k) = max(real(eig(L0 + lams(k)*LJ)));
end
absPoi  = max(real(eig(L0 + lamNom*LJ)));
% critical Poisson rate lambda* by BISECTION ON THE SOS FEASIBILITY: the
% Brownian channel C1 enters the flow as .Ac, the Poisson resets as the
% spontaneous jumps .Jf weighted by the rate .lamf (dwell type 'arbitrary');
% stable for small lambda, destabilized by the expansive reset as lambda grows.
sfeasPoi = @(lam) mvst_stab_sos( ...
    struct('A0',A,'Ac',C1,'Jf',J1,'lamf',lam), struct('type','arbitrary')).slack;
lamHi = 4.0;                                            % bracket the true crossing
lamCrit = sos_bisect(sfeasPoi, 0, lamHi, 'feaslo');     % last feasible lambda
% exact abscissa crossing (analytic bisection over the same bracket; the
% crossing lies beyond the plotted [0,2.5] window, so the on-grid crossing is
% empty, hence the wider bracket here)
lamgrid = linspace(0, lamHi, 401);
lamCritExact = crossing(lamgrid, arrayfun(@(l)max(real(eig(L0+l*LJ))), lamgrid), 0);
fprintf('Stability (Poisson): abscissa %.4f at lambda %.2f, critical rate %.4f (SOS), %.4f (exact), rel. dev. %.2e\n', ...
    absPoi, lamNom, lamCrit, lamCritExact, abs(lamCrit-lamCritExact)/lamCritExact);
assert(abs(lamCrit-lamCritExact)/lamCritExact < 1e-2, ...
    'SOS critical Poisson rate disagrees with the exact abscissa crossing beyond 1e-2 relative');

% ---------- (c) sampled-data: monodromy radius versus period ----------
D = ex_data(); s = D.sd; TbarNom = s.Tbar;
out = mvst_igre_periodic_sd(s.A,s.B,s.Ac,s.Bc,s.Q,s.R,TbarNom);
K = out.K; m = size(s.B,2); na = n + m;
Ab0 = [s.A s.B; zeros(m,na)];
Abl = [s.Ac s.Bc; zeros(m,na)];
Jb  = [eye(n) zeros(n,m); K zeros(m,m)];
Lsd = mvst_liftsym(Ab0,{Abl});
Jlift = mvst_liftsym(zeros(na),{Jb});
Tbars = linspace(0.05, 2.0, 120); rsd = zeros(size(Tbars));   % exact radius curve
for k = 1:numel(Tbars)
    rsd(k) = max(abs(eig(Jlift*expm(Tbars(k)*Lsd))));
end
% maximum stable sampling period Tbar* by BISECTION ON THE SOS FEASIBILITY:
% augmented flow Ab0 is .A0, the Brownian channel Abl is .Ac, the sampling map
% Jb is the scheduled jump .J (dwell type 'periodic', degree-4 certificate).
sfeasSd = @(Tb) mvst_stab_sos( ...
    struct('A0',Ab0,'Ac',Abl,'J',Jb), struct('type','periodic','Tbar',Tb,'deg',4)).slack;
TbarCrit = sos_bisect(sfeasSd, 0.05, 2.0, 'feaslo');    % last feasible period
TbarCritExact = crossing(Tbars, rsd, 1);                % exact radius crossing
fprintf('Stability (sampled-data): nominal Tbar %.2f, critical period %.4f (SOS), %.4f (exact), rel. dev. %.2e\n', ...
    TbarNom, TbarCrit, TbarCritExact, abs(TbarCrit-TbarCritExact)/TbarCritExact);
assert(abs(TbarCrit-TbarCritExact)/TbarCritExact < 1e-2, ...
    'SOS max sampling period disagrees with the exact monodromy crossing beyond 1e-2 relative');
% ranged-dwell certificate (th:stab:RDT): an admissible aperiodic interval
sdTmin = 0.05; sdTmax = 0.85*TbarCrit;
sRangeSd = mvst_stab_sos(struct('A0',Ab0,'Ac',Abl,'J',Jb), ...
    struct('type','range','Tmin',sdTmin,'Tmax',sdTmax,'deg',4)).slack;
fprintf('  ranged-dwell admissible sampling interval [%.4f, %.4f] certified (slack %.2e)\n', ...
    sdTmin, sdTmax, sRangeSd);
assert(sRangeSd > 1e-7, 'ranged-dwell certificate infeasible on the reported interval');

% ---------- (d) switched system: monodromy radius versus dwell ----------
Asw1 = [-1 4; 0 -1];   Asw2 = [-1 0; 4 -1];      % Hurwitz, unstable average
Lsw1 = mvst_liftsym(Asw1,{}); Lsw2 = mvst_liftsym(Asw2,{});
Tds = linspace(0.02, 3.5, 160); rsw = zeros(size(Tds));    % exact radius curve
for k = 1:numel(Tds)
    rsw(k) = max(abs(eig(expm(Tds(k)*Lsw2)*expm(Tds(k)*Lsw1))));
end
% minimum dwell-time T* by BISECTION ON THE MODE-INDEXED SOS FEASIBILITY: each
% mode carries a timer certificate P_i(theta) with the all-pairs switch
% conditions (dwell type 'switched'); infeasible for short dwell, feasible above.
sfeasSw = @(Td) mvst_stab_sos_modes( ...
    struct('A0',cat(3,Asw1,Asw2)), struct('type','switched','Tbar',Td,'deg',4)).slack;
TdCrit = sos_bisect(sfeasSw, 0.1, 3.5, 'feashi');       % first feasible dwell
TdCritExact = crossing(Tds, rsw, 1);                    % exact radius crossing
fprintf('Stability (switched): minimum dwell time %.4f (SOS), %.4f (exact), rel. dev. %.2e\n', ...
    TdCrit, TdCritExact, abs(TdCrit-TdCritExact)/TdCritExact);
assert(abs(TdCrit-TdCritExact)/TdCritExact < 1e-2, ...
    'SOS minimum dwell-time disagrees with the exact monodromy crossing beyond 1e-2 relative');

% ---------- export data macros to LaTeX ----------
w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_stab.tex'));
w.mat('StabGenMjls',Pig);  w.num('StabAbsMjls',absMjls);
w.mat('StabDtMjlsAone',Ad1); w.mat('StabDtMjlsAtwo',Ad2); w.mat('StabDtMjlsP',Pdt);
w.num('StabDtMjlsChan',0.2,2); w.num('StabDtMjlsRhoMode',rhoDt(1));
w.num('StabDtMjlsRad',radDtMjls); w.num('StabDtMjlsRadIid',radDtMjlsIid);
w.num('StabDtMjlsRadNil',radDtMjlsNil,2);
w.num('StabAbsPoi',absPoi); w.num('StabLamCrit',lamCrit);
w.num('StabTbarNom',TbarNom); w.num('StabTbarCrit',TbarCrit);
w.mat('StabSwAone',Asw1); w.mat('StabSwAtwo',Asw2);
w.num('StabTdCrit',TdCrit);
w.close();

% ---------- figure 1: abscissa vs rate + monodromy vs period ----------
fig = figure('Color','w','Units','centimeters','Position',[2 2 14 12]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, lams, absc, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(ax1, 0, ':','LineWidth',1.2,'Color',[0.4 0.4 0.4]);
xline(ax1, lamCrit, '--','LineWidth',1.6,'Color',[0.5 0.5 0.5]);
xlabel(ax1,'jump rate $\lambda_1$','Interpreter','latex');
ylabel(ax1,'lifted-flow abscissa','Interpreter','latex');
lg1 = legend(ax1, {'abscissa','$0$','$\lambda_1^\star$'}, ...
             'Interpreter','latex','Location','northwest'); set(lg1,'Box','off');
xlim(ax1, [lams(1) lams(end)]);   % keep the plotted lambda window (lambda* is off-window)
set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, Tbars, rsd, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(ax2, 1, ':','LineWidth',1.2,'Color',[0.4 0.4 0.4]);
xline(ax2, TbarCrit, '--','LineWidth',1.6,'Color',[0.5 0.5 0.5]);
xline(ax2, TbarNom, ':','LineWidth',1.4,'Color',[0.5 0.5 0.5]);
xlabel(ax2,'sampling period $\bar T$','Interpreter','latex');
ylabel(ax2,'monodromy radius','Interpreter','latex');
lg2 = legend(ax2, {'$\rho(\bar T)$','$1$','$\bar T^\star$','nominal $\bar T$'}, ...
             'Interpreter','latex','Location','northwest'); set(lg2,'Box','off');
set(ax2,'FontSize',11);
mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_stab_rates.pdf'));

% ---------- figure 2: switched monodromy radius versus dwell ----------
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 13 9]);
ax = axes(fig2); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, Tds, rsw, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(ax, 1, ':','LineWidth',1.2,'Color',[0.4 0.4 0.4]);
xline(ax, TdCrit, '--','LineWidth',1.6,'Color',[0.5 0.5 0.5]);
xlabel(ax,'dwell-time $T$','Interpreter','latex');
ylabel(ax,'monodromy radius over one cycle','Interpreter','latex');
lg = legend(ax, {'$\rho(T)$','$1$','$T^\star$'},'Interpreter','latex', ...
            'Location','northeast'); set(lg,'Box','off');
set(ax,'FontSize',11);
mvst_savepdf(fig2, fullfile(mvst_figdir(),'ex_stab_dwell.pdf'));

% ========================================================================
function xc = sos_bisect(sfun, lo, hi, mode)
%SOS_BISECT  Bisection of a monotone SOS-feasibility boundary.  SFUN(x) returns
%  the optimal slack; feasible iff slack > tol.  MODE 'feaslo' brackets a system
%  feasible at LO and infeasible at HI (returns the last feasible x); 'feashi'
%  brackets infeasible at LO and feasible at HI (returns the first feasible x).
tol = 1e-7;
switch mode
    case 'feaslo'
        assert(sfun(lo) >  tol, 'sos_bisect: lo not feasible');
        assert(sfun(hi) <= tol, 'sos_bisect: hi not infeasible');
        for it = 1:30
            mid = 0.5*(lo+hi);
            if sfun(mid) > tol, lo = mid; else, hi = mid; end
        end
    case 'feashi'
        assert(sfun(lo) <= tol, 'sos_bisect: lo not infeasible');
        assert(sfun(hi) >  tol, 'sos_bisect: hi not feasible');
        for it = 1:30
            mid = 0.5*(lo+hi);
            if sfun(mid) > tol, hi = mid; else, lo = mid; end
        end
end
xc = 0.5*(lo+hi);
end

% ========================================================================
function rad = dtmjls_msrad(A, C, P)
%DTMJLS_MSRAD  Coupled DT-MJLS second-moment spectral radius on (S^n)^M.
%  Block (i,j) = P(j,i)*T_j, with T_j the REDUCED lift on the symmetric
%  subspace of X -> A_j X A_j' + C X C' (drift and multiplicative channel).
%  Mean-square stable iff RAD < 1 (mirrors DTMJLS_TRANSITION / EX_STABZ_DTMJLS).
n = size(A,1); M = size(A,3); Us = mvst_svecbasis(n); d = size(Us,2);
T = cell(1,M);
for j = 1:M
    T{j} = Us' * (kron(A(:,:,j),A(:,:,j)) + kron(C,C)) * Us;
end
G = zeros(M*d);
for i = 1:M
    for j = 1:M
        G((i-1)*d+(1:d), (j-1)*d+(1:d)) = P(j,i)*T{j};
    end
end
rad = max(abs(eig(G)));
end

% ========================================================================
function xc = crossing(x, y, level)
%CROSSING  Last downward crossing of Y through LEVEL (linear interpolation).
xc = NaN;
for k = numel(x)-1:-1:1
    if (y(k) - level)*(y(k+1) - level) <= 0
        xc = x(k) + (level - y(k))*(x(k+1) - x(k))/(y(k+1) - y(k));
        return
    end
end
end
