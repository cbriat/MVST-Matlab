%% EX_HIGH_COMPOUND  Second additive compound and 2-positivity.
%  Chapter "Higher-order systems", numerical illustration of
%  Theorem th:high:kpos and Example ex:high:compound.  For the 3x3 instance
%      A = [-3 1 -1; 2 -4 2; -1 1 -2]
%  the script forms the second additive compound A^{[2]} from the explicit
%  formula eq:high:A2, checks the Metzler sign conditions that decide
%  2-positivity (Theorem th:high:kpos), and confirms:
%    * A itself is NOT Metzler (its corner entries a13=a31=-1<0), so the
%      system is NOT 1-positive (does not preserve the nonnegative orthant);
%    * A^{[2]} IS Metzler, so the system IS 2-positive: it preserves the
%      sign-variation cone P^2_- = {s^-(x) <= 1}.
%  Two numerical cross-checks:
%    * the exponential identity eq:high:compexp, (e^{At})^{(2)} = e^{A^{[2]} t},
%      to machine precision, and the entrywise nonnegativity of the
%      multiplicative 2-compound (e^{At})^{(2)} that IS 2-positivity;
%    * a sign-variation-diminishing simulation: trajectories started in
%      P^2_- keep s^-(x(t)) <= 1 for all t.
%  Figure figures/high_compound.pdf, data figures/data/high_compound.tex.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

A   = [-3 1 -1; 2 -4 2; -1 1 -2];
A2c = addcomp2(A);                     % eq:high:A2

% ---- Metzler sign conditions (Theorem th:high:kpos) ----------------------
isMetz1 = metzler(A);                  % 1-positive?  (false here)
isMetz2 = metzler(A2c);                % 2-positive?  (true  here)

% ---- exponential identity eq:high:compexp: (e^{At})^{(2)} = e^{A^{[2]} t} -
tg   = linspace(0.01, 3, 120);
expR = 0; phiMin = inf; g1 = zeros(size(tg)); g2 = zeros(size(tg));
for k = 1:numel(tg)
    Phi  = expm(A*tg(k));
    Phi2 = multcomp2(Phi);
    expR   = max(expR, norm(Phi2 - expm(A2c*tg(k)), 'fro'));
    phiMin = min(phiMin, min(Phi2(:)));         % >= 0  <=> 2-positive
    g1(k)  = min(Phi(:));                         % dips < 0  (not 1-positive)
    g2(k)  = min(Phi2(:));                        % stays >= 0 (2-positive)
end

% ---- sign-variation-diminishing simulation -------------------------------
rng(0); svMax = 0; ntraj = 400;
tsim = linspace(0, 4, 60);
kept = 0;
for r = 1:ntraj
    x0 = randn(3,1);
    if svar(x0) > 1, continue; end                 % start in P^2_- (s^- <= 1)
    kept = kept + 1;
    for k = 1:numel(tsim)
        svMax = max(svMax, svar(expm(A*tsim(k))*x0));
    end
end

% ---- soundness verification ----------------------------------------------
assert(~isMetz1, 'A unexpectedly Metzler (should not be 1-positive).');
assert(isMetz2,  'A^{[2]} not Metzler (should be 2-positive).');
assert(expR < 1e-9, 'exponential compound identity failed.');
assert(phiMin > -1e-9, 'multiplicative 2-compound went negative.');
assert(svMax <= 1, 'sign variation exceeded 1: P^2_- not invariant.');
assert(any(g1 < 0), 'e^{At} stayed nonnegative: instance is not a clean 1-positive counterexample.');

fprintf('EX_HIGH_COMPOUND: Metzler(A)=%d, Metzler(A^[2])=%d\n', isMetz1, isMetz2);
fprintf('  exp-identity resid=%.2e, min entry (e^{At})^{(2)}=%.2e, max s^- =%d (%d trajectories)\n', ...
    expR, phiMin, svMax, kept);

% ---- figure: 1-positivity fails, 2-positivity holds ----------------------
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 8.5]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h1 = plot(ax, tg, g1, '-',  'LineWidth',1.8, 'Color', gc.order(2,:));
h2 = plot(ax, tg, g2, '-',  'LineWidth',1.8, 'Color', gc.order(1,:));
yline(ax, 0, ':', 'LineWidth',1.2, 'Color',[0.4 0.4 0.4]);
xlabel(ax, 'time $t$', 'Interpreter','latex');
ylabel(ax, 'smallest entry', 'Interpreter','latex');
title(ax, 'Compound positivity: $e^{At}$ vs $\left(e^{At}\right)^{(2)}=e^{A^{[2]}t}$', ...
    'Interpreter','latex');
lg = legend(ax, [h1 h2], ...
    {'$\min_{ij}(e^{At})_{ij}$ (dips $<0$: not $1$-positive)', ...
     '$\min_{ij}(e^{At})^{(2)}_{ij}$ ($\ge0$: $2$-positive)'}, ...
    'Interpreter','latex','Location','southeast');
set(lg,'Box','off');
mvst_savepdf(fig, fullfile(mvst_figdir(),'high_compound.pdf'));

% ---- data macros ---------------------------------------------------------
w = mvst_texdata(fullfile(mvst_figdir(),'data','high_compound.tex'));
w.mat('HighCompA',    A,   0);
w.mat('HighCompAtwo', A2c, 0);
w.raw('HighCompExpResid', sprintf('%.0e', expR));
w.raw('HighCompPhiMin',   sprintf('%.1e', phiMin));
w.num('HighCompSvMax',    svMax, 0);
w.close();

% ======================= local functions ==============================
function B = addcomp2(A)
%ADDCOMP2 Second additive compound of a 3x3 matrix (subsets {1,2},{1,3},{2,3}).
B = [A(1,1)+A(2,2), A(2,3),        -A(1,3);
     A(3,2),        A(1,1)+A(3,3),  A(1,2);
    -A(3,1),        A(2,1),         A(2,2)+A(3,3)];
end
function B = multcomp2(A)
%MULTCOMP2 Second multiplicative compound of a 3x3 matrix (all 2x2 minors).
S = [1 2; 1 3; 2 3]; B = zeros(3);
for r=1:3, for c=1:3, B(r,c)=det(A(S(r,:),S(c,:))); end, end
end
function tf = metzler(M)
%METZLER True iff all off-diagonal entries of M are nonnegative.
off = M - diag(diag(M)); tf = all(off(~eye(size(M))) >= -1e-12);
end
function s = svar(z)
%SVAR Number of sign changes in Z after deleting zero entries.
z = z(abs(z) > 1e-12); s = sum(diff(sign(z)) ~= 0);
end
