%% EX_DISS_MJLS  Stabilization by random switching, and its dissipativity.
%  Two-mode Markov-jump diffusion in which mode 1 is OPEN-LOOP UNSTABLE and
%  mode 2 is stable.  Neither the unstable mode in isolation nor a frozen
%  average is mean-square stable a priori, yet the Markov-chain switching
%  between the two modes renders the overall second-moment dynamics
%  mean-square stable (the classic "stabilization by random switching"
%  illustration).  Because mode 1 is not bounded-real on its own, its
%  worst-mode L2-gain is infinite, whereas the COUPLED mean-square gain of
%  the switched system is finite: switching also restores a finite
%  disturbance-to-output gain.
%
%  Two figures are produced:
%    diss_mjls_sample.pdf : Monte-Carlo sample paths of |x(t)| under the
%                           random mode signal (segments grow while the
%                           chain sits in the unstable mode, but the paths
%                           decay overall), with the active mode shaded;
%    diss_mjls_ms.pdf     : the mean-square E|x(t)|^2, computed EXACTLY by
%                           propagating the coupled second moments and
%                           cross-checked by the Monte-Carlo average, shown
%                           decaying to zero on a logarithmic scale, next to
%                           the diverging mean-square of the frozen unstable
%                           mode for contrast.
%  METHOD.  The coupled mean-square L2-gain is certified GENERICALLY on the
%  product cone (H^n)^2 through the SOSTOOLS routine MVST_DISS_SOS_MODES (dwell
%  type 'ctmjls', supply 'l2gain', mode 'gamma'): per-mode storages {P_i} with
%  the flow dissipation operator C_M(P_i) >= 0 augmented by the Markov-generator
%  coupling sum_j Pi(i,j) P_j.  The previous exact coupled bounded-real Riccati
%  value (MVST_L2GAIN_CT_MJLS) is printed as a cross-check.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
rng(7);                                        % reproducible sample paths

% ----- two-mode instance: mode 1 UNSTABLE, mode 2 stable -----------------
A1 = [ 1  2;  0  1];         % unstable: eigenvalues +1, +1
A2 = [-3  0;  0 -3];         % stable:   eigenvalues -3, -3
% canonical performance triple (Bw = [1;0], C = [1 0], Dw = 0.5) with the
% w-actuated diffusion: state-multiplicative channel 0.2 I2 per mode plus the
% disturbance-multiplicative column 0.1[1;0] so w acts in
% the diffusion as well as the drift.
Bw = [1;0];  C1 = [1 0];  C2 = [0.5 1];  Dw = 0.5;
C  = C1;
Cchan = {0.2*eye(2), 0.2*eye(2)};              % per-mode state-mult channel C_{i,1}
Dwc   = {0.1*[1;0]};                           % shared disturbance-mult channel D_{w,1}
Cscale = 0.2;                                  % scalar C_{i,1} = Cscale I2 (for the ms operator)
lam12 = 4; lam21 = 1;                          % 1->2 fast, 2->1 slow
Pig = [-lam12 lam12; lam21 -lam21];            % generator
n = 2; Mn = 2;

% ----- mode spectra and the mean-square spectral abscissa ----------------
ab1 = max(real(eig(A1)));  ab2 = max(real(eig(A2)));
Lms = mjls_ms_operator({A1,A2}, Pig, Cscale);  % coupled 2nd-moment operator (w-diffusion incl.)
absMS = max(real(eig(Lms)));                   % < 0  <=>  mean-square stable
% stationary distribution of the chain
ns = null(Pig.');  pistat = abs(ns(:,1)); pistat = pistat/sum(pistat);

% ----- coupled mean-square L2-gain (finite despite the unstable mode) ----
%  AUTHORITATIVE value: exact coupled bounded-real Riccati carrying BOTH the
%  state-multiplicative channel Cchan and the disturbance-multiplicative
%  channel Dwc (MVST_L2GAIN_CT_MJLS).  The product-cone SOS routine
%  MVST_DISS_SOS_MODES now folds BOTH channels (.Ac per mode and shared .Dwc),
%  so it is a full cross-check on the exact value.
% Long shooting horizon (1000) and tight relative tolerance (1e-6): the
% finite-horizon escape test under-estimates the gain (the
% default horizon 60 / tolerance 1e-3 gave 2.0791 against the LMI 2.0839).
gamC    = mvst_l2gain_ct_mjls({A1,A2}, Bw, {C1,C2}, Dw, Pig, 1000, Cchan, Dwc, 1e-6);  % exact (full)
Ac4 = zeros(n,n,1,Mn); Ac4(:,:,1,1) = Cscale*eye(n); Ac4(:,:,1,2) = Cscale*eye(n);
rC = mvst_diss_sos_modes(struct('A0',cat(3,A1,A2),'Bw',Bw,'C',cat(3,C1,C2), ...
        'Dw',Dw,'Pig',Pig,'Ac',Ac4,'Dwc',Dwc{1}), struct('type','l2gain'), struct('type','ctmjls'), 'gamma');
gamSOS = rC.gamma;                             % SOS (state + disturbance channels)
% mode-1 alone is NOT bounded-real (unstable): worst-mode gain is infinite.
% ISP index of the stable mode 2 as a stand-alone diffusion (output C2, state
% channel 0.2 I2, disturbance column Dwc), in the chapter convention
% He(w'z) - nu|w|^2 (the convention of MVST_PASSIVITY and MVST_PASSIVITY_CT;
% unlike them, this program also handles the diffusion channels).
nu = mvst_isp_margin_ct(A2, Bw, C2, Dw, {Cscale*eye(n)}, Dwc);

fprintf('MJLS: eig(A1) max = %.3f (unstable), eig(A2) max = %.3f (stable)\n',ab1,ab2);
fprintf('MJLS: mean-square spectral abscissa = %.4f (stable iff < 0)\n',absMS);
fprintf('MJLS: stationary dist = [%.3f %.3f]\n',pistat(1),pistat(2));
fprintf('MJLS: coupled mean-square L2-gain (exact, full channels) = %.4f\n',gamC);
fprintf('  [cross-check] exact=%.4f  SOS(full channels)=%.4f  |SOS-exact|=%.2e\n', ...
        gamC, gamSOS, abs(gamSOS-gamC));
fprintf('MJLS: ISP index of the stable mode (with diffusion) nu = %.4f\n', nu);

% ----- EXACT mean-square by propagating the coupled second moments -------
x0 = [1;1];  X0 = x0*x0';
Tend = 8; nt = 801; tg = linspace(0,Tend,nt);
Ems = mjls_ms_curve({A1,A2}, Pig, pistat, X0, tg, Cscale);   % E|x(t)|^2 (exact, w-diffusion incl.)
% frozen unstable mode 1 for contrast: mode-1 second-moment flow WITH the
% state-multiplicative diffusion, d/dt X = A1 X + X A1' + Cscale^2 X.
Efr = zeros(1,nt);
Lfr = kron(eye(n),A1) + kron(conj(A1),eye(n)) + Cscale^2*eye(n*n);
y0fr = X0(:);
for k = 1:nt
    yk = expm(Lfr*tg(k))*y0fr;
    Efr(k) = real(trace(reshape(yk,n,n)));
end

% ----- Monte-Carlo sample paths of the switched diffusion ----------------
npath = 200;  showp = 12;
[Xpaths, modeSig] = mjls_montecarlo({A1,A2}, Pig, pistat, x0, tg, npath, Cscale);
normPaths = squeeze(sqrt(sum(Xpaths.^2,1)));    % nt x npath  = |x(t)|
Emc = mean(sum(Xpaths.^2,1),3);                 % Monte-Carlo E|x|^2 (1 x nt)

% ----- export data macros to LaTeX ---------------------------------------
w = mvst_texdata(fullfile(mvst_figdir(),'data','diss_mjls.tex'));
w.mat('DAmjlsOne',A1); w.mat('DAmjlsTwo',A2); w.mat('DBwmjls',Bw);
w.mat('DCmjls',C1);    w.mat('DCmjlsTwo',C2);  w.num('DDwmjls',Dw);    w.mat('DGenMjls',Pig);
w.num('DGamMjlsC',gamC);  w.num('DNuMjls',nu);
w.num('DAbsMS',absMS);    w.num('DPiStatTwo',pistat(2));
w.num('DAbOne',ab1);      w.num('DAbTwo',ab2);
w.num('DGamMjlsSOS',gamSOS); w.raw('DGamMjlsGap',scifmt(abs(gamSOS-gamC)));
w.close();

% ----- figure 1: sample paths with the active mode shaded ----------------
mjls_sample_fig(tg, normPaths(:,1:showp), modeSig(:,1), ...
    fullfile(mvst_figdir(),'diss_mjls_sample.pdf'), gc);

% ----- figure 2: mean-square, exact + Monte-Carlo, vs frozen unstable ----
mjls_ms_fig(tg, Ems, Emc, Efr, ...
    fullfile(mvst_figdir(),'diss_mjls_ms.pdf'), gc);

% ----- combined 2-up: sample paths (left) + mean square (right) ----------
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tlc = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');
mjls_sample_ax(nexttile(tlc), tg, normPaths(:,1:showp), modeSig(:,1), gc);
mjls_ms_ax(nexttile(tlc), tg, Ems, Efr, gc);
mvst_savepdf(figc, fullfile(mvst_figdir(),'diss_mjls_combined.pdf'));

% ========================================================================
function s = scifmt(x)
%SCIFMT  LaTeX scientific notation m\cdot10^{e} with one decimal.
if x <= 0; s = '0'; return; end
e = floor(log10(x)); mnt = x/10^e;
if round(mnt*10)/10 >= 10; mnt = mnt/10; e = e+1; end
s = sprintf('%.1f\\cdot10^{%d}', mnt, e);
end

% ========================================================================
function L = mjls_ms_operator(Acell, Pig, Cscale)
%MJLS_MS_OPERATOR  Coupled second-moment operator on stacked vec(Xi).
%   The optional CSCALE folds the per-mode state-multiplicative diffusion
%   channel C_{i,1} = CSCALE I into each mode block (adds CSCALE^2 I).
if nargin < 3 || isempty(Cscale); Cscale = 0; end
Mn = numel(Acell); n = size(Acell{1},1); nn = n*n;
L = zeros(Mn*nn);
for i = 1:Mn
    Ai = Acell{i};
    blk = kron(eye(n),Ai) + kron(conj(Ai),eye(n)) + Cscale^2*eye(nn); % + diffusion
    L((i-1)*nn+(1:nn),(i-1)*nn+(1:nn)) = blk;
end
for i = 1:Mn
    for j = 1:Mn
        idx_i = (i-1)*nn+(1:nn); idx_j = (j-1)*nn+(1:nn);
        L(idx_i,idx_j) = L(idx_i,idx_j) + Pig(j,i)*eye(nn);
    end
end
end

% ========================================================================
function E = mjls_ms_curve(Acell, Pig, pi0, X0, tg, Cscale)
%MJLS_MS_CURVE  Exact E|x(t)|^2 by integrating the coupled moments.
if nargin < 6 || isempty(Cscale); Cscale = 0; end
Mn = numel(Acell); n = size(Acell{1},1); nn = n*n;
L = mjls_ms_operator(Acell, Pig, Cscale);
y0 = zeros(Mn*nn,1);
for i = 1:Mn
    Xi0 = pi0(i)*X0; y0((i-1)*nn+(1:nn)) = Xi0(:);
end
E = zeros(1,numel(tg));
for k = 1:numel(tg)
    yk = expm(L*tg(k))*y0;
    s = 0;
    for i = 1:Mn
        Xi = reshape(yk((i-1)*nn+(1:nn)),n,n); s = s + trace(Xi);
    end
    E(k) = real(s);
end
end

% ========================================================================
function [X, modeSig] = mjls_montecarlo(Acell, Pig, pi0, x0, tg, npath, Cscale)
%MJLS_MONTECARLO  Sample paths of dx = A_sigma x dt + Cscale x dW, CTMC sigma.
%   Euler-Maruyama on the fine grid TG with the state-multiplicative diffusion
%   channel C_{sigma,1} = CSCALE I (scalar Brownian dW), and the CTMC mode
%   switched with per-step probability rate*dt.
if nargin < 7 || isempty(Cscale); Cscale = 0; end
Mn = numel(Acell); n = numel(x0); nt = numel(tg);
X = zeros(n,nt,npath); modeSig = zeros(nt,npath);
rates = -diag(Pig);                                   % exit rates
for p = 1:npath
    m = randsample_local(pi0);                        % initial mode ~ pi0
    x = x0;
    for k = 1:nt
        X(:,k,p) = x; modeSig(k,p) = m;
        if k == nt; break; end
        dt = tg(k+1)-tg(k); sq = sqrt(dt);
        x = x + Acell{m}*x*dt + Cscale*x*sq*randn;    % Euler-Maruyama flow + diffusion
        if rates(m) > 0 && rand < rates(m)*dt         % mode switch this step
            pr = Pig(m,:); pr(m) = 0; pr = pr/sum(pr);
            m = randsample_local(pr);
        end
    end
end
end

% ========================================================================
function i = randsample_local(p)
c = cumsum(p(:)); r = rand; i = find(r <= c, 1, 'first');
if isempty(i); i = numel(p); end
end

% ========================================================================
function mjls_sample_fig(tg, paths, mode1, fname, gc)
%MJLS_SAMPLE_FIG  |x(t)| sample paths with the unstable mode shaded.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
mjls_sample_ax(axes(fig), tg, paths, mode1, gc);
mvst_savepdf(fig, fname);
end

% ========================================================================
function mjls_sample_ax(ax, tg, paths, mode1, gc)
%MJLS_SAMPLE_AX  Draw the sample-path panel into an existing axes.
hold(ax,'on'); box(ax,'on');
% shade intervals where the (first shown) path sits in the unstable mode 1
yl = [0, 1.05*max(paths(:))];
inmode1 = (mode1(:).' == 1);
d = diff([0 inmode1 0]);
starts = find(d==1); stops = find(d==-1)-1;
for s = 1:numel(starts)
    x1 = tg(starts(s)); x2 = tg(stops(s));
    patch(ax,[x1 x2 x2 x1],[yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
        'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
end
hp = plot(ax, tg, paths, '-', 'LineWidth',0.7, 'Color',[gc.order(1,:) 0.55]);
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$|x(t)|$','Interpreter','latex');
ylim(ax, yl); set(ax,'FontSize',11);
lg = legend(ax, [hp(1)], {'sample paths ($|x(t)|$)'}, ...
    'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off');
text(ax, 0.98, 0.72, 'shaded: unstable mode active', 'Units','normalized', ...
    'HorizontalAlignment','right','Interpreter','latex','FontSize',9, ...
    'Color',[0.5 0.5 0.5]);
end

% ========================================================================
function mjls_ms_fig(tg, Ems, ~, Efr, fname, gc)
%MJLS_MS_FIG  Exact mean-square E|x|^2 of the switched diffusion against the
%  frozen unstable mode.  Linear axes; the Monte-Carlo average is not shown
%  (it merely reproduces the exact mean square).
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
mjls_ms_ax(axes(fig), tg, Ems, Efr, gc);
mvst_savepdf(fig, fname);
end

% ========================================================================
function mjls_ms_ax(ax, tg, Ems, Efr, gc)
%MJLS_MS_AX  Draw the mean-square panel into an existing axes.
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h1 = plot(ax, tg, Efr, '--', 'LineWidth',1.8, 'Color',gc.order(1,:));
h2 = plot(ax, tg, Ems, '-',  'LineWidth',2.0, 'Color',gc.order(2,:));
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'mean square $\,|x(t)|^2$','Interpreter','latex');
lg = legend(ax, [h1 h2], ...
    {'frozen unstable mode 1', 'switched (exact)'}, ...
    'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off');
set(ax,'FontSize',11);
end
