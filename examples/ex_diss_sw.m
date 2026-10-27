%% EX_DISS_SW  Dissipativity of a switched two-mode system (dwell-time).
%  Two Hurwitz modes sharing the disturbance/output channels, with a minimum
%  dwell-time Td.  Computes the worst-mode L2-gain (a LOWER bound for the
%  switched gain: staying in one mode is dwell-time admissible, and
%  switching generally raises the gain; the dwell-time storage conditions
%  certify upper bounds) and the mode-1 passivity margin, and illustrates the
%  convergence of the finite-horizon available-storage matrix.  Two
%  single-panel figures are produced: the storage-matrix convergence and a
%  sample-path/mean-square plot of the disturbance-driven trajectory under
%  the alternating dwell-time schedule (replacing the former Bode certificate).
%  METHOD.  The per-mode gains flow through the generic SOSTOOLS routine
%  MVST_DISS_SOS (dwell type 'arbitrary'), and the minimum-dwell-time certified
%  UPPER bound flows through the product-cone routine MVST_DISS_SOS_MODES
%  (dwell type 'switched'): mode-indexed timer storages P_i(theta) with per-mode
%  flow C_M >= 0 on [0,Td], the frozen-tail flow C_M^0 >= 0, and the all-pairs
%  switch conditions J^*P_j(0)J - P_i(Td) >= 0 (Theorem th:diss:MDT).  The
%  previous exact/reference values are printed as cross-checks.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- two-mode SISO instance (scalar disturbance w, scalar output z) -----
% A positive feedthrough Dw makes the system input-strictly passive
% (nu = Dw + min_w Re G(jw) > 0) while keeping the L2-gain finite.
A1 = [0 1;-2 -1];   A2 = [0 1;-1 -2];
Bw = [1;0];         C = [1 0];   Dw = 0.5;   Td = 0.5;
% canonical performance triple with the disturbance-actuated diffusion channel
% D_{w,1} = 0.1[1;0] per mode (no state-multiplicative part).
Dwcm = 0.1*[1;0];  Dwc = {Dwcm};
aOU  = 5;          % rate of the Ornstein-Uhlenbeck disturbance of the trajectory figure

% ----- worst-mode L2-gain (lower bound for the switched gain) -----
%  Exact per-mode mean-square gains (with Dwc) are authoritative; SOS per mode
%  (which folds Dwc into its flow block) is the cross-check.
sup = struct('type','l2gain'); dwa = struct('type','arbitrary');
%  Long-horizon (1000), tight-tolerance (1e-6) shooting, with the diffusion column DWC.
g1  = mvst_l2gain_ct_poisson(A1,Bw,C,Dw,{zeros(2)},Dwc,{},[],{},1000,1e-6);   % exact
g2  = mvst_l2gain_ct_poisson(A2,Bw,C,Dw,{zeros(2)},Dwc,{},[],{},1000,1e-6);   % exact
gam = max(g1,g2);
g1SOS = mvst_diss_sos(struct('A0',A1,'Bw',Bw,'C',C,'Dw',Dw,'Ac',zeros(2,2,1),'Dwc',Dwcm), sup, dwa, 'gamma').gamma;
g2SOS = mvst_diss_sos(struct('A0',A2,'Bw',Bw,'C',C,'Dw',Dw,'Ac',zeros(2,2,1),'Dwc',Dwcm), sup, dwa, 'gamma').gamma;

% ----- dwell-time (MDT) certified UPPER bound on the switched gain -----
%  Exact MDT Riccati-shooting bound (with Dwc) is authoritative; the product-
%  cone SOS routine MVST_DISS_SOS_MODES now folds the shared .Dwc channel into
%  its timer flow block, so it is a full cross-check (residual = MDT/SOS
%  timer-degree conservatism).  The true switched gain lies in [gam, gamDT].
%  MVST_L2GAIN_CT_DWELL is NOT exact: it searches the
%  terminal storages of the MDT certificate by multi-start Nelder-Mead, so it
%  returns the smallest level at which its SEARCH found a certificate (an upper
%  bound); the SOS program solves the same MDT conditions (polynomial timer
%  storages of degree 4) as an SDP and returns a smaller certified bound.
%  Both are upper bounds; the SOS one is reported.
gamDTRic = mvst_l2gain_ct_dwell({A1,A2}, Bw, {C,C}, Dw, Td, {{},{}}, Dwc);  % Riccati-search bound
sysm  = struct('A0',cat(3,A1,A2),'Bw',Bw,'C',cat(3,C,C),'Dw',Dw,'Dwc',Dwcm);
gamDTSOS = mvst_diss_sos_modes(sysm, sup, ...
                            struct('type','switched','Tbar',Td,'deg',4), 'gamma').gamma;
gamDT = gamDTSOS;                                % reported certified upper bound

% ----- mode-1 passivity margin (drift) -----
% ISP index in the chapter convention, with the diffusion column
nu = mvst_isp_margin_ct(A1, Bw, C, Dw, {}, Dwc);  isp = nu > 0;

% ----- convergence of the finite-horizon storage matrix (mode 1) -----
[Tg,S0] = diss_pT_dwc(A1,Bw,C,Dw,Dwcm, gam*1.05, 6, 40);   % with Dwc; S = -P_T(0)

fprintf('Switched dissipativity: mode-1 gain %.4f, mode-2 gain %.4f, worst %.4f\n', g1, g2, gam);
fprintf('Switched dissipativity: dwell-time certified upper bound (SOS) %.4f, Riccati-search bound %.4f (interval [%.4f, %.4f])\n', gamDT, gamDTRic, gam, gamDT);
fprintf('Switched passivity (mode 1): isp=%d, nu=%.4f, dwell Td=%.2f\n', isp, nu, Td);
fprintf(['  [cross-check] modes: exact g1=%.4f (SOS %.4f), g2=%.4f (SOS %.4f); ', ...
         'MDT bound: exact=%.4f, SOS(full channels)=%.4f, |d|=%.2e\n'], ...
        g1, g1SOS, g2, g2SOS, gamDTRic, gamDTSOS, abs(gamDTRic-gamDTSOS));

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','diss_sw.tex'));
w.mat('DAswOne',A1); w.mat('DAswTwo',A2); w.mat('DBwsw',Bw); w.mat('DCsw',C);
w.num('DDwsw',Dw); w.num('DTdwellsw',Td);
w.num('DGamSw',gam); w.num('DNuSw',nu);
w.num('DGamSwDT',gamDT);
w.num('DGamSwDTRic',gamDTRic);
w.num('DGamSwOneSOS',g1SOS); w.num('DGamSwTwoSOS',g2SOS);
w.raw('DGamSwGap',scifmt(max(abs(g1SOS-g1),abs(g2SOS-g2))));
w.num('DOUrate',aOU,0);
w.close();

% ----- figure 1: convergence of the available storage -----
diss_conv_fig(Tg, S0, 'Switched: convergence of the available storage', ...
    fullfile(mvst_figdir(),'diss_sw_conv.pdf'), gc);

% ----- figure 2: sample paths + mean square under the switching schedule -----
%  Alternate modes every dwell-time Td (a periodic admissible schedule).
rng(13);
sched = @(t) 1 + mod(floor(t/Td), 2);
diss_traj_fig({A1,A2}, sched, Bw, Dwcm, aOU, 12, ...
    fullfile(mvst_figdir(),'diss_sw_traj.pdf'), gc);

% ----- combined 2-up: storage convergence (left) + trajectory (right) -----
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tlc = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');
diss_conv_ax(nexttile(tlc), Tg, S0, '(a) available-storage convergence', gc);
rng(13);
diss_traj_ax(nexttile(tlc), {A1,A2}, sched, Bw, Dwcm, aOU, 12, '(b) disturbed trajectory', gc);
mvst_savepdf(figc, fullfile(mvst_figdir(),'diss_sw_combined.pdf'));

% ========================================================================
function s = scifmt(x)
%SCIFMT  LaTeX scientific notation m\cdot10^{e} with one decimal.
if x <= 0; s = '0'; return; end
e = floor(log10(x)); mnt = x/10^e;
if round(mnt*10)/10 >= 10; mnt = mnt/10; e = e+1; end
s = sprintf('%.1f\\cdot10^{%d}', mnt, e);
end

% ========================================================================
function [Tgrid, S0] = diss_pT_dwc(A, Bw, C, Dw, D1, gam, Tmax, nT)
%DISS_PT_DWC  Finite-horizon available storage S_T(0) = -P_T(0) versus T for
%   dx = (A x + Bw w)dt + D1 w dW, z = C x + Dw w: the bounded-real Riccati
%     -dS/dt = A'S + S A + C'C + v Theta^{-1} v',
%     v = S Bw + C'Dw,  Theta = gam^2 I - Dw'Dw - D1'S D1,
%   integrated backward from S(T) = 0.
n = size(A,1); m = size(Bw,2);
Tgrid = linspace(0.1, Tmax, nT); S0 = zeros(n,n,nT);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
rhs = @(s,q) bre(q);
for k = 1:nT
    [~,Y] = ode45(rhs, [0 Tgrid(k)], zeros(n*n,1), opts);
    S = reshape(Y(end,:).',n,n); S0(:,:,k) = (S+S.')/2;
end
    function dq = bre(q)
        S = reshape(q,n,n); S = (S+S.')/2;
        v = S*Bw + C.'*Dw; Th = gam^2*eye(m) - Dw.'*Dw - D1.'*S*D1;
        G = A.'*S + S*A + C.'*C + v*(Th\v.'); G = (G+G.')/2; dq = G(:);
    end
end

% ========================================================================
function diss_conv_fig(Tg, S0, ttl, fname, gc) %#ok<INUSL>
%DISS_CONV_FIG  Standalone wrapper around DISS_CONV_AX.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
diss_conv_ax(axes(fig), Tg, S0, '', gc);
mvst_savepdf(fig, fname);
end

% ========================================================================
function diss_conv_ax(ax, Tg, S, ttl, gc)
%DISS_CONV_AX  Entries of the available storage S_T(0) = -P_T(0) versus T.
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h(1) = plot(ax, Tg, reshape(S(1,1,:),1,[]), '-o', 'LineWidth',1.8, ...
            'MarkerSize',3, 'Color',gc.order(1,:));
h(2) = plot(ax, Tg, reshape(S(1,2,:),1,[]), '-s', 'LineWidth',1.8, ...
            'MarkerSize',3, 'Color',gc.order(2,:));
h(3) = plot(ax, Tg, reshape(S(2,2,:),1,[]), '-^', 'LineWidth',1.8, ...
            'MarkerSize',3, 'Color',gc.order(3,:));
xlabel(ax,'$T$','Interpreter','latex');
ylabel(ax,'available storage matrix $-P_T(0)$ entries','Interpreter','latex');
lg = legend(ax, h, {'$-P_{11}$','$-P_{12}$','$-P_{22}$'}, ...
            'Interpreter','latex','Orientation','horizontal', ...
            'Location','southoutside');
set(lg,'Box','off'); set(ax,'FontSize',11);
if ~isempty(ttl); title(ax, ttl, 'Interpreter','latex','FontWeight','normal'); end
end

% ========================================================================
function diss_traj_fig(Acell, sched, Bw, D1, aOU, Tend, fname, gc)
%DISS_TRAJ_FIG  Standalone wrapper around DISS_TRAJ_AX.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
diss_traj_ax(axes(fig), Acell, sched, Bw, D1, aOU, Tend, '', gc);
mvst_savepdf(fig, fname);
end

% ========================================================================
function diss_traj_ax(ax, Acell, sched, Bw, D1, aOU, Tend, ttl, gc)
%DISS_TRAJ_AX  Sample paths and mean square of the stated diffusion
%     dx = (A_sched(t) x + Bw w)dt + D1 w dW,
%   driven by the Ornstein-Uhlenbeck disturbance dw = -aOU w dt + sqrt(2 aOU) dV
%   (V independent of W; unit stationary variance), from rest x(0)=0, w(0)=0.
%   Euler-Maruyama sample paths of |x(t)|^2 (thin), the EXACT mean square from
%   the second moment Sigma of (x,w),
%     dSigma/dt = Ab Sigma + Sigma Ab' + Cb Sigma Cb' + bb bb',
%     Ab = [A Bw; 0 -aOU],  Cb = [0 D1; 0 0],  bb = [0; sqrt(2 aOU)]  (thick),
%   integrated exactly mode by mode (matrix exponential of the vectorized
%   operator, constant mode on each grid step), and the Monte-Carlo average of
%   the sample paths (dotted).  SCHED(t) returns the active mode index at time t.
dt = 2e-3; tg = 0:dt:Tend; nt = numel(tg);
npath = 300; showp = 14; n = size(Acell{1},1); N = n+1;
sq = zeros(nt,npath); sOU = sqrt(2*aOU);
for p = 1:npath
    x = zeros(n,1); w = 0;
    for k = 1:nt
        sq(k,p) = x.'*x;
        if k==nt; break; end
        A = Acell{sched(tg(k))};
        dW = sqrt(dt)*randn; dV = sqrt(dt)*randn;
        x = x + (A*x + Bw*w)*dt + D1*w*dW;
        w = w - aOU*w*dt + sOU*dV;
    end
end
% exact second moment of (x,w): one-step propagators per mode
nm = numel(Acell); Phi = cell(1,nm); bvec = [zeros(n,1); sOU];
for i = 1:nm
    Ab = [Acell{i}, Bw; zeros(1,n), -aOU];
    Cb = [zeros(n), D1; zeros(1,N)];
    L  = kron(eye(N),Ab) + kron(Ab,eye(N)) + kron(Cb,Cb);
    f  = reshape(bvec*bvec.', [], 1);
    Ex = expm([L, f; zeros(1,N*N+1)]*dt);
    Phi{i} = Ex;
end
v = zeros(N*N,1); Ems = zeros(1,nt);
for k = 1:nt
    Sig = reshape(v,N,N); Ems(k) = real(trace(Sig(1:n,1:n)));
    if k==nt; break; end
    y = Phi{sched(tg(k))}*[v; 1]; v = y(1:N*N);
end
Emc = mean(sq,2).';
dd = max(1, round(nt/1000)); id = 1:dd:nt;   % decimate for a light vector PDF
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
hp = plot(ax, tg(id), sq(id,1:showp), '-', 'LineWidth',0.6, 'Color',[gc.order(1,:) 0.4]);
h2 = plot(ax, tg(id), Ems(id), '-', 'LineWidth',2.2, 'Color',gc.order(2,:));
h3 = plot(ax, tg(id), Emc(id), ':', 'LineWidth',1.8, 'Color',gc.order(3,:));
xlabel(ax,'$t$','Interpreter','latex');
ylabel(ax,'$|x(t)|^2$ and its mean','Interpreter','latex');
lg = legend(ax, [hp(1) h2 h3], ...
    {'sample paths $|x(t)|^2$','mean square (exact)','Monte-Carlo average'}, ...
    'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off');
set(ax,'FontSize',11); ylim(ax,[0, 1.05*max(Ems)*4]);
if ~isempty(ttl); title(ax, ttl, 'Interpreter','latex','FontWeight','normal'); end
end
