%% EX_DISS_SD  Dissipativity of a sampled-data (continuous-flow) system.
%  Continuous flow with a scalar disturbance and a dwell/sample period Tbar.
%  Computes the continuous-time L2-gain and the passivity margin, and
%  illustrates the convergence of the finite-horizon available-storage matrix
%  to its infinite-horizon limit.  Two single-panel figures are produced:
%  the storage-matrix convergence and a sample-path/mean-square plot of the
%  disturbance-driven trajectory (replacing the former Bode certificate).
%  METHOD.  Both certified gains flow through the generic SOSTOOLS routine
%  MVST_DISS_SOS: the continuous flow L2-gain via dwell type 'arbitrary'
%  (flow C_M >= 0, the bounded-real LMI eq:diss:bri), and the hybrid
%  sampled-OUTPUT gain via dwell type 'periodic' (flow with pure disturbance
%  credit plus the sampling jump D_M >= 0 carrying the sampled output CJ = C,
%  the periodic dual DLMI th:diss:periodic:LMI).  The previous exact values
%  (MVST_L2GAIN_CT and the periodic I-D-GRE MVST_L2GAIN_SD_HYBRID) are printed
%  as cross-checks.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- SISO instance (scalar disturbance w, scalar output z) -----
% A positive feedthrough Dw makes the system input-strictly passive
% (nu = Dw + min_w Re G(jw) > 0) while keeping the L2-gain finite.
A = [0 1;-1 -1];   Bw = [1;0];   C = [1 0];   Dw = 0.5;   Tbar = 0.3;
% canonical performance triple with the disturbance-actuated diffusion channel
% D_{w,1} = 0.1[1;0] (no state-multiplicative part), so w
% acts in the diffusion as well as the drift.
Dwcm = 0.1*[1;0];  Dwc = {Dwcm};
aOU  = 5;          % rate of the Ornstein-Uhlenbeck disturbance of the trajectory figure

% ----- flow L2-gain: exact mean-square shooting-Riccati (with Dwc) is the
%       authoritative value; the SOS program (which folds Dwc into its flow
%       block) is the cross-check.  Passivity margin from the drift.
%  The exact value is the stochastic bounded-real shooting test with a long
%  horizon (1000) and relative tolerance 1e-6 (horizon 60 under-estimates the gain).
gam      = mvst_l2gain_ct_poisson(A,Bw,C,Dw,{zeros(2)},Dwc,{},[],{},1000,1e-6);  % exact (with Dwc)
gamSOS   = mvst_diss_sos(struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw,'Ac',zeros(2,2,1),'Dwc',Dwcm), ...
                         struct('type','l2gain'), struct('type','arbitrary'), 'gamma').gamma;
% ISP index in the chapter convention He(w'z) - nu|w|^2, with the diffusion
% column Dwc (same convention as MVST_PASSIVITY_CT, which handles the
% drift-only case).
nu = mvst_isp_margin_ct(A, Bw, C, Dw, {}, Dwc);  isp = nu > 0;

% ----- hybrid sampled-output gain z_k = C x(k*Tbar): exact periodic I-D-GRE
%       (with Dwc) is authoritative; periodic deg-4 SOS cross-check.
gamH     = mvst_l2gain_sd_hybrid(A, Bw, C, Tbar, ceil(1000/Tbar), {}, Dwc, 1e-6);   % exact (with Dwc), horizon 1000
gamHSOS  = mvst_diss_sos(struct('A0',A,'Bw',Bw,'J',eye(size(A,1)),'CJ',C,'Ac',zeros(2,2,1),'Dwc',Dwcm), ...
                         struct('type','l2gain'), ...
                         struct('type','periodic','Tbar',Tbar,'deg',4), 'gamma').gamma;

% ----- convergence of the finite-horizon storage matrix -----
%  (with the diffusion column Dwc in the Riccati, and the storage
%   S = -P_T(0) >= 0 returned directly; the horizon is 20, since at T = 6 the
%   storage entries are still rising and they settle near T = 20.)
[Tg,S0] = diss_pT_dwc(A,Bw,C,Dw,Dwcm, gam*1.05, 20, 60);

fprintf('Sampled-data dissipativity: flow L2-gain %.4f, Tbar %.2f\n', gam, Tbar);
fprintf('Sampled-data dissipativity: hybrid (sampled-output) gain %.4f\n', gamH);
fprintf('Sampled-data passivity: isp=%d, nu=%.4f\n', isp, nu);
fprintf('  [cross-check] flow: exact=%.4f SOS=%.4f |d|=%.2e ; hybrid: exact=%.4f SOS=%.4f |d|=%.2e\n', ...
        gam, gamSOS, abs(gam-gamSOS), gamH, gamHSOS, abs(gamH-gamHSOS));

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','diss_sd.tex'));
w.mat('DAsd',A); w.mat('DBwsd',Bw); w.mat('DCsd',C);
w.num('DDwsd',Dw); w.num('DTbarsd',Tbar);
w.num('DGamSd',gam); w.num('DNuSd',nu);
w.num('DGamSdH',gamH);
w.num('DGamSdSOS',gamSOS); w.raw('DGamSdGap',scifmt(abs(gamSOS-gam)));
w.num('DGamSdHSOS',gamHSOS); w.num('DOUrate',aOU,0);
w.close();

% ----- figure 1: convergence of the available storage -----
diss_conv_fig(Tg, S0, 'Sampled-data: convergence of the available storage', ...
    fullfile(mvst_figdir(),'diss_sd_conv.pdf'), gc);

% ----- figure 2: sample paths + mean square of the disturbed trajectory -----
rng(11);
diss_traj_fig({A}, @(t) 1, Bw, Dwcm, aOU, 12, ...
    fullfile(mvst_figdir(),'diss_sd_traj.pdf'), gc);

% ----- combined 2-up: storage convergence (left) + trajectory (right) -----
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tlc = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');
diss_conv_ax(nexttile(tlc), Tg, S0, '(a) available-storage convergence', gc);
rng(11);
diss_traj_ax(nexttile(tlc), {A}, @(t) 1, Bw, Dwcm, aOU, 12, '(b) disturbed trajectory', gc);
mvst_savepdf(figc, fullfile(mvst_figdir(),'diss_sd_combined.pdf'));

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
