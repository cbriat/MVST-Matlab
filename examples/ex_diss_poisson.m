%% EX_DISS_POISSON  Dissipativity of a stochastic system with Poissonian jumps.
%  Computes the mean-square L2-gain and the input-strict-passivity margin.
%  METHOD.  The mean-square L2-gain gamma* is certified GENERICALLY through the
%  SOSTOOLS routine MVST_DISS_SOS (dwell type 'arbitrary', supply 'l2gain',
%  mode 'gamma'): the routine assembles the flow dissipation operator C_M >= 0
%  of the Dissipativity chapter (eq:diss:bri:poisson), carrying the
%  state-multiplicative diffusion sum_l C_l^*P C_l, the disturbance channel
%  D_{w,l} (pivot loading), the Poisson jump-rate term lam (J^*P J - P), and
%  the Poisson jump-disturbance channel Br1 = D_1 (the reset x^+ = J1 x + D1 w,
%  loading the cross and pivot blocks) in the flow, and bisects on gamma.  The
%  previous exact shooting-Riccati value (via the MVST_L2GAIN dispatcher,
%  routing to MVST_L2GAIN_CT_POISSON) is printed as a cross-check.  Figures: convergence
%  of the finite-horizon available storage, and gamma* versus the jump rate
%  lambda (the sweep is recomputed through the same SOS routine).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- instance: stable drift, geometric noise, expansive resets -----
%  The disturbance acts on the drift (Bw), on the diffusion (Dw1) AND on the
%  spontaneous resets (Br1): the noise channel is Ct_{W,1} = [C1 Dw1] acting
%  on (x; w), and the jump channel is Jt_1 = [J1 Br1], so the reset is
%  x^+ = J1 x + Br1 w (a disturbance-driven Poissonian jump).
A  = [0 1; -2 -1];  Bw = [1;0];
C  = [1 0];         Dw = 0.5;
C1 = 0.15*eye(2);   Dw1 = [0.1; 0];
J1 = 1.15*eye(2);   lam = 0.5;   Br1 = [0; 0.2];   % Br1 = jump disturbance D_1
n = 2;

% ----- mean-square L2-gain: exact stochastic shooting-Riccati authoritative,
%       SOS cross-check.  (Exact keeps Ch08/Ch11 cross-consistent, both routed
%       through the same dispatched exact backend MVST_L2GAIN_CT_POISSON.)
sysd = struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw,'Ac',C1,'Dwc',Dw1, ...
              'Jf',J1,'lamf',lam,'Bdf',Br1);
%  The exact value is the shooting
%  test over a long horizon (1000) with relative bisection tolerance 1e-6, as
%  in the other Ch08 examples; the dispatched default (horizon 60, tolerance
%  1e-3) gave 2.0073, 1.2e-2 below the SOS value, and under-estimates the gain.
gam    = mvst_l2gain_ct_poisson(A,Bw,C,Dw,{C1},{Dw1},{J1},lam,{Br1},1000,1e-6);  % exact (long horizon)
gamDef = mvst_l2gain(sysd).gamma;                                             % former dispatched default (diagnostic)
gamSOS = mvst_diss_sos(sysd, struct('type','l2gain'), struct('type','arbitrary'), 'gamma').gamma;
nu  = isp_margin(A, Bw, C, Dw, C1, Dw1, J1, lam, Br1);
fprintf('Poisson dissipativity: gamma* %.4f, ISP margin nu %.4f\n', gam, nu);
fprintf('  [cross-check] exact shooting-Riccati=%.6f  SOS gamma*=%.6f  |diff|=%.2e (rel %.2e)\n', ...
        gam, gamSOS, abs(gam-gamSOS), abs(gam-gamSOS)/gam);
fprintf('  [diagnostic] dispatched default (horizon 60, tol 1e-3) = %.6f\n', gamDef);

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','diss_poisson.tex'));
w.mat('DApoi',A); w.mat('DBwpoi',Bw); w.mat('DCpoi',C); w.num('DDwpoi',Dw);
w.mat('DDwlpoi',Dw1);
w.mat('DJpoi',J1); w.num('DLamPoi',lam); w.mat('DBdpoi',Br1);
w.num('DGamPoi',gam); w.num('DNuPoi',nu);
w.num('DGamPoiSOS',gamSOS); w.raw('DGamPoiGap',scifmt(abs(gamSOS-gam)));
w.close();

% ----- figure 1: convergence of the available storage -----
gplot = 1.05*gam;
Tg = linspace(0.1, 8, 40); P0 = zeros(n,n,numel(Tg));
oo = odeset('RelTol',1e-7,'AbsTol',1e-9);
Rg = gplot^2 - Dw.'*Dw;
for k = 1:numel(Tg)
    [~, Y] = ode45(@(s,p) bre_rhs(p,A,Bw,C,Dw,C1,Dw1,J1,lam,Br1,Rg,n), ...
                   [0 Tg(k)], zeros(n*n,1), oo);
    Pk = reshape(Y(end,:).',n,n); P0(:,:,k) = (Pk+Pk.')/2;
end
S = P0;                                       % storage matrix -P_T(0) = +P here
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h = gobjects(3,1);
h(1) = plot(ax, Tg, reshape(S(1,1,:),1,[]), '-o','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(1,:));
h(2) = plot(ax, Tg, reshape(S(1,2,:),1,[]), '-s','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(2,:));
h(3) = plot(ax, Tg, reshape(S(2,2,:),1,[]), '-^','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(3,:));
xlabel(ax,'$T$','Interpreter','latex');
ylabel(ax,'available storage matrix $-P_T(0)$ entries','Interpreter','latex');
lg = legend(ax, h, {'$-P_{11}$','$-P_{12}$','$-P_{22}$'},'Interpreter','latex', ...
            'Orientation','horizontal','Location','southoutside');
set(lg,'Box','off'); set(ax,'FontSize',11);
mvst_savepdf(fig, fullfile(mvst_figdir(),'diss_poisson_conv.pdf'));

% ----- figure 2: gamma* versus the jump rate lambda -----
lams = linspace(0, 1.6, 22); gams = nan(size(lams));
for k = 1:numel(lams)
    sk = struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw,'Ac',C1,'Dwc',Dw1, ...
                'Jf',J1,'lamf',lams(k),'Bdf',Br1);
    rk = mvst_diss_sos(sk, struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
    if rk.gamma < 1e3; gams(k) = rk.gamma; end   % drop post-blow-up values
end
fig2 = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
ax2 = axes(fig2); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, lams, gams, '-o','LineWidth',2.0,'MarkerSize',3,'Color',gc.order(1,:));
xline(ax2, lam, ':','LineWidth',1.4,'Color',[0.5 0.5 0.5]);
xlabel(ax2,'jump rate $\lambda_1$','Interpreter','latex');
ylabel(ax2,'mean-square gain $\gamma^\star$','Interpreter','latex');
set(ax2,'FontSize',11);
mvst_savepdf(fig2, fullfile(mvst_figdir(),'diss_poisson_rate.pdf'));

% ----- combined 2-up: storage convergence (left) + gain vs rate (right) -----
figc = figure('Color','w','Units','centimeters','Position',[2 2 14 15]);
tlc = tiledlayout(figc,2,1,'TileSpacing','compact','Padding','compact');
axc = nexttile(tlc); hold(axc,'on'); box(axc,'on'); grid(axc,'on');
hc = gobjects(3,1);
hc(1) = plot(axc, Tg, reshape(S(1,1,:),1,[]), '-o','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(1,:));
hc(2) = plot(axc, Tg, reshape(S(1,2,:),1,[]), '-s','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(2,:));
hc(3) = plot(axc, Tg, reshape(S(2,2,:),1,[]), '-^','LineWidth',1.8,'MarkerSize',3,'Color',gc.order(3,:));
xlabel(axc,'$T$','Interpreter','latex');
ylabel(axc,'available storage matrix $-P_T(0)$ entries','Interpreter','latex');
lg = legend(axc, hc, {'$-P_{11}$','$-P_{12}$','$-P_{22}$'},'Interpreter','latex', ...
            'Orientation','horizontal','Location','southoutside');
set(lg,'Box','off'); set(axc,'FontSize',11);
title(axc,'(a) available-storage convergence','Interpreter','latex','FontWeight','normal');
axc = nexttile(tlc); hold(axc,'on'); box(axc,'on'); grid(axc,'on');
plot(axc, lams, gams, '-o','LineWidth',2.0,'MarkerSize',3,'Color',gc.order(1,:));
xline(axc, lam, ':','LineWidth',1.4,'Color',[0.5 0.5 0.5]);
xlabel(axc,'jump rate $\lambda_1$','Interpreter','latex');
ylabel(axc,'mean-square gain $\gamma^\star$','Interpreter','latex');
set(axc,'FontSize',11);
title(axc,'(b) gain vs jump rate','Interpreter','latex','FontWeight','normal');
mvst_savepdf(figc, fullfile(mvst_figdir(),'diss_poisson_combined.pdf'));

% ========================================================================
function dp = bre_rhs(p, A, Bw, C, Dw, C1, Dw1, J1, lam, Br1, Rg, n)
P = reshape(p,n,n); P = (P+P.')/2;
Sg = P*Bw + C1.'*P*Dw1 + lam*(J1.'*P*Br1) + C.'*Dw;
Th = Rg - Dw1.'*P*Dw1 - lam*(Br1.'*P*Br1);   % storage-dependent pivot Theta_S
G = A.'*P + P*A + C1.'*P*C1 + lam*(J1.'*P*J1 - P) + C.'*C + Sg*(Th\Sg.');
G = (G+G.')/2; dp = G(:);
end

function nu = isp_margin(A, Bw, C, Dw, C1, Dw1, J1, lam, Br1)
% Input-strict-passivity margin: the largest nu < He(Dw) such that the
% stochastic positive-real Riccati
%   -Sdot = A'S + S A + C1'S C1 + lam (J1'S J1 - S)
%           + v(S)(He(Dw) - nu I - Dw1'S Dw1 - lam Br1'S Br1)^{-1} v(S)',
%   v(S) = S Bw + C1'S Dw1 + lam J1'S Br1 - C',
% stays bounded when integrated backward from S(T)=0 (shooting test).  The
% jump disturbance Br1 loads the passivity pivot exactly as Dw1 does.
n = size(A,1); hd = Dw + Dw.';
lo = 0; hi = hd - 1e-6; Thor = 40; big = 1e6;
    function ok = feas(v)
        Rv = hd - v;
        if min(real(eig((Rv+Rv.')/2))) <= 0; ok = false; return; end
        % Stop the shooting either on finite escape (norm > big) OR when the
        % storage-dependent pivot loses positive-definiteness (the Riccati is
        % diverging): both mark the trial infeasible.  MaxStep bounds the step
        % count so the integration cannot grind to a halt near the boundary.
        opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',Thor/200, ...
                      'Events',@(s,p) esc_event(p,Rv,big));
        [ss,~,~,~,ie] = ode45(@(s,p) pr_rhs(p,Rv), [0 Thor], zeros(n*n,1), opts);
        ok = (ss(end) >= Thor - 1e-6) && isempty(ie);
    end
    function dp = pr_rhs(p, Rv)
        P = reshape(p,n,n); P = (P+P.')/2;
        Sg = P*Bw + C1.'*P*Dw1 + lam*(J1.'*P*Br1) - C.';
        Piv = Rv - Dw1.'*P*Dw1 - lam*(Br1.'*P*Br1);   % storage-dependent pivot
        pm = min(real(eig((Piv+Piv.')/2)));
        if pm <= 1e-9; dp = zeros(n*n,1); return; end   % guard against blow-up
        G = A.'*P + P*A + C1.'*P*C1 + lam*(J1.'*P*J1 - P) + Sg*(Piv\Sg.');
        G = (G+G.')/2; dp = G(:);
    end
    function [v, tm, dr] = esc_event(p, Rv, cap)
        P = reshape(p,n,n); P = (P+P.')/2;
        Piv = Rv - Dw1.'*P*Dw1 - lam*(Br1.'*P*Br1);
        pm = min(real(eig((Piv+Piv.')/2)));
        v = min(cap - norm(p), pm - 1e-8); tm = 1; dr = -1;
    end
if ~feas(0); nu = 0; return; end
while (hi - lo) > 1e-4
    mid = 0.5*(lo+hi);
    if feas(mid); lo = mid; else; hi = mid; end
end
nu = lo;
end

% ========================================================================
function s = scifmt(x)
%SCIFMT  LaTeX scientific notation m\cdot10^{e} with one decimal.
if x <= 0; s = '0'; return; end
e = floor(log10(x)); mnt = x/10^e;
if round(mnt*10)/10 >= 10; mnt = mnt/10; e = e+1; end
s = sprintf('%.1f\\cdot10^{%d}', mnt, e);
end
