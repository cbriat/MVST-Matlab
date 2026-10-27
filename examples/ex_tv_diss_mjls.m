%% EX_TV_DISS_MJLS  Time-varying dissipativity of a Markov-jump-with-diffusion.
%  Time-varying counterpart of EX_DISS_MJLS (switching-stabilized instance:
%  mode 1 unstable, mode 2 stable).  The stable-mode flow drift decays in time,
%  A_2(t) = A_2 + dA*exp(-t), so the system is genuinely time-varying but
%  settles to the LTI Markov-jump system as t grows.  The reported gain is the
%  coupled mean-square gain and the reported passivity index that of the
%  stable mode, both computed exactly as in EX_DISS_MJLS, so the limits
%  coincide.
%  The figures carry the
%  diffusion channels of EX_DISS_MJLS (state-multiplicative C_{2,1} = 0.2 I2
%  and disturbance-multiplicative column D_{w,1} = 0.1[1;0], same Brownian
%  motion):
%    conv : finite-horizon TIME-VARYING available storage S_T(0) = -P_T(0) >= 0
%           of the stable mode, from the stochastic bounded-real Riccati with
%           both channels;
%    traj : sample paths and mean square of the time-varying stable mode
%           driven by the Ornstein-Uhlenbeck disturbance (eq:diss:ex:ou,
%           rate 5), the mean square from the second moment of (x,w).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- switching-stabilized performance variant (matches EX_DISS_MJLS) -----
A1 = [1 2;0 1];     A2 = [-3 0;0 -3];        % mode 1 unstable, mode 2 stable
Bw = [1;0];         C1 = [1 0];   C2 = [0.5 1];   Dw = 0.5;
Pig = [-4 4; 1 -1];
Cs  = 0.2*eye(2);   Dwcm = 0.1*[1;0];
Cchan = {Cs, Cs};   Dwc = {Dwcm};
aOU = 5;                                      % OU rate of the trajectory figure (as EX_DISS_SD/SW)
dA = [0 0; 1.5 0];                           % decaying drift perturbation
A0 = A2;                                      % stable-mode asymptotic drift (storage converges)
Afun = @(t) A0 + dA*exp(-t);                  % time-varying flow drift (stable mode)

% ----- coupled mean-square gain (exact), coincides with the asymptotic limit -
gam = mvst_l2gain_ct_mjls({A1,A2}, Bw, {C1,C2}, Dw, Pig, 1000, Cchan, Dwc, 1e-6);
nu  = mvst_isp_margin_ct(A2, Bw, C2, Dw, {Cs}, Dwc);  isp = nu > 0;   % stable-mode ISP index

% ----- convergence of the finite-horizon TIME-VARYING storage (both channels) -
g0 = tv_gain0(Afun,Bw,C2,Dw,Cs,Dwcm,200,1e-6);          % gain of the plotted time-varying system from t_0 = 0
gplot = 1.05*max(gam, g0);                               % storage level (5% above both)
[Tg,S0] = tv_diss_pT_dwc(Afun,Bw,C2,Dw,Cs,Dwcm, gplot, 20, 50);
assert(all(arrayfun(@(k) min(eig(S0(:,:,k))), 1:numel(Tg)) >= -1e-10), 'storage must be PSD');

fprintf('TV-MJLS dissipativity: coupled mean-square gain %.4f (switching-stabilized)\n', gam);
fprintf('TV-MJLS passivity (stable mode): isp=%d, nu=%.4f\n', isp, nu);
fprintf('  time-varying gain from t_0=0: %.4f; storage level %.4f; S_T(0) at T=%.1f: %s (change over last step %.1e)\n', ...
    g0, gplot, Tg(end), mat2str(S0(:,:,end),4), norm(S0(:,:,end)-S0(:,:,end-1)));

% ----- export data macros to LaTeX -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_diss_mjls.tex'));
w.mat('TVdAdmjls',dA); w.num('TVDGamMjls',gam); w.num('TVDNuMjls',nu); w.num('TVDGamMjlsTV',g0);
w.close();

% ----- figure 1: convergence of the (time-varying) available storage -----
diss_conv_fig(Tg, S0, fullfile(mvst_figdir(),'ex_tv_diss_mjls_conv.pdf'), gc);

% ----- figure 2: sample paths + mean square of the time-varying stable mode -----
rng(7);
tv_traj_fig({Afun}, @(t) 1, Bw, Cs, Dwcm, aOU, 12, ...
    fullfile(mvst_figdir(),'ex_tv_diss_mjls_traj.pdf'), gc);

% ========================================================================
function [Tgrid, S0] = tv_diss_pT_dwc(Afun, Bw, C, Dw, Ch, D1, gam, Tmax, nT)
%TV_DISS_PT_DWC  Finite-horizon TIME-VARYING available storage S_T(0) = -P_T(0)
%   versus the horizon T, for the diffusion
%     dx = (A(t) x + Bw w)dt + (Ch x + D1 w)dW,   z = C x + Dw w,
%   from the stochastic bounded-real Riccati (the time-varying version of the
%   one of the time-invariant examples)
%     -dS/dt = A(t)'S + S A(t) + Ch'S Ch + C'C + v Theta^{-1} v',
%     v = S Bw + Ch'S D1 + C'Dw,  Theta = gam^2 I - Dw'Dw - D1'S D1,
%   integrated backward in t from S(T) = 0 to t = 0.  S0(:,:,k) >= 0 is the
%   storage for the horizon Tgrid(k) (diffusion column D1 and
%   state channel Ch included; the storage is returned and plotted directly).
n = size(Bw,1); m = size(Bw,2);
Tgrid = linspace(0.1, Tmax, nT); S0 = zeros(n,n,nT);
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
for k = 1:nT
    [~,Y] = ode45(@(t,q) -bre(t,q), [Tgrid(k) 0], zeros(n*n,1), opts);
    S = reshape(Y(end,:).',n,n); S0(:,:,k) = (S+S.')/2;     % last row is t = 0
end
    function dq = bre(t,q)
        A = Afun(t);
        S = reshape(q,n,n); S = (S+S.')/2;
        v = S*Bw + Ch.'*S*D1 + C.'*Dw;
        Th = gam^2*eye(m) - Dw.'*Dw - D1.'*S*D1; Th = (Th+Th.')/2;
        G = A.'*S + S*A + Ch.'*S*Ch + C.'*C + v*(Th\v.'); G = (G+G.')/2; dq = G(:);
    end
end

% ========================================================================
function diss_conv_fig(Tg, S0, fname, gc)
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
function tv_traj_fig(Afuns, sched, Bw, Ch, D1, aOU, Tend, fname, gc)
%TV_TRAJ_FIG  Standalone wrapper around TV_TRAJ_AX.
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
tv_traj_ax(axes(fig), Afuns, sched, Bw, Ch, D1, aOU, Tend, '', gc);
mvst_savepdf(fig, fname);
end

% ========================================================================
function tv_traj_ax(ax, Afuns, sched, Bw, Ch, D1, aOU, Tend, ttl, gc)
%TV_TRAJ_AX  Sample paths and mean square of the time-varying diffusion
%     dx = (A_sched(t)(t) x + Bw w)dt + (Ch x + D1 w)dW,
%   driven by the Ornstein-Uhlenbeck disturbance dw = -aOU w dt + sqrt(2 aOU) dV
%   (V independent of W; unit stationary variance, eq:diss:ex:ou), from rest
%   x(0)=0, w(0)=0.  AFUNS{i}(t) is the (time-varying) drift of mode i and
%   SCHED(t) the active mode.  Euler-Maruyama sample paths of |x(t)|^2 (thin),
%   the mean square from the second moment Sigma of (x,w),
%     dSigma/dt = Ab(t) Sigma + Sigma Ab(t)' + Cb Sigma Cb' + bb bb',
%     Ab(t) = [A(t) Bw; 0 -aOU],  Cb = [Ch D1; 0 0],  bb = [0; sqrt(2 aOU)]  (thick),
%   propagated by the matrix exponential of the vectorized operator on each
%   grid step of length 2e-3 with the drift evaluated at the step midpoint
%   (the mode is constant on each step), and the Monte-Carlo average (dotted).
dt = 2e-3; tg = 0:dt:Tend; nt = numel(tg);
npath = 300; showp = 14; n = size(Bw,1); N = n+1;
sq = zeros(nt,npath); sOU = sqrt(2*aOU);
for p = 1:npath
    x = zeros(n,1); w = 0;
    for k = 1:nt
        sq(k,p) = x.'*x;
        if k==nt; break; end
        A = Afuns{sched(tg(k))}(tg(k));
        dW = sqrt(dt)*randn; dV = sqrt(dt)*randn;
        x = x + (A*x + Bw*w)*dt + (Ch*x + D1*w)*dW;
        w = w - aOU*w*dt + sOU*dV;
    end
end
bvec = [zeros(n,1); sOU]; f = reshape(bvec*bvec.', [], 1);
Cb = [Ch, D1; zeros(1,N)];
v = zeros(N*N,1); Ems = zeros(1,nt);
for k = 1:nt
    Sig = reshape(v,N,N); Ems(k) = real(trace(Sig(1:n,1:n)));
    if k==nt; break; end
    tm = tg(k) + dt/2;
    Ab = [Afuns{sched(tg(k))}(tm), Bw; zeros(1,n), -aOU];
    L  = kron(eye(N),Ab) + kron(Ab,eye(N)) + kron(Cb,Cb);
    y  = expm([L, f; zeros(1,N*N+1)]*dt)*[v; 1]; v = y(1:N*N);
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
    {'sample paths $|x(t)|^2$','mean square (second moment)','Monte-Carlo average'}, ...
    'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off');
set(ax,'FontSize',11); ylim(ax,[0, 1.05*max(Ems)*4]);
if ~isempty(ttl); title(ax, ttl, 'Interpreter','latex','FontWeight','normal'); end
end

% ========================================================================
function g0 = tv_gain0(Afun, Bw, C, Dw, Ch, D1, H, rtol)
%TV_GAIN0  Mean-square L2-gain on [0, H] of the TIME-VARYING diffusion
%     dx = (A(t) x + Bw w)dt + (Ch x + D1 w)dW,   z = C x + Dw w,   x(0) = 0:
%   the smallest gamma for which the time-varying stochastic bounded-real
%   Riccati of TV_DISS_PT_DWC, integrated backward from S(H) = 0, reaches t = 0
%   without finite escape and with Theta > 0 (bisection, relative bracket RTOL;
%   the long horizon H makes it the gain from t_0 = 0 on [0, infinity) once A(t)
%   has settled).  The storage curves are plotted 5% above max(gamma*, g0), so
%   that they converge (at 1.05 gamma* the sampled-data storage
%   escapes for horizons beyond about 8, since A(0) is not Hurwitz).
n = size(Bw,1); m = size(Bw,2);
lo = norm(Dw) + 1e-6; hi = max(2*lo, 1);
while ~feas(hi); lo = hi; hi = 2*hi; end
while hi - lo > rtol*hi
    mid = 0.5*(lo+hi); if feas(mid); hi = mid; else; lo = mid; end
end
g0 = hi;
    function ok = feas(g)
        opts = odeset('RelTol',1e-9,'AbsTol',1e-11,'Events',@(t,q) ev(q,g));
        ws = warning('off','all');                     % escapes are expected below the gain
        try
            [t,Y] = ode45(@(t,q) -bre(t,q,g), [H 0], zeros(n*n,1), opts);
            ok = abs(t(end)) < 1e-9 && all(isfinite(Y(end,:)));
        catch
            ok = false;
        end
        warning(ws);
    end
    function [v,te,d] = ev(q,g)
        S = reshape(q,n,n); S = (S+S.')/2;
        Th = g^2*eye(m) - Dw.'*Dw - D1.'*S*D1;
        v = [1e6 - norm(q); min(eig((Th+Th.')/2)) - 1e-9]; te = [1;1]; d = [0;0];
    end
    function dq = bre(t,q,g)
        A = Afun(t);
        S = reshape(q,n,n); S = (S+S.')/2;
        v = S*Bw + Ch.'*S*D1 + C.'*Dw;
        Th = g^2*eye(m) - Dw.'*Dw - D1.'*S*D1; Th = (Th+Th.')/2;
        G = A.'*S + S*A + Ch.'*S*Ch + C.'*C + v*(Th\v.'); G = (G+G.')/2; dq = G(:);
    end
end
