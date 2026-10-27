%% EX_TV_SD  Time-varying sampled-data optimal control: augmented dual P(tau).
%  Time-varying counterpart of EX02_SAMPLEDDATA_GAINS: the flow drift decays
%  in time, A(t) = A + dA*exp(-t), so the system is genuinely time-varying but
%  settles to the LTI sampled-data system of EX02 as t grows.  The dual P(tau)
%  lives on the AUGMENTED state z = [x;u] (size n+m): between samples we
%  integrate the augmented flow-Riccati of the dual with TIME-VARYING flow
%  Ab0(t) = [A(t) B; 0 0], and at each sample we apply the Schur (jump) update
%  P <- blkdiag(Pxx - Pxu*(Puu\Pux), 0) with the held-input gain K=-(Puu\Pux).
%  This illustrates the main contribution of the paper for the impulsive
%  (sampled-data) case over a time-varying horizon.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
s = D.sd; n = size(s.A,1); m = size(s.B,2); Tbar = s.Tbar;
dA = [0 0; 1.5 0];                          % decaying drift perturbation
Afun = @(t) s.A + dA*exp(-t);               % time-varying flow drift

% ----- backward pass over [0,T] with several sampling periods -----
% Augmented dual Riccati Pdot = -(Ab0(t)'P + P Ab0(t) + Abl'P Abl + diag(Q,R))
% integrated backward between samples, with the Schur jump at each sample.
T = 6*Tbar; d = n+m; ix = 1:n; iu = n+1:d;
Abl = [s.Ac s.Bc; zeros(m,d)];
Zf  = blkdiag((s.Q+s.Q.')/2, (s.R+s.R.')/2);
tk  = Tbar:Tbar:(T-1e-9); N = numel(tk); edges = [0 tk T];
rt  = 1e-8; oo = odeset('RelTol',rt,'AbsTol',rt*1e-2);

% store the within-period dual P(tau) of the first TWO periods near t=0 for
% plotting; also extract the converged (last solved, smallest-t) gain Klast.
P = zeros(d);                               % terminal P(T)=0
Klast = zeros(m,n);
Pfun1 = []; Pfun2 = [];                     % P(tau) over the two periods at t~0
for i = N+1:-1:1
    [ts, Pv] = ode45(@(t,pv) flowback(t,pv,Afun,s.B,Abl,Zf,d), ...
                     [edges(i+1) edges(i)], P(:), oo);
    P = reshape(Pv(end,:),d,d); P = (P+P.')/2;
    if i == 2; Pfun2 = struct('t',ts,'P',Pv); end   % second period from t=0
    if i == 1; Pfun1 = struct('t',ts,'P',Pv); end   % first  period from t=0
    if i >= 2
        Puu = P(iu,iu); Pux = P(iu,ix); Pxx = P(ix,ix); Pxu = P(ix,iu);
        K = -(Puu \ Pux); Klast = K;        % last computed = nearest to t=0
        sc = Pxx - Pxu*(Puu \ Pux); P = blkdiag((sc+sc.')/2, zeros(m));
    end
end
fprintf('TV sampled-data: T=%.2f (%d samples), converged gain K=[% .4f % .4f]\n', ...
        T, N, Klast);

% ----- sample the augmented P(tau) over the first two periods (near t=0) -----
[tau1, Pq1] = sample_period(Pfun1, edges(1), edges(2), d);
[tau2, Pq2] = sample_period(Pfun2, edges(2), edges(3), d);
tt = [tau1, NaN, tau2 + Tbar];             % second period plotted over [Tbar,2Tbar]
Px1x1 = [squeeze(Pq1(1,1,:)).'   NaN squeeze(Pq2(1,1,:)).'];
Px2x2 = [squeeze(Pq1(2,2,:)).'   NaN squeeze(Pq2(2,2,:)).'];
Px1x2 = [squeeze(Pq1(1,2,:)).'   NaN squeeze(Pq2(1,2,:)).'];
Px1u  = [squeeze(Pq1(1,n+1,:)).' NaN squeeze(Pq2(1,n+1,:)).'];
Px2u  = [squeeze(Pq1(2,n+1,:)).' NaN squeeze(Pq2(2,n+1,:)).'];
Puuv  = [squeeze(Pq1(n+1,n+1,:)).' NaN squeeze(Pq2(n+1,n+1,:)).'];

% ----- export data macros -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_sd.tex'));
w.mat('TVdAsd', dA); w.num('TVTbarsd', Tbar); w.mat('TVKsd', Klast);
w.close();

% ----- figure: augmented dual P(tau) over two periods (jump visible) -----
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, tt, Px1x1, '-','LineWidth',2.2,'Color',gc.order(1,:));
plot(ax, tt, Px2x2, '-','LineWidth',2.2,'Color',gc.order(2,:));
plot(ax, tt, Px1x2, '-','LineWidth',2.2,'Color',gc.order(3,:));
plot(ax, tt, Px1u , '-','LineWidth',2.2,'Color',gc.order(4,:));
plot(ax, tt, Px2u , '-','LineWidth',2.2,'Color',gc.order(5,:));
plot(ax, tt, Puuv , '-','LineWidth',2.2,'Color',gc.order(6,:));
xline(ax, Tbar, ':','Color',[0.6 0.6 0.6]);
xlabel(ax,'timer $\tau$ (sampling period $\bar T$)','Interpreter','latex');
ylabel(ax,'entries of augmented $P(\tau)$','Interpreter','latex');
lg = legend(ax, {'$P_{x_1x_1}$','$P_{x_2x_2}$','$P_{x_1x_2}$', ...
      '$P_{x_1u}$','$P_{x_2u}$','$P_{uu}$'}, ...
      'Interpreter','latex','Location','northeastoutside'); set(lg,'Box','off');
set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_tv_sd.pdf'));

% ========================================================================
function dp = flowback(t, pv, Afun, B, Abl, Zf, d)
%FLOWBACK  -RHS of the augmented time-varying dual flow-Riccati.
n = size(B,1); m = size(B,2);
Ab0 = [Afun(t) B; zeros(m,n+m)];            % time-varying augmented flow
P = reshape(pv,d,d); P = (P+P.')/2;
dP = Ab0.'*P + P*Ab0 + Abl.'*P*Abl + Zf;
dp = reshape(-(dP+dP.')/2, [], 1);
end

% ========================================================================
function [tau, Pq] = sample_period(seg, t0, t1, d)
%SAMPLE_PERIOD  Resample a stored ode45 segment onto a uniform timer grid.
%   SEG has fields .t, .P (rows = P(:) at SEG.t).  Returns the within-period
%   timer tau in [0, t1-t0] and Pq(:,:,k) = P at the corresponding time.
tau = linspace(0, t1-t0, 60);
tt  = t0 + tau;                             % absolute times to query
[st, ii] = sort(seg.t); Pv = seg.P(ii,:);
[st, iu] = unique(st);  Pv = Pv(iu,:);
Pq = zeros(d,d,numel(tau));
Vq = interp1(st, Pv, min(max(tt,st(1)),st(end)), 'pchip');
for k = 1:numel(tau)
    Pk = reshape(Vq(k,:),d,d); Pq(:,:,k) = (Pk+Pk.')/2;
end
end
