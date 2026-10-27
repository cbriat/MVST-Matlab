%% EX_STAB_STABZ  Feedback stabilization: CT, DT, and impulsive instances.
%  Chapter "Stability", Section "Examples", subsection "Feedback
%  stabilization".  For each system class a stabilizing state feedback is
%  synthesized and the closed-loop second moment X is compared with the
%  open-loop one.  All second-moment trajectories are propagated EXACTLY
%  through the reduced matrix representation of the (flow / transition /
%  monodromy) operator on the symmetric subspace (mvst_liftsym /
%  mvst_svecbasis), so no ODE solver or toolbox is needed.  Stability is
%  read as: CT lifted-flow abscissa < 0; DT lifted-transition radius < 1;
%  impulsive lifted-monodromy radius over one period < 1.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();
n = 2; U = mvst_svecbasis(n); X0 = eye(n);
sx0 = U' * reshape(X0,[],1);
X0b = diag([2.2, 0.5]);            % 2nd initial condition (show stable convergence)
sx0b = U' * reshape(X0b,[],1);

% ======================================================================
% (a) Continuous-time stabilization: dx = (A0 x + B0 u) dt + A1 x dW1
% ======================================================================
A0 = [0.3 1; 0 0.2];  B0 = eye(n);  A1 = 0.2*eye(n);
Mopen = mvst_liftsym(A0, {A1});
absCtOpen = max(real(eig(Mopen)));

% synthesize the CT gain through the associated (I-)Riccati equation
md = mvst_mode(A0, B0, eye(n), eye(n), 'Ac', A1, 'Bc', zeros(n,n), 'X0', X0);
[~, Kct] = mvst_care_ct(md);
Mcl = mvst_liftsym(A0 + B0*Kct, {A1});
if max(real(eig(Mcl))) >= 0                       % safety fallback
    Kct = B0 \ (-eye(n) - A0);
    Mcl = mvst_liftsym(A0 + B0*Kct, {A1});
end
absCtClosed = max(real(eig(Mcl)));

tc = linspace(0, 6, 200);
trO = zeros(size(tc)); trCa = zeros(size(tc)); trCb = zeros(size(tc));
for k = 1:numel(tc)
    trO(k)  = trace_of(U, expm(tc(k)*Mopen)*sx0,  n);   % open loop  (unstable)
    trCa(k) = trace_of(U, expm(tc(k)*Mcl  )*sx0,  n);   % closed loop (stable), IC 1
    trCb(k) = trace_of(U, expm(tc(k)*Mcl  )*sx0b, n);   % closed loop (stable), IC 2
end
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, tc, trCa, '-', 'LineWidth',2.0,'Color',gc.order(1,:));                          % stable = blue
plot(ax, tc, trCb, '-', 'LineWidth',2.0,'Color',gc.order(1,:),'HandleVisibility','off'); % stable, 2nd IC
plot(ax, tc, trO,  '--','LineWidth',2.0,'Color',gc.order(2,:));                          % unstable = red
set(ax,'YScale','log');
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'$\mathrm{tr}\,X(t)$','Interpreter','latex');
lg = legend(ax,{'closed loop (stable)','open loop (unstable)'},'Interpreter','latex', ...
            'Location','northwest'); set(lg,'Box','off');
mvst_savepdf(fig, fullfile(figdir,'ex_stab_ct_stab.pdf'));
fprintf('CT stabilization: open abscissa %.4f, closed %.4f\n', ...
        absCtOpen, absCtClosed);

% ======================================================================
% (b) Discrete-time stabilization: x+ = A0 x + B0 u + nu A1 x
% ======================================================================
A0d = [1.1 0.3; 0 1.05];  B0d = eye(n);  A1d = 0.15*eye(n);
Topen = mvst_liftsym(zeros(n), {A0d, A1d});
radDtOpen = max(abs(eig(Topen)));
Ades = 0.4*eye(n);  Kdt = B0d \ (Ades - A0d);      % eigenvalue assignment
Tcl = mvst_liftsym(zeros(n), {A0d + B0d*Kdt, A1d});
radDtClosed = max(abs(eig(Tcl)));

Nstep = 15; kk = 0:Nstep;
trOd = zeros(size(kk)); trCda = zeros(size(kk)); trCdb = zeros(size(kk));
sO = sx0; sCa = sx0; sCb = sx0b;
for j = 1:numel(kk)
    trOd(j) = trace_of(U, sO, n);  trCda(j) = trace_of(U, sCa, n);  trCdb(j) = trace_of(U, sCb, n);
    sO = Topen*sO;  sCa = Tcl*sCa;  sCb = Tcl*sCb;
end
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
stairs(ax, kk, trCda, '-', 'LineWidth',2.0,'Color',gc.order(1,:));                          % stable = blue
stairs(ax, kk, trCdb, '-', 'LineWidth',2.0,'Color',gc.order(1,:),'HandleVisibility','off'); % stable, 2nd IC
stairs(ax, kk, trOd,  '--','LineWidth',2.0,'Color',gc.order(2,:));                          % unstable = red
set(ax,'YScale','log');
xlabel(ax,'step $k$','Interpreter','latex');
ylabel(ax,'$\mathrm{tr}\,X(k)$','Interpreter','latex');
lg = legend(ax,{'closed loop (stable)','open loop (unstable)'},'Interpreter','latex', ...
            'Location','northwest'); set(lg,'Box','off');
mvst_savepdf(fig, fullfile(figdir,'ex_stab_dt_stab.pdf'));
fprintf('DT stabilization: open radius %.4f, closed %.4f\n', ...
        radDtOpen, radDtClosed);

% ======================================================================
% (c) Impulsive stabilization under periodic dwell-time:
%     flow dx=(A0 x + B0 u)dt (unstable), jump x+ = J x (contractive)
% ======================================================================
A0i = [0.8 1; 0 0.5];  B0i = eye(n);  Jmat = 0.7*eye(n);  Tbar = 0.6;
MfOpen = mvst_liftsym(A0i, {});
MJ = mvst_liftsym(zeros(n), {Jmat});
monoOpen = max(abs(eig(MJ*expm(Tbar*MfOpen))));
Ades_i = -0.3*eye(n);  Kimp = B0i \ (Ades_i - A0i);   % stabilize the flow
MfCl = mvst_liftsym(A0i + B0i*Kimp, {});
monoClosed = max(abs(eig(MJ*expm(Tbar*MfCl))));

Ncyc = 6; npt = 30;
[tO2,  trO2 ] = saw_traj(U, MfOpen, MJ, Tbar, Ncyc, npt, sx0,  n);   % open loop  (unstable)
[tC2a, trC2a] = saw_traj(U, MfCl,   MJ, Tbar, Ncyc, npt, sx0,  n);   % closed loop (stable), IC 1
[tC2b, trC2b] = saw_traj(U, MfCl,   MJ, Tbar, Ncyc, npt, sx0b, n);   % closed loop (stable), IC 2
fig = figure('Color','w','Units','centimeters','Position',[2 2 12 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, tC2a, trC2a, '-', 'LineWidth',1.8,'Color',gc.order(1,:));                          % stable = blue
plot(ax, tC2b, trC2b, '-', 'LineWidth',1.8,'Color',gc.order(1,:),'HandleVisibility','off'); % stable, 2nd IC
plot(ax, tO2,  trO2,  '--','LineWidth',1.8,'Color',gc.order(2,:));                          % unstable = red
set(ax,'YScale','log');
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'$\mathrm{tr}\,X(t)$','Interpreter','latex');
lg = legend(ax,{'closed loop (stable)','open loop (unstable)'},'Interpreter','latex', ...
            'Location','northwest'); set(lg,'Box','off');
mvst_savepdf(fig, fullfile(figdir,'ex_stab_imp_stab.pdf'));
fprintf('Impulsive stabilization: open monodromy %.4f, closed %.4f\n', ...
        monoOpen, monoClosed);

% ---------- export data macros to LaTeX ----------
w = mvst_texdata(fullfile(figdir,'data','ex_stab_stabz.tex'));
w.num('StabCtOpen',absCtOpen);   w.num('StabCtClosed',absCtClosed);
w.num('StabDtOpen',radDtOpen);   w.num('StabDtClosed',radDtClosed);
w.num('StabImpOpen',monoOpen);   w.num('StabImpClosed',monoClosed);
w.num('StabImpTbar',Tbar);
w.close();

% ========================================================================
function tr = trace_of(U, sx, n)
X = reshape(U*sx, n, n);  tr = trace(X);
end

function [t, tr] = saw_traj(U, Mf, MJ, Tbar, Ncyc, npt, sx0, n)
%SAW_TRAJ  Trace of X(t) over Ncyc dwell-time cycles: flow then jump.
t = []; tr = []; sx = sx0; t0 = 0;
for c = 1:Ncyc
    tt = linspace(0, Tbar, npt);
    for j = 1:npt
        t(end+1)  = t0 + tt(j);                          %#ok<AGROW>
        tr(end+1) = trace_of(U, expm(tt(j)*Mf)*sx, n);   %#ok<AGROW>
    end
    sx = MJ * (expm(Tbar*Mf)*sx);       % apply the jump at end of cycle
    t0 = t0 + Tbar;
end
end
