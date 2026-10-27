function ex_pos_switched()
%EX_POS_SWITCHED  Linear positive switched system under a minimum dwell-time
%  (Chapter "Linear Positive Systems"; the switched analog of the impulsive
%  Example ex_pos, following the switched-system results of the paper
%  "Dwell-time stability and stabilization conditions for linear positive
%  impulsive and switched systems", Corollary 5.5).
%
%  Two Metzler, Hurwitz shear modes with no common linear copositive
%  Lyapunov function:
%     A1 = [-1 0; 3 -1],   A2 = [-1 3; 0 -1]   (both eig = -1, -1).
%  Each mode is stable, but A1+A2 = [-2 3; 3 -2] has an eigenvalue +1, so
%  arbitrarily fast switching is UNSTABLE. Stability is recovered only above
%  a minimum dwell-time. The exact threshold is the worst-case periodic cycle
%  (dwell Tbar in each mode): rho(e^{A2 Tbar} e^{A1 Tbar}) crosses 1 at
%  Tbar*. Above it the mode-dependent linear copositive certificate
%     A_i' p_i < 0,   e^{A_j' Tbar} p_i < p_j   (i != j),
%  a linear program in p1, p2 > 0 (Corollary 5.5(b)), is feasible. The
%  switched state stays nonnegative (both modes Metzler) throughout.
%
%  Figure: pos_switched.pdf, a positive switched trajectory above the
%  threshold next to the worst-case cycle radius versus the dwell-time.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

A1 = [-1 0; 3 -1];                 % Metzler, Hurwitz (eig -1,-1)
A2 = [-1 3; 0 -1];                 % Metzler, Hurwitz (eig -1,-1)
assert(all(A1(~eye(2))>=0) && all(A2(~eye(2))>=0),'A1,A2 not Metzler');
Asum = A1+A2;
fprintf('eig(A1)=%s  eig(A2)=%s  eig(A1+A2)=%s (fast switching unstable)\n', ...
        mat2str(round(eig(A1)',3)), mat2str(round(eig(A2)',3)), ...
        mat2str(round(eig(Asum)',3)));

% ---- exact minimum dwell-time: rho(e^{A2 T} e^{A1 T}) crosses 1 ----
Tg  = linspace(0.01, 3.0, 500);
rho = arrayfun(@(T) max(abs(eig(expm(A2*T)*expm(A1*T)))), Tg);
i0  = find(rho < 1, 1, 'first');
lo = Tg(max(i0-1,1)); hi = Tg(i0);
for it=1:50, m=0.5*(lo+hi);
    if max(abs(eig(expm(A2*m)*expm(A1*m))))<1, hi=m; else, lo=m; end
end
Tstar = hi;
fprintf('exact minimum dwell-time threshold Tbar* = %.4f (rho of cycle = 1)\n', Tstar);

% ---- mode-dependent copositive LP: certificate + feasibility threshold ----
% One call to the generic switched minimum-dwell-time program: the certificate
% p1,p2 at the stable dwell Tstab, plus the bisected LP feasibility threshold.
Tstab = 2.0;                       % above both thresholds
rc = mvst_pos_stability(struct('A',{{A1,A2}}), ...
        struct('class','switched','dwell','min','Tmin',Tstab, ...
               'threshold',true,'Trange',[0.01 4],'Titers',60));
Tlp = rc.Tstar;
fprintf('static copositive-LP feasibility threshold Tbar_LP = %.4f (>= Tbar*)\n', Tlp);
p1 = rc.p{1}; p2 = rc.p{2}; ok = rc.stable;
fprintf(['copositive certificate at Tbar=%.2f (feasible=%d):\n', ...
         '  p1=[%.3f %.3f]>0  A1''p1=[%.3f %.3f]<0  e^{A2''T}p1-p2=[%.3f %.3f]<0\n', ...
         '  p2=[%.3f %.3f]>0  A2''p2=[%.3f %.3f]<0  e^{A1''T}p2-p1=[%.3f %.3f]<0\n'], ...
        Tstab, ok, p1(1),p1(2), A1.'*p1, expm(A2.'*Tstab)*p1-p2, ...
        p2(1),p2(2), A2.'*p2, expm(A1.'*Tstab)*p2-p1);

% ---- simulate the switched trajectory (periodic, dwell Tstab) ----
[ts, xs, sw] = sw_traj(A1, A2, Tstab, 6, [1;0.3]);   % above threshold: stable
fprintf('min state entry over trajectory = %.4f (>= 0: positive)\n', min(xs(:)));

% ================= figure =================
figure('Color','w','Units','centimeters','Position',[2 2 13 13]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
% shade the mode-2 dwell intervals for context
yl = [0 max(xs(:))*1.15];
for s = sw
    if s.mode==2
        patch(ax, [s.t0 s.t1 s.t1 s.t0], [yl(1) yl(1) yl(2) yl(2)], ...
              [.5 .5 .5],'FaceAlpha',0.08,'EdgeColor','none');
    end
end
h1 = plot(ax, ts, xs(1,:), '-','LineWidth',1.8,'Color',gc.order(1,:));
h2 = plot(ax, ts, xs(2,:), '-','LineWidth',1.8,'Color',gc.order(2,:));
yline(ax, 0, '-','Color','k','LineWidth',0.8);
xlabel(ax,'time $t$','Interpreter','latex'); ylabel(ax,'state $x(t)\ge0$','Interpreter','latex');
title(ax,sprintf('$\\bar T=%.1f>\\bar T^\\star$ (stable, positive)',Tstab),'Interpreter','latex','FontWeight','normal');
lg=legend(ax,[h1 h2],{'$x_1$','$x_2$'},'Interpreter','latex','Location','northeast'); set(lg,'Box','off');
set(ax,'FontSize',10); ylim(ax,yl);

ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, Tg, rho, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(ax, 1, '--','Color','k','LineWidth',1.0);
xline(ax, Tstar, ':','Color',[.5 .5 .5],'LineWidth',1.4);
text(ax, Tstar, 1.9, ' $\bar T^\star$','Interpreter','latex','Color',[.5 .5 .5],'FontSize',10);
xlabel(ax,'dwell-time $\bar T$','Interpreter','latex');
ylabel(ax,'cycle radius $\rho(e^{A_2\bar T}e^{A_1\bar T})$','Interpreter','latex');
title(ax,'minimum dwell-time threshold','Interpreter','latex','FontWeight','normal');
set(ax,'FontSize',10); ylim(ax,[0 2.4]); xlim(ax,[0 3]);
mvst_savepdf(gcf, fullfile(figdir,'pos_switched.pdf'));

% ---- data macros ----
try
    w = mvst_texdata(fullfile(figdir,'data','pos_switched.tex'));
    w.num('PosSwTstar',Tstar); w.num('PosSwTlp',Tlp);
    w.num('PosSwPone',p1(1)); w.num('PosSwPtwo',p1(2));
    w.num('PosSwPtre',p2(1)); w.num('PosSwPfou',p2(2));
    w.close();
catch, end
end

% =========================================================================
function [t,x,sw] = sw_traj(A1, A2, Tbar, Ncyc, x0)
%SW_TRAJ  Switched trajectory: alternate modes, each active for Tbar.
t=[]; x=[]; xk=x0; t0=0; npt=60; sw=struct('t0',{},'t1',{},'mode',{});
As={A1,A2};
for c=1:Ncyc
    m = mod(c-1,2)+1;                         % alternate 1,2,1,2,...
    tt=linspace(0,Tbar,npt);
    xx=cell2mat(arrayfun(@(s) expm(As{m}*s)*xk, tt,'uni',0));
    t=[t t0+tt]; x=[x xx];                     %#ok<AGROW>
    sw(end+1)=struct('t0',t0,'t1',t0+Tbar,'mode',m); %#ok<AGROW>
    xk = expm(As{m}*Tbar)*xk;                  % state at switch
    t0=t0+Tbar;
end
end
