function ex_pos()
%EX_POS  Linear positive impulsive system under a maximum dwell-time
%  (Chapter "Linear Positive Systems"). This is the running impulsive
%  example of the paper "Dwell-time stability and stabilization conditions
%  for linear positive impulsive and switched systems" (its Matlab code
%  Maximum_DT_stability.m):
%     A = [0.5 1; 0 0.5]   (Metzler, UNSTABLE, eig +0.5: the flow expands),
%     J = [0.1 0.2; 0 0.1] (>=0, rho(J)=0.1<1: the reset contracts).
%  The constant-dwell-time monodromy J e^{A Tbar} is Schur iff the dwell-time
%  is below a MAXIMUM value Tbar* (the contractive resets must catch the
%  expanding flow). Below Tbar* the paper's linear copositive certificate
%     p>0,   A' p > 0,   (J e^{A Tbar})' p < p   (a linear program in p)
%  is feasible; positivity (x(t) >= 0 for x(0) >= 0) holds throughout.
%
%  Figures: pos_dwell.pdf, the positive trajectory below the threshold next
%  to the monodromy radius versus the dwell-time.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

A = [0.5 1; 0 0.5];               % Metzler, UNSTABLE (eig +0.5, expanding flow)
J = [0.1 0.2; 0 0.1];             % >=0, rho(J)=0.1<1 (contractive reset)
assert(all(A(~eye(2))>=0) && all(J(:)>=0));

% ---- maximum dwell-time: rho(J e^{A Tbar}) crosses 1 from below ----
Tg = linspace(0.05, 8, 400);
rho = arrayfun(@(T) max(abs(eig(J*expm(A*T)))), Tg);
i0 = find(rho >= 1, 1, 'first');
lo = Tg(max(i0-1,1)); hi = Tg(i0);
for it=1:40, m=0.5*(lo+hi); if max(abs(eig(J*expm(A*m))))<1, lo=m; else, hi=m; end; end
Tstar = lo;
fprintf('maximum dwell-time threshold Tbar* = %.4f (rho=1)\n', Tstar);

% ---- copositive certificate (LP) at a stable dwell-time Tbar < Tbar* ----
Tstab = Tstar*0.6;
rc = mvst_pos_stability(struct('A',A,'J',J,'Tbar',Tstab), ...
                        struct('class','impulsive'));
p = rc.p;
fprintf('copositive certificate at Tbar=%.3f: p=[%.3f %.3f]>0, A''p=[%.3f %.3f]>0, (J''e^{A''T}p-p)=[%.3f %.3f]<0\n', ...
        Tstab, p(1),p(2), A.'*p, J.'*expm(A.'*Tstab)*p - p);

% ---- simulate the impulsive trajectory below the threshold (stable, positive) ----
[ts, xs] = imp_traj(A, J, Tstab, 6, [1;0.2]);

% ================= figure =================
figure('Color','w','Units','centimeters','Position',[2 2 13 13]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h1 = plot(ax, ts, xs(1,:), '-','LineWidth',1.8,'Color',gc.order(1,:));
h2 = plot(ax, ts, xs(2,:), '-','LineWidth',1.8,'Color',gc.order(2,:));
xlabel(ax,'time $t$','Interpreter','latex'); ylabel(ax,'state $x(t)\ge0$','Interpreter','latex');
title(ax,sprintf('$\\bar T=%.2f<\\bar T^\\star$ (stable, positive)',Tstab),'Interpreter','latex','FontWeight','normal');
lg=legend(ax,[h1 h2],{'$x_1$','$x_2$'},'Interpreter','latex','Location','northeast'); set(lg,'Box','off');
set(ax,'FontSize',10);

ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, Tg, rho, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(ax, 1, '--','Color','k','LineWidth',1.0);
xline(ax, Tstar, ':','Color',[.5 .5 .5],'LineWidth',1.4);
text(ax, Tstar, 1.18, ' $\bar T^\star$','Interpreter','latex','Color',[.5 .5 .5],'FontSize',10);
xlabel(ax,'dwell-time $\bar T$','Interpreter','latex');
ylabel(ax,'monodromy radius $\rho(J e^{A\bar T})$','Interpreter','latex');
title(ax,'maximum dwell-time threshold','Interpreter','latex','FontWeight','normal');
set(ax,'FontSize',10); ylim(ax,[0 2]);
mvst_savepdf(gcf, fullfile(figdir,'pos_dwell.pdf'));

% ================= L1 gain of a positive system (Example ex:pos:gain) =====
% Compartmental system with a disturbance into compartment 1 and output the
% total mass z=x1+x2. For a positive system the induced L1 gain equals the
% copositive LP value gamma*=lambda_1 with A'lambda=-1, which is also the DC
% gain C(-A)^{-1}E and the area under the (nonnegative) impulse response.
k1=0.4; k2=0.6; k12=0.5;
Ag=[-(k1+k12) 0; k12 -k2]; Eg=[1;0]; Cg=[1 1];
rg = mvst_pos_gain(struct('A',Ag,'E',Eg,'C',Cg), struct('norm','L1'));
gamLP = rg.gamma;                               % copositive LP value gamma*=lambda_1
gamDC = Cg*((-Ag)\Eg);                          % DC gain C(-A)^{-1}E
gamCF = (1/(k1+k12))*(1+k12/k2);                % closed form (chapter)
tt = linspace(0,16,400);
zstep = arrayfun(@(t) Cg*((-Ag)\(eye(2)-expm(Ag*t))*Eg), tt);  % response to unit step
zimp  = arrayfun(@(t) Cg*expm(Ag*t)*Eg, tt);                   % impulse response (>=0)
gamArea = trapz(tt, zimp);                       % area under impulse response = L1 gain
fprintf('L1 gain gamma*: LP=%.4f  DC=%.4f  closed-form=%.4f  impulse-area=%.4f\n', ...
        gamLP, gamDC, gamCF, gamArea);

figure('Color','w','Units','centimeters','Position',[2 2 13 13]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, tt, zstep, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(ax, gamLP, '--','Color','k','LineWidth',1.0);
text(ax, 9, gamLP, ' $\gamma^\star$','Interpreter','latex','VerticalAlignment','bottom','FontSize',11);
xlabel(ax,'time $t$','Interpreter','latex'); ylabel(ax,'output $z(t)$ to a unit step','Interpreter','latex');
title(ax,'step response $\to\gamma^\star$','Interpreter','latex','FontWeight','normal');
set(ax,'FontSize',10);
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
area(ax, tt, zimp, 'FaceColor',gc.order(2,:),'FaceAlpha',0.25,'EdgeColor','none');
plot(ax, tt, zimp, '-','LineWidth',2.0,'Color',gc.order(2,:));
xlabel(ax,'time $t$','Interpreter','latex'); ylabel(ax,'impulse response $Ce^{At}E\ge0$','Interpreter','latex');
title(ax,sprintf('area $=\\gamma^\\star=%.3f$',gamArea),'Interpreter','latex','FontWeight','normal');
set(ax,'FontSize',10); xlim(ax,[0 16]);
mvst_savepdf(gcf, fullfile(figdir,'pos_gain.pdf'));

% ---- data macros ----
try
    w = mvst_texdata(fullfile(figdir,'data','pos.tex'));
    w.num('PosTstar',Tstar); w.num('PosPone',p(1)); w.num('PosPtwo',p(2));
    w.num('PosGain',gamLP);
    w.close();
catch, end
end

% =========================================================================
function [t,x] = imp_traj(A, J, Tbar, Ncyc, x0)
%IMP_TRAJ  Impulsive trajectory: flow e^{At} between jumps, reset J at each.
t=[]; x=[]; xk=x0; t0=0; npt=40;
for c=1:Ncyc
    tt=linspace(0,Tbar,npt);
    xx=cell2mat(arrayfun(@(s) expm(A*s)*xk, tt,'uni',0));
    t=[t t0+tt]; x=[x xx];                 %#ok<AGROW>
    xk = J*(expm(A*Tbar)*xk);              % dose (jump) at end of cycle
    t0=t0+Tbar;
end
end
