function ex_pos_linf()
%EX_POS_LINF  Induced L-infinity (peak-to-peak) gain of a positive system
%  (Chapter "Linear Positive Systems"; the peak-gain counterpart of the L1
%  Example ex_pos, following the L-infinity x l-infinity theory of
%  "Hybrid Linf x linf-performance analysis and control of linear
%  time-varying impulsive and switched positive systems" and the sliced
%  conditions of Corollary cor:pos:Linfslice).
%
%  For a positive, stable LTI system  xdot = A x + E w,  z = C x  (+ F w),
%  with A Metzler-Hurwitz, E,C,F >= 0, the induced L-infinity -> L-infinity
%  (peak-to-peak) gain equals the DC gain of the nonnegative impulse
%  response,
%     gamma_inf = C(-A)^{-1} E ( + F ),
%  and is certified by the ROW-type copositive linear program (primal tube)
%     A xi + E 1 <= 0,   C xi + F 1 <= gamma 1,   xi >= 0,
%  whose least value is attained at xi* = (-A)^{-1} E 1 >= 0. Under a bounded
%  unit disturbance the state is dominated by the constant tube xi*, and
%  under the worst case (constant w = 1) it converges to xi*, so the output
%  reaches exactly gamma_inf: the peak bound is tight.
%
%  Figure: pos_linf.pdf, the output step response saturating at gamma_inf
%  next to the state components rising to (and bounded by) the tube xi*.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---- a 3-state positive (Metzler-Hurwitz) network ----
A = [-3 1 0; 1 -3 1; 0 1 -2];      % Metzler, Hurwitz (diagonally dominant)
E = [1;0;0];                       % disturbance into node 1
C = [1 1 1];                       % total-mass output
F = 0;
n = size(A,1);
assert(all(A(~eye(n))>=0) && all(eig(A)<0),'A not Metzler-Hurwitz');

% ---- three coincident evaluations of the peak-to-peak gain ----
rp      = mvst_pos_gain(struct('A',A,'E',E,'C',C,'F',F), struct('norm','Linf'));
xiStar  = rp.weights;              % tube xi* = (-A)^{-1} E 1  (>= 0)
gamRow  = rp.gammaDC;              % row-type copositive value (closed form C xi*+F)
gamDC   = C*((-A)\E) + F;          % DC gain C(-A)^{-1}E (+F)
gamLP   = rp.gamma;                % row-type LP solved numerically (linprog)
tt = linspace(0,8,500);
zstep = arrayfun(@(t) C*((-A)\(eye(n)-expm(A*t))*E) + F, tt);  % step response
Xstep = cell2mat(arrayfun(@(t) (-A)\(eye(n)-expm(A*t))*E, tt,'uni',0)); % states
gamStep = zstep(end);              % steady-state output under unit disturbance
fprintf('tube xi* = [%.4f %.4f %.4f] >= 0\n', xiStar(1),xiStar(2),xiStar(3));
fprintf('peak-to-peak gain gamma_inf: row-LP(closed)=%.4f  DC-gain=%.4f  YALMIP-LP=%.4f  step-ss=%.4f\n', ...
        gamRow, gamDC, gamLP, gamStep);

% ================= figure =================
figure('Color','w','Units','centimeters','Position',[2 2 13 13]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
plot(ax, tt, zstep, '-','LineWidth',2.0,'Color',gc.order(1,:));
yline(ax, gamRow, '--','Color','k','LineWidth',1.0);
text(ax, 4.5, gamRow, ' $\gamma_\infty$','Interpreter','latex','VerticalAlignment','bottom','FontSize',11);
xlabel(ax,'time $t$','Interpreter','latex'); ylabel(ax,'output $z(t)$ to a unit disturbance','Interpreter','latex');
title(ax,'peak-to-peak: step response $\to\gamma_\infty$','Interpreter','latex','FontWeight','normal');
set(ax,'FontSize',10); ylim(ax,[0 gamRow*1.15]);

ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
cols = [gc.order(1,:); gc.order(2,:); gc.order(3,:)];
h = gobjects(1,n);
for i=1:n
    h(i) = plot(ax, tt, Xstep(i,:), '-','LineWidth',1.8,'Color',cols(i,:));
    yline(ax, xiStar(i), ':','Color',cols(i,:),'LineWidth',1.2);
end
xlabel(ax,'time $t$','Interpreter','latex'); ylabel(ax,'states $x_i(t)\le\xi^\star_i$','Interpreter','latex');
title(ax,'primal tube $\xi^\star=(-A)^{-1}E\mathbf{1}$','Interpreter','latex','FontWeight','normal');
lg=legend(ax,h,{'$x_1$','$x_2$','$x_3$'},'Interpreter','latex','Location','east'); set(lg,'Box','off');
set(ax,'FontSize',10);
mvst_savepdf(gcf, fullfile(figdir,'pos_linf.pdf'));

% ---- data macros ----
try
    w = mvst_texdata(fullfile(figdir,'data','pos_linf.tex'));
    w.num('PosLinfGain',gamRow);
    w.num('PosLinfXione',xiStar(1)); w.num('PosLinfXitwo',xiStar(2));
    w.num('PosLinfXitre',xiStar(3));
    w.close();
catch, end
end
