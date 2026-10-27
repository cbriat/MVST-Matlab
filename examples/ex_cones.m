function ex_cones()
%EX_CONES  Cone-invariant linear systems and linear copositive Lyapunov
%  functions, for the chapter "Generalization to Proper Cones".
%
%  Three cases illustrate the reach and the limits of the conic viewpoint.
%
%  (1) Orthant cone K = R^n_+ (positive systems). A is Metzler (off-diagonal
%      >= 0) and Hurwitz, so the flow keeps R^n_+ invariant, and a LINEAR
%      copositive Lyapunov function V(x)=v'x with v>0, A'v<0 certifies
%      stability. Because A is Metzler-Hurwitz, -A^{-1} >= 0, so v=-(A')^{-1}1
%      solves the copositive LP in closed form.
%
%  (2) A non-orthant polyhedral cone. A is NOT Metzler in the standard basis,
%      yet the simplicial cone Omega=cone(V) is flow-invariant. The result
%      used is: if (A + tau I) Omega subset Omega for some tau, then Omega is
%      invariant, since M Omega subset Omega implies (I+sM)Omega subset Omega
%      implies e^{sM}Omega subset Omega, and e^{At}=e^{-tau t}e^{(A+tau I)t}.
%      In the generator basis, (A+tau I)Omega subset Omega reads
%      V^{-1}(A+tau I)V >= 0, i.e. V^{-1}A V is Metzler; a copositive
%      Lyapunov function on Omega follows exactly as in case (1).
%
%  (3) A rotation. A has purely imaginary eigenvalues; its orbit sweeps every
%      direction, so NO proper (pointed) cone is invariant (Vandergraft):
%      cone invariance needs a real dominant eigenstructure. Stability of the
%      damped focus is thus NOT accessible to a copositive certificate,
%      delimiting the method.
%
%  All trajectories are propagated by the matrix exponential (no solver).

here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ================= (1) orthant: Metzler + copositive LP =================
A1 = [-2 1; 0.5 -3];                       % Metzler (off-diag >= 0), Hurwitz
assert(all(A1(~eye(2))>=0),'A1 not Metzler');
v = -(A1.')\[1;1];                          % solves A1' v = -1  (copositive LP)
fprintf('(1) orthant: v = [%.3f %.3f] (>0: %d),  max(A1''v) = %.3f (<0)\n', ...
        v(1), v(2), all(v>0), max(A1.'*v));

figure('Color','w','Units','centimeters','Position',[2 2 13 7]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
X0s = [3 0.4; 2 2; 0.5 3; 1 1.5].';         % initial conditions in the orthant
tt = linspace(0,4,200);
for c = 1:size(X0s,2)
    xx = cell2mat(arrayfun(@(s) expm(A1*s)*X0s(:,c), tt,'uni',0));
    plot(ax, xx(1,:), xx(2,:), '-','LineWidth',1.6,'Color',gc.order(1,:));
    plot(ax, X0s(1,c), X0s(2,c), 'o','MarkerSize',4, ...
         'MarkerFaceColor',gc.order(1,:),'MarkerEdgeColor','none');
end
plot(ax, 0,0,'k.','MarkerSize',10);
xlabel(ax,'$x_1$','Interpreter','latex'); ylabel(ax,'$x_2$','Interpreter','latex');
title(ax,'orthant $\mathbb{R}^2_+$ invariant','Interpreter','latex','FontWeight','normal');
axis(ax,[0 3.2 0 3.2]); set(ax,'FontSize',10);

ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
for c = 1:size(X0s,2)
    Vt = arrayfun(@(s) v.'*expm(A1*s)*X0s(:,c), tt);
    plot(ax, tt, Vt, '-','LineWidth',1.6,'Color',gc.order(1,:));
end
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'$V(x)=v^\top x$','Interpreter','latex');
title(ax,'linear copositive Lyapunov function','Interpreter','latex','FontWeight','normal');
set(ax,'FontSize',10);
mvst_savepdf(gcf, fullfile(figdir,'cones_orthant.pdf'));

% ============ (2) non-Metzler A, invariant polyhedral cone ============
V  = [1 1; 0 1];                            % generators (1,0) and (1,1) of Omega
D  = diag([-1 -2]);                          % Hurwitz spectrum (real)
A2 = V*D/V;                                  % A2 = [-1 -1; 0 -2], NON-Metzler
Atil = V\A2*V;                               % = D, Metzler (diagonal) in gen. basis
tau = -min(diag(D)) + 0.5;                    % (A2+tau I) maps Omega into Omega
isMetzlerStd = all(A2(~eye(2))>=0);
fprintf('(2) polyhedral: A2 Metzler in std basis? %d ; V^{-1}A2V Metzler? %d ; (Atil+tauI)>=0? %d\n', ...
        isMetzlerStd, all(Atil(~eye(2))>=-1e-9), all(all(Atil+tau*eye(2) >= -1e-9)));
w = -(D.')\[1;1];                             % copositive Lyapunov weights on Omega

figure('Color','w','Units','centimeters','Position',[2 2 8.5 8]);
ax = axes; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
% shade the simplicial cone Omega = cone(V(:,1), V(:,2))
R = 4; g1 = V(:,1)/norm(V(:,1)); g2 = V(:,2)/norm(V(:,2));
patch(ax, R*[0 g1(1) g2(1)], R*[0 g1(2) g2(2)], [0.85 0.85 0.85], ...
      'EdgeColor','none','FaceAlpha',0.6);
plot(ax, R*[0 g1(1)], R*[0 g1(2)], '-','Color',[.5 .5 .5],'LineWidth',1.0);
plot(ax, R*[0 g2(1)], R*[0 g2(2)], '-','Color',[.5 .5 .5],'LineWidth',1.0);
tt = linspace(0,4,200);
X0c = [V*[2.6;0.3], V*[0.3;2.4], V*[1.4;1.4]];   % start inside Omega
for c = 1:size(X0c,2)
    xx = cell2mat(arrayfun(@(s) expm(A2*s)*X0c(:,c), tt,'uni',0));
    plot(ax, xx(1,:), xx(2,:), '-','LineWidth',1.6,'Color',gc.order(1,:));
    plot(ax, X0c(1,c), X0c(2,c), 'o','MarkerSize',4, ...
         'MarkerFaceColor',gc.order(1,:),'MarkerEdgeColor','none');
end
plot(ax, 0,0,'k.','MarkerSize',10);
xlabel(ax,'$x_1$','Interpreter','latex'); ylabel(ax,'$x_2$','Interpreter','latex');
title(ax,'invariant cone $\Omega=\mathrm{cone}(V)$ (non-Metzler $A$)', ...
      'Interpreter','latex','FontWeight','normal');
axis(ax,'equal'); axis(ax,[-0.4 3.6 -0.4 3]); set(ax,'FontSize',10);
mvst_savepdf(gcf, fullfile(figdir,'cones_polyhedral.pdf'));

% ==================== (3) rotation: no invariant cone ====================
Arot = [0 1; -1 0];                           % eigenvalues +-i (pure rotation)
Afoc = [-0.15 1; -1 -0.15];                   % stable focus (eigenvalues -0.15+-i)
tt = linspace(0,2*pi,400); ts = linspace(0,10,600);
orb = cell2mat(arrayfun(@(s) expm(Arot*s)*[1;0], tt,'uni',0));
spi = cell2mat(arrayfun(@(s) expm(Afoc*s)*[1;0], ts,'uni',0));
fprintf('(3) rotation eig(Arot) = %s ; no proper invariant cone.\n', ...
        mat2str(round(eig(Arot),3)));
figure('Color','w','Units','centimeters','Position',[2 2 8.5 8]);
ax = axes; hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h1 = plot(ax, orb(1,:), orb(2,:), '-','LineWidth',1.8,'Color',gc.order(1,:));
h2 = plot(ax, spi(1,:), spi(2,:), '-','LineWidth',1.6,'Color',gc.order(2,:));
plot(ax, 0,0,'k.','MarkerSize',10);
xlabel(ax,'$x_1$','Interpreter','latex'); ylabel(ax,'$x_2$','Interpreter','latex');
title(ax,'rotation: no proper invariant cone','Interpreter','latex','FontWeight','normal');
lg = legend(ax,[h1 h2],{'orbit ($\pm i$)','stable focus'}, ...
      'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off'); axis(ax,'equal'); axis(ax,[-1.4 1.4 -1.4 1.4]); set(ax,'FontSize',10);
mvst_savepdf(gcf, fullfile(figdir,'cones_rotation.pdf'));

% ============ native combined 2-up: polyhedral (a) + rotation (b) =========
% Single-figure merge replacing the book's minipage tiling of the two phase
% portraits that delimit the cone-invariant class.
figc = figure('Color','w','Units','centimeters','Position',[2 2 17 8]);
tlc = tiledlayout(figc,1,2,'TileSpacing','compact','Padding','compact');

% (a) invariant non-orthant polyhedral cone (non-Metzler A)
ax = nexttile(tlc); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
R = 4; g1 = V(:,1)/norm(V(:,1)); g2 = V(:,2)/norm(V(:,2));
patch(ax, R*[0 g1(1) g2(1)], R*[0 g1(2) g2(2)], [0.85 0.85 0.85], ...
      'EdgeColor','none','FaceAlpha',0.6);
plot(ax, R*[0 g1(1)], R*[0 g1(2)], '-','Color',[.5 .5 .5],'LineWidth',1.0);
plot(ax, R*[0 g2(1)], R*[0 g2(2)], '-','Color',[.5 .5 .5],'LineWidth',1.0);
ttc = linspace(0,4,200);
X0c = [V*[2.6;0.3], V*[0.3;2.4], V*[1.4;1.4]];
for c = 1:size(X0c,2)
    xx = cell2mat(arrayfun(@(s) expm(A2*s)*X0c(:,c), ttc,'uni',0));
    plot(ax, xx(1,:), xx(2,:), '-','LineWidth',1.6,'Color',gc.order(1,:));
    plot(ax, X0c(1,c), X0c(2,c), 'o','MarkerSize',4, ...
         'MarkerFaceColor',gc.order(1,:),'MarkerEdgeColor','none');
end
plot(ax, 0,0,'k.','MarkerSize',10);
xlabel(ax,'$x_1$','Interpreter','latex'); ylabel(ax,'$x_2$','Interpreter','latex');
title(ax,'(a) invariant cone $\Omega=\mathrm{cone}(V)$','Interpreter','latex','FontWeight','normal');
axis(ax,'equal'); axis(ax,[-0.4 3.6 -0.4 3]); set(ax,'FontSize',10);

% (b) rotation: no invariant cone
ax = nexttile(tlc); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
h1 = plot(ax, orb(1,:), orb(2,:), '-','LineWidth',1.8,'Color',gc.order(1,:));
h2 = plot(ax, spi(1,:), spi(2,:), '-','LineWidth',1.6,'Color',gc.order(2,:));
plot(ax, 0,0,'k.','MarkerSize',10);
xlabel(ax,'$x_1$','Interpreter','latex'); ylabel(ax,'$x_2$','Interpreter','latex');
title(ax,'(b) rotation: no invariant cone','Interpreter','latex','FontWeight','normal');
lg = legend(ax,[h1 h2],{'orbit ($\pm i$)','stable focus'}, ...
      'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off'); axis(ax,'equal'); axis(ax,[-1.4 1.4 -1.4 1.4]); set(ax,'FontSize',10);
mvst_savepdf(figc, fullfile(figdir,'cones_polyhedral_rotation.pdf'));

% ---- data macros ----
try
    wt = mvst_texdata(fullfile(figdir,'data','cones.tex'));
    wt.num('ConeOrthVone',v(1)); wt.num('ConeOrthVtwo',v(2));
    wt.num('ConeTau',tau); wt.mat('ConeA2',A2); wt.mat('ConeV',V);
    wt.close();
catch, end
end
