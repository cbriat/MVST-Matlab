function ex_cones_soc()
%EX_CONES_SOC  Stability on the second-order (Lorentz, ice cream) cone L^3,
%  for the chapter "Generalization to Proper Cones" (case (iv) of the
%  numerical illustrations).
%
%  AXIS-FIRST convention (matches Ch14 and the corrected journal paper):
%  x = (x_0, x_1, ..., x_{n-1}) with x_0 the cone axis, and
%  L^n = { x in R^n : ||(x_1,...,x_{n-1})||_2 <= x_0 }, written here as
%  L^n = { x : x'Q x >= 0, x_0 >= 0 } with the indefinite form
%  Q = diag(1, -I_{n-1}). It is proper and SELF-DUAL, and neither
%  polyhedral (case (ii)) nor a copy of the PSD cone. Its stability
%  certificate is therefore a SECOND-ORDER-CONE PROGRAM, not an LP
%  (orthant) nor an SDP (PSD cone).
%
%  In MATLAB's 1-based arrays the axis is x(1)=x_0 and the space part is
%  x(2:3)=(x_1,x_2).
%
%  System.  dot x = A x with A = -alpha I + W, where W is in the Lie
%  algebra so(2,1) = { W : Q W + W'Q = 0 } of the Lorentz group, so that
%  e^{tW} preserves the form x'Q x and hence L^3. W combines an x1-x2
%  rotation with a genuine x0-x2 boost, so the flow is a true Lorentz
%  transformation (not a spatial rotation) and, with the radial factor
%  e^{-alpha t}, interior trajectories move to the origin along shrinking
%  hyperboloids x'Q x = const > 0, while boundary starts stay on the
%  surface of the ice cream cone.
%
%  Verifications performed:
%   (a) W in so(2,1);  A Hurwitz (spectral abscissa < 0).
%   (b) Stern-Wolkowicz / Finsler cross-positivity condition
%       exists nu:  Q A + A'Q - nu Q  >=  0  (POSITIVE semidefinite).
%       Here Q A + A'Q = -2 alpha Q, so it holds with equality at nu=-2alpha.
%   (c) e^{tA} maps a fine sampling of the boundary of L^3 into L^3 for a
%       range of t >= 0 (direct invariance check).
%   (d) SOCP dual test (Prop. 13.16): p in int L^3 with -A'p in int L^3,
%       closed-form certificate p = e_0 (verified).
%
%  All trajectories are propagated by the matrix exponential (no solver).

here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% =================== system on the Lorentz cone L^3 ====================
Q = diag([1 -1 -1]);                        % L^3 = { x'Qx>=0, x0>=0 }, x0 axis
a = 1.0; b = 0.5; alpha = 1.5;              % rotation, boost, radial rate
W = [0 0 b; 0 0 -1; b 1 0];                 % element of so(2,1), axis-first
A = -alpha*eye(3) + W;
assert(norm(Q*W + W.'*Q,'fro') < 1e-12, 'W not in so(2,1)');
absc = max(real(eig(A)));
fprintf('(iv) Lorentz: eig(A) = %s\n', mat2str(round(eig(A),4)));
fprintf('     spectral abscissa = %.4f (Hurwitz: %d)\n', absc, absc<0);

% ---- (b) Stern-Wolkowicz / Finsler inequality  QA+A'Q - nu Q >= 0 ----
S  = Q*A + A.'*Q;                            % = -2 alpha Q for this A
nu = -2*alpha;                              % candidate multiplier
eigSW = sort(real(eig(S - nu*Q)));
fprintf('     QA+A''Q = -2*alpha*Q ? %d ; nu=%.2f ; min eig(S-nuQ)=%.2e (>=0)\n', ...
        norm(S + 2*alpha*Q,'fro')<1e-12, nu, min(eigSW));

% ---- (c) direct invariance: e^{tA} maps boundary of L^3 into L^3 ----
th = linspace(0,2*pi,721);
bd = [ones(size(th)); cos(th); sin(th)];    % boundary rays (x0 = 1, axis first)
tv = linspace(0,30,300);
worst = inf; worstQAx = inf;
for t = tv
    Y = expm(A*t)*bd;
    worst = min(worst, min(Y(1,:) - hypot(Y(2,:),Y(3,:))));
end
for k = 1:size(bd,2)
    x = bd(:,k); worstQAx = min(worstQAx, x.'*Q*A*x);
end
fprintf('     min invariance margin over dL^3, t in[0,30] = %.2e (>=0)\n', worst);
fprintf('     min x''QAx on dL^3 = %.2e (>=0, cross-positive)\n', worstQAx);

% ---- (d) SOCP dual test: p in int L^3, -A'p in int L^3 (Prop 13.16) ----
p = solve_socp(A);                          % max-margin certificate, p0=1
g = -A.'*p;
mp = p(1) - hypot(p(2),p(3));
mg = g(1) - hypot(g(2),g(3));
fprintf('     SOCP cert p = [%.3f %.3f %.3f], margin %.3f (in intL^3: %d)\n', ...
        p, mp, mp>0);
fprintf('     -A''p       = [%.3f %.3f %.3f], margin %.3f (in intL^3: %d)\n', ...
        g, mg, mg>0);

% ========================= trajectories in L^3 =========================
tt  = linspace(0,7,400);
% interior initial conditions strictly inside L^3 (space radius < height)
phis = [0 2*pi/5 4*pi/5 6*pi/5 8*pi/5];
hgt  = 2.4; rad = 1.6;                       % axis height, space radius
X0int = [hgt*ones(size(phis)); rad*cos(phis); rad*sin(phis)];
% boundary initial conditions on the surface (space radius = height)
phib = [pi/5 pi];
rb   = 2.4;
X0bd = [rb*ones(size(phib)); rb*cos(phib); rb*sin(phib)];

figure('Color','w','Units','centimeters','Position',[2 2 15 7]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ---- left panel: 3D ice cream cone with trajectories (x0 vertical) ----
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'ColorOrder',gc.order);                       % gem palette, reset
Hc = hgt*1.15;
[thc,rc] = meshgrid(linspace(0,2*pi,72), linspace(0,Hc,2));
% cone surface x0 = ||(x1,x2)||: plot (x1,x2) in the horizontal plane, x0 up
surf(ax, rc.*cos(thc), rc.*sin(thc), rc, 'FaceColor',[0.82 0.82 0.85], ...
     'FaceAlpha',0.35,'EdgeColor','none');
% a few rulings of the cone boundary for readability
for pth = linspace(0,2*pi,13)
    plot3(ax,[0 Hc*cos(pth)],[0 Hc*sin(pth)],[0 Hc], ...
          '-','Color',[0.6 0.6 0.65],'LineWidth',0.4);
end
% interior trajectories (move to origin along shrinking hyperboloids)
for c = 1:size(X0int,2)
    xx = cell2mat(arrayfun(@(s) expm(A*s)*X0int(:,c), tt,'uni',0));
    plot3(ax, xx(2,:), xx(3,:), xx(1,:), '-','LineWidth',1.5, ...
          'Color',gc.order(1,:));
    plot3(ax, X0int(2,c),X0int(3,c),X0int(1,c),'o','MarkerSize',4, ...
          'MarkerFaceColor',gc.order(1,:),'MarkerEdgeColor','none');
end
% boundary trajectories (stay on the surface, spiral down to the origin)
for c = 1:size(X0bd,2)
    xx = cell2mat(arrayfun(@(s) expm(A*s)*X0bd(:,c), tt,'uni',0));
    plot3(ax, xx(2,:), xx(3,:), xx(1,:), '-','LineWidth',1.5, ...
          'Color',gc.order(2,:));
    plot3(ax, X0bd(2,c),X0bd(3,c),X0bd(1,c),'s','MarkerSize',4, ...
          'MarkerFaceColor',gc.order(2,:),'MarkerEdgeColor','none');
end
plot3(ax, 0,0,0,'k.','MarkerSize',10);
xlabel(ax,'$x_1$','Interpreter','latex');
ylabel(ax,'$x_2$','Interpreter','latex');
zlabel(ax,'$x_0$','Interpreter','latex');
title(ax,'ice cream cone $\mathcal{L}^3$ invariant', ...
      'Interpreter','latex','FontWeight','normal');
view(ax,[-37 22]); axis(ax,'tight'); set(ax,'FontSize',10);

% ---- right panel: height x_0 = <p,x(t)> decreasing along trajectories ----
ax = nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
set(ax,'ColorOrder',gc.order);                       % gem palette, reset
for c = 1:size(X0int,2)
    Vt = arrayfun(@(s) p.'*expm(A*s)*X0int(:,c), tt);
    plot(ax, tt, Vt, '-','LineWidth',1.6,'Color',gc.order(1,:));
end
for c = 1:size(X0bd,2)
    Vt = arrayfun(@(s) p.'*expm(A*s)*X0bd(:,c), tt);
    plot(ax, tt, Vt, '-','LineWidth',1.6,'Color',gc.order(2,:));
end
xlabel(ax,'time $t$','Interpreter','latex');
ylabel(ax,'$V(x)=\langle p,x\rangle=x_0$','Interpreter','latex');
title(ax,'linear copositive Lyapunov function (SOCP certificate)', ...
      'Interpreter','latex','FontWeight','normal');
set(ax,'FontSize',10);
mvst_savepdf(gcf, fullfile(figdir,'cones_lorentz.pdf'));

% --------------------------- data macros ---------------------------
try
    wt = mvst_texdata(fullfile(figdir,'data','cones_soc.tex'));
    wt.mat('ConeSocA', A);
    wt.mat('ConeSocW', W);
    wt.num('ConeSocAlpha', alpha);
    wt.num('ConeSocNu', nu);
    wt.num('ConeSocAbsc', absc);
    wt.mat('ConeSocP', p);
    wt.mat('ConeSocAdp', g);
    wt.num('ConeSocMargP', mp);
    wt.num('ConeSocMargG', mg);
    wt.close();
catch, end
end

% ======================================================================
function p = solve_socp(A)
%SOLVE_SOCP  Max-margin certificate p in int L^3 with -A'p in int L^3.
%   max t  s.t.  ||p_{2:3}|| <= p_1 - t,  ||(-A'p)_{2:3}|| <= (-A'p)_1 - t,
%                p_1 = 1.   Closed-form solution p = e_0 (verified).
%   (axis-first: the axis coordinate is index 1, the space part indices 2:3.)
% closed-form certificate: the gauge e_0 is strictly interior and, since
% -A'e_0 = (3/2,0,-1/2) in int L^3 for this A, is a valid certificate.
p = [1;0;0];
g = -A.'*p;
assert(g(1) > norm(g(2:3)), 'ex_cones_soc:socp', 'p = e_0 is not a valid certificate for this A.');
p = p / p(1);                                % normalize p_0 = 1
end
