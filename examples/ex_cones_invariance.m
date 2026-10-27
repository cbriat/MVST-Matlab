function ex_cones_invariance()
%EX_CONES_INVARIANCE  Cone invariance on a genuinely NON-SIMPLICIAL polyhedral
%   cone (Chapter "Generalization to Proper Cones", Section sec:cones:LTI:num):
%   the Metzler generator factorization A V = V H and the copositive-Lyapunov LP
%   certificate, for a cone with more generators than the dimension (m>n).
%
%   The cone is the square pyramid cone in R^3
%     Omega = cone(V),  V = [v1 v2 v3 v4],  v_k the four unit-height corners,
%   which is proper (pointed, solid) but NON-simplicial (m=4 generators, n=3).
%   A dynamics preserving Omega is built from a Metzler generator H0 by
%   A = V H0 V^+ (V^+ the pseudoinverse); H0 cyclically contracts the generators,
%   A v_k = -v_k + 0.3 v_{k+1} (indices mod 4), which rotates the corners so that
%   A is NOT Metzler (the orthant is not invariant).  The 'cone' class of
%   MVST_POS_STABILITY recovers a Metzler H with A V = V H (invariance,
%   Vandergraft / Berman-Plemmons) and certifies Hurwitz stability on Omega by
%   the LP exists p>0 with H' p < 0 (the pulled-back certificate V(x)=<p,xi>).
%
%   It also computes the least shift tau with (A+tau I)Omega in Omega by the LP
%   (A+tau I)V = V N, N >= 0 (Prop. prop:cones:shift), and a certificate y of
%   the generator-coordinate dual test V'y > 0, V'A'y < 0 (Cor.
%   cor:cones:copLP(d)), whose image V'y is the certificate p of H'p < 0.
%
%   Outputs: figures/cones_nonsimplicial.pdf (the cone, its four generators, and
%   a trajectory started inside Omega staying inside and decaying), macros in
%   figures/data/cones_invariance.tex.
%
%   See also MVST_POS_STABILITY.

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---- non-simplicial square-pyramid cone in R^3 ----------------------------
v1=[1;1;1]; v2=[1;-1;1]; v3=[-1;-1;1]; v4=[-1;1;1];
V = [v1 v2 v3 v4];                              % 3 x 4 generators (m=4 > n=3)
P4 = [0 0 0 1;1 0 0 0;0 1 0 0;0 0 1 0];         % cyclic shift v_k -> v_{k+1}
H0 = -eye(4) + 0.3*P4;                           % Metzler, Hurwitz (cyclic)
A  = V*H0*pinv(V);                              % dynamics preserving Omega

r = mvst_pos_stability(struct('A',A,'V',V), struct('class','cone'));
mu = max(real(eig(A)));
Hp = r.H.'*r.p;                                 % must be < 0 entrywise

% ---- least shift tau with (A+tau I)Omega in Omega (Prop. prop:cones:shift) --
% Polyhedral case: (A+tau I_n)V = V N for some entrywise nonnegative N.  LP in
% (vec(N), tau):  min tau  s.t.  (I_m kron V) vec(N) - vec(V) tau = vec(A V),
% N >= 0, tau >= 0.
[n,m] = size(V);
lpo = optimoptions('linprog','Display','off');
fS  = [zeros(m*m,1); 1];
AeS = [kron(eye(m),V), -V(:)];
beS = reshape(A*V,[],1);
lbS = zeros(m*m+1,1);
[zS,~,flS] = linprog(fS,[],[],AeS,beS,lbS,[],lpo);
assert(flS==1,'least-shift LP failed');
tau = zS(end);

% ---- generator-coordinate dual test, Cor. cor:cones:copLP(d) ---------------
% exists y: V'y > 0, V'A'y < 0.  Homogeneous; normalize 1'V'y = m and
% minimize t with V'A'y <= t 1, V'y >= 0 (t < 0 and V'y > 0 certify).
fY  = [zeros(n,1); 1];
AiY = [ (A*V).', -ones(m,1); -V.', zeros(m,1) ];
biY = zeros(2*m,1);
AeY = [ones(1,m)*V.', 0];  beY = m;
[zY,~,flY] = linprog(fY,AiY,biY,AeY,beY,[],[],lpo);
assert(flY==1,'dual-test LP failed');
y = zY(1:n); y(abs(y)<1e-10) = 0;
assert(all(V.'*y > 0) && all(V.'*(A.'*y) < 0), 'condition (d) not certified');
fprintf('   least shift tau=%.6f ; y=[%.3f %.3f %.3f] ; V''y=[%s] ; V''A''y=[%s]\n', ...
        tau, y, sprintf('%.3f ',V.'*y), sprintf('%.3f ',V.'*A.'*y));

fprintf('[ex_cones_invariance] coneInvariant=%d stable=%d mu(A)=%.4f\n', ...
        r.coneInvariant, r.stable, mu);
fprintf('   ||A V - V H||=%.2e ; max(H''p)=%.4f ; p=[%.3f %.3f %.3f %.3f]\n', ...
        norm(A*V-V*r.H,'fro'), max(Hp), r.p);

% ---- figure: the cone, its generators, a decaying interior trajectory -----
figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
L = 3;                                           % ray length for drawing
hold on
% cone faces (four triangles from the origin along consecutive generators)
Vs = V./vecnorm(V);                              % unit generators for drawing
order = [1 2 3 4 1];
for k = 1:4
    a = L*Vs(:,order(k)); b = L*Vs(:,order(k+1));
    patch([0 a(1) b(1)],[0 a(2) b(2)],[0 a(3) b(3)], gc.order(1,:), ...
          'FaceAlpha',0.12,'EdgeColor','none');
end
for k = 1:4
    a = L*Vs(:,k);
    plot3([0 a(1)],[0 a(2)],[0 a(3)],'-','Color',gc.order(1,:),'LineWidth',1.6);
end
% a trajectory started strictly inside Omega
x0 = V*[0.6;0.9;0.5;0.8];
tt = linspace(0,6,400); Xtraj = zeros(3,numel(tt));
for i=1:numel(tt); Xtraj(:,i) = expm(A*tt(i))*x0; end
plot3(Xtraj(1,:),Xtraj(2,:),Xtraj(3,:),'-','Color',gc.order(2,:),'LineWidth',1.8);
plot3(x0(1),x0(2),x0(3),'o','Color',gc.order(2,:),'MarkerFaceColor',gc.order(2,:),'MarkerSize',5);
grid on; box on; view(35,22); axis equal
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');
zlabel('$x_3$','Interpreter','latex');
mvst_savepdf(gcf, fullfile(figdir,'cones_nonsimplicial.pdf'));

% ---- data macros ----------------------------------------------------------
try
    w = mvst_texdata(fullfile(figdir,'data','cones_invariance.tex'));
    w.mat('ConeInvA', A);  w.mat('ConeInvV', V);  w.mat('ConeInvH', r.H);
    w.mat('ConeInvP', r.p.');
    w.num('ConeInvMu', mu);
    w.raw('ConeInvN', sprintf('%d',size(V,1)));
    w.raw('ConeInvM', sprintf('%d',size(V,2)));
    w.num('ConeInvTau', tau);
    w.mat('ConeInvY', y);
    w.close();
catch, end
end
