%% EX_STEER_POISSON  Finite-time covariance steering with Poissonian jumps.
%  Steers the state covariance from X0 to X1f over [0,T] for the GBM +
%  Poissonian-jump system: backward pass = generalized differential Riccati
%  from P(T)=M, forward pass = closed-loop moment equation with the
%  diffusion and jump-rate terms, multiplier M found by Newton.
%  Thin caller: the per-class gap closure and forward propagation now live in
%  the generic wrapper MVST_STEER (class 'poisson').
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();

A  = [0 1; 2 -0.5];  B  = [0;1];
C1 = 0.2*eye(2);     D1 = [0;0];
J1 = 1.25*eye(2);    lam = 0.5;
Q  = eye(2);         R  = 1;
T = 2.0;
X0 = [2 0.3; 0.3 1]; Xf = [0.8 0; 0 0.8];

sys = struct('A',A,'B',B,'Cc',{{C1}},'Dc',{{D1}},'Jc',{{J1}}, ...
             'lam',lam,'Q',Q,'R',R);
% The isotropic default guess M0=0 leaves the damped Newton stuck at the free
% evolution (gap ~5): the jump-inflated covariance needs an aggressive control
% whose multiplier sits in a basin separated from M0=0 by a Riccati-escape
% barrier.  Seeding M0=I lands inside that basin and Newton converges.
res = mvst_steer(sys, Xf, T, struct('class','poisson','X0',X0,'M0',eye(2)));
M = res.M; it = res.iters; tg = res.t; Sig = res.Sigma;
fprintf('Poisson steering: %d Newton iterations, terminal gap %.2e\n', it, res.res);

w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_steer_poisson.tex'));
w.mat('XzeroPoi', X0); w.mat('XfPoi', Xf); w.mat('MstarPoi', M);
w.num('TPoiS', T); w.num('ItPoiS', it, 0); w.close();

mvst_steer_ellipse(Sig, X0, Xf, 'Poissonian-jump covariance steering', ...
   fullfile(mvst_figdir(),'ex_steer_poisson.pdf'));
mvst_steer_components(tg, Sig, 'Poissonian-jump covariance components', ...
   fullfile(mvst_figdir(),'ex_steer_poisson_comp.pdf'));
