%% EX_STEER_SD  Finite-time covariance steering by sampled-data control.
%  Steers the state covariance from X0 to X1f over [0,T] using a zero-order-
%  hold (sampled) feedback.  Backward pass: the impulsive I-GRE on the
%  augmented state z=[x;u] (control-free flow + jump Riccati at each sample,
%  yielding the sample gains K_k) from terminal P(T)=blkdiag(M,0); forward
%  pass: the augmented covariance flowing and being reset by sample-and-hold;
%  the multiplier M is found by Newton.
%  Thin caller: the per-class gap closure and forward propagation now live in
%  the generic wrapper MVST_STEER (class 'sd').
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data();
s = D.sd; T = 3.0;
X0 = [3 0;0 0.5]; Xf = [1.0 0;0 1.6];

% The isotropic default guess M0=0 strands the damped Newton at a terminal gap
% of order 1 (the far trial step escapes the backward I-GRE, then the line
% search stalls).  Seeding it in the solution basin restores convergence.
res = mvst_steer(s, Xf, T, struct('class','sd','X0',X0,'M0',-eye(2)));
M = res.M; it = res.iters; tg = res.t; Sig = res.Sigma;
fprintf('Sampled-data steering: %d Newton iterations, terminal gap %.2e\n', it, res.res);

w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_steer_sd.tex'));
w.mat('XzeroSd', X0); w.mat('XfSd', Xf); w.mat('MstarSd', M);
w.num('TSd', T); w.num('TbarSdS', s.Tbar); w.close();

mvst_steer_ellipse(Sig, X0, Xf, 'Sampled-data covariance steering', ...
   fullfile(mvst_figdir(),'ex_steer_sd.pdf'));
tk = s.Tbar:s.Tbar:(T-1e-9);                          % sampling instants
mvst_steer_components(tg, Sig, 'Sampled-data covariance components', ...
   fullfile(mvst_figdir(),'ex_steer_sd_comp.pdf'), tk);
