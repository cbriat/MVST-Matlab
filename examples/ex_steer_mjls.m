%% EX_STEER_MJLS  Finite-time covariance steering of a Markov jump diffusion.
%  Steers the aggregate second moment Sigma = sum_i X_i from X0 to X1f over
%  [0,T] for the two-mode Markov-jump-with-diffusion system, via the
%  forward-backward algorithm with a common terminal multiplier M
%  (P_i(T)=M).  Backward pass: the coupled finite-horizon I-GRE (gains
%  K_i(t)); forward pass: the coupled moment equation; M is found by Newton.
%  Thin caller: the per-class gap closure and forward propagation now live in
%  the generic wrapper MVST_STEER (class 'mjls').
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data();
sys = D.mjls; T = 2.0; X0 = sys.modes{1}.X0; Xf = [1.0 0;0 1.6];

res = mvst_steer(sys, Xf, T, struct('class','mjls','X0',X0));
M = res.M; it = res.iters; tg = res.t; Sig = res.Sigma;
fprintf('MJLS steering: %d Newton iterations, terminal gap %.2e\n', it, res.res);

w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_steer_mjls.tex'));
w.mat('XzeroSm', X0); w.mat('XfSm', Xf); w.mat('MstarSm', M); w.num('TSm', T);
w.num('ItSm', it, 0); w.close();   % integer format (iteration count)

mvst_steer_ellipse(Sig, X0, Xf, 'Markov-jump covariance steering', ...
   fullfile(mvst_figdir(),'ex_steer_mjls.pdf'));
mvst_steer_components(tg, Sig, 'Markov-jump covariance components', ...
   fullfile(mvst_figdir(),'ex_steer_mjls_comp.pdf'));
