%% EX_TV_STEER_SD  Time-varying covariance steering by sampled-data control.
%  Time-varying counterpart of EX_STEER_SD: the flow drift decays in time,
%  A(t) = A + dA*exp(-t), so the system is genuinely time-varying but settles
%  to the LTI sampled-data system as t grows.  Steers the state covariance
%  from X0 to Xf over [0,T] by zero-order-hold (sampled) feedback.  Backward
%  pass: the impulsive I-GRE on the augmented state z=[x;u] with TIME-VARYING
%  flow Ab0(t)=[A(t) B; 0 0]; forward pass: the augmented covariance flowing
%  with the SAME time-varying Ab0(t) and being reset by sample-and-hold; the
%  multiplier M is found by Newton.  Thin caller: the per-class gap closure and
%  forward propagation now live in the generic wrapper MVST_STEER (class 'sd',
%  time-varying flow supplied as a function handle in SYS.A).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data();
s = D.sd; T = 3.0;
X0 = [3 0;0 0.5]; Xf = [1.0 0;0 1.6];
dA = [0 0; 1.5 0];                          % decaying drift perturbation
s.A = @(t) s.A + dA*exp(-t);                % time-varying flow drift

res = mvst_steer(s, Xf, T, struct('class','sd','X0',X0));
tg = res.t; Sig = res.Sigma;
fprintf('TV sampled-data steering: %d Newton iterations, terminal gap %.2e\n', res.iters, res.res);

w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_steer_sd.tex'));
w.mat('TVXzeroSd', X0); w.mat('TVXfSd', Xf); w.num('TVTSd', T);
w.mat('TVdAsd', dA); w.close();

mvst_steer_ellipse(Sig, X0, Xf, 'Time-varying sampled-data covariance steering', ...
   fullfile(mvst_figdir(),'ex_tv_steer_sd.pdf'));
tk = s.Tbar:s.Tbar:(T-1e-9);                          % sampling instants
mvst_steer_components(tg, Sig, 'Time-varying sampled-data covariance components', ...
   fullfile(mvst_figdir(),'ex_tv_steer_sd_comp.pdf'), tk);
