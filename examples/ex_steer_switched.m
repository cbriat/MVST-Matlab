%% EX_STEER_SWITCHED  Finite-time covariance steering of a switched system.
%  Steers the state covariance from X0 to X1f over [0,T] under a dwell-time
%  switching schedule.  Backward pass: the mode-dependent I-GRE integrated
%  backward through the schedule from P(T)=M (the state is continuous at
%  switches, so P is too); forward pass: the switched closed-loop moment
%  equation; the multiplier M is found by Newton.
%  Thin caller: the per-class gap closure and forward propagation now live in
%  the generic wrapper MVST_STEER (class 'switched').
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data();
Td = D.sw.Tdwell; T = 4.0;
X0 = [3 0;0 0.5]; Xf = [1.0 0;0 1.6];

res = mvst_steer(D.sw, Xf, T, struct('class','switched','X0',X0));
M = res.M; it = res.iters; tg = res.t; Sig = res.Sigma;
fprintf('Switched steering: %d Newton iterations, terminal gap %.2e\n', it, res.res);

w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_steer_switched.tex'));
w.mat('XzeroSw', X0); w.mat('XfSw', Xf); w.mat('MstarSw', M); w.num('TSw', T);
w.num('TdwellSwS', Td); w.close();

mvst_steer_ellipse(Sig, X0, Xf, 'Switched covariance steering', ...
   fullfile(mvst_figdir(),'ex_steer_switched.pdf'));
sw = Td:Td:(T-1e-9);                                  % switching instants
mvst_steer_components(tg, Sig, 'Switched covariance components', ...
   fullfile(mvst_figdir(),'ex_steer_switched_comp.pdf'), sw);
