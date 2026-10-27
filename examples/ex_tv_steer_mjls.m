%% EX_TV_STEER_MJLS  Time-varying covariance steering of a Markov jump diffusion.
%  Time-varying counterpart of EX_STEER_MJLS: each mode's flow drift decays in
%  time, A_i(t) = A_i + dA*exp(-t), so the system is genuinely time-varying but
%  settles to the LTI Markov-jump system as t grows.  We steer the aggregate
%  second moment Sigma = sum_i X_i from X0 to Xf over [0,T] via the
%  forward-backward algorithm with a common terminal multiplier M (P_i(T)=M).
%  Backward pass: the coupled TIME-VARYING I-GRE (gains K_i(t)); forward pass:
%  the coupled moment equation with the time-varying closed loop; M is found by
%  Newton.  Thin caller: the per-class gap closure and forward propagation now
%  live in the generic wrapper MVST_STEER (class 'mjls', time-varying path).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data();
base = D.mjls; n = base.n; m = base.m; Mn = base.M; T = 2.0;
X0 = base.modes{1}.X0; Xf = [1.0 0;0 1.6];
dA = [0 0; 1.5 0];                          % decaying drift perturbation

% ----- build the time-varying modes A_i(t) = A_i + dA*exp(-t) -----
md = cell(1,Mn);
for i = 1:Mn
    b = base.modes{i};
    md{i} = mvst_mode_tv(n, m, ...
        'A', @(t) b.A + dA*exp(-t), 'B', b.B, 'Q', b.Q, 'R', b.R, ...
        'Ac', b.Ac, 'Bc', b.Bc, 'W', b.W, 'X0', b.X0);
end
sys = mvst_system(md, 'Pi', base.Pi, 'pi0', base.pi0);

res = mvst_steer(sys, Xf, T, struct('class','mjls','X0',X0));
tg = res.t; Sig = res.Sigma;
fprintf('TV-MJLS steering: %d Newton iterations, terminal gap %.2e\n', res.iters, res.res);

w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_steer_mjls.tex'));
w.mat('TVXzeroSm', X0); w.mat('TVXfSm', Xf); w.num('TVTSm', T);
w.mat('TVdAsm', dA); w.close();

mvst_steer_ellipse(Sig, X0, Xf, 'Time-varying Markov-jump covariance steering', ...
   fullfile(mvst_figdir(),'ex_tv_steer_mjls.pdf'));
mvst_steer_components(tg, Sig, 'Time-varying Markov-jump covariance components', ...
   fullfile(mvst_figdir(),'ex_tv_steer_mjls_comp.pdf'));
