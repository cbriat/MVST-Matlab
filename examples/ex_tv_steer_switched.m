%% EX_TV_STEER_SWITCHED  Time-varying covariance steering of a switched system.
%  Time-varying counterpart of EX_STEER_SWITCHED: each mode's flow drift
%  decays in time, A_i(t) = A_i + dA*exp(-t), so the system is genuinely
%  time-varying but settles to the LTI switched system as t grows.  Steers the
%  state covariance from X0 to Xf over [0,T] under a dwell-time switching
%  schedule.  Backward pass: the mode-dependent TIME-VARYING I-GRE integrated
%  backward through the schedule from P(T)=M; forward pass: the switched
%  closed-loop moment equation with the time-varying closed loop; M is found by
%  Newton.  Thin caller: the per-class gap closure and forward propagation now
%  live in the generic wrapper MVST_STEER (class 'switched', time-varying path).
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data();
modes = D.sw.modes; Td = D.sw.Tdwell; n = modes{1}.n; m = modes{1}.m; T = 4.0;
Mn = numel(modes);
X0 = [3 0;0 0.5]; Xf = [1.0 0;0 1.6];
dA = [0 0; 1.5 0];                          % decaying drift perturbation

% ----- build the time-varying modes A_i(t) = A_i + dA*exp(-t) -----
md = cell(1,Mn);
for i = 1:Mn
    b = modes{i};
    md{i} = mvst_mode_tv(n, m, ...
        'A', @(t) b.A + dA*exp(-t), 'B', b.B, 'Q', b.Q, 'R', b.R, ...
        'Ac', b.Ac, 'Bc', b.Bc, 'W', b.W, 'X0', b.X0);
end
sys = struct('Tdwell', Td); sys.modes = md;

res = mvst_steer(sys, Xf, T, struct('class','switched','X0',X0));
tg = res.t; Sig = res.Sigma;
fprintf('TV switched steering: %d Newton iterations, terminal gap %.2e\n', res.iters, res.res);

w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_steer_switched.tex'));
w.mat('TVXzeroSw', X0); w.mat('TVXfSw', Xf); w.num('TVTSw', T);
w.mat('TVdAsw', dA); w.close();

mvst_steer_ellipse(Sig, X0, Xf, 'Time-varying switched covariance steering', ...
   fullfile(mvst_figdir(),'ex_tv_steer_switched.pdf'));
sw = Td:Td:(T-1e-9);                                  % switching instants
mvst_steer_components(tg, Sig, 'Time-varying switched covariance components', ...
   fullfile(mvst_figdir(),'ex_tv_steer_switched_comp.pdf'), sw);
