%% EX_STEER  Minimum-energy covariance steering of a linear diffusion (warm-up).
%  Single continuous-time linear system with multiplicative (Ito) noise, steered
%  in second moment from X0 to Xf over [0,T] by the forward-backward algorithm:
%  backward pass = finite-horizon I-GRE gain K(t) from P(T)=M; forward pass =
%  closed-loop covariance; the multiplier M is found by Newton so Sigma(T;M)=Xf.
%  Minimum-energy cost (Q=0, R=I).  This simplest instance already prices the
%  terminal equality with an INDEFINITE multiplier M*, which is the point of the
%  warm-up.  Thin caller over the generic wrapper MVST_STEER (class 'ct').
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); D = ex_data();
st = D.steer; md = st.md; T = st.T; X0 = md.X0; Xf = st.X1f;

res = mvst_steer(md, Xf, T, struct('class','ct','X0',X0));
M = res.M; it = res.iters; Sig = res.Sigma;
fprintf('CT steering: %d Newton iterations, terminal gap %.2e\n', it, res.res);
ev = eig((M+M.')/2);
fprintf('M* eigenvalues: [% .3f, % .3f]\n', ev(1), ev(2));

w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_steer.tex'));
w.mat('AstE', md.A); w.mat('BstE', md.B);
w.mat('XzeroSt', X0); w.mat('XfSt', Xf);
w.mat('Mstar', M); w.num('Tst', T); w.num('ItSt', it, 0);
w.close();

mvst_steer_ellipse(Sig, X0, Xf, 'Continuous-time covariance steering', ...
   fullfile(mvst_figdir(),'ex_steer.pdf'));
