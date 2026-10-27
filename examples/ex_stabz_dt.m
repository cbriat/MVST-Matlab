function ex_stabz_dt()
%EX_STABZ_DT  Discrete-time stochastic system with random parameters:
%  mean-square state-feedback stabilization.  Synthesis mirror of the analysis
%  example ex:stab:dtstoch.  The controlled system
%     x(k+1) = A0 x + B0 u + nu(k) (A1 x + B1 u),   u = K x,
%  with zero-mean unit-variance nu, is open-loop mean-square unstable; K is
%  obtained by solving the discrete-time synthesis condition (congruence of the
%  DT mean-square Lyapunov inequality) as a matrix sum-of-squares program
%  (discrete dwell-time) through mvst_stabz_sos.
%
%  Two figures + one data file:
%    figures/stabz_dt_paths.pdf : open-loop paths (top), closed-loop (middle),
%                                 control u (bottom);
%    figures/stabz_dt_cov.pdf   : tr X(k), open loop RED vs closed loop BLUE;
%    figures/data/ex_stabz_dt.tex : \StabzDtOpen, \StabzDtClosed (lifted
%                                   transition spectral radius, open / closed).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

% ------------------- system: unstable DT stochastic system ------------------
n = 2; m = 1;
A0 = [1.1 0.3; 0 1.05];    % unstable mean transition
B0 = [0; 1];
A1 = 0.2*eye(n);           % state-multiplicative random channel
B1 = 0.1*[0; 1];           % input selector actuating the random channel (eq:stabz:dt:data)

Topen   = mvst_liftsym(zeros(n), {A0, A1});      % lifted open-loop transition
radOpen = max(abs(eig(Topen)));                  % > 1 : mean-square unstable

% ------------------- chapter synthesis (matrix SOS, SOSTOOLS + SeDuMi) ------
sys = struct('A0', A0, 'B0', B0, 'Ac', A1, 'Bc', B1);
res = mvst_stabz_sos(sys, struct('type','discrete'));
K   = res.Kfun(0);
Tcl       = mvst_liftsym(zeros(n), {A0 + B0*K, A1 + B1*K});
radClosed = max(abs(eig(Tcl)));                  % < 1 : mean-square stable

% ------------------- sample paths: open (K=0) and closed --------------------
Nstep = 16; nreal = 16; x0 = [1; -1]; rng(7);
kk = 0:Nstep;
Xo = zeros(n, Nstep+1, nreal); Xc = zeros(n, Nstep+1, nreal);
Uc = zeros(m, Nstep+1, nreal);
for r = 1:nreal
    xo = x0; xc = x0;
    for j = 1:Nstep+1
        Xo(:,j,r) = xo; Xc(:,j,r) = xc; Uc(:,j,r) = K*xc;
        nu = randn;                                  % same draw for both loops
        xo = A0*xo + nu*(A1*xo);
        uc = K*xc; xc = A0*xc + B0*uc + nu*(A1*xc + B1*uc);
    end
end

% ------------------- covariance trace (exact) -------------------------------
Us = mvst_svecbasis(n); s0 = Us' * reshape(eye(n), [], 1);
trO = zeros(1,Nstep+1); trC = zeros(1,Nstep+1); sO = s0; sC = s0;
for j = 1:Nstep+1
    trO(j) = trace(reshape(Us*sO, n, n));
    trC(j) = trace(reshape(Us*sC, n, n));
    sO = Topen*sO; sC = Tcl*sC;
end

% ------------------- the two standard figures -------------------------------
mvst_stabfig(figdir, 'dt', kk, Xo, Xc, Uc, kk, trO, trC, true, gc);
mvst_stabtable(figdir, 'dt');   %% generated summary-table row

% ------------------- data macros --------------------------------------------
w = mvst_texdata(fullfile(figdir,'data','ex_stabz_dt.tex'));
w.num('StabzDtOpen',   radOpen);
w.num('StabzDtClosed', radClosed);
w.mat('StabzDtK', K);
w.close();
fprintf('DT stabilization: open radius %.4f, closed %.4f\n', radOpen, radClosed);
end
