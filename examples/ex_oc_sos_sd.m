function ex_oc_sos_sd()
%EX_OC_SOS_SD  Sampled-data (timer-dependent) optimal control: SOS certificate
%  of the periodic generalized-Riccati storage function over one period.
%
%  Chapter "Linear-Cost Optimal Control", periodic sampled-data section. For
%  the zero-order-hold control of the diffusion
%     dx = (Ax+Bu)dt + sum_l (Ac_l x + Bc_l u) dW_l,   u held on [t_k,t_{k+1}),
%  the augmented state z=[x;u] flows with no control (du=0), Ab0=[A B;0 0],
%  Abl=[Ac_l Bc_l;0 0], running cost z'Zf z with Zf=blkdiag(Q,R); the gain is
%  the jump decision, obtained by the Schur complement of the u-block at each
%  sample. mvst_igre_periodic_sd returns the within-period value P_gre(tau) on
%  [0,Tbar], the timer-dependent dual variable of the storage functional
%  V(tau,z)=z'P(tau)z. Along the flow it satisfies the LINEAR timer-dependent
%  storage relation (the multiplicative-noise Lyapunov flow, no Riccati term
%  because the input is frozen)
%
%     -dP/dtau = He[Ab0' P(tau)] + sum_l Abl' P(tau) Abl + Zf,   tau in [0,Tbar].
%
%  The exact solution is a matrix exponential of the augmented flow, hence
%  transcendental: a finite SOS *synthesis* cannot reproduce it exactly, so we
%  deliver the VERIFICATION variant. We fit a matrix polynomial P(tau) to the
%  periodic-GRE certificate and certify, by SOS on [0,Tbar] with NO gridding
%  (matrix Positivstellensatz, tau(Tbar-tau) SOS-matrix multipliers, pvar mode
%  as in ex_dof_sos.m), that
%     (i)  P(tau) >= 0  (storage-function positivity), and
%     (ii) the timer-dependent storage DLMI holds to a certified tolerance eps,
%          i.e. -eps I <= R(tau) <= eps I on [0,Tbar] with
%          R(tau) = dP/dtau + He[Ab0'P] + sum_l Abl'P Abl + Zf,
%          so P(tau) satisfies the storage DLMI R(tau)<=0 tightly (equality).
%  The certified eps is minimized as the SOS objective; the interior-point
%  status (pinf/dinf/numerr) and the max pointwise gap to P_gre are printed.
%  Solver: SOSTOOLS + SeDuMi (shipped under matlab/).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();

% ---- periodic sampled-data instance (comp_sos.m data) ----
A=[0 1;0 0]; B=[0;1]; Ac=0.10*[0 0;1 0]; Bc=[0;0.10];
Q=eye(2); R=1.0; Tbar=0.6;
n=size(A,1); m=size(B,2); nm=n+m;
Ab0=[A B; zeros(m,nm)];
Abl=[Ac Bc; zeros(m,nm)];
Zf =blkdiag(Q,R);

% ---- periodic GRE certificate P_gre(tau) over [0,Tbar] ----
out = mvst_igre_periodic_sd(A,B,Ac,Bc,Q,R,Tbar);
fprintf('\n==== sampled-data SOS verification of the periodic-GRE certificate ====\n');
fprintf('  periodic I-GRE fixed point: iters=%d  residual=%.2e\n', out.iters, out.residual);

% ---- fit a matrix polynomial P(tau) to the GRE certificate ----
d = 6;                                   % polynomial degree in the timer
Ns = 80; ts = linspace(0,Tbar,Ns);
Pg = zeros(nm,nm,Ns);
for k=1:Ns, Pg(:,:,k)=out.Pfun(ts(k)); end
C = fit_polymat(ts, Pg, d);              % C{k+1} = coeff of tau^k, nm x nm
% pointwise fit error (fine grid)
tf = linspace(0,Tbar,400); efit=0;
for t=tf, efit=max(efit, norm(evalmat(C,t)-out.Pfun(t),'fro')); end
fprintf('  polynomial fit (deg %d): max_tau ||P_poly-P_gre|| = %.2e\n', d, efit);

% ---- SOS verification program (pvar mode, no gridding, feasibility) ----
% Certify -eps I <= R(tau) <= eps I and P(tau) >= 0 on [0,Tbar] for a fixed
% numeric eps; bisect eps to the smallest SOS-certifiable tolerance (each solve
% is a matrix Positivstellensatz, no gridding). This matches ex_dof_sos.m
% (feasibility, no dpvar objective).
solve_feas = @(e) sd_feas(e, C, Ab0, Abl, Zf, Tbar, nm);
% bracket: grow eps_hi until cleanly feasible, keep eps_lo infeasible
eps_hi = 1e-9; okinfo = solve_feas(eps_hi);
while ~okinfo.solved && eps_hi < 1e-1
    eps_hi = eps_hi*3.16; okinfo = solve_feas(eps_hi);
end
info = okinfo.info; solved = okinfo.solved;
eps_lo = eps_hi/3.16;
for it=1:18                                   % bisection on log-eps
    em = sqrt(max(eps_lo,1e-14)*eps_hi);
    oi = solve_feas(em);
    if oi.solved, eps_hi=em; info=oi.info; solved=true; else, eps_lo=em; end
end
eps_cert = eps_hi;

% ---- pointwise cross-check of the DLMI residual and positivity ----
maxR = 0; minP = inf;
for t=tf
    maxR = max(maxR, norm(resmat(C,t,Ab0,Abl,Zf),'fro'));
    minP = min(minP, min(real(eig(evalmat(C,t)))));
end
fprintf('  SOS solve: solved=%d  (pinf=%d dinf=%d numerr=%d)\n', ...
    solved, info.pinf, info.dinf, info.numerr);
fprintf('  certified DLMI tolerance  eps = %.2e\n', eps_cert);
fprintf('  pointwise max ||R(tau)||  (storage DLMI residual) = %.2e\n', maxR);
fprintf('  pointwise min eig P(tau)  (positivity)            = %.2e\n', minP);

% ================== figure ==================
lamP = arrayfun(@(t) min(real(eig(evalmat(C,t)))), tf);
normR= arrayfun(@(t) norm(resmat(C,t,Ab0,Abl,Zf),'fro'), tf);
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
yyaxis(ax,'left');
h1 = plot(ax, tf, lamP, '-','LineWidth',2.0,'Color',gc.order(1,:));
ylabel(ax,'$\lambda_{\min}(P(\tau))$','Interpreter','latex');
set(ax,'YColor',gc.order(1,:));
yyaxis(ax,'right');
h2 = semilogy(ax, tf, max(normR,1e-16), '-','LineWidth',2.0,'Color',gc.order(2,:));
ylabel(ax,'$\|R(\tau)\|$ (storage DLMI residual)','Interpreter','latex');
set(ax,'YColor',gc.order(2,:));
yline(ax, eps_cert, '--','Color','k','LineWidth',1.2);
xlabel(ax,'timer $\tau$','Interpreter','latex');
lg = legend(ax,[h1 h2],{'$\lambda_{\min}(P(\tau))\ge0$', ...
    '$\|R(\tau)\|\le\epsilon$'},'Interpreter','latex','Location','best');
set(lg,'Box','off'); set(ax,'FontSize',10);
title(ax,'Sampled-data: SOS certifies the periodic-GRE storage function', ...
    'Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(mvst_figdir(),'oc_sos_sd.pdf'));
end

% =========================================================================
function out = sd_feas(epsval, C, Ab0, Abl, Zf, Tbar, nm)
%SD_FEAS  SOS feasibility of the storage DLMI at a fixed tolerance epsval:
%  certify  -eps I <= R(tau) <= eps I  and  P(tau) >= 0  on [0,Tbar] by a
%  matrix Positivstellensatz (tau(Tbar-tau) SOS-matrix multipliers, no grid).
pvar tau
prog = sosprogram(tau);
P = polymat_pvar(C, tau); Pdot = polymat_pvar_der(C, tau);
R = Pdot + Ab0.'*P + P*Ab0 + Zf + Abl.'*P*Abl;   % storage DLMI residual (~0)
[prog, Sp] = sossosmatrixvar(prog, monomials(tau,0:2), [nm nm]);
[prog, Sm] = sossosmatrixvar(prog, monomials(tau,0:2), [nm nm]);
[prog, Sz] = sossosmatrixvar(prog, monomials(tau,0:2), [nm nm]);
prog = sosmatrixineq(prog, epsval*eye(nm) - R - Sp*tau*(Tbar-tau));  % R <=  eps I
prog = sosmatrixineq(prog, epsval*eye(nm) + R - Sm*tau*(Tbar-tau));  % R >= -eps I
prog = sosmatrixineq(prog, P              - Sz*tau*(Tbar-tau));      % P(tau) >= 0
opt.solver='sedumi'; opt.params.fid=0;
[~, info] = sossolve(prog, opt);
out.info = info;
% strict acceptance: a trustworthy SOS certificate needs a clean solve, else
% SeDuMi returns pinf=dinf=0 with numerr=1 near the feasibility boundary and
% the "solution" may violate the constraints (false positive in the bisection).
out.solved = (info.pinf==0 && info.dinf==0 && info.numerr==0);
end

% =========================================================================
function C = fit_polymat(ts, Pg, d)
%FIT_POLYMAT  Least-squares matrix-polynomial fit: C{k+1} = coeff of tau^k.
nm = size(Pg,1); C = cell(1,d+1); for k=1:d+1, C{k}=zeros(nm); end
for a=1:nm
    for b=a:nm
        y = squeeze(Pg(a,b,:));
        p = polyfit(ts(:), y(:), d);      % highest degree first
        p = fliplr(p);                    % now p(k) = coeff of tau^(k-1)
        for k=1:d+1
            C{k}(a,b)=p(k); C{k}(b,a)=p(k);
        end
    end
end
end

function M = evalmat(C,t)
%EVALMAT  Numeric value of the matrix polynomial C at tau=t.
M = zeros(size(C{1}));
for k=1:numel(C), M = M + C{k}*t^(k-1); end
M=(M+M.')/2;
end

function R = resmat(C,t,Ab0,Abl,Zf)
%RESMAT  Numeric storage-DLMI residual R(t)=dP/dt+He[Ab0'P]+Abl'P Abl+Zf.
P = evalmat(C,t);
Pd = zeros(size(P));
for k=2:numel(C), Pd = Pd + (k-1)*C{k}*t^(k-2); end
Pd=(Pd+Pd.')/2;
R = Pd + Ab0.'*P + P*Ab0 + Zf + Abl.'*P*Abl; R=(R+R.')/2;
end

function P = polymat_pvar(C, tau)
%POLYMAT_PVAR  Build the pvar matrix polynomial from numeric coefficients.
P = C{1};
for k=2:numel(C), P = P + C{k}*tau^(k-1); end
end

function Pd = polymat_pvar_der(C, tau)
%POLYMAT_PVAR_DER  Derivative pvar matrix polynomial.
Pd = zeros(size(C{1}));
for k=2:numel(C), Pd = Pd + (k-1)*C{k}*tau^(k-2); end
end
