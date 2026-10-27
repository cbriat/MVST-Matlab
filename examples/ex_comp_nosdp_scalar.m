function R = ex_comp_nosdp_scalar()
%EX_COMP_NOSDP_SCALAR  Scalar interval minimization: Gram SDP (SOSTOOLS +
%  SeDuMi) versus the weighted SOS cone solved directly (SOSNOSDP).
%
%  Chapter "Semidefinite Programming and Sums of Squares",
%  closing comparison of Section 3.7, Example 1.
%
%  Instance (degree d even, K = (d-2)/2):
%     p_d(tau) = (tau - 1/3)^2 * prod_{k=1}^{K} (1 + tau^2/k) + 1.
%  Every factor 1 + tau^2/k is >= 1, so p_d >= 1 on [0,1] with equality at
%  tau = 1/3 only: the minimum of p_d on [0,1] is EXACTLY 1 (checked below
%  on a 10^6-point grid).  The normalization 1 + tau^2/k = (tau^2 + k)/k
%  divides (p_d - 1) by K!; without it the data reach 1e23 at d = 48 and
%  the value 1 is lost in double precision.
%
%  Program (same SOS program, two solvers):
%     lambda*_d = max lambda  s.t.  p_d - lambda = s0 + tau(1-tau) s1,
%                 s0 SOS of degree d, s1 SOS of degree d-2,
%  i.e. p_d - lambda in the weighted SOS cone of [0,1] with weights 1 and
%  tau(1-tau).  The Lukacs-Markov theorem makes lambda*_d = 1 exactly.
%   (i)  SOSTOOLS: explicit multiplier s1 (sossosvar), Gram SDP, SeDuMi.
%   (ii) SOSNOSDP: sosineq_interval, interpolant basis, no Gram matrix.
%
%  Reported per degree and route: lambda, |lambda - 1|, size of the conic
%  program actually solved (read from the solver structures), iterations,
%  wall-clock time of the solve call only (tic/toc around sossolve or
%  sossolve_nosdp; single computational thread; one warm-up solve first;
%  median of 5 repeats when a solve takes less than 2 s).
%
%  R = EX_COMP_NOSDP_SCALAR() returns the results; the figure and the data
%  macros (both examples) are written by EX_COMP_NOSDP_MATRIX, which calls
%  this function.
%
%  Requires SOSNOSDP (https://github.com/cbriat/sosnosdp, bundled in
%  thirdparty/sosnosdp).  Solvers: SOSTOOLS + SeDuMi (bundled in thirdparty/).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
if isempty(which('sossolve_nosdp'))
    error('mvst:example:nosdp', ['SOSNOSDP is not on the path (sossolve_nosdp not found). ' ...
          'It ships in thirdparty/sosnosdp: run mvst_setup, or add that folder and run its setup_nosdp.']);
end

% Single computational thread for both routes: on the reference machine the
% multithreaded LAPACK of MATLAB R2025a (Apple Accelerate BLAS + NAG LAPACK)
% makes chol of a dense n-by-n matrix ~150x slower for n >= ~170 than for
% n = 160 (10 ms vs 70 us), a platform artifact that would dominate the
% SOSNOSDP timings (dense Hessian factorizations); SeDuMi uses its own
% sparse Cholesky.  One thread removes it and makes the timings comparable.
nth = maxNumCompThreads(1); restoreThreads = onCleanup(@() maxNumCompThreads(nth));

degs   = [4 8 16 24 32 48 64 96 128 192 256];   % 4..48 requested, beyond: scaling
tcap   = 120;     % a route stops once one solve exceeds tcap seconds
truth  = 1;

% ---- exact-truth verification on a fine grid ----
tg = linspace(0,1,1e6+1);
R.grid = struct('n',numel(tg),'min',[],'argmin',[],'pAtThird',[]);
for d = degs
    pv = pd_eval(tg,d);
    [m,i] = min(pv);
    R.grid.min(end+1) = m; R.grid.argmin(end+1) = tg(i);
    R.grid.pAtThird(end+1) = pd_eval(1/3,d);
    fprintf('[truth] d=%2d: grid min %.15f at tau=%.6f, p_d(1/3)=%.15f, max %.3f\n', ...
        d, m, tg(i), pd_eval(1/3,d), max(pv));
end

% ---- warm-up (JIT, path caches), discarded ----
solve_sedumi(4); solve_nosdp(4);

% ---- sweep ----
R.degs = degs; R.truth = truth;
R.sdp   = sweep(@solve_sedumi, degs, tcap, 'SOSTOOLS+SeDuMi');
R.nosdp = sweep(@solve_nosdp,  degs, tcap, 'SOSNOSDP');

fprintf('\n==== Example 1: max lambda, p_d - lambda in WSOS([0,1]), truth 1 ====\n');
print_table(R, 'd');
end

% ======================================================================
function S = sweep(fun, degs, tcap, name)
S = struct('val',nan(size(degs)),'err',nan(size(degs)),'nvar',nan(size(degs)), ...
    'nsym',nan(size(degs)),'neq',nan(size(degs)),'iter',nan(size(degs)), ...
    'time',nan(size(degs)),'tsolver',nan(size(degs)),'status',{cell(size(degs))}, ...
    'blocks',{cell(size(degs))},'msg',{cell(size(degs))});
stop = false;
for j = 1:numel(degs)
    if stop, S.status{j} = 'skipped'; continue; end
    try
        o = fun(degs(j));
        if o.time < 2
            tt = o.time; ts = o.tsolver;
            for r = 1:4
                o2 = fun(degs(j)); tt(end+1) = o2.time; ts(end+1) = o2.tsolver; %#ok<AGROW>
            end
            o.time = median(tt); o.tsolver = median(ts);
        end
    catch err
        S.status{j} = 'error'; S.msg{j} = getReport(err,'extended','hyperlinks','off');
        fprintf(2,'[%s] d=%d ERROR:\n%s\n', name, degs(j), S.msg{j});
        stop = true; continue
    end
    f = fieldnames(o);
    for k = 1:numel(f)
        if iscell(S.(f{k})), S.(f{k}){j} = o.(f{k}); else, S.(f{k})(j) = o.(f{k}); end
    end
    S.err(j) = abs(o.val - 1);
    fprintf('[%s] d=%2d: lambda=%.12f err=%.2e vars=%d (distinct %d) eq=%d it=%d t=%.3fs (%s)\n', ...
        name, degs(j), o.val, S.err(j), o.nvar, o.nsym, o.neq, o.iter, o.time, o.status);
    if o.time > tcap, stop = true; end
end
end

function o = solve_sedumi(d)
pvar tau; dpvar lam          % SOSTOOLS 4: decision variables as dpvar
p = pd_poly(tau, d);
prog = sosprogram(tau, lam);
[prog, s1] = sossosvar(prog, monomials(tau, 0:(d-2)/2));
prog = sosineq(prog, p - lam - tau*(1-tau)*s1);
prog = sossetobj(prog, -lam);
opts.solver = 'sedumi';
opts.params = struct('tol',1e-9,'alg',2,'fid',0);   %#ok<STRNU> SOSTOOLS defaults, silent
t0 = tic;
[~, prog] = evalc('sossolve(prog, opts)');   %#ok<NODEF> silences SOSTOOLS printout
o.time = toc(t0);
o.val = double(sosgetsol(prog, lam));
o = sdp_size(prog, o);
end

function o = solve_nosdp(d)
pvar tau lam
p = pd_poly(tau, d);
prog = sosprogram(tau, lam);
prog = sosineq_interval(prog, p - lam, tau, [0 1]);
prog = sossetobj(prog, -lam);
ws = warning('off','MATLAB:nearlySingularMatrix');
w2 = warning('off','MATLAB:singularMatrix');
t0 = tic;
[prog, info] = sossolve_nosdp(prog);
o.time = toc(t0);
warning(ws); warning(w2);
o.val = double(sosgetsol(prog, lam));
o = nosdp_size(info, o);
end

% ---- p_d as a pvar polynomial and as a numerical function ----
function p = pd_poly(tau, d)
p = (tau - 1/3)^2;
for k = 1:(d-2)/2, p = p*(1 + tau^2/k); end
p = p + 1;
end
function v = pd_eval(t, d)
v = (t - 1/3).^2;
for k = 1:(d-2)/2, v = v.*(1 + t.^2/k); end
v = v + 1;
end

function print_table(R, lbl)
fprintf('%4s | %-8s | %16s %9s %7s %7s %6s %4s %9s %9s | %s\n', lbl,'route', ...
    'value','error','vars','distinct','eq','it','time[s]','solver[s]','status');
rt = {'sdp','nosdp'}; nm = {'SeDuMi','SOSNOSDP'};
for j = 1:numel(R.degs)
    for r = 1:2
        S = R.(rt{r});
        fprintf('%4d | %-8s | %16.12f %9.2e %7d %7d %6d %4d %9.4f %9.4f | %s\n', R.degs(j), nm{r}, ...
            S.val(j), S.err(j), S.nvar(j), S.nsym(j), S.neq(j), S.iter(j), S.time(j), ...
            S.tsolver(j), S.status{j});
    end
end
end

% ---- size of the conic program actually solved ----
function o = sdp_size(prog, o)
% SeDuMi primal: x = [free; vec(Gram blocks)] with A x = b (one row per
% matched coefficient).  Sizes read back from the solved SOSTOOLS program.
si = prog.solinfo;
s = [];
for i = 1:numel(si.extravar.primal), s(end+1) = size(si.extravar.primal{i},1); end %#ok<AGROW>
for i = 1:prog.var.num
    if strcmp(prog.var.type{i},'sos'), s(end+1) = size(si.var.primal{i},1); end %#ok<AGROW>
end
o.nvar = numel(si.x);                          % length of SeDuMi's x
nf = o.nvar - sum(s.^2);                       % free variables (K.f)
o.nsym = nf + sum(s.*(s+1)/2);                 % distinct scalar unknowns
o.neq = numel(si.y);                           % equality constraints
o.iter = si.info.iter;
o.tsolver = si.info.wallsec;
o.blocks = sprintf('K.f=%d, K.s=[%s]', nf, num2str(s));
if si.info.pinf || si.info.dinf, o.status = 'infeasible';
elseif si.info.numerr == 1, o.status = 'inaccurate';
elseif si.info.numerr == 2, o.status = 'failed';
else, o.status = 'optimal'; end
end

function o = nosdp_size(info, o)
% SOSNOSDP primal (nsd_ipm): x in K* of dimension cone_dim, Q' x = b with
% rank(L) equality rows; dual: s = c + L z in K (n_free coordinates z).
o.nvar = info.cone_dim;
o.nsym = info.cone_dim;
o.neq = info.rank;
o.iter = info.iters;
o.tsolver = info.time;
o.blocks = sprintf('n_free=%d, nu=%d, %s', info.n_free, info.barrier_parameter, ...
    strjoin(info.cones, '; '));
o.status = info.status;
end
