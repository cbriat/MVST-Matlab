function ex_comp_nosdp_matrix()
%EX_COMP_NOSDP_MATRIX  Matrix peak eigenvalue on an interval: Gram SDP
%  (SOSTOOLS + SeDuMi) versus the matrix weighted SOS cone solved directly
%  (SOSNOSDP); also writes the figure and the data macros of BOTH
%  comparison examples (it calls EX_COMP_NOSDP_SCALAR for Example 1).
%
%  Chapter "Semidefinite Programming and Sums of Squares",
%  closing comparison of Section 3.7, Example 2.
%
%  Instance (fixed 2x2 symmetric polynomial matrix on [0,1]):
%     M(tau) = [ tau^2 + 1     tau(1-tau)        ]
%              [ tau(1-tau)    2 + tau - 2 tau^2 ].
%  (The (2,2) entry 2 - tau first proposed puts the maximum of
%  lambda_max(M) at both endpoints, with value 2; the concave entry
%  2 + tau - 2 tau^2 makes the maximizer interior and unique.)
%  Truth: gamma* = max_{tau in [0,1]} lambda_max(M(tau)), from the explicit
%  2x2 eigenvalue formula and Newton's method on its derivative.
%
%  Program at certificate degree 2d (same SOS program, two solvers):
%     min gamma  s.t.  gamma I_2 - M = S0 + tau(1-tau) S1,
%                      S0, S1 2x2 SOS matrices of degrees 2d and 2d-2.
%  By the matrix Lukacs theorem the value is gamma* for every d >= 1, so
%  the sweep shows size/time scaling and accuracy saturation.
%   (i)  SOSTOOLS: S1 by sosposmatrvar (Gram), sosmatrixineq, SeDuMi.
%   (ii) SOSNOSDP: sosmatrixineq_interval in the matrix WSOS cone of
%        degree 2d.  The cone degree is the degree of the constraint, so
%        for 2d > 2 the constraint is padded with cpad tau^(2d) I_2 and the
%        equality cpad = 0 (soseq), which the bridge eliminates.
%
%  Metrics as in EX_COMP_NOSDP_SCALAR.  Requires SOSNOSDP
%  (bundled in thirdparty/sosnosdp).  Solvers: SOSTOOLS + SeDuMi (bundled in thirdparty/).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
if isempty(which('sossolve_nosdp'))
    error('mvst:example:nosdp', ['SOSNOSDP is not on the path (sossolve_nosdp not found). ' ...
          'It ships in thirdparty/sosnosdp: run mvst_setup, or add that folder and run its setup_nosdp.']);
end
gc = mvst_gemcolors(); mvst_latexdefaults();

% Single computational thread for both routes: on the reference machine the
% multithreaded LAPACK of MATLAB R2025a (Apple Accelerate BLAS + NAG LAPACK)
% makes chol of a dense n-by-n matrix ~150x slower for n >= ~170 than for
% n = 160 (10 ms vs 70 us), a platform artifact that would dominate the
% SOSNOSDP timings (dense Hessian factorizations); SeDuMi uses its own
% sparse Cholesky.  One thread removes it and makes the timings comparable.
nth = maxNumCompThreads(1); restoreThreads = onCleanup(@() maxNumCompThreads(nth));

degs = [2 4 8 16 32 64 128];     % 2d; 2..16 requested, beyond: scaling
tcap = 120;                      % a route stops once one solve exceeds tcap seconds

% ---- exact truth: gamma* and its maximizer ----
T = truth();
fprintf('[truth] gamma* = %.16f at tau* = %.16f (f''(tau*) = %.1e, f''''(tau*) = %.4f)\n', ...
    T.gamma, T.tau, T.fp, T.fpp);
fprintf('[truth] endpoint values %.6f, %.6f; eigenvalue gap at tau*: %.4f\n', ...
    T.f0, T.f1, T.gap);
fprintf('[truth] grid (10^6+1 pts): max %.15f at %.6f, sign changes of f'': %d\n', ...
    T.gridmax, T.gridarg, T.nsign);
if ~isnan(T.vpa)
    fprintf('[truth] 40-digit check: gamma* = %s, |double - vpa| = %.1e\n', ...
        T.vpastr, abs(T.gamma - T.vpa));
end

% ---- warm-up (JIT, path caches), discarded ----
solve_sedumi(2); solve_nosdp(2);

% ---- sweep ----
R2.degs = degs; R2.truth = T.gamma; R2.T = T;
R2.sdp   = sweep(@solve_sedumi, degs, tcap, 'SOSTOOLS+SeDuMi', T.gamma);
R2.nosdp = sweep(@solve_nosdp,  degs, tcap, 'SOSNOSDP', T.gamma);
fprintf('\n==== Example 2: min gamma, gamma I - M in matrix WSOS([0,1]), truth %.12f ====\n', T.gamma);
print_table(R2, '2d');

% ---- Example 1 (scalar), for the joint figure and macros ----
R1 = ex_comp_nosdp_scalar();

% ---- figure: solve time vs degree, two stacked panels, log time axis ----
fig = figure('Color','w','Units','centimeters','Position',[2 2 15 16]);
tl = tiledlayout(fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');
panel(nexttile(tl), R1, gc, 'degree $d$ of $p_d$', ...
    '(a) scalar: $\max\lambda$, $p_d-\lambda$ weighted SOS on $[0,1]$');
panel(nexttile(tl), R2, gc, 'certificate degree $2d$', ...
    '(b) matrix: $\min\gamma$, $\gamma I_2-M$ weighted SOS on $[0,1]$');
mvst_savepdf(fig, fullfile(mvst_figdir(),'comp_nosdp.pdf'));

% ---- data macros ----
w = mvst_texdata(fullfile(mvst_figdir(),'data','comp_nosdp.tex'));
w.raw('CompNosdpMatlab', version('-release'));
w.raw('CompNosdpCpu', cpu_name());
w.raw('CompNosdpMatGamma', sprintf('%.12f', T.gamma));
w.raw('CompNosdpMatTau', sprintf('%.12f', T.tau));
w.raw('CompNosdpScalMaxDeg', sprintf('%d', max(R1.degs(isfinite(R1.sdp.val) & isfinite(R1.nosdp.val)))));
w.raw('CompNosdpMatMaxDeg', sprintf('%d', max(R2.degs(isfinite(R2.sdp.val) & isfinite(R2.nosdp.val)))));
emit_runs(w, 'CompNosdpScal', R1);
emit_runs(w, 'CompNosdpMat',  R2);
w.close();
end

% ======================================================================
function panel(ax, R, gc, xl, ttl)
hold(ax,'on'); box(ax,'on'); grid(ax,'on');
k1 = isfinite(R.sdp.time); k2 = isfinite(R.nosdp.time);
loglog(ax, R.degs(k1), R.sdp.time(k1), '-o', 'LineWidth',2, 'MarkerSize',6, ...
    'Color',gc.order(1,:), 'MarkerFaceColor',gc.order(1,:));
loglog(ax, R.degs(k2), R.nosdp.time(k2), '-s', 'LineWidth',2, 'MarkerSize',6, ...
    'Color',gc.order(2,:), 'MarkerFaceColor',gc.order(2,:));
set(ax,'XScale','log','YScale','log','XTick',R.degs,'FontSize',12, ...
    'XLim',[R.degs(1)/1.2 R.degs(end)*1.2]);
xlabel(ax, xl); ylabel(ax, 'solve time [s]');
title(ax, ttl, 'FontWeight','normal', 'FontSize',12);
lg = legend(ax, {'SOSTOOLS + SeDuMi (Gram SDP)','SOSNOSDP (weighted SOS cone)'}, ...
    'Location','northwest'); set(lg,'Box','off','FontSize',11);
end

function emit_runs(w, pre, R)
rt = {'sdp','nosdp'}; nm = {'Sdp','Nosdp'};
for j = 1:numel(R.degs)
    dw = numword(R.degs(j));
    for r = 1:2
        S = R.(rt{r}); p = [pre dw nm{r}];
        if ~isfinite(S.val(j))
            w.raw([p 'Val'], '--'); w.raw([p 'Err'], '--'); w.raw([p 'Time'], '--');
            w.raw([p 'Tsolver'], '--');
            w.raw([p 'Vars'], '--'); w.raw([p 'Distinct'], '--'); w.raw([p 'Eq'], '--');
            w.raw([p 'Iter'], '--'); continue
        end
        w.raw([p 'Val'],  sprintf('%.10f', S.val(j)));
        w.raw([p 'Err'],  scinum(S.err(j)));
        w.raw([p 'Time'], sprintf('%.3g', S.time(j)));
        w.raw([p 'Tsolver'], sprintf('%.3g', S.tsolver(j)));
        w.raw([p 'Vars'], sprintf('%d', S.nvar(j)));
        w.raw([p 'Distinct'], sprintf('%d', S.nsym(j)));
        w.raw([p 'Eq'],   sprintf('%d', S.neq(j)));
        w.raw([p 'Iter'], sprintf('%d', S.iter(j)));
    end
end
end

function s = scinum(x)
if x == 0, s = '0'; return; end
e = floor(log10(x)); m = x/10^e;
if round(m*10)/10 >= 10, m = m/10; e = e+1; end
s = sprintf('\\ensuremath{%.1f\\cdot10^{%d}}', m, e);
end

function s = numword(n)
words = containers.Map('KeyType','double','ValueType','char');
words(2)='Two'; words(4)='Four'; words(8)='Eight'; words(16)='Sixteen';
words(24)='TwentyFour'; words(32)='ThirtyTwo'; words(48)='FortyEight';
words(64)='SixtyFour'; words(96)='NinetySix'; words(128)='OneTwentyEight';
words(192)='OneNinetyTwo'; words(256)='TwoFiftySix';
s = words(n);
end

function s = cpu_name()
s = computer('arch');
if ismac
    [st, out] = system('sysctl -n machdep.cpu.brand_string');
    if st == 0, s = strtrim(out); end
end
end

% ---- exact truth ----
function T = truth()
% M = [a b; b c]; lambda_max = m + r, m = (a+c)/2, h = (a-c)/2, r = sqrt(h^2+b^2)
a  = @(t) t.^2 + 1;   b = @(t) t.*(1-t);   c = @(t) 2 + t - 2*t.^2;
f  = @(t) (a(t)+c(t))/2 + sqrt(((a(t)-c(t))/2).^2 + b(t).^2);
h  = @(t) (3*t.^2 - t - 1)/2;  h1 = @(t) (6*t - 1)/2;  h2 = 3;
m1 = @(t) (1 - 2*t)/2;         m2 = -1;
b1 = @(t) 1 - 2*t;             b2 = -2;
r  = @(t) sqrt(h(t).^2 + b(t).^2);
g  = @(t) h(t).*h1(t) + b(t).*b1(t);
fp = @(t) m1(t) + g(t)./r(t);
fpp = @(t) m2 + (h1(t).^2 + h(t)*h2 + b1(t).^2 + b(t)*b2)./r(t) - g(t).^2./r(t).^3;
tg = linspace(0,1,1e6+1); fg = f(tg);
[T.gridmax, i] = max(fg); T.gridarg = tg(i);
T.nsign = sum(abs(diff(sign(fp(tg)))) > 0);
t = T.gridarg;
for it = 1:50
    dt = fp(t)/fpp(t); t = t - dt;
    if abs(dt) < 1e-17, break; end
end
T.tau = t; T.gamma = f(t); T.fp = fp(t); T.fpp = fpp(t);
T.f0 = f(0); T.f1 = f(1); T.gap = 2*r(t);
T.vpa = NaN; T.vpastr = '';
if exist('vpa','file') && license('test','Symbolic_Toolbox')
    try
        syms s real
        fs = (s^2+1 + 2+s-2*s^2)/2 + sqrt(((s^2+1-2-s+2*s^2)/2)^2 + (s*(1-s))^2);
        ts = vpasolve(diff(fs,s) == 0, s, [T.tau-0.05, T.tau+0.05]);
        gs = vpa(subs(fs, s, ts), 40);
        T.vpa = double(gs); T.vpastr = char(vpa(gs, 20));
    catch
    end
end
end

% ---- sweep and solves ----
function S = sweep(fun, degs, tcap, name, truth)
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
            for q = 1:4
                o2 = fun(degs(j)); tt(end+1) = o2.time; ts(end+1) = o2.tsolver; %#ok<AGROW>
            end
            o.time = median(tt); o.tsolver = median(ts);
        end
    catch err
        S.status{j} = 'error'; S.msg{j} = getReport(err,'extended','hyperlinks','off');
        fprintf(2,'[%s] 2d=%d ERROR:\n%s\n', name, degs(j), S.msg{j});
        stop = true; continue
    end
    f = fieldnames(o);
    for k = 1:numel(f)
        if iscell(S.(f{k})), S.(f{k}){j} = o.(f{k}); else, S.(f{k})(j) = o.(f{k}); end
    end
    S.err(j) = abs(o.val - truth);
    fprintf('[%s] 2d=%3d: gamma=%.12f err=%.2e vars=%d (distinct %d) eq=%d it=%d t=%.3fs (%s) %s\n', ...
        name, degs(j), o.val, S.err(j), o.nvar, o.nsym, o.neq, o.iter, o.time, o.status, o.blocks);
    if o.time > tcap, stop = true; end
end
end

function M = Mpoly(tau)
M = [tau^2 + 1, tau*(1-tau); tau*(1-tau), 2 + tau - 2*tau^2];
end

function o = solve_sedumi(D)
pvar tau; dpvar gam                    % SOSTOOLS 4: decision variables as dpvar
prog = sosprogram(tau, gam);
[prog, S1] = sosposmatrvar(prog, 2, D-2, tau);
prog = sosmatrixineq(prog, gam*eye(2) - Mpoly(tau) - tau*(1-tau)*S1);
prog = sossetobj(prog, gam);
opts.solver = 'sedumi';
opts.params = struct('tol',1e-9,'alg',2,'fid',0);   %#ok<STRNU> SOSTOOLS defaults, silent
t0 = tic;
[~, prog] = evalc('sossolve(prog, opts)');   %#ok<NODEF> silences SOSTOOLS printout
o.time = toc(t0);
o.val = double(sosgetsol(prog, gam));
o = sdp_size(prog, o);
end

function o = solve_nosdp(D)
pvar tau; dpvar gam cpad
P = gam*eye(2) - Mpoly(tau);
if D > 2
    prog = sosprogram(tau, [gam; cpad]);
    P = P + cpad*tau^D*eye(2);         % raises the cone degree to D ...
    prog = soseq(prog, cpad);          % ... with cpad = 0 (eliminated)
else
    prog = sosprogram(tau, gam);
end
prog = sosmatrixineq_interval(prog, P, tau, [0 1]);
prog = sossetobj(prog, gam);
ws = warning('off','MATLAB:nearlySingularMatrix');
w2 = warning('off','MATLAB:singularMatrix');
t0 = tic;
[prog, info] = sossolve_nosdp(prog);
o.time = toc(t0);
warning(ws); warning(w2);
o.val = double(sosgetsol(prog, gam));
o = nosdp_size(info, o);
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
