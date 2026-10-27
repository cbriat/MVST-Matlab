function res = mvst_passivity(sys, opts)
%MVST_PASSIVITY  Routed passivity analysis with indices (exact CT test by default).
%   RES = MVST_PASSIVITY(SYS, OPTS) tests passivity of the SQUARE-port system
%        z = C x + Dw w,  w, z in R^m,
%   whose dynamics are described by the struct SYS (the MVST_L2GAIN naming:
%   .A0 or .A, .Bw, .C, .Dw, plus class fields .Ac/.Dwc, .Pig, .Td, .Tbar,
%   .Jf/.lamf/.Bdf, or .domain='dt'), and returns
%     RES.passive      true when the system is passive (nu_isp >= -tol);
%     RES.nu_isp       input-strict-passivity index nu, the largest nu with
%                      supply  He(w'z) - nu |w|^2  dissipative;
%     RES.rho_osp      output-strict-passivity index rho >= 0, the largest rho
%                      with supply  He(w'z) - rho |z|^2  dissipative (NaN when
%                      not passive, or when not computed, see below);
%     RES.certificate  storage struct .S (V(x)=x'S x, S >= 0 for the ISP
%                      supply at nu_cert) or [] when unavailable;
%     RES.method       the test that produced the numbers (routing table);
%     RES.class        the resolved class tag;
%     RES.diagnostics  struct (.tol, .exact, .note, ...).
%
%   NU CONVENTION (single, the book's).  Dissipativity chapter, Section
%   sec:diss:passivity, eq:diss:M:pass: the ISP supply is He(w'z) - nu|w|^2,
%   i.e. M22 = He(Dw) - nu I, so for an LTI system
%        nu = min_w lambda_min( G(jw) + G(jw)' ),   G + G' >= nu I.
%   MVST_ISP_MARGIN_CT and MVST_PASSIVITY_CT (changed 2026-10-27, it used to
%   return nu/2) follow the same convention.  rho is the book's OSP index
%   (He(w'z) - rho|z|^2).
%
%   ROUTING TABLE
%     class                        condition           method (exact?)
%     ct, no noise, A Hurwitz      Dw+Dw' nonsing.         'hamiltonian' : imaginary
%                                                      axis test of the positive-
%                                                      real Hamiltonian, bisection
%                                                      in nu and rho (exact)
%     dt, A Schur                  nonsingular after   'hamiltonian-bilinear':
%                                  the map             bilinear map to CT (exact)
%     ct / dt                      Dw+Dw' singular     'isp-sos' : MVST_ISP_MARGIN_CT
%                                                      SOS/LMI program (nu only)
%     ct with .Ac/.Dwc (Ito noise) any                 'isp-sos' (nu only)
%     sd (flow port)               any                 'isp-sos' on the flow
%     mjls                         any                 'diss-sos-modes' bisection
%                                                      (MVST_DISS_SOS_MODES)
%     poisson                      any                 'diss-sos' bisection
%                                                      (MVST_DISS_SOS)
%     dwell (switched)             any                 'diss-sos-modes' (minimum
%                                                      dwell-time Td)
%   OPTS fields (all optional; struct or MVST.OPTIONS):
%     .method 'auto' (default) | 'hamiltonian' | 'grid' (MVST_PASSIVITY_CT
%             frequency sweep, approximate) | 'sos' (force the SOS/LMI route);
%     .rho    logical: also compute rho_osp on SOS routes (default false; it
%             is always computed on the exact route);
%     .class, .Td, .Tbar, .solver as in MVST_L2GAIN; .tol passivity threshold
%             (default 1e-6).
%
%   Example:
%       r = mvst_passivity(struct('A0',[0 1;-1 -1],'Bw',[0;1],'C',[1 0],'Dw',0.5));
%       r.passive, r.nu_isp, r.rho_osp, r.method
%
%   See also MVST_PASSIVITY_CT, MVST_ISP_MARGIN_CT, MVST_HINFNORM, MVST_KYP,
%   MVST_DISS_SOS.
if nargin < 2 || isempty(opts); opts = struct(); end
oo = mvst_opts(opts, false);
meth = lower(getopt(oo, 'method', 'auto'));
tolP = getopt(oo, 'tol', 1e-6);
wantRho = isequal(getopt(oo, 'rho', false), true);

A  = getfld(sys, {'A0','A'}, []);
Bw = getfld(sys, {'Bw'}, []);
C  = getfld(sys, {'C'}, []);
if isempty(A);  error('mvst:passivity:noA',  'SYS must carry a drift .A0 (or .A).'); end
if isempty(Bw); error('mvst:passivity:noBw', 'SYS must carry an input .Bw.'); end
if isempty(C);  error('mvst:passivity:noC',  'SYS must carry an output .C.'); end
m = size(Bw,2);
if maxdim(C,1) ~= m
    error('mvst:passivity:square', ...
        'Passivity needs a square port (size(C,1) == size(Bw,2)); got %d outputs, %d inputs.', ...
        maxdim(C,1), m);
end
Dw = getfld(sys, {'Dw'}, zeros(m));
cls = resolve_class(sys, oo);
noisy = hasnoise(sys);
diag = struct('tol', NaN, 'exact', false, 'note', '');

if strcmp(meth,'grid') && ~any(strcmp(cls,{'ct','dt'}))
    error('mvst:passivity:grid', 'The grid method covers deterministic LTI (ct, dt) only.');
end

switch cls
    case {'ct','dt'}
        if noisy && strcmp(cls,'dt')
            error('mvst:passivity:dtnoise', 'Discrete-time multiplicative noise is not supported.');
        end
        if noisy
            res = sos_flow(sys, A, Bw, C, Dw, oo, tolP, cls, 'isp-sos', ...
                'Ito noise: flow ISP program (nu only)');
            return;
        end
        A1 = A(:,:,1); C1 = C(:,:,1);
        if strcmp(cls,'dt')
            [A1, Bw, C1, Dw] = bilinear(A1, Bw, C1, Dw);
            diag.note = 'DT mapped to CT by the bilinear map (exact; gains/indices preserved)';
        end
        if max(real(eig(A1))) >= 0
            error('mvst:passivity:unstable', ...
                'The system must be asymptotically stable (A Hurwitz / Schur) for passivity indices.');
        end
        smin = min(eig((Dw+Dw.')/2))*2;           % lambda_min(Dw+Dw')
        tolD = 1e-8*max(1, norm(Dw));
        switch meth
            case 'grid'
                [pa, nu] = mvst_passivity_ct(A1, Bw, C1, Dw);
                diag.tol = 1e-7; diag.note = 'frequency grid, approximate';
                res = pack(pa, nu, NaN, [], 'grid', cls, diag);
            case 'sos'
                res = sos_flow(sys, A1, Bw, C1, Dw, oo, tolP, cls, 'isp-sos', 'forced SOS route');
            otherwise
                if abs(smin) > tolD
                    [nu, rho, S, nuc] = exact_indices(A1, Bw, C1, Dw, smin);
                    diag.tol = 1e-9; diag.exact = true;
                    diag.note = 'positive-real Hamiltonian imaginary-axis test, bisection';
                    cert = [];
                    if ~isempty(S)
                        % DT route: S is the storage of the bilinear-mapped CT realization;
                        % the DT storage satisfying the DT KYP is S/2 (measure 2 in the map).
                        if strcmp(cls,'dt'); S = S/2; end
                        cert = struct('S', S, 'nu', nuc);
                    end
                    pa = nu >= -tolP;
                    if ~pa; rho = NaN; end
                    res = pack(pa, nu, rho, cert, 'hamiltonian', cls, diag);
                    if strcmp(cls,'dt'); res.method = 'hamiltonian-bilinear'; end
                else
                    if exist('sosprogram','file') == 2
                        res = sos_flow(sys, A1, Bw, C1, Dw, oo, tolP, cls, 'isp-sos', ...
                            'Dw+Dw'' singular: exact test not applicable, ISP program');
                    else
                        warning('mvst:passivity:gridfallback', ...
                            ['Dw+Dw'' is singular and SOSTOOLS is not on the path: falling ' ...
                             'back to the approximate frequency grid.']);
                        [pa, nu] = mvst_passivity_ct(A1, Bw, C1, Dw);
                        diag.tol = 1e-7; diag.note = 'fallback: frequency grid, approximate';
                        res = pack(pa, nu, NaN, [], 'grid', cls, diag);
                    end
                end
        end

    case 'sd'
        res = sos_flow(sys, A, Bw, C, Dw, oo, tolP, cls, 'isp-sos', 'flow port of the sampled-data system');

    case {'mjls','dwell','poisson'}
        res = sos_dissipativity(sys, cls, oo, tolP, wantRho);

    otherwise
        error('mvst:passivity:class', 'Unknown class tag ''%s''.', cls);
end
end

% ========================================================================
function res = pack(pa, nu, rho, cert, method, cls, diag)
res = struct('passive', logical(pa), 'nu_isp', nu, 'rho_osp', rho, ...
             'certificate', cert, 'method', method, 'class', cls, ...
             'diagnostics', diag);
end

% ------------------------------------------------------------------------
function [A, B, C, D] = bilinear(A, B, C, D)
%BILINEAR  z=(1+s)/(1-s): unit circle -> imaginary axis, G(e^{jt}) = Gc(jw).
n = size(A,1);
if max(abs(eig(A))) >= 1
    error('mvst:passivity:unstable', 'The DT system must be Schur stable.');
end
IA = eye(n) + A;
Ac = IA \ (A - eye(n)); Bc = IA \ B; Cc = 2*(C / IA); Dc = D - C*(IA \ B);
A = Ac; B = Bc; C = Cc; D = Dc;
end

% ------------------------------------------------------------------------
function [nu, rho, S, nuc] = exact_indices(A, B, C, D, smin)
%EXACT_INDICES  nu_isp and rho_osp by bisection on the imaginary-axis test.
scl = max(1, norm(D) + norm(C,'fro')*norm(B,'fro')/max(1,norm(A,'fro')));
% ---- nu: feasible iff nu < smin and the Hamiltonian has no imaginary eig
hiN = smin;
if testvsp(A,B,C,D,0,0)
    loN = 0; stepN = max(1, smin);
else
    stepN = 1; loN = -stepN;
    while ~testvsp(A,B,C,D,loN,0)
        stepN = 2*stepN; loN = -stepN;
        if stepN > 1e12; break; end
    end
end
nu = bisect(@(v) testvsp(A,B,C,D,v,0), loN, hiN, 1e-9*scl);
% ---- rho (only meaningful for a passive system)
rho = NaN;
if nu >= -1e-9
    loR = 0; hiR = 1;
    while testvsp(A,B,C,D,0,hiR)
        loR = hiR; hiR = 2*hiR;
        if hiR > 1e12; rho = Inf; break; end
    end
    if isnan(rho); rho = bisect(@(v) testvsp(A,B,C,D,0,v), loR, hiR, 1e-9*max(1,hiR)); end
end
% ---- certificate: stabilizing ARE solution just below nu*
nuc = nu - 1e-6*scl;
S = storage(A,B,C,D,nuc);
end

% ------------------------------------------------------------------------
function x = bisect(f, lo, hi, tol)
%BISECT  Largest x in [lo,hi] with f true (f true at lo, false at hi).
while (hi - lo) > tol
    mid = 0.5*(lo + hi);
    if f(mid); lo = mid; else; hi = mid; end
end
x = lo;
end

% ------------------------------------------------------------------------
function ok = testvsp(A, B, C, D, nu, rho)
%TESTVSP  True iff G+G* - nu I - rho G*G > 0 on the imaginary axis (incl. inf).
%   Realization of Phi = G + G~ - nu I - rho G~ G on the state (x, xt):
%     x' = A x + B w,   xt' = rho C'C x - A' xt - C'(I - rho D) w,
%     out = (C - rho D'C) x + B' xt + (D + D' - rho D'D - nu I) w.
%   Phi(jw) is singular iff jw is an eigenvalue of Aphi - Bphi Dphi^{-1} Cphi,
%   and Phi(inf) = Dphi must be positive definite.
n = size(A,1); m = size(B,2);
Dp = D + D.' - rho*(D.'*D) - nu*eye(m); Dp = (Dp + Dp.')/2;
[~, pc] = chol(Dp);
if pc ~= 0; ok = false; return; end
Ap = [A, zeros(n); rho*(C.'*C), -A.'];
Bp = [B; -C.'*(eye(m) - rho*D)];
Cp = [C - rho*(D.'*C), B.'];
Ai = Ap - Bp*(Dp \ Cp);
ev = eig(Ai);
scl = max(1, norm(Ai,1));
ok = all(abs(real(ev)) > 1e-9*scl);
end

% ------------------------------------------------------------------------
function S = storage(A, B, C, D, nu)
%STORAGE  Stabilizing solution of  A'S+SA+(SB-C')R^{-1}(B'S-C)=0, R=D+D'-nu I.
%   V(x)=x'S x then satisfies dV/dt <= He(w'z) - nu|w|^2 (S >= 0).
S = [];
try
    n = size(A,1); m = size(B,2);
    R = D + D.' - nu*eye(m);
    F = A - B*(R\C); G = B*(R\B.'); Q = C.'*(R\C);
    H = [F, G; -Q, -F.'];
    [U, T] = schur(H, 'real'); [U, ~] = ordschur(U, T, real(ordeig(T)) < 0);
    S = U(n+1:end,1:n) / U(1:n,1:n); S = (S + S.')/2;
    L = [A.'*S + S*A, S*B - C.'; B.'*S - C, -R];
    if min(eig(S)) < -1e-7*(1+norm(S)) || max(eig((L+L.')/2)) > 1e-6*(1+norm(L))
        S = [];
    end
catch
    S = [];
end
end

% ------------------------------------------------------------------------
function res = sos_flow(sys, A, Bw, C, Dw, oo, tolP, cls, method, note)
%SOS_FLOW  ISP index of the (Ito) flow by the SOS/LMI program.
oo2 = mvst_opts(oo, true, 'sos');
n = size(A,1);
Cc = tocell3(getfld(sys, {'Ac'}, {})); Dwc = tocell3(getfld(sys, {'Dwc'}, {}));
o = oo2;
[nu, P] = mvst_isp_margin_ct(A(:,:,1), Bw, C(:,:,1), Dw, Cc, Dwc, o);
diag = struct('tol', 1e-6, 'exact', false, 'note', note);
cert = struct('S', -P, 'nu', nu);
rho = NaN;
res = pack(nu >= -tolP, nu, rho, cert, method, cls, diag);
if n == 0; res.certificate = []; end
end

% ------------------------------------------------------------------------
function res = sos_dissipativity(sys, cls, oo, tolP, wantRho)
%SOS_DISSIPATIVITY  nu (and rho) by bisection on the SOS passivity feasibility.
oo2 = mvst_opts(oo, true, 'sos');
sol = oo2;   % normalized options (reltol/verbosity) travel to the solves
if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
Dw = getfld(sys, {'Dw'}, []);
m = size(sys.Bw, 2);
if isempty(Dw); Dw = zeros(m); end
hiN = min(eig((Dw+Dw.')/2))*2;
feas = @(nu, rho) feasible(sys, cls, nu, rho, sol);
if feas(0, 0); loN = 0; else; loN = -1; while ~feas(loN,0) && loN > -1e6; loN = 2*loN; end; end
nu = bisect(@(v) feas(v,0), loN, hiN, 1e-4*max(1, abs(hiN)));
rho = NaN;
if wantRho && nu >= -tolP
    loR = 0; hiR = 1; while feas(0,hiR) && hiR < 1e6; loR = hiR; hiR = 2*hiR; end
    rho = bisect(@(v) feas(0,v), loR, hiR, 1e-4*max(1, hiR));
end
diag = struct('tol', 1e-4, 'exact', false, ...
    'note', 'bisection on the SOS passivity feasibility (certified lower bound on nu)');
if strcmp(cls,'mjls'); method = 'diss-sos-modes'; elseif strcmp(cls,'dwell'); method = 'diss-sos-modes';
else; method = 'diss-sos'; end
res = pack(nu >= -tolP, nu, rho, [], method, cls, diag);
end

function tf = feasible(sys, cls, nu, rho, sol)
sup = struct('type','passivity','nu',nu,'rho',rho);
switch cls
    case 'mjls'
        r = mvst_diss_sos_modes(sys, sup, struct('type','ctmjls'), 'feas', sol);
    case 'dwell'
        Td = getfld(sys, {'Td'}, 1);
        r = mvst_diss_sos_modes(sys, sup, struct('type','switched','Tbar',Td), 'feas', sol);
    otherwise
        r = mvst_diss_sos(sys, sup, struct('type','arbitrary'), 'feas', sol);
end
tf = logical(r.feas);
end

% ------------------------------------------------------------------------
function cls = resolve_class(sys, opts)
if isfield(opts, 'class') && ~isempty(opts.class)
    cls = lower(opts.class);
    if strcmp(cls,'switched'); cls = 'dwell'; end
    if strcmp(cls,'hybrid');   cls = 'sd';    end
    return;
end
A = getfld(sys, {'A0','A'}, []);
multi = (iscell(A) && numel(A) > 1) || (~isempty(A) && size(A,3) > 1) || isfield(sys,'modes');
isdt = (isfield(sys,'domain') && strcmpi(sys.domain,'dt')) || (isfield(sys,'dt') && isequal(sys.dt,true));
if isdt && ~multi
    cls = 'dt';
elseif isfield(sys,'Jf') || isfield(sys,'lamf')
    cls = 'poisson';
elseif isfield(sys,'Tbar') && ~multi
    cls = 'sd';
elseif multi && (isfield(sys,'Pig') || isfield(sys,'Pi'))
    cls = 'mjls';
elseif multi && (isfield(sys,'Td') || (isfield(opts,'Td') && ~isempty(opts.Td)))
    cls = 'dwell';
else
    cls = 'ct';
end
end

function tf = hasnoise(sys)
tf = false;
for f = {'Ac','Dwc'}
    if isfield(sys, f{1}) && ~isempty(sys.(f{1}))
        x = sys.(f{1});
        if iscell(x); tf = tf || any(cellfun(@(y) nnz(y) > 0, x)); else; tf = tf || nnz(x) > 0; end
    end
end
end

function v = getfld(s, names, default)
v = default;
for k = 1:numel(names)
    if isfield(s, names{k}) && ~isempty(s.(names{k})); v = s.(names{k}); return; end
end
end

function v = getopt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name)); v = opts.(name); else; v = default; end
end

function d = maxdim(X, k)
if isempty(X); d = 0; elseif iscell(X); d = size(X{1}, k); else; d = size(X, k); end
end

function c = tocell3(X)
if iscell(X); c = X; elseif isempty(X); c = {};
else; p = size(X,3); c = cell(1,p); for k = 1:p; c{k} = X(:,:,k); end; end
end
