function res = mvst_perf(sys, arg2, opts)
%MVST_PERF  Generic performance-analysis dispatcher (Chapter "Further
%  Performance Measures", Ch. 11).  One entry point for the performance
%  measures the book states conditions for, routed per system class to the
%  matching exact engine (Riccati / lifted-eigenvalue for the induced gains,
%  matrix-exponential / Lyapunov Gramian for the covariance measures).
%
%  ------------------------------------------------------------------------
%  TWO CALL FORMS
%  ------------------------------------------------------------------------
%
%  (1) DISPATCHER (measure by name):
%
%      RES = MVST_PERF(SYS, MEASURE)
%      RES = MVST_PERF(SYS, MEASURE, OPTS)
%
%  SYS is a thin/rich system struct (or an mvst.System) carrying the drift
%  .A (or .A0), the disturbance input .Bw and the performance output .C, plus
%  the per-class data below.  MEASURE is one of the strings:
%
%    'l2'        hybrid mean-square L2 x l2-gain (H-infinity type), the
%                worst-case energy amplification of Chapter "Dissipativity"
%                (def:diss:msgain).  Delegates verbatim to MVST_L2GAIN, so it
%                inherits every class of that router.
%    'tracegain' hybrid trace gain of def:perf:gains / th:perf:tracegain:
%                the SQUARE of the hybrid L2 x l2-gain.  Computed as
%                MVST_L2GAIN(...)^2 (same class coverage as 'l2').
%    'e2p'       hybrid energy-to-peak (generalized H2) gain, the induced
%                L2 x l2 -> Linf x linf gain of def:perf:e2p.  Exactly
%                sqrt( sup_t lmax( C Pi(t) C' ) ) on the controllability
%                Gramian (th:perf:e2p:seq, fixed sequence; th:perf:e2p:periodic,
%                constant dwell-time).
%    'peak'      hybrid peak gain of def:perf:gains / th:perf:peakgain: the
%                least gamma with H(t,E* Xibar E) <= gamma I on the
%                reachability envelope Xibar, which coincides with the
%                controllability Gramian Pi.  Equal to sup_t lmax( C Pi(t) C' )
%                (= the SQUARE of the energy-to-peak value).
%    'h2'        hybrid H2 cost of def:perf:h2 / th:perf:h2NS: the average
%                output power sqrt( (1/Tbar) int_0^Tbar tr( C Pi(tau) C' ) ).
%                Reduces to the classical H2 norm in the jump-free case
%                (prop:perf:h2classic); periodic form th:perf:h2:periodic.
%    'p2p'       deterministic peak-to-peak (star-norm, induced
%                Linf x linf -> Linf x linf) of Section sec:perf:ptpstar /
%                th:perf:ptp:seq / th:perf:ptp:dwell.  The decay-rate-certified
%                invariant-set LMIs, minimized over (P,gamma^2) for fixed
%                (lambda,nu) and BISECTED over (lambda,nu) (backend
%                MVST_PERF_DWELL).  Carries the DC-gain lower bound companion
%                res.gammaLB = ||C A^{-1} Bw||_2 (prop:perf:ptplower).
%    'peakfull'  amplitude-bounded peak, full-block dual certificate of
%                prop:perf:peakgainfull: the same (lambda,nu) LMIs with a linear
%                output term (backend MVST_PERF_DWELL).
%
%   DWELL-TIME FAMILIES.  For 'e2p','peak','h2', OPTS.dwell in {'min','range'}
%   selects the genuinely timer-dependent supersolution tube / budget
%   certificates (prop:perf:peakgaindwell, th:perf:e2p:MDT/RDT,
%   th:perf:h2:MDT/RDT, prop:perf:ptpint:dwell), solved by the SOS-timer backend
%   MVST_PERF_DWELL; OPTS.dwell absent (or 'constant') keeps the exact steady /
%   periodic Gramian engine.
%
%  CLASS COVERAGE (which class-measure pairs the book states conditions for
%  and this dispatcher computes; every other pair raises
%  'mvst:perf:unsupported' naming the nearest supported pair):
%
%      measure   | ct  mjls  switched  impulsive/sd  MDT/RDT dwell
%      ----------+---------------------------------------------------------
%      l2        |  S   S      S          S           (MVST_L2GAIN)
%      tracegain |  S   S      S          S           (= l2^2)
%      e2p       |  S   -      -          S            S  (Gramian; SOS tube)
%      peak      |  S   -      -          S            S  (Gramian; SOS tube)
%      h2        |  S   -      -          S            S  (Gramian; SOS budget)
%      p2p       |  S   -      -          S            S  ((lambda,nu) LMI + LB)
%      peakfull  |  S   -      -          S            S  (full-block (lambda,nu))
%
%  ('S' under the MDT/RDT column marks the timer-dependent minimum/range
%  dwell-time tube (peak/e2p), budget (h2), or invariant-set (p2p/peakfull)
%  certificate, backend MVST_PERF_DWELL; the multi-mode (mjls, switched) and
%  stochastic (poisson) covariance measures remain unsupported for the Gramian
%  measures, routing to 'l2'/'tracegain' via MVST_L2GAIN.)
%
%  The Gramian measures ('e2p','peak','h2') are book-characterized for the
%  time-dependent / continuous-time system (ct) and its constant dwell-time
%  impulsive/sampled instance; the multi-mode (mjls, switched) and stochastic
%  (poisson) covariance measures, and the minimum/range dwell-time tube
%  certificates (th:perf:e2p:MDT/RDT, th:perf:h2:MDT/RDT), are SDP/DLMI
%  problems not provided by the exact engines and are reported unsupported.
%
%  RES (uniform):
%    .gamma    the measure value (l2: gain; tracegain: gain^2; e2p: gain;
%              peak: peak gain = lmax; h2: cost).
%    .measure  the resolved measure name.
%    .class    the resolved class tag ('ct'|'impulsive' for the Gramian
%              measures, or the MVST_L2GAIN tag for 'l2'/'tracegain').
%    plus measure-specific extras: 'l2'/'tracegain' carry nothing further
%    ('tracegain' also .gammaL2); the Gramian measures carry .Pi (steady /
%    periodic-reset Gramian) and .H2sq / .e2psq / .p2psq.
%
%  OPTS (dispatcher): .class forces the class; .Tbar / .J supply the period
%  and reset for the impulsive Gramian measures; .npts / .nprofper / .profpts
%  tune the Gramian quadrature (see below); any field is also forwarded to
%  MVST_L2GAIN for the induced gains.
%
%  ------------------------------------------------------------------------
%  (2) LEGACY GRAMIAN BATCH (backward compatible):
%
%      RES = MVST_PERF(SYS, OUT, OPTS)
%
%  When the second argument is a STRUCT (not a measure name) the historical
%  batch interface is used unchanged: it returns ALL Gramian-based measures
%  of the constant-dwell-time impulsive system at once, with an optional
%  white-noise Monte-Carlo cross-check.  The system is
%       dx = (A x + Bw w) dt   on (0,Tbar),      x+ = J x   at each period,
%       z  = C x,
%  driven by unit-intensity white disturbance w, with periodic controllability
%  Gramian Pi(tau) and output Gramian S(tau) = C Pi(tau) C'.  The measures are
%       H2^2           = (1/Tbar) int_0^Tbar tr S(tau) dtau      (avg power),
%       energy2peak^2  = max_tau lambda_max( S(tau) )            (gen. H2),
%       peak2peak^2    = max_tau tr S(tau)                       (peak power).
%  The peak2peak returned is the INTENSITY peak-to-peak (peak output second
%  moment, the reachability envelope), which is Gramian-computable; the
%  DETERMINISTIC peak-to-peak (star-norm) is not (see the 'p2p' note above).
%
%  SYS fields (legacy):  .A (or .A0), .J (default I), .Tbar (required).
%  OUT fields (legacy):  .Bw (required), .C (required).
%  OPTS fields (legacy): .npts (200), .nprofper (1), .profpts (120), .mc
%                        (false; logical or struct with .seed/.dt/.ncyc/.nreal).
%  RES fields (legacy):  .H2 .energy2peak .peak2peak .H2sq .e2psq .p2psq
%                        .Pi .Pifun .Wc .M .tau .power .lambda .profile .mc.
%
%  See also MVST_L2GAIN, MVST_POS_STABILITY, MVST_POS_GAIN, MVST_DISS_SOS.

if nargin < 3 || isempty(opts); opts = struct(); end
if nargin < 2
    error('mvst:perf:nargin', ...
        'MVST_PERF needs at least a system and a measure name (or OUT struct).');
end

% ---- dispatch on the shape of the second argument ----
if ischar(arg2) || (isstring(arg2) && isscalar(arg2))
    res = perf_dispatch(sys, lower(char(arg2)), opts);
else
    res = perf_gramian_batch(sys, arg2, opts);   % legacy interface (OUT struct)
end
end

% =========================================================================
%  NEW DISPATCHER
% =========================================================================
function res = perf_dispatch(sys, measure, opts)
sys = normalize_sys(sys);

switch measure
    case {'l2', 'tracegain'}
        res = perf_induced(sys, measure, opts);

    case {'e2p', 'peak', 'h2'}
        % Genuinely timer-dependent dwell families (minimum / range) route to the
        % SOS-timer tube / budget backends (prop:perf:peakgaindwell,
        % th:perf:e2p:MDT/RDT, th:perf:h2:MDT/RDT); the constant-dwell / ct
        % instances stay on the exact Gramian engine.
        dwl = lower(mvst_getf(opts,'dwell',''));
        if any(strcmp(dwl, {'min','range'}))
            res = mvst_perf_dwell(sys, measure, opts);
        else
            res = perf_gramian_measure(sys, measure, opts);
        end

    case 'p2p'
        % Deterministic peak-to-peak (star-norm) via the (lambda,nu)-bisected
        % invariant-set LMIs (th:perf:ptp:seq / th:perf:ptp:dwell), with the
        % DC-gain lower bound companion (prop:perf:ptplower).
        res = mvst_perf_dwell(sys, 'p2p', opts);

    case 'peakfull'
        % Amplitude-bounded peak, full-block dual certificate
        % (prop:perf:peakgainfull), same (lambda,nu) LMIs, linear output.
        res = mvst_perf_dwell(sys, 'peakfull', opts);

    otherwise
        error('mvst:perf:measure', ...
            ['Unknown measure ''%s''.  Valid measures: ' ...
             '''l2'', ''tracegain'', ''e2p'', ''peak'', ''h2'', ''p2p'', ' ...
             '''peakfull''.'], ...
            measure);
end
end

% -------------------------------------------------------------------------
%  Induced gains: 'l2' (delegate) and 'tracegain' (= l2^2)
% -------------------------------------------------------------------------
function res = perf_induced(sys, measure, opts)
% Strip the Gramian-only tuning fields before forwarding to MVST_L2GAIN.
lopts = rmfields(opts, {'npts', 'nprofper', 'profpts', 'mc'});
r = mvst_l2gain(sys, lopts);          % RES.gamma, RES.method

res = struct();
res.measure = measure;
res.class   = r.method;
if strcmp(measure, 'l2')
    res.gamma = r.gamma;              % the hybrid L2 x l2-gain
else                                  % 'tracegain'
    res.gammaL2 = r.gamma;
    res.gamma   = r.gamma^2;          % th:perf:tracegain: trace gain = L2-gain^2
end
end

% -------------------------------------------------------------------------
%  Covariance (Gramian) measures: 'e2p', 'peak', 'h2'
% -------------------------------------------------------------------------
function res = perf_gramian_measure(sys, measure, opts)
cls = resolve_gramian_class(sys, opts, measure);

A  = mvst_getf(sys, 'A', mvst_getf(sys, 'A0', []));
if isempty(A)
    error('mvst:perf:A', 'SYS must carry a drift .A (or .A0).');
end
A  = A(:,:,1);
n  = size(A, 1);
Bw = mvst_getf(sys, 'Bw', []);
C  = mvst_getf(sys, 'C',  []);
if isempty(Bw) || isempty(C)
    error('mvst:perf:io', ...
        'The Gramian measure ''%s'' needs the disturbance .Bw and output .C.', ...
        measure);
end
C  = C(:,:,1);

switch cls
    case 'ct'
        % Jump-free steady Gramian: He(A Pi) + Bw Bw' = 0 (prop:perf:h2classic).
        if max(real(eig(A))) >= 0
            error('mvst:perf:unstable', ...
                'The ct Gramian measure needs a Hurwitz drift .A.');
        end
        Pi = sylvester(A, A.', -Bw*Bw.');   % A Pi + Pi A' + Bw Bw' = 0 (base MATLAB, no Control System Toolbox)
        Pi = (Pi + Pi.')/2;
        S  = C*Pi*C.';  S = (S + S.')/2;
        H2sq  = trace(S);              % constant in time -> its own average
        e2psq = max(real(eig(S)));
        p2psq = trace(S);
        extra = struct('Pi', Pi);

    case 'impulsive'
        Tbar = mvst_getf(sys, 'Tbar', mvst_getf(opts, 'Tbar', []));
        if isempty(Tbar)
            error('mvst:perf:Tbar', ...
                'The impulsive Gramian measure needs the period .Tbar.');
        end
        J = mvst_getf(sys, 'J', mvst_getf(opts, 'J', eye(n)));
        g = gramian_core(A, J, Tbar, Bw, C, opts);
        H2sq = g.H2sq;  e2psq = g.e2psq;  p2psq = g.p2psq;
        extra = struct('Pi', g.Pi, 'Pifun', g.Pifun, 'tau', g.tau, ...
                       'power', g.power, 'lambda', g.lambda);
end

res = struct();
res.measure = measure;
res.class   = cls;
res.H2sq = H2sq;  res.e2psq = e2psq;  res.p2psq = p2psq;
switch measure
    case 'h2'
        res.gamma = sqrt(H2sq);        % avg-power H2 cost
    case 'e2p'
        res.gamma = sqrt(e2psq);       % energy-to-peak gain (th:perf:e2p:seq)
    case 'peak'
        res.gamma = e2psq;             % hybrid peak gain (th:perf:peakgain)
end
fn = fieldnames(extra);
for k = 1:numel(fn); res.(fn{k}) = extra.(fn{k}); end
end

% -------------------------------------------------------------------------
%  Class resolution for the Gramian measures.  Only 'ct' and 'impulsive' are
%  book-characterized (th:perf:e2p:*, th:perf:h2:*).  Everything else routes
%  to 'mvst:perf:unsupported' naming the nearest supported pair.
% -------------------------------------------------------------------------
function cls = resolve_gramian_class(sys, opts, measure)
if isfield(opts, 'class') && ~isempty(opts.class)
    cls = lower(char(opts.class));
    if any(strcmp(cls, {'sd', 'sampled', 'sampleddata'})); cls = 'impulsive'; end
    if ~any(strcmp(cls, {'ct', 'impulsive'}))
        unsupported(cls, measure);
    end
    return;
end

A = mvst_getf(sys, 'A', mvst_getf(sys, 'A0', []));
multi = (iscell(A) && numel(A) > 1) || (~isempty(A) && size(A,3) > 1) ...
        || isfield(sys, 'modes');

if isfield(sys, 'Jf') || isfield(sys, 'lamf')
    unsupported('poisson', measure);
elseif multi && (isfield(sys, 'Pig') || isfield(sys, 'Pi'))
    unsupported('mjls', measure);
elseif multi || isfield(sys, 'Td')
    unsupported('switched', measure);
elseif isfield(sys, 'J') || isfield(sys, 'Tbar') ...
        || (isfield(opts, 'Tbar') && ~isempty(opts.Tbar))
    cls = 'impulsive';
else
    cls = 'ct';
end
end

% -------------------------------------------------------------------------
function unsupported(cls, measure)
error('mvst:perf:unsupported', ...
    ['Measure ''%s'' is not covered for the ''%s'' class: the book states ' ...
     'no Gramian (covariance) condition for it (the induced-gain theorems ' ...
     'cover this class instead).  Nearest supported pair: ''tracegain'' or ' ...
     '''l2'' on the ''%s'' class (both delegate to MVST_L2GAIN).'], ...
    measure, cls, cls);
end

% =========================================================================
%  LEGACY BATCH INTERFACE (returns all measures + optional Monte-Carlo)
% =========================================================================
function res = perf_gramian_batch(sys, out, opts)
A  = mvst_getf(sys,'A', mvst_getf(sys,'A0',[]));
if isempty(A); error('mvst:perf:A','SYS must carry .A (or .A0).'); end
n  = size(A,1);
J  = mvst_getf(sys,'J', eye(n));
Tbar = mvst_getf(sys,'Tbar',[]);
if isempty(Tbar); error('mvst:perf:Tbar','SYS.Tbar (period) is required.'); end
Bw = mvst_getf(out,'Bw',[]);
C  = mvst_getf(out,'C',[]);
if isempty(Bw) || isempty(C)
    error('mvst:perf:io','OUT must carry .Bw and .C.');
end

g = gramian_core(A, J, Tbar, Bw, C, opts);

res = struct();
res.H2          = sqrt(g.H2sq);
res.energy2peak = sqrt(g.e2psq);
res.peak2peak   = sqrt(g.p2psq);
res.H2sq  = g.H2sq;   res.e2psq = g.e2psq;   res.p2psq = g.p2psq;
res.Pi    = g.Pi;     res.Pifun = g.Pifun;
res.Wc    = g.Wc;     res.M = g.M;
res.tau   = g.tau;    res.power = g.power;    res.lambda = g.lambda;
res.profile = g.profile;

% ---- optional Monte-Carlo cross-check of the H2 cost ----
mc = mvst_getf(opts,'mc',false);
res.mc = [];
if ~isequal(mc,false)
    if ~isstruct(mc); mc = struct(); end
    seed  = mvst_getf(mc,'seed',1);
    dt    = mvst_getf(mc,'dt',5e-3);
    ncyc  = mvst_getf(mc,'ncyc',400);
    nreal = mvst_getf(mc,'nreal',60);
    h2 = montecarlo_power(A, Bw, C, J, Tbar, seed, dt, ncyc, nreal);
    res.mc = struct('h2',h2,'seed',seed,'dt',dt,'ncyc',ncyc,'nreal',nreal);
end
end

% =========================================================================
%  SHARED GRAMIAN CORE (identical numerics to the historical implementation)
% =========================================================================
function g = gramian_core(A, J, Tbar, Bw, C, opts)
npts     = mvst_getf(opts,'npts',200);
nprofper = mvst_getf(opts,'nprofper',1);
profpts  = mvst_getf(opts,'profpts',120);

Q = Bw*Bw.';                                     % disturbance intensity

% ---- one-period flow Gramian and monodromy ----
Wc = vanloan_gram(A, Q, Tbar);
M  = J*expm(A*Tbar);

% ---- steady periodic Gramian at the reset instant: Pi0 - M Pi0 M' = J Wc J' ---
Pi0 = dlyap_kron(M, J*Wc*J.');

% ---- propagate Pi(tau) over one period ----
Pifun = @(tau) expm(A*tau)*Pi0*expm(A.'*tau) + vanloan_gram(A, Q, tau);

tau = linspace(0, Tbar, npts);
power  = zeros(size(tau));
lambda = zeros(size(tau));
for i = 1:numel(tau)
    S = C*Pifun(tau(i))*C.';
    power(i)  = trace(S);
    lambda(i) = max(real(eig(S)));
end

g = struct();
g.H2sq  = trapz(tau, power)/Tbar;    % average output power
g.e2psq = max(lambda);               % peak eigenvalue along the Gramian
g.p2psq = max(power);                % peak output power (intensity ptp)
g.Pi    = Pi0;  g.Pifun = Pifun;  g.Wc = Wc;  g.M = M;
g.tau   = tau;  g.power = power;   g.lambda = lambda;

% ---- output-power profile over several periods (for plotting) ----
tt = []; pw = []; Pcur = Pi0;
for c = 1:nprofper
    ts = linspace(0, Tbar, profpts);
    for i = 1:numel(ts)
        Pit = expm(A*ts(i))*Pcur*expm(A.'*ts(i)) + vanloan_gram(A, Q, ts(i));
        tt(end+1) = (c-1)*Tbar + ts(i);   %#ok<AGROW>
        pw(end+1) = trace(C*Pit*C.');      %#ok<AGROW>
    end
    PT = expm(A*Tbar)*Pcur*expm(A.'*Tbar) + Wc;  Pcur = J*PT*J.';   % jump
end
g.profile = struct('t', tt, 'power', pw);
end

% =========================================================================
function s = normalize_sys(sys)
%NORMALIZE_SYS  Accept an mvst.System or a thin/rich struct, returning a
%  struct carrying the drift and the I/O blocks at the top level.
if isa(sys, 'mvst.System')
    r = sys.rich();                       % .A0 naming, may drop the I/O blocks
    t = sys.thin();
    if isstruct(t) && isfield(t, 'modes'); m1 = t.modes{1};
    elseif iscell(t);                      m1 = t{1};
    else;                                  m1 = t;
    end
    for f = {'A','A0','Bw','C','Dw','J','Tbar','Td'}
        if isstruct(m1) && isfield(m1, f{1}) && ~isfield(r, f{1})
            r.(f{1}) = m1.(f{1});
        end
    end
    s = r;
elseif isstruct(sys)
    s = sys;
else
    error('mvst:perf:input', ...
        'SYS must be a system struct or an mvst.System object.');
end
end

% =========================================================================
function s = rmfields(s, names)
%RMFIELDS  Remove any of NAMES present in the struct S (no error if absent).
for k = 1:numel(names)
    if isfield(s, names{k}); s = rmfield(s, names{k}); end
end
end

% =========================================================================
function W = vanloan_gram(A, Q, T)
%VANLOAN_GRAM  Finite-horizon Gramian W = int_0^T e^{As} Q e^{A's} ds by van
%  Loan's method: exponentiate the 2n x 2n block matrix [-A Q; 0 A'] and read
%  the Gramian off the upper-right and lower-right blocks of the result.
n = size(A,1);
F = [ -A, Q; zeros(n), A.' ];
G = expm(F*T);
W = G(n+1:end, n+1:end).' * G(1:n, n+1:end);
W = (W + W.')/2;
end

% =========================================================================
function X = dlyap_kron(M, Q)
%DLYAP_KRON  Solve the discrete Lyapunov equation X - M X M' = Q by the
%  Kronecker (vectorized) form vec(X) = (I - conj(M) kron M)^{-1} vec(Q).
n = size(M,1);
X = reshape((eye(n*n) - kron(conj(M),M)) \ Q(:), n, n);
X = (X + X.')/2;
end

% =========================================================================
function h2 = montecarlo_power(A, Bw, C, J, Tbar, seed, dt, ncyc, nreal)
%MONTECARLO_POWER  Average output second moment E|z|^2 under discrete-time
%  white noise, by Euler-Maruyama with the EXACT per-step disturbance
%  covariance (van Loan Gramian over one step), averaged over the second half
%  of each realization to discard the transient.  The RNG is seeded here so
%  the cross-check is deterministic and reproducible.
rng(seed);
n = size(A,1);
Ad = expm(A*dt);
Bd = vanloan_gram(A, Bw*Bw.', dt);
Lc = chol(Bd + 1e-12*eye(n),'lower');
acc = 0; cnt = 0; steps = round(Tbar/dt);
for r = 1:nreal
    x = zeros(n,1);
    for c = 1:ncyc
        for s = 1:steps
            x = Ad*x + Lc*randn(n,1);       % exact per-step covariance
            if c > ncyc/2
                z = C*x; acc = acc + z.'*z; cnt = cnt + 1;
            end
        end
        x = J*x;                            % jump
    end
end
h2 = acc/cnt;
end
