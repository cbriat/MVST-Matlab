function res = mvst_stab_critical(sys, param, range, opts)
%MVST_STAB_CRITICAL  Critical parameter value on the mean-square stability margin.
%   RES = MVST_STAB_CRITICAL(SYS,PARAM,RANGE) finds, by bisection, the value of
%   a scalar PARAM in the interval RANGE=[lo hi] at which the exact mean-square
%   stability margin of SYS (computed by MVST_STABILITY) crosses the stability
%   threshold, i.e. the CRITICAL parameter separating the stable and unstable
%   regimes.  It is the pointwise counterpart of the stability-vs-parameter
%   curves MVST_STABILITY returns over a grid: where the curve returns the whole
%   margin sweep, this routine pins the single zero crossing.
%
%   PARAM is the char name of the swept scalar.  By default it is forwarded as
%   an OPTS field to MVST_STABILITY, so the natural use is the dwell-time /
%   sampling period 'T' of the dwell-dependent classes ('switched','sampled',
%   'impulsive'):
%
%       res = mvst_stab_critical(S, 'T', [0.05 3]);   % critical dwell-time
%
%   For a genuinely arbitrary parameter (a noise scale, a coupling rate, ...)
%   supply a builder OPTS.setter = @(sys,x) -> sys that returns the system
%   rebuilt at parameter value X; PARAM is then only a label and the margin is
%   read on setter(SYS,X):
%
%       o.setter = @(S,c) mvst.System.ct([0 1;-2 -1], zeros(2,0), 'Ac', c*eye(2));
%       res = mvst_stab_critical(S, 'noise', [0 2], o);
%
%   RES = MVST_STAB_CRITICAL(SYS,PARAM,RANGE,OPTS) accepts:
%     OPTS.tol     stability threshold on the margin (default 0): the crossing
%                  solves MARGIN(x) = OPTS.tol.  Forwarded to MVST_STABILITY too.
%     OPTS.class   force the stability class (forwarded to MVST_STABILITY).
%     OPTS.setter  builder @(sys,x) -> sys (see above); when present PARAM is a
%                  label only and is NOT written into the MVST_STABILITY options.
%     OPTS.maxit   maximum bisection iterations (default 80).
%     OPTS.xtol    absolute tolerance on the bracket width (default 1e-8).
%     any other field is forwarded to MVST_STABILITY unchanged.
%
%   RES has fields
%     .crit        the critical parameter value (bracket midpoint at convergence);
%     .margin_at   the stability margin MVST_STABILITY reports at .crit;
%     .iterations  the number of bisection steps taken;
%     .bracket     the final [lo hi] bracket;
%     .param       the swept parameter name;
%     .tol         the margin threshold used.
%
%   A missing sign change over RANGE (both endpoints on the same side of the
%   threshold, so no crossing is bracketed) raises
%   'mvst:stab_critical:nobracket'.
%
%   ASSUMPTION.  The margin is assumed to cross the threshold exactly ONCE
%   (monotonically) on RANGE, as it does for the standard dwell-time / noise-
%   scale sweeps this wrapper targets.  Only the two endpoints are tested for
%   the initial bracket, so if the margin has SEVERAL crossings inside RANGE
%   the routine converges to an ARBITRARY one of them (whichever the bisection
%   happens to bracket), and if it has an EVEN number of crossings the endpoints
%   share a sign and 'mvst:stab_critical:nobracket' is raised even though
%   crossings exist.  For a non-monotone or multi-crossing margin, first locate
%   the crossings with an MVST_STABILITY sweep and pass a RANGE that isolates the
%   one of interest.
%
%   See also MVST_STABILITY, MVST_STABILIZE, MVST.SYSTEM.

if nargin < 4 || isempty(opts); opts = struct(); end
if ~(ischar(param) || (isstring(param) && isscalar(param)))
    error('mvst:stab_critical:input', 'PARAM must be a parameter name (char).');
end
param = char(param);
if ~(isnumeric(range) && numel(range) == 2 && all(isfinite(range)))
    error('mvst:stab_critical:input', 'RANGE must be a finite [lo hi] pair.');
end
if ~isstruct(opts)
    error('mvst:stab_critical:input', 'OPTS must be a struct.');
end
lo = double(range(1)); hi = double(range(2));
if ~(lo < hi)
    error('mvst:stab_critical:input', 'RANGE must satisfy lo < hi.');
end

tol   = getf(opts, 'tol',   0);
maxit = getf(opts, 'maxit', 80);
xtol  = getf(opts, 'xtol',  1e-8);
setter = [];
if isfield(opts, 'setter') && ~isempty(opts.setter)
    setter = opts.setter;
    if ~isa(setter, 'function_handle')
        error('mvst:stab_critical:input', 'OPTS.setter must be a function handle.');
    end
end

% Options forwarded to MVST_STABILITY: strip this routine's control fields; the
% swept parameter is written per evaluation (unless a setter rebuilds SYS).
stabopts = rmfields(opts, {'setter', 'maxit', 'xtol'});

    function mrg = margin_at(x)
        o = stabopts;
        if isempty(setter)
            o.(param) = x;
            mrg = mvst_stability(sys, o).margin;
        else
            mrg = mvst_stability(setter(sys, x), o).margin;
        end
        if ~isscalar(mrg)
            error('mvst:stab_critical:input', ...
                ['MVST_STABILITY returned a non-scalar margin for PARAM=%s; ', ...
                 'sweep a SCALAR parameter (RANGE endpoints are scalars).'], param);
        end
    end

flo = margin_at(lo) - tol;
fhi = margin_at(hi) - tol;

% A genuine bracket needs opposite signs (a zero exactly at an endpoint counts).
if flo * fhi > 0
    error('mvst:stab_critical:nobracket', ...
        ['No stability-margin crossing is bracketed on [%.6g, %.6g] for ', ...
         'PARAM ''%s'': the margin is %.6g at lo and %.6g at hi (both %s the ', ...
         'threshold %.6g).  Widen RANGE.'], lo, hi, param, flo + tol, fhi + tol, ...
         ternary(flo > 0, 'above', 'below'), tol);
end

% Bisection maintaining opposite signs at the bracket ends.
iters = 0;
while (hi - lo) > xtol && iters < maxit
    mid  = 0.5 * (lo + hi);
    fmid = margin_at(mid) - tol;
    if fmid == 0
        lo = mid; hi = mid; break;
    elseif flo * fmid <= 0
        hi = mid;                 % keep [lo,hi] straddling; flo stays the anchor
    else
        lo = mid; flo = fmid;
    end
    iters = iters + 1;
end

crit = 0.5 * (lo + hi);
res  = struct('crit', crit, 'margin_at', margin_at(crit), ...
              'iterations', iters, 'bracket', [lo hi], ...
              'param', param, 'tol', tol);
end

% ------------------------------------------------------------------------
function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end

% ------------------------------------------------------------------------
function s = rmfields(s, fields)
for k = 1:numel(fields)
    if isfield(s, fields{k}); s = rmfield(s, fields{k}); end
end
end

% ------------------------------------------------------------------------
function v = ternary(c, a, b)
if c; v = a; else; v = b; end
end
