function o = mvst_opts(opts, probe, requires)
%MVST_OPTS  Normalize solver options and verify the SDP/SOS solver is present.
%   O = MVST_OPTS(OPTS) accepts OPTS as [] (or omitted), a plain struct, or an
%   MVST.OPTIONS object, and returns a plain struct O carrying at least the
%   field O.solver (default 'sedumi').  Any extra fields of a struct OPTS are
%   preserved (for example .fid, .psd, .riccati, .class, .Td), so the helper is
%   a safe drop-in normalizer for every engine and dispatcher in the package.
%
%   By default the requested solver is validated against the MATLAB path: if it
%   cannot be resolved the function errors with identifier
%   'mvst:opts:solverMissing' rather than letting a backend fall back silently
%   to a different solver.  Recognized solvers and their probe entry points:
%
%     'sedumi' -> sedumi ,   'sdpt3' -> sdpt3 ,   'mosek' -> mosekopt ,
%     'nosdp'  -> sossolve_nosdp (SOSNOSDP; SOS engines only, needs SOSTOOLS).
%
%   A 'sedumi' resolving to the placeholder sedumi.m that SOSNOSDP's setup adds
%   when no SDP solver is installed is REFUSED ('mvst:opts:sedumiPlaceholder'):
%   it is not a real SeDuMi.
%
%   O = MVST_OPTS(OPTS, false) normalizes without probing the path.  Use this in
%   routers whose backends are purely numerical (no SDP/SOS solve), so that a
%   missing SDP solver does not spuriously block a computation that never needs
%   one; the solver NAME is still validated.
%
%   O = MVST_OPTS(OPTS, PROBE, REQUIRES) additionally checks the dependency of
%   the calling engine (only when PROBE is true): REQUIRES = 'sos' demands
%   SOSTOOLS (sosprogram) and accepts every solver incl. 'nosdp'; REQUIRES =
%   'yalmip' demands YALMIP (sdpvar) and rejects 'nosdp' (error
%   'mvst:opts:nosdpYalmip': the SOSNOSDP backend solves SOS programs only);
%   '' (default) checks nothing beyond the solver.  Errors:
%   'mvst:opts:sostoolsMissing', 'mvst:opts:yalmipMissing'.
%
%   Example:
%       o = mvst_opts(struct('solver','mosek'));      % errors if MOSEK absent
%       o = mvst_opts(mvst.Options('solver','sedumi'));
%       o = mvst_opts(struct('Td',1.0), false);       % accept solver-free opts
%
%   See also MVST.OPTIONS, MVST.SOLVERCHECK.
arguments
    opts        = struct()
    probe (1,1) logical = true
    requires (1,:) char {mustBeMember(requires, {'', 'sos', 'yalmip'})} = ''
end

% ---- reject char/string inputs outright (an empty '' is NOT a valid ------
% ---- "defaults" request: only [] or an omitted argument means defaults) ---
if ischar(opts) || isstring(opts)
    error('mvst:opts:input', ...
          ['OPTS must be a struct, an mvst.Options, or []/omitted for ' ...
           'defaults (got a char/string).']);
end

% ---- normalize to a plain struct with a .solver field --------------------
if isa(opts, 'mvst.Options')
    o = struct('solver', opts.solver, 'reltol', opts.reltol, ...
               'verbosity', opts.verbosity);
elseif isempty(opts)
    o = struct();
elseif isstruct(opts)
    o = opts;
else
    error('mvst:opts:type', ...
          'OPTS must be a struct, an mvst.Options, or empty (got %s).', class(opts));
end
if ~isfield(o, 'solver') || isempty(o.solver)
    o.solver = 'sedumi';
end
o.solver = lower(char(o.solver));

% ---- validate the solver name --------------------------------------------
entries = struct('sedumi', 'sedumi', 'sdpt3', 'sdpt3', 'mosek', 'mosekopt', ...
                 'nosdp', 'sossolve_nosdp');
if ~isfield(entries, o.solver)
    error('mvst:opts:solver', ...
          ['Unknown solver ''%s'' (expected ''sedumi'', ''sdpt3'', ''mosek'', ' ...
           'or ''nosdp'').'], o.solver);
end

% ---- probe the path (no silent fallback) ---------------------------------
if probe
    entry = entries.(o.solver);
    if strcmp(requires, 'sos') && ~mvst_opts_onpath('sosprogram')
        error('mvst:opts:sostoolsMissing', ...
              'This engine needs SOSTOOLS (sosprogram) on the MATLAB path: run mvst_setup.');
    end
    if strcmp(requires, 'yalmip')
        if strcmp(o.solver, 'nosdp')
            error('mvst:opts:nosdpYalmip', ...
                  ['Solver ''nosdp'' (SOSNOSDP) solves SOS programs only; this ' ...
                   'engine uses YALMIP: choose ''sedumi'', ''sdpt3'' or ''mosek''.']);
        end
        if ~mvst_opts_onpath('sdpvar')
            error('mvst:opts:yalmipMissing', ...
                  ['This (legacy) engine needs YALMIP (sdpvar) on the MATLAB path; YALMIP is ' ...
                   'an optional external dependency, not bundled with MVST.']);
        end
    end
    if strcmp(o.solver, 'nosdp') && ~mvst_opts_onpath('sosprogram')
        error('mvst:opts:solverMissing', ...
              'Solver ''nosdp'' needs SOSTOOLS (sosprogram) on the path as well.');
    end
    if strcmp(o.solver, 'sedumi') && mvst_opts_onpath(entry) && mvst_opts_isplaceholder()
        error('mvst:opts:sedumiPlaceholder', ...
              ['The sedumi on the path is the SOSNOSDP placeholder (added when no ' ...
               'SDP solver is installed), not SeDuMi.  Install SeDuMi ' ...
               '(install_sedumi) or select solver ''nosdp''.']);
    end
    if ~mvst_opts_onpath(entry)
        error('mvst:opts:solverMissing', ...
              ['Requested solver ''%s'' is not on the MATLAB path (probe ' ...
               '''%s'' failed).  Install it or choose another solver; no ' ...
               'silent fallback is performed.'], o.solver, entry);
    end
end
end

% -------------------------------------------------------------------------
function tf = mvst_opts_onpath(name)
%MVST_OPTS_ONPATH  True if NAME resolves as a function, script, or MEX file.
tf = (exist(name, 'file') ~= 0) || ~isempty(which(name));
end

% -------------------------------------------------------------------------
function tf = mvst_opts_isplaceholder()
%MVST_OPTS_ISPLACEHOLDER  True if the sedumi on the path is SOSNOSDP's stub.
f = which('sedumi');
tf = false;
if isempty(f) || ~endsWith(f, '.m'), return; end
try
    fid = fopen(f, 'r'); t = fread(fid, 600, '*char')'; fclose(fid);
    tf = contains(t, 'sosnosdp:placeholder');
catch
end
end
