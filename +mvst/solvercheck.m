function report = solvercheck(opts)
%SOLVERCHECK  Detect and probe the SDP / SOS toolchain on the MATLAB path.
%   MVST.SOLVERCHECK() probes for SeDuMi, SDPT3, MOSEK, SOSTOOLS, SOSNOSDP
%   (sossolve_nosdp) and, as an OPTIONAL tool, YALMIP (needed only by the legacy
%   opts.engine = 'yalmip' paths; not bundled), wherever they are on the MATLAB path: the bundled
%   copies under thirdparty/ (added by MVST_SETUP) or your own installation
%   (an installation already on the path wins, MVST_SETUP does not shadow it).
%   Prints a setup report, and (with an output) returns a struct array with
%   one row per tool: fields .name, .found (logical), .version (char),
%   .location (the resolved path of the probed entry point, or ''), .isSolver
%   (logical, true for the three SDP solvers), .probed (logical, true when a
%   trivial-SDP capability probe was attempted) and .capable (logical, true
%   when that probe actually solved).
%
%   Beyond mere path presence, each detected SDP solver is PROBED with a
%   trivial 2x2 semidefinite program (minimize trace(P) subject to P >= I,
%   whose optimum is P = I).  The probe is written as a SOSTOOLS program
%   (degree-0 matrix SOS constraint) and solved with the backend under test, so
%   it exercises the real backend; .capable is set only if the solve returns
%   success and recovers P ~ I.  This reports TRUE capability, not just that the
%   files are on the path.  (The probe needs SOSTOOLS on the path; without it the
%   SDP solvers are marked .probed = false.  YALMIP is not used.)
%
%   R = MVST.SOLVERCHECK() also returns the report struct silently-friendly:
%   the table is still printed unless 'quiet' is set.
%
%   Name/value options:
%     'quiet'  (logical, default false) suppress the printed report.
%
%   This is an onboarding aid for the release repository: run it once after
%   cloning to confirm the numerical backends resolve and actually solve.
%
%   Example:
%       mvst.solvercheck();
%       r = mvst.solvercheck('quiet', true);
%
%   See also MVST.OPTIONS, MVST_OPTS.
arguments
    opts.quiet (1,1) logical = false
end

% name, probe entry point, version callback, isSolver, SOSTOOLS solver tag
tools = {
    'SeDuMi',   'sedumi',    @ver_sedumi,   true,  'sedumi'
    'SDPT3',    'sdpt3',     @ver_sdpt3,    true,  'sdpt3'
    'MOSEK',    'mosekopt',  @ver_mosek,    true,  'mosek'
    'SOSTOOLS', 'sosprogram',@ver_sostools, false, ''
    'SOSNOSDP', 'sossolve_nosdp', @ver_sosnosdp, false, ''
    'YALMIP',   'yalmip',    @ver_yalmip,   false, ''
};

nT = size(tools, 1);
report = repmat(struct('name', '', 'found', false, 'version', '', ...
    'location', '', 'isSolver', false, 'probed', false, ...
    'capable', false, 'probeInfo', '', 'hint', ''), nT, 1);

hasSos = ~isempty(which('sosprogram'));

for k = 1:nT
    name  = tools{k, 1};
    probe = tools{k, 2};
    vfun  = tools{k, 3};
    loc   = which(probe);
    found = ~isempty(loc);
    if strcmp(name, 'SeDuMi') && found && is_placeholder(loc)
        found = false;                      % SOSNOSDP stub, not a real SeDuMi
        loc = '';
    end
    vstr  = '';
    if found
        try
            vstr = vfun();
        catch
            vstr = '';
        end
    end
    report(k).name     = name;
    report(k).found    = found;
    report(k).version  = vstr;
    report(k).location = loc;
    report(k).hint     = '';
    report(k).isSolver = tools{k, 4};

    % capability probe for the SDP solvers that are present
    if report(k).isSolver && found
        if hasSos
            [cap, info] = probe_sdp(tools{k, 5});
            report(k).probed    = true;
            report(k).capable   = cap;
            report(k).probeInfo = info;
            if ~cap
                report(k).hint = install_hint(name);
            end
        else
            report(k).probeInfo = 'not probed (needs SOSTOOLS on path)';
        end
    end
end

if ~opts.quiet
    print_report(report);
end
end

% ------------------------------------------------------------------------
function [capable, info] = probe_sdp(tag)
%PROBE_SDP  Solve the trivial SDP  min trace(P) s.t. P >= I (2x2)  with the backend
%   TAG, posed as a SOSTOOLS program.  CAPABLE is true iff the solve succeeds and
%   the optimizer is P ~ I (trace ~ 2).  INFO carries a short human-readable note.
capable = false; info = '';
try
    pvar x;
    prog = sosprogram(x);
    [prog, P] = sospolymatrixvar(prog, monomials(x, 0), [2 2], 'symmetric');
    prog = sosmatrixineq(prog, P - eye(2));
    prog = sossetobj(prog, trace(P));
    so = struct('solver', tag, 'params', struct('fid', 0));
    [prog, sinfo] = sossolve(prog, so);
    Pv = double(sosgetsol(prog, P));
    okval = ~any(isnan(Pv(:))) && norm(Pv - eye(2), 'fro') < 1e-3;
    capable = (sinfo.pinf == 0) && (sinfo.dinf == 0) && (sinfo.numerr < 2) && okval;
    if capable, info = 'solved'; else, info = sprintf('pinf %g dinf %g numerr %g', sinfo.pinf, sinfo.dinf, sinfo.numerr); end
catch ME
    info = ME.message;
end
end

% ------------------------------------------------------------------------
function h = install_hint(name)
%INSTALL_HINT  Remedy for a solver that is on the path but fails the probe.
switch name
    case 'SeDuMi', h = 'run install_sedumi once (compiles/locates its MEX files).';
    case 'SDPT3',  h = 'run install_sdpt3 once (compiles its MEX files).';
    otherwise,     h = 'check the license and installation.';
end
end

function tf = is_placeholder(f)
%IS_PLACEHOLDER  True if F is SOSNOSDP's stand-in sedumi.m (not a real SeDuMi).
tf = false;
if ~endsWith(f, '.m'), return; end
try
    fid = fopen(f, 'r'); t = fread(fid, 600, '*char')'; fclose(fid);
    tf = contains(t, 'sosnosdp:placeholder');
catch
end
end

% ------------------------------------------------------------------------
function print_report(report)
fprintf('\nMVST solver / SOS toolchain check\n');
fprintf('%s\n', repmat('-', 1, 72));
fprintf('  %-10s %-7s %-8s %-12s %s\n', 'tool', 'found', 'solves', 'version', 'location');
fprintf('%s\n', repmat('-', 1, 72));
for k = 1:numel(report)
    r = report(k);
    mark = 'no'; if r.found; mark = 'yes'; end
    if ~r.isSolver
        smark = '-';
    elseif ~r.probed
        smark = '?';
    elseif r.capable
        smark = 'yes';
    else
        smark = 'no';
    end
    ver = r.version; if isempty(ver); ver = '-'; end
    loc = r.location; if isempty(loc); loc = '(not on path)'; end
    if strcmp(r.name, 'YALMIP'), if r.found, loc = [loc '  [optional]']; else, loc = '(not on path; optional, legacy engines only)'; end, end
    fprintf('  %-10s %-7s %-8s %-12s %s\n', r.name, mark, smark, ver, loc);
end
fprintf('%s\n', repmat('-', 1, 72));
nfound = sum([report.found]);
ncap   = sum([report.capable]);
fprintf('  %d of %d tools detected; %d SDP solver(s) verified to solve.\n', ...
        nfound, numel(report), ncap);
hints = {report.hint};
for k = 1:numel(report)
    if ~isempty(hints{k})
        fprintf('  hint: %s found but does not solve: %s\n', report(k).name, hints{k});
    end
end
fprintf('  solves column:  yes = trivial 2x2 SDP solved,  no = failed,\n');
fprintf('                  ? = not probed (no SOSTOOLS),  - = not an SDP solver.\n\n');
end

% ------------------------------------------------------------------------
%  Per-tool version probes (best-effort; return '' when unavailable).
% ------------------------------------------------------------------------
function v = ver_sedumi()
v = '';
p = fileparts(which('sedumi'));
f = fullfile(p, 'Changelog.md');
if exist(f, 'file')
    t = regexp(fileread(f), '##\s*SeDuMi\s+([\d.]+)', 'tokens', 'once');
    if ~isempty(t), v = t{1}; end
end
if isempty(v) && exist('sedumi_version', 'file')
    try, v = num2str(sedumi_version()); catch, end
end
end

function v = ver_sdpt3()
v = '4.0';   % SDPT3 4.x; the sources carry no programmatic version string
end

function v = ver_sosnosdp()
v = '';
p = fileparts(fileparts(which('sossolve_nosdp')));   % core/ -> package root
f = fullfile(p, 'CITATION.cff');
if exist(f, 'file')
    t = regexp(fileread(f), '(?m)^version:\s*([\d.]+)', 'tokens', 'once');
    if ~isempty(t), v = t{1}; end
end
end

function v = ver_mosek()
v = '';
[~, res] = mosekopt('version echo(0)');
if isfield(res, 'version')
    mv = res.version;
    v = sprintf('%d.%d.%d', mv.major, mv.minor, mv.revision);
end
end

function v = ver_sostools()
v = '';
f = fullfile(fileparts(which('sosprogram')), 'Contents.m');
if exist(f, 'file')
    t = regexp(fileread(f), 'Version\s+([\d.]+)', 'tokens', 'once');
    if ~isempty(t), v = t{1}; end
end
end

function v = ver_yalmip()
v = yalmip('version');
if isnumeric(v); v = num2str(v); end
end
