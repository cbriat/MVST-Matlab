function res = mvst_stabilize_selftest()
%MVST_STABILIZE_SELFTEST  Validation of the MVST_STABILIZE SOS delegation.
%   Confirms that the discrete-time and dwell-time-constrained impulsive classes,
%   which have no closed-form Riccati route, are now DELEGATED to MVST_STABZ_SOS
%   (rather than raising 'mvst:stabilize:unsupported') and that the result is
%   flagged with .method = 'sos' and carries a usable timer gain.  Needs an SDP
%   backend (SeDuMi + SOSTOOLS), so it is skipped (with a message) when the solvers
%   are missing; the test class TSELFTESTS wraps it.  Returns a struct array (.name/.ok) and
%   prints a table.
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
addpath(root); addpath(fullfile(root,'lib'));
if exist('sossolve','file') ~= 2 || exist('sedumi','file') == 0
    fprintf('%s: SKIPPED (SOSTOOLS + SeDuMi required; run mvst_setup).\n', mfilename);
    res = [];
    return
end

R = struct('name',{},'ok',{},'detail',{});

% ---- (1) impulsive: minimum dwell-time SOS synthesis --------------------
S = mvst.System.impulsive([0 1; -1 0.2], [0; 1], 0.6*eye(2), 'Q', eye(2), 'R', 1);
try
    r = S.stabilize(struct('T', 0.5));
    ok = strcmp(r.method,'sos') && isfield(r,'Kfun') && isa(r.Kfun,'function_handle');
    R(end+1) = row('impulsive -> sos delegation', ok, ...
                   sprintf('method=%s feas=%d', r.method, logical(r.feas)));
catch ME
    R(end+1) = row('impulsive -> sos delegation', false, ME.message);
end

% ---- (2) discrete-time single mode: forced class 'dt' -------------------
Sdt = mvst.System.fromstruct(mvst_mode([0.9 0.2; 0 0.8], [0;1], eye(2), 1));
try
    r = mvst_stabilize(Sdt, struct('class','dt'));
    ok = strcmp(r.method,'sos') && isfield(r,'Kfun');
    R(end+1) = row('dt -> sos delegation', ok, sprintf('method=%s', r.method));
catch ME
    R(end+1) = row('dt -> sos delegation', false, ME.message);
end

% ---- (3) the old unsupported error must be gone -------------------------
try
    mvst_stabilize(S, struct('class','impulsive','T',0.5));
    R(end+1) = row('no mvst:stabilize:unsupported', true, 'delegated');
catch ME
    ok = ~strcmp(ME.identifier,'mvst:stabilize:unsupported');
    R(end+1) = row('no mvst:stabilize:unsupported', ok, ME.identifier);
end

% ---- print --------------------------------------------------------------
fprintf('\n==== mvst_stabilize self-test ====\n');
fprintf('%-34s %5s   %s\n','case','ok','detail');
for k = 1:numel(R)
    fprintf('%-34s %5d   %s\n', R(k).name, R(k).ok, R(k).detail);
end
allok = all([R.ok]);
fprintf('---- %s ----\n', ternary(allok,'ALL OK','FAILURES'));
if nargout == 0; clear res; else; res = R; end
end

function r = row(name, ok, detail)
r = struct('name',name,'ok',logical(ok),'detail',detail);
end

function v = ternary(c,a,b); if c; v=a; else; v=b; end; end
