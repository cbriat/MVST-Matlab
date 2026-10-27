function mvst_runscripts(scripts)
%MVST_RUNSCRIPTS  Run a list of example scripts with a fresh workspace each.
%   MVST_RUNSCRIPTS(SCRIPTS) adds the library and examples folders to the
%   path and runs each script in SCRIPTS (a cell array of script names,
%   without extension) in this function's workspace, clearing leftover
%   variables between scripts so that no script inherits stale graphics
%   handles or mismatched variables from the previous one.  Shared core of
%   the per-chapter runners (run_stability, run_oc, ...) and of run_all.
%
%   See also RUN_ALL, RUN_STABILITY, RUN_OC.
libdir = fileparts(mfilename('fullpath'));
addpath(libdir);
addpath(fullfile(libdir,'..','examples'));
for s = scripts
    fprintf('=== %s ===\n', s{1});
    clearvars -except s scripts libdir
    run(fullfile(libdir,'..','examples',[s{1} '.m']));
end
libdir = fileparts(mfilename('fullpath'));
fprintf('Done. Figures are in %s\n', ...
    mvst_figdir());
end
