function run_all_tests
%RUN_ALL_TESTS  Consolidated test entry point for the MVST companion library.
%   Runs every test suite in tests/ (the classdef matlab.unittest classes and
%   the script-based tests) plus the golden regression harness, prints one
%   consolidated PASS/FAIL/SKIPPED summary.  Tests skipped through assumptions
%   (a missing solver) are reported as SKIPPED and are not failures; the
%   function errors (nonzero -batch exit) only if a test truly FAILS.  Classdef suites are discovered automatically, so a new
%   t*.m TestCase is picked up without editing this file.
%
%   From the matlab/ directory:
%       /Applications/MATLAB_R2025a.app/bin/matlab -batch "cd tests; run_all_tests"
%
%   See also RUN_REGRESS, MVST_REGRESS, MVST_SETUP.

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);                     % matlab/
    addpath(root);                              % +mvst package parent
    addpath(fullfile(root, 'lib'));             % engine layer
    addpath(fullfile(root, 'bookutil'));        % figure/style helpers
    addpath(here);                              % expose the test files

    fprintf('\n==================== MVST test battery ====================\n');

    totalPass = 0; totalFail = 0; totalSkip = 0;
    rows = cell(0, 4);                           % {name, npass, nfail, status}

    % ---- 1. classdef matlab.unittest suites (auto-discovered) ---------- %
    classFiles = dir(fullfile(here, 't*.m'));
    scriptTests = {};                            % script-based ones deferred
    for k = 1:numel(classFiles)
        name = classFiles(k).name(1:end-2);      % strip the .m
        txt  = fileread(fullfile(here, classFiles(k).name));
        if ~contains(txt, 'classdef')
            scriptTests{end+1} = name;           %#ok<AGROW>
            continue;
        end
        fprintf('\n>>> running suite %s ...\n', name);
        res = runtests(name);
        np  = sum([res.Passed]);
        nf  = sum([res.Failed]);
        ni  = sum([res.Incomplete]);             % assumption-skips (missing solver etc.)
        totalPass = totalPass + np; totalFail = totalFail + nf; totalSkip = totalSkip + ni;
        if nf > 0
            st = 'FAIL';
        elseif ni > 0
            st = sprintf('PASS (%d SKIPPED)', ni);
        else
            st = 'PASS';
        end
        rows(end+1, :) = {name, np, nf, st}; %#ok<AGROW>
        if ni > 0
            for q = find([res.Incomplete])
                fprintf('  SKIPPED: %s\n', res(q).Name);
            end
        end
    end

    % ---- 2. script-based tests (run top-to-bottom, error on failure) --- %
    for k = 1:numel(scriptTests)
        nm = scriptTests{k};
        fprintf('\n>>> running script test %s ...\n', nm);
        try
            runScript(fullfile(here, [nm '.m']));
            totalPass = totalPass + 1;
            rows(end+1, :) = {nm, 1, 0, 'PASS'}; %#ok<AGROW>
        catch ME
            totalFail = totalFail + 1;
            rows(end+1, :) = {nm, 0, 1, 'FAIL'}; %#ok<AGROW>
            fprintf(2, '  %s FAILED: %s\n', nm, ME.message);
        end
    end

    % ---- 3. golden regression harness --------------------------------- %
    fprintf('\n>>> running golden regression run_regress ...\n');
    %   Excluded (with a message) when the companion-document data macros are
    %   absent; regenerate them with mvst_regress('Regenerate', true).
    if exist('mvst_regress','file') ~= 2
        rep = struct('skipped', true);
    else
        liveDir = fullfile(fileparts(root), 'figures', 'data');   % book figures/data, independent of pwd
        if isfolder(liveDir)
            rep = mvst_regress('Report', true, 'SkipIfMissing', true, 'LiveDir', liveDir);
        else
            rep = mvst_regress('Report', true, 'SkipIfMissing', true);
        end
    end
    if rep.skipped
        rows(end+1, :) = {'run_regress', 0, 0, 'SKIPPED (no live data)'}; %#ok<AGROW>
    else
        totalPass = totalPass + rep.nPass;
        totalFail = totalFail + rep.nFail;
        rows(end+1, :) = {sprintf('run_regress (%d files)', rep.nFiles), ...
                          rep.nPass, rep.nFail, ternary(rep.pass, 'PASS', 'FAIL')}; %#ok<AGROW>
    end

    % ---- consolidated summary ----------------------------------------- %
    fprintf('\n---------------------------- summary ----------------------------\n');
    fprintf('  %-26s %6s %6s   %s\n', 'suite', 'pass', 'fail', 'status');
    fprintf('  %s\n', repmat('-', 1, 58));
    for i = 1:size(rows, 1)
        fprintf('  %-26s %6d %6d   %s\n', rows{i,1}, rows{i,2}, rows{i,3}, rows{i,4});
    end
    fprintf('  %s\n', repmat('-', 1, 58));
    fprintf('  TOTAL: %d passed, %d failed, %d skipped\n', totalPass, totalFail, totalSkip);
    fprintf('-----------------------------------------------------------------\n');

    if totalFail > 0
        error('mvst:tests:failed', ...
              '%d test(s)/suite(s) failed; see the summary above.', totalFail);
    end
    fprintf('  ALL GREEN.\n\n');
end

% ===================================================================== %
%  Helpers                                                               %
% ===================================================================== %
function runScript(f)
    %RUNSCRIPT  Execute a script test in an ISOLATED workspace so its local
    %   variables cannot clobber the caller's bookkeeping.
    run(f);
end

function v = ternary(c, a, b)
    if c; v = a; else; v = b; end
end
