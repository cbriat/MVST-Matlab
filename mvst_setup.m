function mvst_setup(opts)
%MVST_SETUP  Put the MVST companion library on the MATLAB path and run the
%   solver doctor.  Run ONCE after cloning, from the matlab/ directory:
%       run('mvst_setup.m')          % or simply:  mvst_setup
%
%   MVST_SETUP adds
%     * this folder (the parent of the +mvst package, so MVST.SYSTEM etc.
%       resolve),
%     * the engine layer lib/,
%     * the plotting/figure helpers bookutil/ (figure style, PDF export),
%     * the test battery tests/,
%     * every EXISTING bundled solver folder under thirdparty/ (sostools,
%       sedumi, sdpt3, sosnosdp),
%   to the MATLAB path, then calls MVST.SOLVERCHECK to print the third-party
%   toolchain report (which SDP / SOS backends resolve and actually solve).
%
%   The bundled solvers are independent, unmodified snapshots under their own
%   licenses (see THIRD-PARTY.md and thirdparty/VERSIONS.md).  A missing
%   subfolder is silently skipped, so any of them can be deleted or replaced
%   by a newer upstream version.  A copy of a solver ALREADY on the path (your
%   own installation) is respected: the bundled copy of that solver is not
%   added.  SDPT3 ships without compiled MEX files: run install_sdpt3 once.
%   SeDuMi needs its MEX files for fast solves: run install_sedumi once (the
%   doctor prints this hint when a solver is found but does not solve).
%
%   Name/value options:
%     'examples'   (logical, default false) also add examples/ to the path.
%     'quiet'      (logical, default false) suppress the solver report.
%     'thirdparty' (logical, default true)  add the bundled thirdparty/ solvers.
%
%   Example:
%       mvst_setup('examples', true);        % also expose the example scripts
%       mvst_setup('quiet', true);           % path only, no doctor report
%       mvst_setup('thirdparty', false);     % use your own solver installation
%
%   See also MVST.SOLVERCHECK, MVST.OPTIONS, MVST.SYSTEM, RUN_ALL_TESTS.
    arguments
        opts.examples   (1,1) logical = false
        opts.quiet      (1,1) logical = false
        opts.thirdparty (1,1) logical = true
    end

    root = fileparts(mfilename('fullpath'));
    addpath(root);                              % +mvst package parent
    addpath(fullfile(root, 'lib'));             % engine layer (mvst_*)
    addpath(fullfile(root, 'bookutil'));        % plotting style / figure helpers
    addpath(fullfile(root, 'tests'));           % test battery
    if opts.examples
        addpath(fullfile(root, 'examples'));    % companion example scripts
    end

    fprintf('MVST companion library added to the path.\n');
    fprintf('  root : %s\n', root);
    fprintf('  lib  : %s\n', fullfile(root, 'lib'));
    fprintf('  bookutil : %s\n', fullfile(root, 'bookutil'));
    if opts.examples
        fprintf('  examples : %s\n', fullfile(root, 'examples'));
    end

    if opts.thirdparty
        tp = fullfile(root, 'thirdparty');
        % folder, entry point used to detect an installation already on the path
        entries = {
            'sostools',  'sosprogram'
            'sedumi',    'sedumi'
            'sdpt3',     'sdpt3'
            'sosnosdp',  'sossolve_nosdp'
        };
        for k = 1:size(entries, 1)
            d = fullfile(tp, entries{k, 1});
            if ~isfolder(d), continue; end
            if ~isempty(which(entries{k, 2}))
                fprintf('  %-9s: already on the path, bundled copy not added\n', entries{k, 1});
                continue;
            end
            switch entries{k, 1}
                case 'sosnosdp'
                    % mirrors sosnosdp/setup_nosdp.m (without its placeholder logic)
                    sub = {'core', 'bridge', 'examples'};
                    for j = 1:numel(sub)
                        if isfolder(fullfile(d, sub{j})), addpath(fullfile(d, sub{j})); end
                    end
                case 'sdpt3'
                    addpath(d);
                    for sub = {'Solver', 'HSDSolver'}
                        s = fullfile(d, sub{1});
                        if isfolder(s), addpath(genpath(s)); end
                    end
                otherwise
                    addpath(genpath(d));
            end
            fprintf('  %-9s: added from %s\n', entries{k, 1}, d);
        end
        % sosnosdp's placeholder sedumi.m is only for machines with no SDP solver
        % (needed so SOSTOOLS can build programs); never shadow a real solver.
        pl = fullfile(tp, 'sosnosdp', 'placeholder');
        if isfolder(pl) && ~isempty(which('sossolve_nosdp')) && ...
                isempty(which('sedumi')) && isempty(which('sqlp')) && isempty(which('mosekopt'))
            addpath(pl);
            fprintf('  no SDP solver found: sosnosdp placeholder added so SOSTOOLS can build programs\n');
        end
    end

    if ~opts.quiet
        mvst.solvercheck();
    end
end
