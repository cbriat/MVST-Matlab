classdef tselftests < matlab.unittest.TestCase
%TSELFTESTS  matlab.unittest wrapper around the four stand-alone SOS selftests
%   (MVST_DISS_SOS_SELFTEST, MVST_DISSYN_SOS_SELFTEST, MVST_GAME_SOS_SELFTEST,
%   MVST_STABILIZE_SELFTEST).  Each selftest stays runnable on its own from
%   the MATLAB prompt; here every row of its result table must pass.
%   Requires SOSTOOLS + SeDuMi (on the path after MVST_SETUP); the suite is
%   skipped when they are missing.

    methods (TestClassSetup)
        function addpaths(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);
            addpath(root); addpath(fullfile(root, 'lib')); addpath(here);
            tc.assumeTrue(exist('sossolve','file')==2, ...
                'SOSTOOLS (sossolve) not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, ...
                'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods (Test)
        function dissSos(tc),    tc.checkAll(mvst_diss_sos_selftest()); end
        function dissynSos(tc),  tc.checkAll(mvst_dissyn_sos_selftest()); end
        function gameSos(tc),    tc.checkAll(mvst_game_sos_selftest()); end
        function stabilize(tc),  tc.checkAll(mvst_stabilize_selftest()); end
    end

    methods (Access = private)
        function checkAll(tc, R)
            tc.verifyNotEmpty(R);
            for k = 1:numel(R)
                tc.verifyTrue(logical(R(k).ok), sprintf('selftest row failed: %s', R(k).name));
            end
        end
    end
end
