classdef tnosdp < matlab.unittest.TestCase
%TNOSDP  The 'nosdp' (SOSNOSDP) backend agrees with SeDuMi through MVST_SOSSOLVE.
%   Runs MVST_DISS_SOS with OPTS.solver = 'nosdp' and 'sedumi' on a CT L2-gain
%   case (the gammas agree with each other and with the exact H-infinity norm),
%   a feasible passivity case and an infeasible (non-passive) case (the verdicts
%   agree), plus the status mapping of MVST_SOSSOLVE and the MVST_OPTS guards.
%   Skipped when SOSTOOLS, sossolve_nosdp or SeDuMi is missing.

    methods (TestClassSetup)
        function addpaths(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);
            addpath(root); addpath(fullfile(root, 'lib'));
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS not on the path: run mvst_setup');
            tc.assumeTrue(exist('sossolve_nosdp','file')==2, 'sosnosdp not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods (Test)
        function l2gainAgrees(tc)
            sys = struct('A0',[-3 0;0 -3],'Bw',[0;1],'C',[0.5 1],'Dw',0.5);
            ex  = mvst_l2gain_ct(sys.A0,sys.Bw,sys.C,sys.Dw);
            rs  = mvst_diss_sos(sys, struct('type','l2gain'), struct('type','arbitrary'), 'gamma', struct('solver','sedumi'));
            rn  = mvst_diss_sos(sys, struct('type','l2gain'), struct('type','arbitrary'), 'gamma', struct('solver','nosdp'));
            tc.verifyEqual(rn.gamma, rs.gamma, 'RelTol', 2e-2);
            tc.verifyEqual(rn.gamma, ex, 'RelTol', 2e-2);
        end
        function passivityFeasible(tc)
            sys = struct('A0',[0 1;-1 -1],'Bw',[0;1],'C',[1 0],'Dw',0.5);
            for s = {'sedumi','nosdp'}
                r = mvst_diss_sos(sys, struct('type','passivity'), struct('type','arbitrary'), 'feas', struct('solver',s{1}));
                tc.verifyTrue(logical(r.feas), s{1});
            end
        end
        function passivityInfeasible(tc)
            sys = struct('A0',[0 1;-1 -1],'Bw',[0;1],'C',[0 1],'Dw',0);
            for s = {'sedumi','nosdp'}
                r = mvst_diss_sos(sys, struct('type','passivity'), struct('type','arbitrary'), 'feas', struct('solver',s{1}));
                tc.verifyFalse(logical(r.feas), s{1});
            end
        end
        function statusMapping(tc)
            % infeasible SOS feasibility program: -1 - x^2 >= 0 is not SOS
            pvar x; prog = sosprogram(x); prog = sosineq(prog, -1 - x^2);
            [~, info] = mvst_sossolve(prog, struct('solver','nosdp'));
            tc.verifyEqual([info.pinf info.dinf], [0 1]);
            tc.verifyEqual(info.feasratio, -1);
            tc.verifyGreaterThanOrEqual(info.numerr, 0);
            prog = sosprogram(x); prog = sosineq(prog, 1 + x^2);
            [~, info] = mvst_sossolve(prog, struct('solver','nosdp'));
            tc.verifyEqual([info.pinf info.dinf info.numerr], [0 0 0]);
        end
        function reltolReachesBothBackends(tc)
            % normalized opts (reltol) must reach the solver: a loose tolerance changes the answer
            sys = struct('A0',[-3 0;0 -3],'Bw',[0;1],'C',[0.5 1],'Dw',0.5);
            ex  = mvst_l2gain_ct(sys.A0,sys.Bw,sys.C,sys.Dw);
            for s = {'sedumi','nosdp'}
                gl = mvst_diss_sos(sys, struct('type','l2gain'), struct('type','arbitrary'), 'gamma', struct('solver',s{1},'reltol',1e-1)).gamma;
                gt = mvst_diss_sos(sys, struct('type','l2gain'), struct('type','arbitrary'), 'gamma', struct('solver',s{1},'reltol',1e-9)).gamma;
                tc.verifyEqual(gt, ex, 'RelTol', 1e-4, s{1});
                tc.verifyGreaterThan(abs(gl-gt), 1e-3, [s{1} ': loose reltol ignored']);
            end
        end
        function optsGuards(tc)
            o = mvst_opts(struct('solver','nosdp'), true, 'sos');
            tc.verifyEqual(o.solver, 'nosdp');
            tc.verifyError(@() mvst_opts(struct('solver','nosdp'), true, 'yalmip'), 'mvst:opts:nosdpYalmip');
            tc.verifyEqual(mvst.Options('solver','nosdp').solver, 'nosdp');
        end
    end
end
