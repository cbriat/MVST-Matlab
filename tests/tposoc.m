classdef tposoc < matlab.unittest.TestCase
%TPOSOC  Unit tests for the positive-systems linear-cost optimal-control engine
%   MVST_POS_OC (Chapter "Linear Positive Systems", Section sec:pos:oc): the
%   copositive costate / adjoint programs prop:pos:ocfinite, prop:pos:ocinf,
%   prop:pos:ocdwell, prop:pos:ocbangbang, prop:pos:ocstateconstr.  These
%   modes are ODE / closed-form / LP (no SOS), so the suite runs fast.
%
%   Validation anchors:
%     * infinite-horizon LTI, no penalty: closed form p = -(A')^{-1} q, exact;
%     * finite-horizon costate converges to the infinite-horizon value as T->inf;
%     * bang-bang optimal cost = zero-control cost when no control helps
%       (sigma_c >= 0), and strictly below it when r_c makes sigma_c < 0
%       (eq:pos:ocvalueopt);
%     * state-constrained saturating feedback obeys |K_c| <= E_c
%       (eq:pos:ocKbound, eq:pos:ocKrow).

    properties (Constant)
        tol = 1e-6;
    end

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
        end
    end

    methods (Test)

        function infClosedForm(tc)
            A=[-1 0;1 -2]; q=[1;1]; x0=[1;1];
            r = mvst_pos_oc(struct('A',A), struct('mode','inf','q',q,'x0',x0));
            tc.verifyTrue(r.feas);
            tc.verifyEqual(r.p0, -(A.')\q, 'AbsTol', tc.tol);
            tc.verifyEqual(r.Jopt, (-(A.')\q).'*x0, 'AbsTol', tc.tol);
        end

        function finiteConvergesToInf(tc)
            A=[-1 0;1 -2]; q=[1;1]; x0=[1;1]; pinf = -(A.')\q;
            r = mvst_pos_oc(struct('A',A), ...
                    struct('mode','finite','q',q,'qT',[0;0],'T',30,'x0',x0));
            tc.verifyEqual(r.p0, pinf, 'AbsTol', 1e-4);
        end

        function finiteWithJumps(tc)
            A=[-1 0;1 -2];
            r = mvst_pos_oc(struct('A',A,'J',0.5*eye(2)), ...
                    struct('mode','finite','q',[1;1],'qd',[0.1;0.1],'qT',[0;0], ...
                           'T',3,'tk',[1 2],'x0',[1;1]));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThanOrEqual(min(r.p0), -1e-9);
        end

        function dwellPeriodic(tc)
            A=[-1 0;1 -2];
            r = mvst_pos_oc(struct('A',A,'J',0.8*eye(2)), ...
                    struct('mode','dwell','q',[1;1],'qd',[0.1;0.1],'Tp',1.5,'x0',[1;1]));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThanOrEqual(min(r.p0), -1e-9);
        end

        function bangbangNoBenefit(tc)
            % r_c = 0, B_c >= 0, p >= 0  => sigma_c >= 0 => zero control optimal.
            A=[-1 0;1 -2]; Bc=[1;0.5];
            r = mvst_pos_oc(struct('A',A,'Bc',Bc), ...
                    struct('mode','bangbang','q',[1;1],'rc',0,'T',5,'ubc',1,'x0',[1;1]));
            tc.verifyEqual(r.Jopt, r.p0.'*[1;1], 'AbsTol', 1e-4);
        end

        function bangbangBenefit(tc)
            A=[-1 0;1 -2]; Bc=[1;0.5];
            r = mvst_pos_oc(struct('A',A,'Bc',Bc), ...
                    struct('mode','bangbang','q',[1;1],'rc',-2,'T',5,'ubc',1,'x0',[1;1]));
            tc.verifyLessThan(r.Jopt, r.p0.'*[1;1]);      % control strictly helps
        end

        function stateconstrGainBounded(tc)
            A=[-1 0;1 -2]; Bc=[1;0.5]; Ec=[0.5 0.5];
            r = mvst_pos_oc(struct('A',A,'Bc',Bc,'Ec',Ec), ...
                    struct('mode','stateconstr','q',[1;1],'qT',[0;0],'rc',-1,'T',5,'x0',[1;1]));
            for t = linspace(0,5,10)
                K = r.Kc(t);
                tc.verifyLessThanOrEqual(max(abs(K(:))), max(Ec(:)) + 1e-9);
            end
        end

        function systemPosOcHook(tc)
            % S.pos('oc',...) delegates to the direct MVST_POS_OC call.
            A=[-1 0;1 -2]; q=[1;1]; x0=[1;1];
            S = mvst.System.ct(A, [0;1]);
            rHook   = S.pos('oc', struct('mode','inf','q',q,'x0',x0));
            rDirect = mvst_pos_oc(struct('A',A), struct('mode','inf','q',q,'x0',x0));
            tc.verifyEqual(rHook.Jopt, rDirect.Jopt, 'AbsTol', tc.tol);
        end

    end
end
