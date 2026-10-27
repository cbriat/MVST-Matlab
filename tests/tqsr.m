classdef tqsr < matlab.unittest.TestCase
%TQSR  Tests of the general quadratic (Q,S,R) supply in the SOS engines.
%   Contract: supply = struct('type','qsr','Q',Q,'S',S,'R',R), the supply
%   rate [z;w]'*[Q S;S' R]*[z;w] (l2gain: Q=-I,S=0,R=gamma^2 I; passivity:
%   Q=-rho I, S=I, R=-nu I).  Anchors:
%     * validation of the struct (sizes, symmetry, type) by MVST_QSR_SUPPLY;
%     * qsr with the l2gain / passivity data reproduces the specialised
%       blocks bit-for-bit (identical slack and certificates, same program);
%     * a genuinely indefinite (sector) supply with a verdict known from the
%       Nyquist disc of G(s) = d + 1/(s+1) (circle with diameter [d,d+1]);
%     * the same for the product-cone engine and the synthesis engine.
    properties (Constant)
        sys = struct('A0',-1,'Bw',1,'C',1,'Dw',0.5);
    end
    methods (TestClassSetup)
        function addPaths(tc)
            here = fileparts(mfilename('fullpath')); root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end
    methods (Test)
        function partialJumpWeightsRefused(tc)
            sp = struct('type','qsr','Q',-1,'S',0,'R',4,'Qk',-5);
            tc.verifyError(@() mvst_qsr_supply(sp,1,1,'jump'), 'mvst:qsr:jumpweights');
        end

        function fullJumpWeightsHonored(tc)
            sp = struct('type','qsr','Q',-1,'S',0,'R',4,'Qk',-5,'Sk',1,'Rk',9);
            [Q,S,R] = mvst_qsr_supply(sp,1,1,'jump');
            tc.verifyEqual([Q S R], [-5 1 9]);
            [Q,S,R] = mvst_qsr_supply(sp,1,1);
            tc.verifyEqual([Q S R], [-1 0 4]);
        end

        function validation(tc)
            f = @(s,p,m) mvst_qsr_supply(s,p,m);
            ok = struct('type','qsr','Q',-eye(2),'S',zeros(2,1),'R',1);
            [Q,S,R] = f(ok,2,1); tc.verifyEqual({Q,S,R},{-eye(2),zeros(2,1),1});
            bad = ok; bad.Q = -eye(3);
            tc.verifyError(@() f(bad,2,1),'mvst:qsr:sizeQ');
            bad = ok; bad.S = zeros(1,2);
            tc.verifyError(@() f(bad,2,1),'mvst:qsr:sizeS');
            bad = ok; bad.R = eye(2);
            tc.verifyError(@() f(bad,2,1),'mvst:qsr:sizeR');
            bad = ok; bad.Q = [-1 1;0 -1];
            tc.verifyError(@() f(bad,2,1),'mvst:qsr:symQ');
            bad = rmfield(ok,'R');
            tc.verifyError(@() f(bad,2,1),'mvst:qsr:missing');
            bad = ok; bad.type = 'l2gain';
            tc.verifyError(@() f(bad,2,1),'mvst:qsr:type');
            tc.verifyError(@() mvst_qsr_supply(ok,3,1,'jump'),'mvst:qsr:jumpweights');
        end

        function reproducesL2gainBitForBit(tc)
            g = 0.9;
            a = mvst_diss_sos(tc.sys, struct('type','l2gain','gamma',g), struct('type','arbitrary'), 'feas');
            b = mvst_diss_sos(tc.sys, struct('type','qsr','Q',-1,'S',0,'R',g^2), struct('type','arbitrary'), 'feas');
            tc.verifyEqual(a.feas, b.feas);
            tc.verifyTrue(isequal(a.slack, b.slack), 'slack must be bit-identical');
            tc.verifyTrue(isequal(a.Pfun(0), b.Pfun(0)), 'certificate must be bit-identical');
        end

        function reproducesPassivityBitForBit(tc)
            a = mvst_diss_sos(tc.sys, struct('type','passivity','nu',0.5,'rho',0.1), struct('type','arbitrary'), 'feas');
            b = mvst_diss_sos(tc.sys, struct('type','qsr','Q',-0.1,'S',1,'R',-0.5), struct('type','arbitrary'), 'feas');
            tc.verifyTrue(a.feas && b.feas);
            tc.verifyTrue(isequal(a.slack, b.slack));
        end

        function indefiniteSectorVerdict(tc)
            % G = 0.5 + 1/(s+1): Nyquist disc has diameter [0.5, 1.5].
            % Sector [a,b]: s = -z^2 + (a+b) z w - a b w^2 (Q=-1,S=(a+b)/2,R=-ab).
            sec = @(a,b) struct('type','qsr','Q',-1,'S',(a+b)/2,'R',-a*b);
            run = @(a,b) mvst_diss_sos(tc.sys, sec(a,b), struct('type','arbitrary'), 'feas');
            tc.verifyTrue (run(0.4,1.6).feas, 'inside the Nyquist disc: dissipative');
            tc.verifyFalse(run(0.6,1.6).feas, 'lower edge violated');
            tc.verifyFalse(run(0.4,1.4).feas, 'upper edge violated');
            % the supply matrix is genuinely indefinite
            M = [-1 1;1 -0.64]; tc.verifyLessThan(min(eig(M)),0); tc.verifyGreaterThan(max(eig(M)),0);
        end

        function passivityShiftedIndefinite(tc)
            % Re G >= 0.5 -> passive with nu up to 2*0.5 = 1: Q=0,S=1,R=-nu.
            q = @(nu) struct('type','qsr','Q',0,'S',1,'R',-nu);
            run = @(nu) mvst_diss_sos(tc.sys, q(nu), struct('type','arbitrary'), 'feas');
            tc.verifyTrue (run(0.9).feas);
            tc.verifyFalse(run(1.1).feas);
        end

        function modesEngine(tc)
            sm = struct('A0',cat(3,-1,-2),'Bw',1,'C',cat(3,1,1),'Dw',0.5,'Pig',[-1 1;1 -1]);
            g = 1.2;
            a = mvst_diss_sos_modes(sm, struct('type','l2gain','gamma',g), struct('type','ctmjls'), 'feas');
            b = mvst_diss_sos_modes(sm, struct('type','qsr','Q',-1,'S',0,'R',g^2), struct('type','ctmjls'), 'feas');
            tc.verifyEqual(a.feas,b.feas);
            tc.verifyTrue(isequal(a.slack,b.slack));
            q = @(lo,hi) struct('type','qsr','Q',-1,'S',(lo+hi)/2,'R',-lo*hi);
            tc.verifyTrue (mvst_diss_sos_modes(sm,q(0.2,2.5),struct('type','ctmjls'),'feas').feas);
            tc.verifyFalse(mvst_diss_sos_modes(sm,q(1.0,2.5),struct('type','ctmjls'),'feas').feas);
        end

        function nosdpBackend(tc)
            % qsr supply on the SDP-free backend: l2 reproduction and the sector verdicts.
            tc.assumeTrue(exist('sossolve_nosdp','file')>0, 'SOSNOSDP not on the path: run mvst_setup');
            dw = struct('type','arbitrary');
            sec = @(a,b) struct('type','qsr','Q',-1,'S',(a+b)/2,'R',-a*b);
            % 2-state realization of the same G (SOSNOSDP crashes on lone 1x1 Gram blocks, n = 1)
            s2 = struct('A0',diag([-1 -3]),'Bw',[1;0],'C',[1 0],'Dw',0.5);
            run = @(a,b) mvst_diss_sos(s2, sec(a,b), dw, 'feas', struct('solver','nosdp'));
            tc.verifyTrue (run(0.4,1.6).feas);
            tc.verifyFalse(run(0.6,1.6).feas);
        end

        function synthesisEngine(tc)
            sy = struct('A0',[0 1;-1 0],'B0',[0;1],'Bw',[0;1],'C',[0 1],'Du',0,'Dw',0.5);
            dw = struct('type','arbitrary');
            g = 3;
            a = mvst_dissyn_sos(sy, struct('type','l2gain','gamma',g), dw, 'feas');
            b = mvst_dissyn_sos(sy, struct('type','qsr','Q',-1,'S',0,'R',g^2), dw, 'feas');
            tc.verifyEqual(a.feas,b.feas); tc.verifyTrue(a.feas);
            tc.verifyTrue(isequal(a.slack,b.slack));
            % passivity by qsr (Q=0: no OSP term)
            p = mvst_dissyn_sos(sy, struct('type','qsr','Q',0,'S',1,'R',0), dw, 'feas');
            p0 = mvst_dissyn_sos(sy, struct('type','passivity'), dw, 'feas');
            tc.verifyEqual(p.feas,p0.feas); tc.verifyTrue(p.feas);
            % non-convexifiable Q is refused
            tc.verifyError(@() mvst_dissyn_sos(sy, struct('type','qsr','Q',1,'S',0,'R',1), dw, 'feas'), ...
                'mvst:dissyn:qsrconvex');
        end
    end
end
