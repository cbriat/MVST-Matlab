classdef tocdwell < matlab.unittest.TestCase
%TOCDWELL  Unit tests for the causal minimum / range dwell-time optimal-control
%   engine MVST_OC_DWELL (Chapter "Linear Cost", Section sec:OC:dwell): the
%   timer-dependent forward-Riccati / jump-GRE programs Theorem th:OC:MDT and
%   Theorem th:OC:RDT, and the average-cost dwell bounds th:IMP:AC:dwell:MDT /
%   th:IMP:AC:dwell:RDT.  These modes are ODE / Riccati (no SOS), so the suite
%   runs on the core library alone.
%
%   Validation anchors:
%     * CARE limit: as Tmin -> inf (constant data) the certificate P_s(0)
%       approaches the algebraic I-GRE solution MVST_CARE_CT;
%     * monotonicity: the certified cost bound <P_s(0),X10> is non-increasing
%       in the minimum dwell-time;
%     * feasibility pattern: the end-of-dwell LMI (iii) holds for the classic
%       MDT scenario (expansive reset stabilized by dwell) and fails for a
%       contractive reset whose frozen tail cannot dissipate;
%     * SOUNDNESS (the theorem's direction): the recovered causal policy
%       realizes a cost J <= <P_s(0),X10> over admissible dwell sequences, for
%       BOTH minimum and range dwell-time;
%     * average mode returns a finite occupation-weighted rate bound.

    properties (Constant)
        tol = 1e-4;
    end

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
        end
    end

    methods (Test)

        function careLimit(tc)
            [md,jump] = sys();
            Pinf = mvst_care_ct(md);
            r = mvst_oc_dwell(md, jump, struct('dwell','min','Tmin',30));
            tc.verifyLessThan(norm(r.P0 - Pinf,'fro'), 1e-5);
        end

        function mdtMonotoneAndFeasible(tc)
            [md,jump] = sys(); Ts = [0.3 0.6 1 2 4]; Jb = zeros(size(Ts));
            for i = 1:numel(Ts)
                r = mvst_oc_dwell(md, jump, struct('dwell','min','Tmin',Ts(i),'X10',eye(2)));
                tc.verifyTrue(r.feas);
                Jb(i) = r.Jbound;
            end
            tc.verifyTrue(all(diff(Jb) <= 1e-6));       % non-increasing in Tmin
        end

        function contractiveJumpInfeasible(tc)
            % A contractive reset's frozen tail cannot dissipate: (iii) fails.
            [md,~] = sys();
            jc = struct('J',0.6*eye(2),'Bd',[0;1],'Qd',0.05*eye(2),'Rd',1);
            r = mvst_oc_dwell(md, jc, struct('dwell','min','Tmin',1));
            tc.verifyFalse(r.feas);
            tc.verifyGreaterThan(r.lmiMDT, 1e-6);       % end-of-dwell LMI violated
        end

        function mdtSoundness(tc)
            [md,jump] = sys();
            r = mvst_oc_dwell(md, jump, struct('dwell','min','Tmin',1,'X10',eye(2)));
            x0 = [1;-0.5]; bound = x0.'*r.P0*x0;
            worst = worst_cost(md, jump, r, 1, 1, x0, 40);
            tc.verifyLessThanOrEqual(worst, bound + 1e-3);   % J <= <P_s(0),X10>
        end

        function rdtSoundness(tc)
            [md,jump] = sys();
            r = mvst_oc_dwell(md, jump, struct('dwell','range','Tmin',1,'Tmax',2,'X10',eye(2)));
            tc.verifyTrue(r.feas);
            x0 = [1;-0.5]; bound = x0.'*r.P0*x0;
            worst = worst_cost(md, jump, r, 1, 2, x0, 40);
            tc.verifyLessThanOrEqual(worst, bound + 1e-3);
        end

        function rdtFamilyCheckHolds(tc)
            % Passing case: the jump-GRE family condition holds on [Tmin,Tmax]
            % (P_s(theta) >= Pm), so the certified bound is sound.
            [md,jump] = sys();
            r = mvst_oc_dwell(md, jump, struct('dwell','range','Tmin',1,'Tmax',2,'X10',eye(2)));
            tc.verifyTrue(r.feas);
            tc.verifyLessThanOrEqual(r.lmiRDT, 1e-6);
        end

        function rdtNonMonotoneCaught(tc)
            % Directed non-monotone case: a strong time-varying cost SPIKE on
            % (Tmin,Tmax) drives the forward certificate P_s(theta) below Pm at an
            % interior theta, violating the jump-GRE family (th:OC:RDT (ii)).  The
            % engine must VERIFY the family on a grid and return feas=false, not
            % assume the forward flow is monotone.  (Confirmed genuinely unsound:
            % the recovered policy realizes a cost far exceeding <P_s(0),X10>.)
            A = [-1 0.5;0 -0.8]; B = [0;1]; R = 1;
            qf = @(t) 1 + 40*(t > 1.3 & t < 1.65);      % cost spike mid-window
            md = @(t) mvst_mode(A, B, qf(t)*eye(2), R);
            jump = struct('J',1.2*eye(2),'Bd',[0;1],'Qd',0.05*eye(2),'Rd',1);
            r = mvst_oc_dwell(md, jump, struct('dwell','range','Tmin',1,'Tmax',2,'X10',eye(2)));
            tc.verifyFalse(r.feas);
            tc.verifyGreaterThan(r.lmiRDT, 1e-6);       % interior violation caught
        end

        function averageMode(tc)
            [md,jump] = sys();
            r = mvst_oc_dwell(md, jump, struct('dwell','min','Tmin',1,'mode','average', ...
                    'W',eye(2),'Winf',eye(2),'Wjp',0.1*eye(2),'nu',0.8));
            tc.verifyTrue(r.feas);
            tc.verifyTrue(isfinite(r.rho) && r.rho > 0);
        end

    end
end

% ======================================================================= %
function [md,jump] = sys()
% Classic MDT scenario: expansive reset stabilized by a minimum dwell-time.
    md = mvst_mode([-1 0.5;0 -0.8], [0;1], eye(2), 1);
    jump = struct('J',1.2*eye(2),'Bd',[0;1],'Qd',0.05*eye(2),'Rd',1);
end

function worst = worst_cost(md, jump, r, Tmin, Tmax, x0, ntrials)
% Worst realized cost of the recovered causal policy over random admissible
% dwell sequences; the theorem asserts this stays <= <P_s(0),X10>.
    A = md.A; B = md.B; Q = md.Q; R = md.R;
    J = jump.J; Bd = jump.Bd; Qd = jump.Qd; Rd = jump.Rd;
    isRange = Tmax > Tmin + 1e-12;
    worst = 0; rng(7);
    oo = odeset('RelTol',1e-8,'AbsTol',1e-10);
    for t = 1:ntrials
        x = x0; Jc = 0;
        for k = 1:6
            D = Tmin + (Tmax-Tmin)*rand;
            [~,Y] = ode45(@(tau,y) flow(tau,y,A,B,Q,R,r,Tmin), [0 D], [x;0], oo);
            x = Y(end,1:2).'; Jc = Jc + Y(end,3);
            if isRange; Kd = r.Kd(D); else; Kd = r.Kd; end
            ud = Kd*x; Jc = Jc + x.'*Qd*x + ud.'*Rd*ud; x = J*x + Bd*ud;
        end
        worst = max(worst, Jc);
    end
end

function dy = flow(tau, y, A, B, Q, R, r, Tmin)
    x = y(1:2); Kc = r.Kc(min(tau,Tmin)); u = Kc*x;
    dy = [(A + B*Kc)*x; x.'*Q*x + u.'*R*u];
end
