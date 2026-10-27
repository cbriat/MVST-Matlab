classdef tpos < matlab.unittest.TestCase
%TPOS  Unit tests for the positive-systems copositive-LP solver family
%   MVST_POS_STABILITY, MVST_POS_GAIN,
%   MVST_POS_STABILIZE, MVST_POS_SUPPLY (shared solver MVST_POS_LP).
%
%   Run from the matlab/tests directory:
%       results = runtests('tpos');
%
%   The gain tests reproduce, to tol 1e-6, numbers currently computed inline
%   by the example scripts ex_pos.m (induced L1 gain and its dual weight) and
%   ex_pos_linf.m (induced peak-to-peak gain and its primal tube).  The
%   dwell-time tests reproduce the maximum/minimum dwell-time thresholds of
%   ex_pos.m / ex_pos_switched.m.  A stabilize test certifies closed-loop
%   positivity and stability, and a supply test cross-checks dissipativity
%   against the L1 gain.

    properties (Constant)
        tol = 1e-6;
    end

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root);
            addpath(fullfile(root, 'lib'));
        end
    end

    methods (Test)

        % =============================================================== %
        %  Gains: reproduce the inline numbers of ex_pos.m / ex_pos_linf.m %
        % =============================================================== %
        function l1GainMatchesExample(tc)
            % ex_pos.m compartmental instance: k1=0.4, k2=0.6, k12=0.5.
            k1 = 0.4; k2 = 0.6; k12 = 0.5;
            A  = [-(k1+k12) 0; k12 -k2];
            E  = [1; 0]; C = [1 1];
            ref     = C*((-A)\E);                    % DC gain = closed form
            refClf  = (1/(k1+k12))*(1+k12/k2);       % chapter closed form
            r = mvst_pos_gain(struct('A',A,'E',E,'C',C), struct('norm','L1'));
            tc.verifyTrue(r.feas);
            tc.verifyEqual(r.gamma,   ref,    'AbsTol', tc.tol);
            tc.verifyEqual(r.gamma,   refClf, 'AbsTol', tc.tol);
            tc.verifyEqual(r.gamma,   r.gammaDC, 'AbsTol', tc.tol);
            % dual weight lambda* = -(A')^{-1} C'  (vertex, exact)
            lamref = -(A.')\(C.');
            tc.verifyEqual(r.weights, lamref, 'AbsTol', tc.tol);
        end

        function linfGainMatchesExample(tc)
            % ex_pos_linf.m 3-state network.
            A = [-3 1 0; 1 -3 1; 0 1 -2];
            E = [1;0;0]; C = [1 1 1];
            ref  = C*((-A)\E);
            xref = (-A)\E;                           % primal tube xi* (exact)
            r = mvst_pos_gain(struct('A',A,'E',E,'C',C), struct('norm','Linf'));
            tc.verifyTrue(r.feas);
            tc.verifyEqual(r.gamma,   ref, 'AbsTol', tc.tol);
            tc.verifyEqual(r.gamma,   r.gammaDC, 'AbsTol', tc.tol);
            tc.verifyEqual(r.weights, xref, 'AbsTol', tc.tol);
        end

        % =============================================================== %
        %  Dwell-time thresholds: ex_pos.m and ex_pos_switched.m           %
        % =============================================================== %
        function maxDwellThreshold(tc)
            % ex_pos.m: unstable Metzler flow, contractive reset.
            % Chapter: Tbar* = 2 ln 10 (monodromy radius crosses 1).
            A = [0.5 1; 0 0.5]; J = [0.1 0.2; 0 0.1];
            r = mvst_pos_stability(struct('A',A,'J',J,'Tbar',2.76), ...
                    struct('class','impulsive','threshold',true, ...
                           'Trange',[0.05 8],'Titers',60));
            tc.verifyTrue(r.stable);                 % feasible at Tbar = 2.76
            tc.verifyEqual(r.Tstar, 2*log(10), 'AbsTol', 1e-3);
        end

        function minDwellThreshold(tc)
            % ex_pos_switched.m: two shear modes, symmetric pair.
            A1 = [-1 0; 3 -1]; A2 = [-1 3; 0 -1];
            % exact spectral threshold (as computed inline in the example)
            f = @(T) max(abs(eig(expm(A2*T)*expm(A1*T))));
            lo = 0.01; hi = 3.0;
            for it = 1:60; m=0.5*(lo+hi); if f(m)<1; hi=m; else; lo=m; end; end
            Tspec = hi;
            r = mvst_pos_stability(struct('A',{{A1,A2}}), ...
                    struct('class','switched','dwell','min','Tmin',2.0, ...
                           'threshold',true,'Trange',[0.01 4],'Titers',60));
            tc.verifyTrue(r.stable);                 % feasible at Tmin = 2
            % LP threshold coincides with the spectral one for this pair
            tc.verifyEqual(r.Tstar, Tspec, 'AbsTol', 5e-3);
        end

        function switchedMinCertificateValid(tc)
            A1 = [-1 0; 3 -1]; A2 = [-1 3; 0 -1];
            Tmin = 2.0; mg = 1e-6;
            r = mvst_pos_stability(struct('A',{{A1,A2}}), ...
                    struct('class','switched','dwell','min','Tmin',Tmin, ...
                           'margin',mg));
            tc.verifyTrue(r.stable);
            p1 = r.p{1}; p2 = r.p{2};
            tc.verifyGreaterThan(min([p1;p2]), 0);
            tc.verifyLessThanOrEqual(A1.'*p1, -mg*ones(2,1) + 1e-9);
            tc.verifyLessThanOrEqual(A2.'*p2, -mg*ones(2,1) + 1e-9);
            % transition inequalities e^{A_i' Tmin} p_j < p_i
            tc.verifyLessThanOrEqual(expm(A1.'*Tmin)*p2 - p1, 1e-9);
            tc.verifyLessThanOrEqual(expm(A2.'*Tmin)*p1 - p2, 1e-9);
            tc.verifyGreaterThan(r.margin, 0);
        end

        % =============================================================== %
        %  ct / dt / mjls stability certificates                          %
        % =============================================================== %
        function ctMetzlerHurwitz(tc)
            A = [-1 0; 3 -1];                        % Metzler, Hurwitz
            r = mvst_pos_stability(struct('A',A), struct('class','ct'));
            tc.verifyTrue(r.stable);
            tc.verifyGreaterThan(min(r.p), 0);
            tc.verifyLessThan(max(A.'*r.p), 0);
        end

        function ctUnstableInfeasible(tc)
            A = [1 0; 0 -1];                         % Metzler, NOT Hurwitz
            r = mvst_pos_stability(struct('A',A), struct('class','ct'));
            tc.verifyFalse(r.stable);
        end

        function dtNonnegativeSchur(tc)
            M = [0.3 0.2; 0.1 0.4];                  % nonnegative, Schur
            r = mvst_pos_stability(struct('A',M), struct('class','dt'));
            tc.verifyTrue(r.stable);
            tc.verifyLessThan(max(M.'*r.p - r.p), 0);
        end

        function mjlsCoupledLP(tc)
            % Two Metzler-Hurwitz modes with a valid generator.
            A1 = [-1 0.5; 0 -2]; A2 = [-2 0; 0.5 -1];
            Pi = [-1 1; 2 -2];
            r = mvst_pos_stability(struct('A',{{A1,A2}},'Pi',Pi), ...
                                   struct('class','mjls'));
            tc.verifyTrue(r.stable);
            p1 = r.p{1}; p2 = r.p{2};
            tc.verifyGreaterThan(min([p1;p2]), 0);
            tc.verifyLessThan(max(A1.'*p1 + Pi(1,1)*p1 + Pi(1,2)*p2), 0);
            tc.verifyLessThan(max(A2.'*p2 + Pi(2,1)*p1 + Pi(2,2)*p2), 0);
        end

        % =============================================================== %
        %  Stabilization: design a gain and certify closure               %
        % =============================================================== %
        function ctStabilizeCertifiesClosure(tc)
            % Metzler but NOT Hurwitz open loop (eig = 1,-2); feedback must
            % render A + Bc K Metzler AND Hurwitz.
            A  = [1 1; 0 -2]; Bc = [1; 0];
            r  = mvst_pos_stabilize(struct('A',A,'Bc',Bc), ...
                                    struct('class','ct'));
            tc.verifyTrue(r.feas);
            Acl = A + Bc*r.K;
            tc.verifyGreaterThanOrEqual(min(Acl(~eye(2))), -1e-9); % Metzler
            tc.verifyLessThan(max(real(eig(Acl))), 0);             % Hurwitz
            tc.verifyTrue(r.cert.stable);            % re-certified by the LP
        end

        % =============================================================== %
        %  General linear supply: dissipativity cross-checks the L1 gain   %
        % =============================================================== %
        function supplyMassBalanceInfeasibleBelowGain(tc)
            % Compartmental system with L1 gain 2.037 > 1, so it is NOT
            % dissipative w.r.t. the mass-balance (gamma = 1) supply.
            A = [-0.9 0; 0.5 -0.6]; E = [1;0]; C = [1 1];
            sys = struct('A',A,'Ec',E,'Cc',C);
            rBad = mvst_pos_supply(sys, ...
                     struct('supply',struct('qc',-1,'rc',1)));
            tc.verifyFalse(rBad.feas);
            % with r = gamma*1 for gamma above the L1 gain it is dissipative
            rOk = mvst_pos_supply(sys, ...
                     struct('supply',struct('qc',-1,'rc',2.1)));
            tc.verifyTrue(rOk.feas);
            tc.verifyGreaterThanOrEqual(min(rOk.lambda), 0);
        end

        % =============================================================== %
        %  System.pos method path == direct MVST_POS_* function path       %
        %  (single- and multi-mode)                                        %
        % =============================================================== %
        function methodPathEqualsDirectCt(tc)
            % Single-mode ct positive system through the factory: the pos
            % method must give the same certificate as the direct call on
            % the equivalent flat struct.
            A = [-1 0; 3 -1];                        % Metzler, Hurwitz
            rDirect = mvst_pos_stability(struct('A',A), struct('class','ct'));
            S = mvst.System.ct(A, zeros(2,0));
            rMethod = S.pos('stability', struct('class','ct'));
            tc.verifyEqual(rMethod.stable, rDirect.stable);
            tc.verifyTrue(rMethod.stable);
            tc.verifyEqual(rMethod.p,      rDirect.p,      'AbsTol', tc.tol);
            tc.verifyEqual(rMethod.margin, rDirect.margin, 'AbsTol', tc.tol);
        end

        function methodPathEqualsDirectSwitched(tc)
            % Multi-mode switched positive system through the factory: the
            % per-mode drifts must be gathered into the cell arrangement the
            % switched copositive LP expects, matching the direct call.
            A1 = [-1 0; 3 -1]; A2 = [-1 3; 0 -1];
            optSw = struct('class','switched','dwell','min','Tmin',2.0);
            rDirect = mvst_pos_stability(struct('A',{{A1,A2}}), optSw);
            S = mvst.System.switched({A1,A2}, {zeros(2,0), zeros(2,0)});
            rMethod = S.pos('stability', optSw);
            tc.verifyEqual(rMethod.stable, rDirect.stable);
            tc.verifyTrue(rMethod.stable);
            tc.verifyEqual(rMethod.margin, rDirect.margin, 'AbsTol', tc.tol);
            tc.verifyEqual(rMethod.p{1}, rDirect.p{1}, 'AbsTol', tc.tol);
            tc.verifyEqual(rMethod.p{2}, rDirect.p{2}, 'AbsTol', tc.tol);
        end

        % =============================================================== %
        %  Solver plumbing                                                 %
        % =============================================================== %
        function unknownSolverErrors(tc)
            lp = struct('nvar',1,'Aub',1,'bub',1,'lb',0);
            tc.verifyError(@() mvst_pos_lp(lp, struct('solver','sedumi')), ...
                           'mvst:pos:solver');
        end
    end
end
