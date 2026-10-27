classdef tposms < matlab.unittest.TestCase
%TPOSMS  Unit tests for the stochastic (second-moment) positive-systems layer
%   MVST_POS_MS (Chapter "Linear Positive Systems", stochastic layer): the
%   covariance-level SDP certificates that see the Brownian / Poisson noise --
%   SPN/DNN cone certificates (th:pos:cop, prop:pos:dnntrace, prop:pos:dnnpeak)
%   and the diagonal quadratic certificates (prop:pos:diag / prop:pos:dnndiag).
%   Requires YALMIP + SeDuMi (skipped when missing; see MVST_SETUP).
%
%   Validation anchors:
%     * prop:pos:dnndiag(b) SMALL-NOISE CLOSED FORM: the canonical diagonal
%       certificate D=diag(e_i/d_i) certifies iff d_i Sigma_i e_i < d_i + e_i for
%       all i; the verdict flips exactly at the threshold (both sides swept);
%     * the general diagonal search is broader than the canonical certificate
%       (feasible on a larger noise range);
%     * zero-noise consistency: the SPN and diagonal stability verdicts agree
%       with the mean-level copositive-LP layer (MVST_POS_STABILITY, class ct);
%     * cone tightness (n<=4, SPN=COP): the SPN stability verdict matches the
%       copositive-LP verdict;
%     * trace/peak gains and a Poissonian-jump stability instance solve.

    properties (Constant)
        tol = 1e-6;
    end

    methods (TestClassSetup)
        function addpaths(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS (sossolve) not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods (Test)

        function smallNoiseThresholdCanonical(tc)
            % Canonical D=diag(e_i/d_i) flips exactly at s* = min_i (d_i+e_i)/(d_i e_i).
            A0 = [-1 0;1 -2];
            d = -A0\ones(2,1); e = -A0.'\ones(2,1);
            sstar = min((d+e)./(d.*e));                 % = 1.6667
            for s = [1.0 1.5 1.62]                       % below threshold: feasible
                r = mvst_pos_ms(struct('A0',A0,'Ac',sqrt(s)*eye(2)), ...
                        struct('what','stability','cert','diag','canonical',true));
                tc.verifyTrue(r.feas, sprintf('s=%.2f should certify (< %.3f)',s,sstar));
            end
            for s = [1.72 1.9 2.5]                       % above threshold: infeasible
                r = mvst_pos_ms(struct('A0',A0,'Ac',sqrt(s)*eye(2)), ...
                        struct('what','stability','cert','diag','canonical',true));
                tc.verifyFalse(r.feas, sprintf('s=%.2f should NOT certify (> %.3f)',s,sstar));
            end
        end

        function generalDiagBroaderThanCanonical(tc)
            % At a noise level just past the canonical threshold the canonical D
            % fails but the general diagonal search still certifies.
            A0 = [-1 0;1 -2]; s = 1.75;
            rc = mvst_pos_ms(struct('A0',A0,'Ac',sqrt(s)*eye(2)), ...
                    struct('what','stability','cert','diag','canonical',true));
            rg = mvst_pos_ms(struct('A0',A0,'Ac',sqrt(s)*eye(2)), ...
                    struct('what','stability','cert','diag'));
            tc.verifyFalse(rc.feas);
            tc.verifyTrue(rg.feas);
        end

        function zeroNoiseConsistencyWithLP(tc)
            % No Brownian/Poisson noise: SPN and diagonal second-moment stability
            % must agree with the mean-level copositive-LP verdict.
            A0 = [-1 0;3 -1];
            lp = mvst_pos_stability(struct('A',A0), struct('class','ct'));
            rs = mvst_pos_ms(struct('A0',A0), struct('what','stability','cert','spn'));
            rd = mvst_pos_ms(struct('A0',A0), struct('what','stability','cert','diag'));
            tc.verifyEqual(rs.feas, lp.stable);
            tc.verifyEqual(rd.feas, lp.stable);
        end

        function unstableInfeasible(tc)
            Au = [1 0;1 -2];                            % eig +1: not Hurwitz
            rs = mvst_pos_ms(struct('A0',Au), struct('what','stability','cert','spn'));
            rd = mvst_pos_ms(struct('A0',Au), struct('what','stability','cert','diag'));
            tc.verifyFalse(rs.feas);
            tc.verifyFalse(rd.feas);
        end

        function coneTightnessN2(tc)
            % n=2 <= 4: SPN = COP, so the SPN verdict matches the copositive LP.
            A0 = [-2 0.5;1 -1.5];
            lp = mvst_pos_stability(struct('A',A0), struct('class','ct'));
            rs = mvst_pos_ms(struct('A0',A0), struct('what','stability','cert','spn'));
            tc.verifyEqual(rs.feas, lp.stable);
        end

        function traceGainSPN(tc)
            r = mvst_pos_ms(struct('A0',[-1 0;1 -2],'Ec',[1;0],'Cc',[1 1],'Fc',0), ...
                    struct('what','tracegain','cert','spn'));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThan(r.gamma, 0);
        end

        function peakGainDNN(tc)
            r = mvst_pos_ms(struct('A0',[-1 0;1 -2],'Ec',[1;0],'Cc',[1 1],'Fc',0), ...
                    struct('what','peakgain','cert','spn'));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThan(r.gamma, 0);
        end

        function poissonStabilityRateThreshold(tc)
            % Metzler-Hurwitz flow with an EXPANSIVE Poisson reset: mean-square
            % stability holds at a low rate and fails at a high rate (the
            % rate-weighted CP term sum_r lam_r (Jf_r X Jf_r' - X)).
            A0 = [-2 0;1 -2]; Jf = 1.6*eye(2);
            rlo = mvst_pos_ms(struct('A0',A0,'Jf',Jf,'lamf',0.2), ...
                    struct('what','stability','cert','diag'));
            rhi = mvst_pos_ms(struct('A0',A0,'Jf',Jf,'lamf',5.0), ...
                    struct('what','stability','cert','diag'));
            tc.verifyTrue(rlo.feas);
            tc.verifyFalse(rhi.feas);
        end

        function dwellDiagConstantSlice(tc)
            % The deg-0 minimum-dwell diagonal certificate is the
            % arbitrary-dwell (constant) certificate slice -> identical verdict,
            % and the deg-2 timer dominates the deg-0 slack.
            s = dwellSys();
            rC = mvst_pos_ms(s, struct('what','stability','cert','diag'));
            r0 = mvst_pos_ms(s, struct('what','stability','cert','diag', ...
                                       'dwell','min','Tmin',0.5,'deg',0));
            r2 = mvst_pos_ms(s, struct('what','stability','cert','diag', ...
                                       'dwell','min','Tmin',0.5,'deg',2));
            tc.verifyEqual(r0.feas, rC.feas);                 % constant slice
            tc.verifyTrue(r2.feas);
            tc.verifyGreaterThanOrEqual(r2.info.slack, r0.info.slack - 1e-6); % dominance
        end

        function dwellDiagRangeFeasible(tc)
            % Range dwell-time (theta in [Tmin,Tmax]) diagonal certificate.
            s = dwellSys();
            rR = mvst_pos_ms(s, struct('what','stability','cert','diag', ...
                                       'dwell','range','Tmin',0.3,'Tmax',0.8,'deg',2));
            tc.verifyTrue(rR.feas);
        end

        function dwellTraceSPNConstantSlice(tc)
            % SPN trace gain (prop:pos:dnntrace tau/theta rows).
            % deg-0 min-dwell reproduces the constant SPN trace gain; the deg-2
            % timer certificate dominates it (never worse).
            s = dwellSys();
            tC = mvst_pos_ms(s, struct('what','tracegain','cert','spn'));
            t0 = mvst_pos_ms(s, struct('what','tracegain','cert','spn', ...
                                       'dwell','min','Tmin',0.5,'deg',0));
            t2 = mvst_pos_ms(s, struct('what','tracegain','cert','spn', ...
                                       'dwell','min','Tmin',0.5,'deg',2));
            tc.verifyEqual(t0.gamma, tC.gamma, 'AbsTol', 5e-3);   % constant slice
            tc.verifyLessThanOrEqual(t2.gamma, tC.gamma + 5e-3);  % dominance
        end

        function dwellPeakSPN(tc)
            % (peak): DNN peak dwell certificates (prop:pos:dnnpeak +
            % th:pos:peakgainthm, eq:pos:peakcond).  deg-0 min-dwell reproduces the
            % constant DNN peak; range-dwell solves and dominates; zero-noise
            % long-dwell upper-bounds the steady flow-Gramian peak.
            s = dwellSys();
            pC = mvst_pos_ms(s, struct('what','peakgain','cert','spn'));
            p0 = mvst_pos_ms(s, struct('what','peakgain','cert','spn', ...
                                       'dwell','min','Tmin',0.5,'deg',0));
            pr = mvst_pos_ms(s, struct('what','peakgain','cert','spn', ...
                                       'dwell','range','Tmin',0.3,'Tmax',0.8,'deg',2));
            tc.verifyEqual(p0.gamma, pC.gamma, 'AbsTol', 5e-3);   % constant slice
            tc.verifyTrue(pr.feas);
            tc.verifyLessThanOrEqual(pr.gamma, pC.gamma + 5e-3);  % dominance
            % zero-noise consistency vs steady Gramian peak
            A0=[-1 0.5;0.2 -1.5]; Ec=[0.3;0.2]; Cc=[1 0]; Fc=0.05;
            Pi=sylvester(A0, A0.', -Ec*Ec.'); peakG=max(real(eig(Cc*Pi*Cc.')))+Fc.'*Fc;
            pz=mvst_pos_ms(struct('A0',A0,'Ec',Ec,'Cc',Cc,'Fc',Fc), ...
                    struct('what','peakgain','cert','spn','dwell','min','Tmin',8,'deg',4));
            tc.verifyGreaterThanOrEqual(pz.gamma, peakG - 1e-2);  % tube upper-bounds
        end

        function dwellZeroNoiseMatchesLP(tc)
            % At zero noise the second-moment dwell certificate must agree in
            % FEASIBILITY with the mean-level copositive-LP layer.
            As = [-1 0.5;0.2 -1.5];    Au = [0.5 0.5;0.2 0.1];
            lps = mvst_pos_stability(struct('A',As));
            lpu = mvst_pos_stability(struct('A',Au));
            ms  = mvst_pos_ms(struct('A0',As), struct('what','stability', ...
                              'cert','diag','dwell','min','Tmin',1.0,'deg',0));
            mu  = mvst_pos_ms(struct('A0',Au), struct('what','stability', ...
                              'cert','diag','dwell','min','Tmin',1.0,'deg',0));
            tc.verifyEqual(ms.feas, lps.stable);   tc.verifyTrue(lps.stable);
            tc.verifyEqual(mu.feas, lpu.stable);   tc.verifyFalse(lpu.stable);
        end

        function enginesAgreeStabilityAndGain(tc)
            % Default engine ('sostools') vs legacy 'yalmip': same verdicts, same gains.
            tc.assumeTrue(exist('sdpvar','file')==2, 'YALMIP (sdpvar) not on the path: run mvst_setup');
            A0 = [-1 0;1 -2]; Am = [-1 0.5;0.2 -1.5];
            C = {
              struct('A0',A0,'Ac',sqrt(1.5)*eye(2)),  struct('what','stability','cert','spn')
              struct('A0',A0,'Ac',sqrt(2.5)*eye(2)),  struct('what','stability','cert','spn')
              struct('A0',A0,'Ac',sqrt(1.75)*eye(2)), struct('what','stability','cert','diag')
              struct('A0',A0,'Ac',sqrt(2.5)*eye(2)),  struct('what','stability','cert','diag')
              struct('A0',[-2 0;1 -2],'Jf',1.6*eye(2),'lamf',0.2), struct('what','stability','cert','diag')
              struct('A0',[-2 0;1 -2],'Jf',1.6*eye(2),'lamf',5),   struct('what','stability','cert','diag')
              struct('A0',Am,'J',0.5*eye(2)), struct('what','stability','cert','spn')
              struct('A0',Am,'J',{{1.3*eye(2)}}), struct('what','stability','cert','diag')
              struct('A0',A0,'Ec',[1;0],'Cc',[1 1],'Fc',0), struct('what','tracegain','cert','spn')
              struct('A0',A0,'Ec',[1;0],'Cc',[1 1],'Fc',0), struct('what','peakgain','cert','spn')
              struct('A0',Am,'Ac',0.3*eye(2),'Jf',1.2*eye(2),'lamf',0.3,'Ec',[0.3;0.2],'Cc',[1 0],'Fc',0.05), ...
                  struct('what','peakgain','cert','spn')
              struct('A0',Am,'Ec',[0.3;0.2],'Cc',[1 0],'Fc',0.05,'J',{{0.6*eye(2)}}), struct('what','peakgain','cert','spn')
              };
            for k = 1:size(C,1)
                a = mvst_pos_ms(C{k,1}, setfield(C{k,2},'engine','sostools')); %#ok<SFLD>
                b = mvst_pos_ms(C{k,1}, setfield(C{k,2},'engine','yalmip'));   %#ok<SFLD>
                tc.verifyEqual(a.feas, b.feas, sprintf('case %d verdict',k));
                if a.feas && ~isempty(b.gamma)
                    tc.verifyEqual(a.gamma, b.gamma, 'AbsTol', 1e-4*max(1,abs(b.gamma)), sprintf('case %d gamma',k));
                end
            end
        end

        function scalarBlocksBothBackends(tc)
            % n = 1: all Gram blocks are 1x1 (SeDuMi pretransfo needs them merged);
            % SeDuMi and SOSNOSDP agree (stable scalar Brownian system, both certificates).
            tc.assumeTrue(exist('sossolve_nosdp','file')>0, 'SOSNOSDP not on the path: run mvst_setup');
            for s = {'sedumi','nosdp'}
                for cert = {'spn','diag'}
                    r = mvst_pos_ms(struct('A0',-2,'Ac',0.5), ...
                            struct('what','stability','cert',cert{1},'solver',s{1}));
                    tc.verifyTrue(r.feas, [s{1} '/' cert{1}]);
                end
            end
        end

        function nosdpBackend(tc)
            % The SOS engine runs without an SDP solver on clear-cut cases.
            tc.assumeTrue(exist('sossolve_nosdp','file')>0, 'SOSNOSDP not on the path: run mvst_setup');
            A0 = [-1 0;1 -2];
            ro = struct('what','stability','cert','spn','solver','nosdp');
            tc.verifyTrue (mvst_pos_ms(struct('A0',A0,'Ac',sqrt(1.0)*eye(2)), ro).feas);
            tc.verifyFalse(mvst_pos_ms(struct('A0',[1 0;1 -2]), ro).feas);
            rg = mvst_pos_ms(struct('A0',A0,'Ec',[1;0],'Cc',[1 1],'Fc',0), ...
                   struct('what','tracegain','cert','spn','solver','nosdp'));
            rs = mvst_pos_ms(struct('A0',A0,'Ec',[1;0],'Cc',[1 1],'Fc',0), ...
                   struct('what','tracegain','cert','spn'));
            tc.verifyTrue(rg.feas);
            tc.verifyEqual(rg.gamma, rs.gamma, 'RelTol', 5e-3);
        end

        function diagStorageDecreases(tc)
            % Soundness: the certified diagonal storage V=<D,X> decreases along
            % the simulated mean-square (covariance) flow with a Brownian channel.
            A0 = [-1 0;1 -2]; Ac = 0.5*eye(2);
            r = mvst_pos_ms(struct('A0',A0,'Ac',Ac), struct('what','stability','cert','diag'));
            tc.verifyTrue(r.feas);
            D = r.D; X = [2 1;1 1];                      % initial covariance (DNN)
            Fop = @(X) A0*X + X*A0.' + Ac*X*Ac.';        % covariance flow F(X)
            V0 = trace(D*X); h = 1e-3;
            for k = 1:200; X = X + h*Fop(X); X = (X+X.')/2; end
            tc.verifyLessThan(trace(D*X), V0);           % V decreased
        end

    end
end

% ======================================================================= %
function s = dwellSys()
% Stochastic positive impulsive system: Metzler-Hurwitz flow, Brownian
% channel, contracting reset, trace-gain channels (used by the dwell-time tests).
    s = struct('A0',[-1 0.5;0.2 -1.5],'Ac',{{0.2*eye(2)}}, ...
               'J',{{[0.6 0;0.1 0.5]}}, ...
               'Ec',[0.3;0.2],'Cc',[1 0],'Fc',0.05, ...
               'Cd',[0 1],'Ed',[0.2;0.1],'Fd',0.05);
end
