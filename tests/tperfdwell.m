classdef tperfdwell < matlab.unittest.TestCase
%TPERFDWELL  Unit tests for the dwell-time and deterministic-peak measures of MVST_PERF
%   via MVST_PERF_DWELL (Chapter "Further Performance Measures"): the dwell-time
%   supersolution tubes (prop:perf:peakgaindwell, th:perf:e2p:MDT/RDT,
%   th:perf:h2:MDT/RDT, prop:perf:ptpint:dwell), the deterministic peak-to-peak
%   star-norm (th:perf:ptp:seq / th:perf:ptp:dwell) with its DC-gain lower bound
%   (prop:perf:ptplower), and the amplitude-bounded full-block peak
%   (prop:perf:peakgainfull).  Requires SOSTOOLS + SeDuMi (skipped when missing; see MVST_SETUP); the
%   legacy YALMIP engine is compared in the engine test (skipped without YALMIP).
%
%   Validation anchors:
%     * consistency: the peak/e2p/h2 tube at the constant-dwell (J=I, long-dwell)
%       instance matches the exact steady Gramian engine;
%     * ordering: range-dwell peak <= minimum-dwell peak (a subset of sequences);
%     * sandwich: gammaLB = ||C A^{-1} Bw|| <= certified p2p gamma, and for the
%       SISO CT slice the certificate is tight to the true star-norm (= L1 norm
%       of the impulse response);
%     * soundness: the true star-norm <= certified p2p gamma (direction asserted);
%     * the full-block peak is feasible and dominates the DC lower bound.

    properties (Constant)
        tol = 1e-2;
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

        function peakTubeMatchesGramian(tc)
            % Consistency: with J=I and a long minimum dwell-time the peak tube
            % settles to the continuous-time steady Gramian peak.
            s = sysdata();
            Pi = sylvester(s.A, s.A.', -s.Bw*s.Bw.'); peakG = max(real(eig(s.C*Pi*s.C.')));
            r = mvst_perf(s, 'peak', struct('dwell','min','Tmin',8,'J',eye(2),'deg',4));
            tc.verifyEqual(r.gamma, peakG, 'AbsTol', tc.tol);
        end

        function e2pTubeMatchesGramian(tc)
            s = sysdata();
            Pi = sylvester(s.A, s.A.', -s.Bw*s.Bw.'); peakG = max(real(eig(s.C*Pi*s.C.')));
            r = mvst_perf(s, 'e2p', struct('dwell','min','Tmin',8,'J',eye(2),'deg',4));
            tc.verifyEqual(r.gamma, sqrt(peakG), 'AbsTol', tc.tol);
        end

        function h2TubeUpperBoundsGramian(tc)
            s = sysdata();
            Pi = sylvester(s.A, s.A.', -s.Bw*s.Bw.'); h2G = sqrt(trace(s.C*Pi*s.C.'));
            r = mvst_perf(s, 'h2', struct('dwell','min','Tmin',8,'J',eye(2),'deg',4));
            tc.verifyGreaterThanOrEqual(r.gamma, h2G - tc.tol);   % tube is an upper bound
        end

        function rangeLeqMinPeak(tc)
            % Range dwell-time restricts sequences to [Tmin,Tmax], a subset of the
            % minimum-dwell family, so its certified peak is no larger.
            s = sysdata(); J = 0.7*eye(2);
            rmin = mvst_perf(s, 'peak', struct('dwell','min','Tmin',1,'J',J,'deg',4));
            rrng = mvst_perf(s, 'peak', struct('dwell','range','Tmin',1,'Tmax',2,'J',J,'deg',4));
            tc.verifyLessThanOrEqual(rrng.gamma, rmin.gamma + tc.tol);
        end

        function p2pSandwichTightSISO(tc)
            % gammaLB <= certified gamma; and for the SISO CT slice the star-norm
            % certificate is tight to the true L1 norm of the impulse response.
            s = sysdata();
            r = mvst_perf(s, 'p2p', struct());
            tg = linspace(0,30,20000); g = zeros(size(tg));
            for i = 1:numel(tg); g(i) = s.C*expm(s.A*tg(i))*s.Bw; end
            L1 = trapz(tg, abs(g));                       % true star-norm (SISO CT)
            tc.verifyLessThanOrEqual(r.gammaLB, r.gamma + 1e-6);       % sandwich lower side
            tc.verifyLessThanOrEqual(L1, r.gamma + tc.tol);           % SOUNDNESS: true <= certified
            tc.verifyEqual(r.gamma, L1, 'AbsTol', 5e-2);              % tight for this SISO slice
        end

        function p2pScaleInvariantInBw(tc)
            % gamma is linear in Bw: no clamp from the output-weight mu cap (gate case Bw*1e-3)
            s = sysdata(); r1 = mvst_perf(s, 'p2p', struct());
            s.Bw = 1e-3*s.Bw; r2 = mvst_perf(s, 'p2p', struct());
            tc.verifyEqual(r2.gamma, 1e-3*r1.gamma, 'RelTol', 1e-3);
        end

        function dcGainLowerBound(tc)
            % prop:perf:ptplower closed form.
            s = sysdata();
            r = mvst_perf(s, 'p2p', struct());
            tc.verifyEqual(r.gammaLB, norm(s.C*(s.A\s.Bw)), 'AbsTol', 1e-9);
        end

        function peakFullFeasibleDominatesLB(tc)
            s = sysdata();
            r = mvst_perf(s, 'peakfull', struct());
            tc.verifyTrue(isfinite(r.gamma) && r.gamma > 0);
            tc.verifyGreaterThanOrEqual(r.gamma, r.gammaLB - tc.tol);
        end

        function scalarMixedBlocksBothBackends(tc)
            % n = 1 (1x1 P blocks next to the 2x2 flow block, with a jump and a
            % jump output Cd): SeDuMi and SOSNOSDP agree and are finite.
            tc.assumeTrue(exist('sossolve_nosdp','file')>0, 'SOSNOSDP not on the path: run mvst_setup');
            sys = struct('A',-1,'Bw',1,'C',1,'J',0.8,'Cd',1);
            g = zeros(1,2); k = 0;
            for s = {'sedumi','nosdp'}
                k = k + 1;
                r = mvst_perf_dwell(sys, 'p2p', struct('dwell','min','Tmin',1,'deg',2,'solver',s{1}));
                g(k) = r.gamma;
            end
            tc.verifyTrue(all(isfinite(g)) && all(g > 0));
            tc.verifyEqual(g(2), g(1), 'RelTol', 5e-2);
        end

        function enginesAgreeP2PAndPeakfull(tc)
            % Default engine ('sostools') vs legacy 'yalmip' for the (lambda,nu)
            % bisection LMIs of p2p and peakfull (same grids, same bisection).
            tc.assumeTrue(exist('sdpvar','file')==2, 'YALMIP (sdpvar) not on the path: run mvst_setup');
            C = { sysdata(),'p2p',struct(); ...
                  impSys(),'p2p',struct('dwell','min','Tmin',0.8); ...
                  sysdata(),'peakfull',struct(); ...
                  impSys(),'peakfull',struct('dwell','min','Tmin',0.8) };
            for k = 1:size(C,1)
                a = mvst_perf_dwell(C{k,1}, C{k,2}, setfield(C{k,3},'engine','sostools')); %#ok<SFLD>
                b = mvst_perf_dwell(C{k,1}, C{k,2}, setfield(C{k,3},'engine','yalmip'));   %#ok<SFLD>
                tc.verifyEqual(a.gamma, b.gamma, 'RelTol', 5e-3, sprintf('case %d (%s)',k,C{k,2}));
            end
        end

        function p2pTimerConstantSlice(tc)
            % The deg-0 timer-varying p2p certificate reproduces
            % the constant-P slice, and the deg-2 timer certificate dominates it.
            s = impSys(); Tmin = 0.8;
            rc = mvst_perf(s,'p2p', struct('dwell','min','Tmin',Tmin));
            r0 = mvst_perf(s,'p2p', struct('dwell','min','Tmin',Tmin,'timer',true,'deg',0));
            r2 = mvst_perf(s,'p2p', struct('dwell','min','Tmin',Tmin,'timer',true,'deg',2));
            tc.verifyEqual(r0.gamma, rc.gamma, 'AbsTol', 2e-3);   % deg0 == constant slice
            tc.verifyLessThanOrEqual(r2.gamma, rc.gamma + 2e-3);  % dominance
        end

        function p2pTimerRangeDominates(tc)
            % Range-dwell timer-varying p2p solves and STRICTLY dominates its
            % constant-P slice (the frozen tail is absent, so P(tau) is free).
            s = impSys();
            rc = mvst_perf(s,'p2p', struct('dwell','range','Tmin',0.5,'Tmax',1.2));
            rt = mvst_perf(s,'p2p', struct('dwell','range','Tmin',0.5,'Tmax',1.2, ...
                                           'timer',true,'deg',2));
            tc.verifyTrue(isfinite(rt.gamma));
            tc.verifyLessThan(rt.gamma, rc.gamma);                % strictly better
        end

        function p2pTimerRangeBackendsAgree(tc)
            % nu is probed on the whole grid, so the non-unique constant-slice nu0
            % does not make SOSNOSDP conservative (was 1.877 vs 1.478 for SeDuMi).
            tc.assumeTrue(exist('sossolve_nosdp','file')>0, 'SOSNOSDP not on the path: run mvst_setup');
            g = zeros(1,2); k = 0;
            for sv = {'sedumi','nosdp'}
                k = k + 1;
                r = mvst_perf(impSys(),'p2p', struct('dwell','range','Tmin',0.5,'Tmax',1.2, ...
                                  'timer',true,'deg',2,'solver',sv{1}));
                g(k) = r.gamma;
            end
            tc.verifyEqual(g(2), g(1), 'RelTol', 5e-2);
        end

        function p2pTimerSoundness(tc)
            % Soundness: the realized peak-to-peak gain (bang-bang disturbance,
            % admissible dwell) does not exceed the certified timer gamma.
            s = impSys(); Tmin = 0.8;
            r = mvst_perf(s,'p2p', struct('dwell','min','Tmin',Tmin,'timer',true,'deg',2));
            rp = realizedP2P(s.A, s.Bw, s.C, s.J, Tmin);
            tc.verifyLessThanOrEqual(rp, r.gamma + 1e-2);
        end

    end
end

% ======================================================================= %
function s = sysdata()
    s = struct('A',[-1 0.3;0 -2], 'Bw',[1;1], 'C',[1 0.5]);
end
function s = impSys()
    s = struct('A',[-0.2 1;-1 -0.2], 'Bw',[0.5;0.3], 'C',[1 0], 'J',0.6*eye(2));
end
function rp = realizedP2P(A,Bw,C,J,Tmin)
    rng(1); rp=0; ns=40;
    for tr=1:300
        x=zeros(2,1); pk=0;
        for kk=1:8
            dwell=Tmin+rand*0.8; h=dwell/ns; wsg=sign(randn);
            for st=1:ns
                w=wsg; if rand<0.15; wsg=sign(randn); w=wsg; end
                M=expm([A Bw*w; zeros(1,3)]*h); x=M(1:2,1:2)*x+M(1:2,3);
                pk=max(pk, abs(C*x));
            end
            x=J*x;
        end
        rp=max(rp,pk);
    end
end
