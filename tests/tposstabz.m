classdef tposstabz < matlab.unittest.TestCase
%TPOSSTABZ  Unit tests for the timer-dependent positive-systems state-feedback
%   SYNTHESIS engine MVST_POS_STABILIZE_SOS and its router dispatch through
%   MVST_POS_STABILIZE (Chapter "Linear Positive Systems").  Covers the
%   copositive change of variables (lem:pos:cov) with the minimum / range /
%   constant dwell-time synthesis of th:pos:stabzfix / prop:pos:stabzdwell, in
%   both the pure-stabilization and the gain-optimal (Linf) objectives.
%
%   Requires SOSTOOLS + SeDuMi (skipped when missing; see MVST_SETUP).
%
%   Validation anchors (book form = ship form):
%     * min-DT gain-optimal synthesis returns a certified Linf gamma that is
%       DOMINATED by the NAHS reference 0.84401 (the book form omits the
%       script-side Kmax cap / differentiability pin, hence is tighter);
%     * range-DT synthesis reproduces the NAHS reference 0.69199 closely;
%     * the recovered gains render the closed loop positive (Metzler flow,
%       nonnegative jump) and realize a peak gain <= the certified gamma;
%     * the optional Kmax cap (script-side regularizer, default OFF) binds the
%       recovered flow gain when requested.

    properties (Constant)
        tol = 5e-3;
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

        function minDTGainDominance(tc)
            s = minSys();
            r = mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                                             'Tmin',0.2,'deg',2,'norm','Linf'));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThan(r.gamma, 0);
            tc.verifyLessThanOrEqual(r.gamma, 0.84401 + tc.tol);   % dominance
        end

        function rangeDTGainMatchesReference(tc)
            s = rangeSys();
            r = mvst_pos_stabilize(s, struct('class','impulsive','dwell','range', ...
                                             'Tmin',0.1,'Tmax',0.3,'deg',2,'norm','Linf'));
            tc.verifyTrue(r.feas);
            tc.verifyEqual(r.gamma, 0.69199, 'AbsTol', 5e-2);
            tc.verifyLessThanOrEqual(r.gamma, 0.69199 + tc.tol);   % dominance
        end

        function certificationFlag(tc)
            % SeDuMi: certified. nosdp: min-DT certified; range-DT optimum is
            % 'uncertified_optimal' (numerr = 2) -> design returned, certified = false.
            s = minSys();
            r = mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                                             'Tmin',0.2,'deg',2,'norm','Linf'));
            tc.verifyTrue(r.info.certified);
            tc.assumeTrue(exist('sossolve_nosdp','file')>0, 'SOSNOSDP not on the path: run mvst_setup');
            tc.verifyTrue(mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                'Tmin',0.2,'deg',2,'norm','Linf','solver','nosdp')).info.certified);
            ws = warning('off','mvst:uncertified');
            c = onCleanup(@() warning(ws));
            rr = mvst_pos_stabilize(rangeSys(), struct('class','impulsive','dwell','range', ...
                'Tmin',0.1,'Tmax',0.3,'deg',2,'norm','Linf','solver','nosdp'));
            tc.verifyTrue(rr.feas);
            tc.verifyFalse(rr.info.certified);
            tc.verifyEqual(rr.gamma, 0.69199, 'AbsTol', 5e-2);
        end

        function closedLoopPositive(tc)
            s = minSys(); Tmin = 0.2;
            r = mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                                             'Tmin',Tmin,'deg',2,'norm','Linf'));
            metz = true;
            for tau = linspace(0,Tmin,15)
                Acl = s.A + s.Bc*r.K(tau); off = Acl - diag(diag(Acl));
                if any(off(:) < -1e-6); metz = false; end
            end
            Jcl = s.J + s.Bd*r.Kd;
            tc.verifyTrue(metz, 'closed-loop flow must be Metzler');
            tc.verifyGreaterThanOrEqual(min(Jcl(:)), -1e-6);   % nonneg jump
        end

        function soundnessRealizedLeqGamma(tc)
            s = minSys(); Tmin = 0.2;
            r = mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                                             'Tmin',Tmin,'deg',2,'norm','Linf'));
            rp = realizedPeak(s, r.K, r.Kd, Tmin);
            tc.verifyLessThanOrEqual(rp, r.gamma + 1e-2);
        end

        function pureStabilizationFeasible(tc)
            s = minSys();
            r = mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                                             'Tmin',0.2,'deg',2));
            tc.verifyTrue(r.feas);        % no norm -> stabilization slack > 0
        end

        function linfStateFeedbackGainOptimal(tc)
            % L-infinity-optimal state feedback: the certified synthesis gamma
            % equals the analysis L-infinity gain of the designed closed loop,
            % which is Metzler-Hurwitz.
            A=[-2 1;0 -1]; Bc=[1;0]; Ec=[0.3;0.2]; Cc=[0 1];
            r = mvst_pos_stabilize(struct('A',A,'Bc',Bc,'Ec',Ec,'C',Cc), ...
                    struct('class','ct','norm','Linf'));
            tc.verifyTrue(r.feas);
            ra = mvst_pos_gain(struct('A',r.Acl,'E',Ec,'C',Cc), struct('norm','Linf'));
            tc.verifyEqual(r.gamma, ra.gamma, 'AbsTol', 1e-4);        % synth == analysis
            off = r.Acl - diag(diag(r.Acl));
            tc.verifyGreaterThanOrEqual(min(off(:)), -1e-9);          % Metzler
            tc.verifyLessThan(max(real(eig(r.Acl))), 0);             % Hurwitz
        end

        function l1DesignDuality(tc)
            % prop:pos:Linfdesign: L-infinity-optimal state feedback for
            % (A,B,E,C) equals L1-optimal output injection for the transposed
            % data (A',B',C',E'); both routes give the same gain.
            A=[-2 1;0 -1]; Bc=[1;0]; Ec=[0.3;0.2]; Cc=[0 1];
            rL = mvst_pos_stabilize(struct('A',A,'Bc',Bc,'Ec',Ec,'C',Cc), ...
                    struct('class','ct','norm','Linf'));
            rO = mvst_pos_stabilize(struct('A',A.','Cinj',Bc.','Bw',Cc.','Cz',Ec.'), ...
                    struct('class','ct','norm','L1'));
            tc.verifyTrue(rO.feas);
            tc.verifyEqual(rO.gamma, rL.gamma, 'AbsTol', 1e-4);       % duality
            % analysis cross-check: L1 gain of the output-injection loop == synth
            rb = mvst_pos_gain(struct('A',rO.Acl,'E',Cc.','C',Ec.'), struct('norm','L1'));
            tc.verifyEqual(rb.gamma, rO.gamma, 'AbsTol', 1e-4);
        end

        function linfSoundnessSimulation(tc)
            % Realized peak output under unit-amplitude nonnegative disturbances
            % must not exceed the certified L-infinity gain.
            A=[-2 1;0 -1]; Bc=[1;0]; Ec=[0.3;0.2]; Cc=[0 1];
            r = mvst_pos_stabilize(struct('A',A,'Bc',Bc,'Ec',Ec,'C',Cc), ...
                    struct('class','ct','norm','Linf'));
            worst = 0; rng(1);
            for t = 1:2000
                u = rand(size(Ec,2),1);                 % in [0,1]
                xss = -r.Acl\(Ec*u);                    % steady state (stable positive)
                worst = max(worst, abs(Cc*xss));
            end
            tc.verifyLessThanOrEqual(worst, r.gamma + 1e-2);
        end

        function kmax20NahsCrossCheck(tc)
            % NAHS min-DT reference (0.84401) was obtained WITH the script cap
            % Kmax=20 essentially active (max|Kc|=19.985).  The book form
            % omits the cap AND the script's eps=1e-4 gain-row margins, so it
            % stays strictly below the reference even with Kmax=20 set; the cap
            % binds near 20 and dominance (book <= NAHS) holds.  The exact
            % reference is reproduced only by a faithful port of the NAHS script.
            s = minSys();
            r = mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                                'Tmin',0.2,'deg',2,'norm','Linf','Kmax',20));
            tc.verifyTrue(r.feas);
            km = 0; for tau = linspace(0,0.2,21); km = max(km, max(abs(r.K(tau)))); end
            tc.verifyLessThanOrEqual(km, 20 + 1e-2);          % cap binds
            tc.verifyGreaterThan(km, 15);                     % cap essentially active
            tc.verifyLessThanOrEqual(r.gamma, 0.84401 + tc.tol);  % dominance
        end

        function kmaxCapBindsWhenRequested(tc)
            s = minSys();
            r = mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                                'Tmin',0.2,'deg',2,'norm','Linf','Kmax',20));
            tc.verifyTrue(r.feas);
            kmax = 0;
            for tau = linspace(0,0.2,15); kmax = max(kmax, max(abs(r.K(tau)))); end
            tc.verifyLessThanOrEqual(kmax, 20 + 1e-2);  % capped (Kmax=5 is infeasible for this system)
            rinf = mvst_pos_stabilize(s, struct('class','impulsive','dwell','min', ...
                                'Tmin',0.2,'deg',2,'norm','Linf','Kmax',5));
            tc.verifyFalse(rinf.feas);                  % infeasible cap is reported as infeasible
        end

    end
end

% ======================================================================= %
function s = minSys()
    s = struct('A',[3 -1;2 -1],'Ec',[0.5;0.1],'Cc',[0 1],'Fc',0.1, ...
               'J',[2 1;0 1.2],'Ed',[0.8;0.3],'Cd',[0 1],'Fd',0.1, ...
               'Bc',[1;0],'Bd',[1;1]);
end
function s = rangeSys()
    s = struct('A',[1 1;0 1],'Ec',[0.2;0.3],'Cc',[0 1],'Fc',0.1, ...
               'J',[1.2 0;1 0.1],'Ed',[0.3;0.3],'Cd',[0 1],'Fd',0.1, ...
               'Bc',[1;0],'Bd',[1;1]);
end
function rp = realizedPeak(s, Kc, Kd, Tmin)
% Fast analytic realized peak of |z| under unit-box nonnegative inputs and
% admissible dwell >= Tmin, timer-varying flow gain Kc(tau) frozen beyond Tmin.
    A=s.A;Bc=s.Bc;J=s.J;Bd=s.Bd;Ec=s.Ec;Cc=s.Cc;Fc=s.Fc;Ed=s.Ed;Cd=s.Cd;Fd=s.Fd;
    Jcl=J+Bd*Kd; rng(3); rp=0; ns=8;
    for tr=1:200
        x=zeros(2,1); pk=0;
        for k=1:6
            uc=rand; ud=rand; dwell=Tmin+rand*0.6; h=dwell/ns;
            for st=1:ns
                tau=min((st-1)*h,Tmin); Acl=A+Bc*Kc(tau);
                M=expm([Acl Ec*uc; zeros(1,3)]*h); x=M(1:2,1:2)*x+M(1:2,3);
                pk=max(pk,abs(Cc*x+Fc*uc));
            end
            pk=max(pk,abs(Cd*x+Fd*ud)); x=Jcl*x+Ed*ud;
        end
        rp=max(rp,pk);
    end
end
