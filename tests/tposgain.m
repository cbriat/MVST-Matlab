classdef tposgain < matlab.unittest.TestCase
%TPOSGAIN  Unit tests for the class-dispatched positive-systems induced-gain
%   router MVST_POS_GAIN and its timer-dependent copositive-SOS backend
%   MVST_POS_GAIN_SOS (Chapter "Linear Positive Systems").  Covers the LINPROG
%   branches (ct, dt, impulsive-arbitrary, switched, mjls, poisson) and the
%   SOSTOOLS timer-polynomial branches (impulsive/switched minimum and range
%   dwell-time), for both the L1 (trace) and Linf (peak) gains.
%
%   The SOS tests require SOSTOOLS + SeDuMi; both are bundled in thirdparty/
%   (run MVST_SETUP); the suite is skipped when they are missing.
%
%   Validation anchors (book form = ship form):
%     * ct branch on a compartmental flow returns the DC gain C(-A)^{-1}E (+F);
%     * dt branch returns the closed-form Cd(I-J)^{-1}Ed + Fd = 0.45254;
%     * impulsive minimum-dwell Linf on the NAHS sweep data returns the pure
%       book-form gamma (forward timer, frozen tail eq:pos:pktail, no pin):
%       gamma(T=1) = 1.13564, gamma(T=2) = 1.10811 (deg 4), each DOMINATED by
%       the NAHS reference values 1.1455 / 1.1129 (tighter program);
%     * realized peak gain (deterministic worst-ish simulation) <= book gamma;
%     * gamma monotone non-increasing in the minimum dwell-time at deg 6.

    properties (Constant)
        tolLP  = 1e-6;
        tolSOS = 5e-3;
    end

    methods (TestClassSetup)
        function addpaths(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                     % matlab/
            addpath(root); addpath(fullfile(root,'lib'));
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS (sossolve) not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods (Test)

        % ---------------- LINPROG branches ------------------------------- %
        function ctDCgain(tc)
            A=[-0.9 0;0.5 -0.6]; E=[1;0]; C=[1 1];
            for nm = ["L1","Linf"]
                r = mvst_pos_gain(struct('A',A,'E',E,'C',C), struct('norm',char(nm)));
                tc.verifyTrue(r.feas);
                tc.verifyEqual(r.gamma, C*((-A)\E), 'AbsTol', tc.tolLP);
                tc.verifyEqual(r.gamma, r.gammaDC,  'AbsTol', tc.tolLP);
            end
        end

        function dtClosedForm(tc)
            J=[0.1 1;0.1 0.1]; Ed=[0.3;0.3]; Cd=[0 1]; Fd=0.03;
            ref = Cd*((eye(2)-J)\Ed) + Fd;              % = 0.45254...
            for nm = ["L1","Linf"]
                r = mvst_pos_gain(struct('J',J,'Ed',Ed,'Cd',Cd,'Fd',Fd), ...
                                  struct('norm',char(nm),'class','dt'));
                tc.verifyTrue(r.feas);
                tc.verifyEqual(r.gamma, ref, 'AbsTol', tc.tolLP);
            end
        end

        function impArbitraryStatic(tc)
            s = sweepSys();
            r = mvst_pos_gain(s, struct('norm','Linf','class','impulsive','dwell','arbitrary'));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThan(r.gamma, 1.13564);     % arbitrary >= min-dwell
        end

        function switchedCommon(tc)
            sw = struct('A',{{[-1 0;1 -2],[-1 1;1 -3]}}, ...
                        'Ec',{{[0.1;0.1],[0.5;0]}}, ...
                        'Cc',{{[0 1],[0 1]}}, 'Fc',{{0.1,0.1}});
            r = mvst_pos_gain(sw, struct('norm','Linf','class','switched','dwell','arbitrary'));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThan(r.gamma, 0);
        end

        function mjlsCoupled(tc)
            sm = struct('A',{{[-1 0;1 -2],[-1 1;1 -3]}}, 'Pi',[-1 1;2 -2], ...
                        'Ec',{{[0.1;0.1],[0.5;0]}}, 'Cc',{{[0 1],[0 1]}}, 'Fc',{{0.1,0.1}});
            r = mvst_pos_gain(sm, struct('norm','L1','class','mjls'));
            tc.verifyTrue(r.feas);
            tc.verifyEqual(numel(r.weights), 2);        % one weight per mode
        end

        function mjlsForwardVsBlockAssembly(tc)
            % Directed test (asymmetric generator): the mode-indexed mjls gain
            % must equal the gain of the INDEPENDENTLY assembled block system.
            %   Linf forward tube: block(i,i)=A_i+pi_ii I, block(i,j)=pi_ji I
            %   L1  backward adj : same block matrix, transposed by the ct L1 row
            % builder.  Both cross-checks use Abig = blkdiag(A_i) + kron(Pi',I).
            A1=[-2 0;1 -3]; A2=[-1 1;0 -4];
            Pi=[-2 2; 0.5 -0.5];                    % asymmetric generator
            E1=[0.2;0.1]; E2=[0.4;0]; C1=[1 0]; C2=[0 1]; F1=0.05; F2=0.1;
            n=2;
            Abig=blkdiag(A1,A2)+kron(Pi.',eye(n));
            Ebig=blkdiag(E1,E2); Cbig=blkdiag(C1,C2); Fbig=blkdiag(F1,F2);
            sm=struct('A',{{A1,A2}},'Pi',Pi,'Ec',{{E1,E2}},'Cc',{{C1,C2}},'Fc',{{F1,F2}});
            for nm = ["Linf","L1"]
                rBlk = mvst_pos_gain(struct('A',Abig,'E',Ebig,'C',Cbig,'F',Fbig), ...
                                     struct('norm',char(nm),'class','ct'));
                rMj  = mvst_pos_gain(sm, struct('norm',char(nm),'class','mjls'));
                tc.verifyTrue(rBlk.feas && rMj.feas);
                tc.verifyEqual(rMj.gamma, rBlk.gamma, 'AbsTol', 1e-6, ...
                    sprintf('mjls %s must match block assembly', char(nm)));
            end
        end

        function poissonRateTerm(tc)
            sp = struct('A',[-1 0;1 -2],'Ec',[0.1;0.1],'Cc',[0 1],'Fc',0.1, ...
                        'Jf',{{[0.5 0;0 0.5]}},'lamf',0.3);
            r = mvst_pos_gain(sp, struct('norm','Linf','class','poisson'));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThan(r.gamma, 0);
        end

        % ---------------- SOS timer-dependent branches ------------------- %
        function impMinDTLinfBookValue(tc)
            s = sweepSys();
            r1 = mvst_pos_gain(s, struct('norm','Linf','class','impulsive', ...
                                         'dwell','min','Tmin',1,'deg',4));
            r2 = mvst_pos_gain(s, struct('norm','Linf','class','impulsive', ...
                                         'dwell','min','Tmin',2,'deg',4));
            tc.verifyTrue(r1.feas && r2.feas);
            tc.verifyEqual(r1.gamma, 1.13564, 'AbsTol', tc.tolSOS);
            tc.verifyEqual(r2.gamma, 1.10811, 'AbsTol', tc.tolSOS);
            % dominance: book form is tighter than the NAHS reference program.
            tc.verifyLessThanOrEqual(r1.gamma, 1.1455 + tc.tolSOS);
            tc.verifyLessThanOrEqual(r2.gamma, 1.1129 + tc.tolSOS);
        end

        function impMinDTSoundness(tc)
            % Realized peak gain of a worst-ish admissible run must not exceed
            % the certified book-form gamma.
            s = sweepSys();
            r = mvst_pos_gain(s, struct('norm','Linf','class','impulsive', ...
                                        'dwell','min','Tmin',1,'deg',4));
            realized = realizedPeak(s, 1);
            tc.verifyLessThanOrEqual(realized, r.gamma + 1e-3);
        end

        function impMinDTMonotoneDeg6(tc)
            s = sweepSys(); Ts = [1 2 3]; g = zeros(1,3);
            for i = 1:3
                r = mvst_pos_gain(s, struct('norm','Linf','class','impulsive', ...
                                            'dwell','min','Tmin',Ts(i),'deg',6));
                g(i) = r.gamma;
            end
            tc.verifyTrue(all(diff(g) <= 1e-3));         % non-increasing at deg 6
        end

        function impRangeLinfFeasible(tc)
            s = sweepSys();
            r = mvst_pos_gain(s, struct('norm','Linf','class','impulsive', ...
                                        'dwell','range','Tmin',0.3,'Tmax',0.5,'deg',4));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThan(r.gamma, 0);
        end

        function impMinDTL1Feasible(tc)
            s = sweepSys();
            r = mvst_pos_gain(s, struct('norm','L1','class','impulsive', ...
                                        'dwell','min','Tmin',1,'deg',4));
            tc.verifyTrue(r.feas);
            tc.verifyGreaterThan(r.gamma, 0);
        end

        % ---------------- ported MatrixTools helpers --------------------- %
        function helperPortsBuild(tc)
            pvar th;                                    %#ok<NASGU>
            prog = sosprogram(th);
            [prog, S] = mvst_sosposmatrvar(prog, 1, 4, th);   % scalar SOS deg 4
            tc.verifyClass(S, 'dpvar');
            [prog, P] = mvst_sosmatrvar(prog, monomials(th,0:2), [2 1]); %#ok<ASGLU>
            tc.verifyEqual(size(P), [2 1]);
        end

    end
end

% ======================================================================= %
function s = sweepSys()
% NAHS min-DT L-infinity sweep data (constant, arbitrary-DT block).
    s = struct('A',[-1 0;1 -2],'Ec',[0.1;1.1],'Cc',[0 1],'Fc',0.3, ...
               'J',[0.1 1;0.1 0.1],'Ed',[0.3;0.3],'Cd',[0 1],'Fd',0.03);
end

function pk = realizedPeak(s, Tmin)
% Deterministic worst-ish realized peak of |z|_inf under unit-box nonnegative
% inputs, dwell exactly Tmin, several jumps (peak of w is 1, so ratio = peak z).
    A=s.A; Ec=s.Ec; Cc=s.Cc; Fc=s.Fc; J=s.J; Ed=s.Ed; Cd=s.Cd; Fd=s.Fd;
    n=size(A,1); x=zeros(n,1); pk=0; uc=ones(size(Ec,2),1); ud=ones(size(Ed,2),1);
    for k=1:8
        for tau=linspace(0,Tmin,40)
            xt=expm(A*tau)*x + (A\(expm(A*tau)-eye(n)))*Ec*uc;
            pk=max(pk,max(abs(Cc*xt+Fc*uc)));
        end
        xm=expm(A*Tmin)*x + (A\(expm(A*Tmin)-eye(n)))*Ec*uc;
        pk=max(pk,max(abs(Cd*xm+Fd*ud)));
        x=J*xm+Ed*ud;
    end
end
