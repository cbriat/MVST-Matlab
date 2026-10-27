classdef thinfnorm < matlab.unittest.TestCase
%THINFNORM  Tests for MVST_HINFNORM, the DT / DT-MJLS l2-gain backends and the
%   sampled-data convergence fix of MVST_L2GAIN_SD_HYBRID.
%
%   Anchors: CT exact vs closed-form norms (first order; resonant second order
%   1/(2 zeta sqrt(1-zeta^2))); DT vs hand-computed norms (1/(z-0.5) -> 2,
%   1+1/(z-0.5) -> 3) and a dense unit-circle grid; DT-MJLS reduces to DT;
%   certificate LMIs; the sampled-data example data of MVST_DISS_SOS (hybrid gain
%   Nper=60 -> 2.053, converged 2.108).  Base MATLAB only (no solver needed).

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
        end
    end

    methods (Test)
        function ctFirstOrder(tc)
            r = mvst_hinfnorm(struct('A0',-1,'Bw',1,'C',1));
            tc.verifyEqual(r.gamma, 1, 'AbsTol', 2e-4);
            tc.verifyEqual(r.class, 'ct');
            tc.verifyTrue(r.diagnostics.exact);
        end

        function ctResonant(tc)
            z = 0.05;
            s = struct('A0',[0 1;-1 -2*z],'Bw',[0;1],'C',[1 0]);
            r = mvst_hinfnorm(s);
            tc.verifyEqual(r.gamma, 1/(2*z*sqrt(1-z^2)), 'RelTol', 5e-4);
        end

        function ctCertificateIsStorage(tc)
            s = struct('A0',[0 1;-1 -0.5],'Bw',[0;1],'C',[1 0],'Dw',0.2);
            r = mvst_hinfnorm(s);
            tc.assertNotEmpty(r.certificate);
            S = r.certificate.S; g = r.certificate.gamma;
            A = s.A0; B = s.Bw; C = s.C; D = s.Dw;
            L = [A'*S+S*A+C'*C, S*B+C'*D; B'*S+D'*C, D'*D-g^2];
            tc.verifyLessThanOrEqual(max(eig((L+L')/2)), 1e-7);
            tc.verifyGreaterThanOrEqual(min(eig(S)), -1e-9);
        end

        function dtHandComputed(tc)
            r = mvst_hinfnorm(struct('A0',0.5,'Bw',1,'C',1,'domain','dt'));
            tc.verifyEqual(r.class, 'dt');
            tc.verifyEqual(r.gamma, 2, 'AbsTol', 2e-4);          % 1/(z-0.5), peak at z=1
            r = mvst_hinfnorm(struct('A0',0.5,'Bw',1,'C',1,'Dw',1,'domain','dt'));
            tc.verifyEqual(r.gamma, 3, 'AbsTol', 2e-4);          % 1 + 1/(z-0.5)
            r = mvst_hinfnorm(struct('A0',-0.5,'Bw',1,'C',1,'domain','dt'));
            tc.verifyEqual(r.gamma, 2, 'AbsTol', 2e-4);          % pole at -0.5, peak at z=-1
        end

        function dtAgainstGrid(tc)
            A = [0.9 0.3;-0.3 0.9]; B = [1;0.5]; C = [1 0]; D = 0.2;
            g = 0;
            for t = linspace(0, pi, 100001)
                z = exp(1i*t); g = max(g, norm(C*((z*eye(2)-A)\B)+D));
            end
            tc.verifyEqual(mvst_l2gain_dt(A,B,C,D), g, 'RelTol', 1e-4);
            tc.verifyError(@() mvst_l2gain_dt(1.1,1,1,0), 'mvst:l2gain_dt:unstable');
        end

        function dtCertificateIsStorage(tc)
            s = struct('A0',[0.5 0.2;0 0.3],'Bw',[1;0],'C',[1 0],'Dw',0.1,'domain','dt');
            r = mvst_hinfnorm(s);
            tc.assertNotEmpty(r.certificate);
            S = r.certificate.S; g = r.certificate.gamma;
            A = s.A0; B = s.Bw; C = s.C; D = s.Dw;
            L = [A'*S*A-S+C'*C, A'*S*B+C'*D; B'*S*A+D'*C, B'*S*B+D'*D-g^2];
            tc.verifyLessThanOrEqual(max(eig((L+L')/2)), 1e-6);
        end

        function dtMjlsReducesToDt(tc)
            A = [0.9 0.3;-0.3 0.9]; B = [1;0.5]; C = [1 0]; D = 0.2;
            gdt = mvst_l2gain_dt(A,B,C,D);
            tc.verifyEqual(mvst_l2gain_dt_mjls({A},B,{C},D,1), gdt, 'RelTol', 5e-4);
            P = [.3 .7;.5 .5];
            tc.verifyEqual(mvst_l2gain_dt_mjls({A,A},B,{C,C},D,P), gdt, 'RelTol', 5e-4);
            tc.verifyError(@() mvst_l2gain_dt_mjls({A,A},B,{C,C},D,[.5 .6;.5 .5]), ...
                'mvst:l2gain_dt_mjls:transition');
        end

        function dtMjlsCouplingMatters(tc)
            % a second, higher-gain mode visited with positive probability raises the gain
            A1 = 0.3; A2 = 0.9; B = 1; C = 1;
            g1 = mvst_l2gain_dt_mjls({A1,A1},B,{C,C},0,[.5 .5;.5 .5]);
            gc = mvst_l2gain_dt_mjls({A1,A2},B,{C,C},0,[.5 .5;.5 .5]);
            tc.verifyGreaterThan(gc, g1 + 0.1);
            r = mvst_hinfnorm(struct('A0',cat(3,A1,A2),'Bw',B,'C',C,'domain','dt', ...
                'Ptr',[.5 .5;.5 .5]));
            tc.verifyEqual(r.class, 'dtmjls');
            tc.verifyEqual(r.gamma, gc, 'RelTol', 1e-3);
        end

        function ctMjlsSmoke(tc)
            s = struct('A0',cat(3,[-1 0.2;0 -2],[-2 0;0.3 -1]),'Bw',[1;0.5], ...
                'C',[1 0],'Pig',[-1 1;2 -2]);
            r = mvst_hinfnorm(s);
            tc.verifyEqual(r.class, 'mjls');
            tc.verifyTrue(isfinite(r.gamma) && r.gamma > 0.2 && r.gamma < 20);
            tc.verifyFalse(r.diagnostics.exact);
        end

        function sampledDataConvergence(tc)
            A = [0 1;-1 -1]; Bw = [0;1]; C = [1 0]; Tbar = 0.3;
            % former fixed default (Nper = 60) under-reports ...
            g60 = mvst_l2gain_sd_hybrid(A,Bw,C,Tbar,60);
            tc.verifyEqual(g60, 2.053, 'AbsTol', 2e-3);
            % ... the adaptive default converges to the SOS/periodic value 2.108
            [g, info] = mvst_l2gain_sd_hybrid(A,Bw,C,Tbar);
            tc.verifyTrue(info.adaptive && info.converged);
            tc.verifyEqual(g, 2.108, 'AbsTol', 3e-3);
            tc.verifyGreaterThan(g, g60 + 0.04);
            tc.verifyGreaterThan(info.Nper, 60);
            r = mvst_hinfnorm(struct('A0',A,'Bw',Bw,'C',C,'Tbar',Tbar));
            tc.verifyEqual(r.class, 'sd');
            tc.verifyEqual(r.gamma, g, 'AbsTol', 1e-9);
            tc.verifyTrue(r.diagnostics.converged);
        end

        function sampledDataExplicitNperUnchanged(tc)
            A = [0 1;-1 -1]; Bw = [0;1]; C = [1 0];
            tc.verifyEqual(mvst_l2gain_sd_hybrid(A,Bw,C,0.3,60), 2.053, 'AbsTol', 2e-3);
            tc.verifyEqual(mvst_l2gain_sd_hybrid(A,Bw,C,0.3,200), 2.102, 'AbsTol', 3e-3);
        end

        function systemMethodAndOptions(tc)
            S = mvst.System.ct([-1 0.5; 0 -2], zeros(2,0), 'Bw', [1;0], 'C', [1 0]);
            r = S.hinfnorm();
            ref = mvst_l2gain_ct([-1 0.5;0 -2],[1;0],[1 0],0);
            tc.verifyEqual(r.gamma, ref, 'AbsTol', 1e-9);
            r2 = S.hinfnorm(mvst.Options('solver','nosdp'));      % Options accepted
            tc.verifyEqual(r2.gamma, ref, 'AbsTol', 1e-9);
            r3 = S.l2gain(mvst.Options());
            tc.verifyEqual(r3.gamma, ref, 'AbsTol', 1e-9);
            tc.verifyError(@() S.hinfnorm(42), 'mvst:System:opts');
        end
    end
end
