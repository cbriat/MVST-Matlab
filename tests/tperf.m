classdef tperf < matlab.unittest.TestCase
%TPERF  Unit tests for the generic performance dispatcher MVST_PERF.
%   A matlab.unittest.TestCase suite (auto-discovered by RUN_ALL_TESTS).
%   Run standalone from the matlab/ directory:
%       /Applications/MATLAB_R2025a.app/bin/matlab -batch "cd tests; runtests('tperf')"
%
%   Coverage (one test method per historical group, so the harness reports
%   per-method pass/fail counts):
%     A1  'l2' equals MVST_L2GAIN on the ct class;
%     A2  'l2' equals MVST_L2GAIN on the sampled-data class;
%     A3  'tracegain' equals the L2-gain squared (and carries .gammaL2);
%     B1  the Gramian measures ('h2','e2p','peak') reproduce the ex_perf.m
%         batch numbers to machine precision;
%     B2  those two exported macros (H2, e2p) match an INDEPENDENT fixed-point
%         ODE Gramian re-derivation to tolerance 1e-6;
%     B3  the scalar closed forms (sec:perf:exscalar) hold;
%     C1  unsupported class-measure pairs raise 'mvst:perf:unsupported' and an
%         unknown measure raises 'mvst:perf:measure'.
%
%   See also MVST_PERF, MVST_L2GAIN, RUN_ALL_TESTS.

    properties
        ctsys
        sdsys
        sysB
        outB
        dispB
        tolx = 1e-6;      % headline reproduction tolerance
        tole = 1e-9;      % engine-identity tolerance (same backend)
    end

    methods (TestClassSetup)
        function addPaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root); addpath(fullfile(root,'lib'));
        end

        function buildFixtures(tc)
            tc.ctsys = struct('A',[-0.5 1; -1 -0.5], 'Bw',[0;1], 'C',[1 0.5]);
            tc.sdsys = tc.ctsys; tc.sdsys.Tbar = 1.2;   % sampled-output class
            tc.sysB  = struct('A',[-0.5 1; -1 -0.5], 'J',0.8*eye(2), 'Tbar',1.2);
            tc.outB  = struct('Bw',[0;1], 'C',[1 0.5]);
            tc.dispB = struct('A',tc.sysB.A, 'J',tc.sysB.J, 'Tbar',tc.sysB.Tbar, ...
                              'Bw',tc.outB.Bw, 'C',tc.outB.C);
        end
    end

    methods (Test)
        % ================================================================= %
        % A. 'l2' delegates to MVST_L2GAIN; 'tracegain' = l2^2
        % ================================================================= %
        function A1_l2_ct(tc)
            r_ct = mvst_l2gain(tc.ctsys);
            p_ct = mvst_perf(tc.ctsys, 'l2');
            tc.verifyEqual(p_ct.class, 'ct', 'l2 ct class tag');
            tc.verifyEqual(p_ct.gamma, r_ct.gamma, 'AbsTol', tc.tole, ...
                'l2 (ct) must equal mvst_l2gain');
        end

        function A2_l2_sd(tc)
            r_sd = mvst_l2gain(tc.sdsys);
            p_sd = mvst_perf(tc.sdsys, 'l2');
            tc.verifyEqual(p_sd.class, 'sd', 'l2 sd class tag');
            tc.verifyEqual(p_sd.gamma, r_sd.gamma, 'AbsTol', tc.tole, ...
                'l2 (sd) must equal mvst_l2gain');
        end

        function A3_tracegain(tc)
            r_ct = mvst_l2gain(tc.ctsys);
            p_tg = mvst_perf(tc.ctsys, 'tracegain');
            tc.verifyEqual(p_tg.gamma, r_ct.gamma^2, 'AbsTol', tc.tole, ...
                'tracegain must equal l2^2');
            tc.verifyEqual(p_tg.gammaL2, r_ct.gamma, 'AbsTol', tc.tole, ...
                'tracegain must carry .gammaL2');
        end

        % ================================================================= %
        % B. Gramian measures on the ex_perf.m impulsive system
        % ================================================================= %
        function B1_dispatcher_equals_batch(tc)
            L    = mvst_perf(tc.sysB, tc.outB, struct('npts',200));
            h2d  = mvst_perf(tc.dispB, 'h2');
            e2pd = mvst_perf(tc.dispB, 'e2p');
            pkd  = mvst_perf(tc.dispB, 'peak');
            tc.verifyEqual(h2d.class, 'impulsive', 'h2 impulsive class tag');
            tc.verifyEqual(h2d.gamma,  L.H2,          'AbsTol', tc.tole, 'h2 dispatcher vs batch');
            tc.verifyEqual(e2pd.gamma, L.energy2peak, 'AbsTol', tc.tole, 'e2p dispatcher vs batch');
            tc.verifyEqual(pkd.gamma,  L.e2psq,       'AbsTol', tc.tole, 'peak = lmax = e2p^2');
        end

        function B2_independent_rederivation(tc)
            % Recompute the batch and dispatcher results this method needs.
            L    = mvst_perf(tc.sysB, tc.outB, struct('npts',200));
            h2d  = mvst_perf(tc.dispB, 'h2');
            e2pd = mvst_perf(tc.dispB, 'e2p');

            % Independent periodic controllability Gramian by fixed-point ODE
            % integration (route independent of the van-Loan/Kronecker engine).
            A = tc.sysB.A; J = tc.sysB.J; Tbar = tc.sysB.Tbar;
            Bw = tc.outB.Bw; C = tc.outB.C;
            n = size(A,1); Q = Bw*Bw.';
            oo  = odeset('RelTol',1e-11,'AbsTol',1e-13);
            rhs = @(~,p) reshape(A*reshape(p,n,n) + reshape(p,n,n)*A.' + Q, [], 1);

            Pi0 = zeros(n);
            for it = 1:5000
                [~,Y] = ode45(rhs, [0 Tbar], Pi0(:), oo);
                PiT = reshape(Y(end,:).', n, n);
                Pn  = J*PiT*J.';  Pn = (Pn + Pn.')/2;
                if norm(Pn - Pi0, 'fro') < 1e-13; Pi0 = Pn; break; end
                Pi0 = Pn;
            end
            tc.verifyLessThanOrEqual(norm(Pi0 - L.Pi, 'fro'), tc.tolx, ...
                'independent Pi0 vs engine Pi0');

            % measures on the same 200-point period grid
            tau = linspace(0, Tbar, 200);
            [~,Yg] = ode45(rhs, tau, Pi0(:), oo);
            pw = zeros(size(tau)); lam = zeros(size(tau));
            for i = 1:numel(tau)
                P = reshape(Yg(i,:).', n, n);
                S = C*P*C.';  S = (S + S.')/2;
                pw(i)  = trace(S);
                lam(i) = max(real(eig(S)));
            end
            H2_ref  = sqrt(trapz(tau, pw)/Tbar);
            e2p_ref = sqrt(max(lam));

            tc.verifyEqual(H2_ref,  h2d.gamma,  'AbsTol', tc.tolx, ...
                'independent H2 vs dispatcher');
            tc.verifyEqual(e2p_ref, e2pd.gamma, 'AbsTol', tc.tolx, ...
                'independent e2p vs dispatcher');
        end

        function B3_scalar_closed_form(tc)
            % ct closed forms (prop:perf:h2classic, sec:perf:exscalar):
            % a=-0.7,b=1.3,c=0.9 -> J_H2 = |cb|/sqrt(2|a|), l2 = |cb|/|a|.
            a = -0.7; b = 1.3; c = 0.9;
            scal = struct('A',a, 'Bw',b, 'C',c);
            h2s  = mvst_perf(scal, 'h2');
            e2ps = mvst_perf(scal, 'e2p');
            l2s  = mvst_perf(scal, 'l2');
            H2_closed = abs(c*b)/sqrt(2*abs(a));
            L2_closed = abs(c*b)/abs(a);
            tc.verifyEqual(h2s.gamma,  H2_closed, 'AbsTol', tc.tolx, 'scalar H2 closed form');
            tc.verifyEqual(e2ps.gamma, H2_closed, 'AbsTol', tc.tolx, 'scalar e2p = H2 closed form');
            % L2-gain is computed by bisection (tol ~1e-4), cross-check looser
            tc.verifyEqual(l2s.gamma, L2_closed, 'AbsTol', 1e-3, 'scalar L2 closed form');
        end

        % ================================================================= %
        % C. Unsupported pairs and unknown measures
        % ================================================================= %
        % ================================================================= %
        % D. Multi-output L2-gain (regression: R block sized by disturbance
        %    count, not output count).  Historically MVST_L2GAIN_CT sized
        %    R = gamma^2 I - Dw'Dw by rows(C); a multi-row C with a column Dw
        %    then errored at R\(Dw'*C).  Cross-checked against a fine
        %    frequency sweep of sigma_max(C(jwI-A)^-1 Bw + Dw).
        % ================================================================= %
        function D1_multi_output_l2gain_ct(tc)
            A  = [-2 1 0; 0 -3 1; 1 0 -4];
            C  = [1 0 0.5; 0 1 -0.5];              % 2 outputs, 3 states
            n  = size(A,1);
            w  = logspace(-3, 3, 200000);

            % (i) 2-output, 1-disturbance with a COLUMN Dw (the failing shape)
            Bw = [1; 0.5; -1];  Dw = [0.2; -0.3];
            g1 = mvst_l2gain_ct(A, Bw, C, Dw);     % must not error
            pk1 = tc.sweepPeak(A, Bw, C, Dw, w, n);
            tc.verifyEqual(g1, pk1, 'AbsTol', 1e-4, ...
                '2-output/1-disturbance gain vs frequency sweep');

            % (ii) genuine 2x2 Dw (2 outputs, 2 disturbances)
            Bw2 = [1 0; 0.5 1; -1 0.5];  Dw2 = [0.2 0.1; -0.3 0.4];
            g2  = mvst_l2gain_ct(A, Bw2, C, Dw2);
            pk2 = tc.sweepPeak(A, Bw2, C, Dw2, w, n);
            tc.verifyEqual(g2, pk2, 'AbsTol', 1e-4, ...
                '2-output/2-disturbance gain vs frequency sweep');

            % (iii) scalar-0 Dw shorthand still accepted on a multi-output C
            g3 = mvst_l2gain_ct(A, Bw, C, 0);      % expands to zeros(2,1)
            pk3 = tc.sweepPeak(A, Bw, C, zeros(2,1), w, n);
            tc.verifyEqual(g3, pk3, 'AbsTol', 1e-4, ...
                'scalar-0 Dw shorthand, multi-output');

            % (iv) a genuinely wrong Dw size is rejected
            tc.verifyError(@() mvst_l2gain_ct(A, Bw, C, [0.2 0.1]), ...
                'mvst:l2gain_ct:dim');
        end

        function C1_unsupported_and_unknown(tc)
            % (Deterministic peak-to-peak 'p2p' and the dwell tubes are now
            % SUPPORTED via the SDP/SOS backend MVST_PERF_DWELL; they are covered
            % in TPERFDWELL, which adds the YALMIP/SOSTOOLS path.  This suite
            % stays on the exact engines.)
            % Gramian measure on a multi-mode (mjls) system: no covariance stated
            mjls = struct('A0', cat(3,[-1 0;0 -2],[-3 1;0 -1]), 'Bw',[0;1], ...
                          'C',[1 0.5], 'Pig',[-1 1; 2 -2]);
            tc.verifyError(@() mvst_perf(mjls,'e2p'), 'mvst:perf:unsupported');
            % forcing an unsupported class
            tc.verifyError(@() mvst_perf(tc.dispB,'h2',struct('class','poisson')), ...
                'mvst:perf:unsupported');
            % unknown measure name
            tc.verifyError(@() mvst_perf(tc.ctsys,'bogus'), 'mvst:perf:measure');
        end
    end

    methods (Static)
        % Independent H-infinity norm by a fine frequency sweep of the
        % largest singular value of G(jw) = C(jwI-A)^{-1} Bw + Dw.
        function pk = sweepPeak(A, Bw, C, Dw, w, n)
            pk = 0;
            for k = 1:numel(w)
                G = C*((1i*w(k)*eye(n) - A)\Bw) + Dw;
                s = norm(G, 2);
                if s > pk, pk = s; end
            end
        end
    end
end
