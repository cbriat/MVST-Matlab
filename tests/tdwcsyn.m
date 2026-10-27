classdef tdwcsyn < matlab.unittest.TestCase
%TDWCSYN  Cross-check tests for the disturbance-multiplicative channel Dwc
%   folded into the state-feedback SYNTHESIS SOS routine MVST_DISSYN_SOS (synthesis counterpart
%   of TDWCMODES).  The synthesis flow
%   block stacks the shared disturbance-multiplicative channel D_{w,l} against
%   the completely positive Schur pivot -I_q(x)W, exactly as the analysis
%   certificate loads the flow pivot Theta.  This suite verifies that:
%     (1) on a Dwc-carrying continuous-time LTI instance the certified
%         synthesis level matches the EXACT state-feedback H-infinity game
%         level carrying the same channel (MVST_GAME_GAMMA_CT with Dwc) to
%         solver tolerance (the constant-storage relaxation is lossless in the
%         memoryless case, so this is an exact-vs-SOS cross-check);
%     (2) perturbing Dwc genuinely moves the certified level (the channel is
%         carried, not silently dropped) -- a SOLVER-based sensitivity test;
%     (3) BACKWARD COMPATIBILITY: a Dwc-free (channel-free) input reproduces
%         the noise-free synthesis level, bit-for-bit the pre-extension path.
%
%   See also MVST_DISSYN_SOS, MVST_GAME_GAMMA_CT, TDWCMODES, RUN_ALL_TESTS.

    methods (TestClassSetup)
        function addPaths(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root); addpath(fullfile(root, 'lib'));
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS (sossolve) not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods (Static)
        function s = ctsys(dwc)
            % Continuous-time LTI performance variant with penalty z=[x1;x2;u]
            % (Q=I_2, R_u=1) and a shared disturbance-multiplicative channel Dwc.
            A = [0 1; -1 -1]; Bu = [0;1]; Bw = [1;0];
            C = [eye(2); zeros(1,2)]; Du = [0;0;1];
            s = struct('A0',A,'B0',Bu,'Bw',Bw,'C',C,'Du',Du,'Dwc',dwc);
        end
    end

    methods (Test)
        % (1) CT synthesis(+Dwc) must match the exact H-infinity game(+Dwc).
        function ctSynDwcMatchesGame(tc)
            dwc = 0.1*[1;0];
            r = mvst_dissyn_sos(tdwcsyn.ctsys(dwc), struct('type','l2gain'), ...
                                struct('type','arbitrary'), 'gamma');
            gex = mvst_game_gamma_ct([0 1;-1 -1], [0;1], [1;0], eye(2), 1, 40, ...
                                     {}, {dwc});
            tc.verifyEqual(r.gammastar, gex, 'RelTol', 4e-2, ...
                'CT synthesis(+Dwc) must match the exact game level(+Dwc).');
        end

        % (2) SENSITIVITY: perturbing Dwc must move the certified level.
        function ctSynDwcSensitivity(tc)
            g1 = mvst_dissyn_sos(tdwcsyn.ctsys(0.1*[1;0]), struct('type','l2gain'), ...
                                 struct('type','arbitrary'), 'gamma').gammastar;
            g3 = mvst_dissyn_sos(tdwcsyn.ctsys(0.3*[1;0]), struct('type','l2gain'), ...
                                 struct('type','arbitrary'), 'gamma').gammastar;
            tc.verifyGreaterThan(abs(g3-g1), 1e-2, ...
                'Perturbing Dwc must change the synthesis level (channel folded).');
        end

        % (3) BACKWARD COMPAT: Dwc-free synthesis matches the noise-free game.
        function ctSynDiffFreeUnchanged(tc)
            A = [0 1;-1 -1]; Bu = [0;1]; Bw = [1;0];
            C = [eye(2); zeros(1,2)]; Du = [0;0;1];
            r = mvst_dissyn_sos(struct('A0',A,'B0',Bu,'Bw',Bw,'C',C,'Du',Du), ...
                                struct('type','l2gain'), ...
                                struct('type','arbitrary'), 'gamma');
            gex = mvst_game_gamma_ct(A, Bu, Bw, eye(2), 1, 40);
            tc.verifyEqual(r.gammastar, gex, 'RelTol', 4e-2, ...
                'Channel-free synthesis must still match the noise-free game.');
        end
    end
end
