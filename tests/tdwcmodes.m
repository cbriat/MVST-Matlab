classdef tdwcmodes < matlab.unittest.TestCase
%TDWCMODES  Cross-check tests for the disturbance-multiplicative channel Dwc
%   folded into the PRODUCT-CONE SOS routine MVST_DISS_SOS_MODES.
%   The exact coupled/dwell backends (MVST_L2GAIN_CT_MJLS, MVST_L2GAIN_CT_DWELL)
%   carry both the state channel .Ac and the shared disturbance channel .Dwc;
%   this suite verifies that the SOS product-cone certificate now agrees with
%   the exact value on a Dwc-carrying MJLS instance (within SOS bisection
%   tolerance), and proves the channel is genuinely carried by a perturbation
%   sensitivity test (SOLVER-based: SOSTOOLS + SeDuMi).
%
%   See also MVST_DISS_SOS_MODES, MVST_L2GAIN_CT_MJLS, MVST_L2GAIN_CT_DWELL,
%   TCHANNELS, RUN_ALL_TESTS.

    methods (TestClassSetup)
        function addPaths(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root); addpath(fullfile(root, 'lib'));
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS (sossolve) not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods (Test)
        % MJLS coupled: SOS(+Dwc) must match the exact coupled BRE (+Dwc).
        function mjlsCoupledDwcMatchesExact(tc)
            A1 = [1 2; 0 1]; A2 = [-3 0; 0 -3]; Bw = [1; 0];
            C1 = [1 0]; C2 = [0.5 1]; Dw = 0.5; Pig = [-4 4; 1 -1];
            Chan = {0.2*eye(2), 0.2*eye(2)}; Dwc = {0.1*[1;0]};
            gex = mvst_l2gain_ct_mjls({A1,A2}, Bw, {C1,C2}, Dw, Pig, [], Chan, Dwc);
            Ac4 = zeros(2,2,1,2); Ac4(:,:,1,1) = 0.2*eye(2); Ac4(:,:,1,2) = 0.2*eye(2);
            gsos = mvst_diss_sos_modes(struct('A0',cat(3,A1,A2),'Bw',Bw, ...
                'C',cat(3,C1,C2),'Dw',Dw,'Pig',Pig,'Ac',Ac4,'Dwc',0.1*[1;0]), ...
                struct('type','l2gain'), struct('type','ctmjls'), 'gamma').gamma;
            tc.verifyEqual(gsos, gex, 'RelTol', 6e-3, ...
                'MJLS coupled SOS(+Dwc) must match the exact coupled BRE.');
        end

        % SENSITIVITY: perturbing Dwc must move the SOS value (channel carried).
        function mjlsDwcSensitivity(tc)
            A1 = [1 2; 0 1]; A2 = [-3 0; 0 -3]; Bw = [1; 0];
            C1 = [1 0]; C2 = [0.5 1]; Dw = 0.5; Pig = [-4 4; 1 -1];
            Ac4 = zeros(2,2,1,2); Ac4(:,:,1,1) = 0.2*eye(2); Ac4(:,:,1,2) = 0.2*eye(2);
            base = @(dw) mvst_diss_sos_modes(struct('A0',cat(3,A1,A2),'Bw',Bw, ...
                'C',cat(3,C1,C2),'Dw',Dw,'Pig',Pig,'Ac',Ac4,'Dwc',dw), ...
                struct('type','l2gain'), struct('type','ctmjls'), 'gamma').gamma;
            g1 = base(0.1*[1;0]); g3 = base(0.3*[1;0]);
            tc.verifyGreaterThan(abs(g3-g1), 1e-2, ...
                'Perturbing Dwc must change the SOS gain (Dwc genuinely folded).');
        end

        % BACKWARD COMPAT: diffusion-free modes SOS unchanged vs exact.
        function mjlsDiffFreeUnchanged(tc)
            A1 = [1 2; 0 1]; A2 = [-3 0; 0 -3]; Bw = [1; 0];
            C1 = [1 0]; C2 = [0.5 1]; Dw = 0.5; Pig = [-4 4; 1 -1];
            gex = mvst_l2gain_ct_mjls({A1,A2}, Bw, {C1,C2}, Dw, Pig);
            gsos = mvst_diss_sos_modes(struct('A0',cat(3,A1,A2),'Bw',Bw, ...
                'C',cat(3,C1,C2),'Dw',Dw,'Pig',Pig), ...
                struct('type','l2gain'), struct('type','ctmjls'), 'gamma').gamma;
            tc.verifyEqual(gsos, gex, 'RelTol', 6e-3, ...
                'Diffusion-free coupled SOS must still match the exact BRE.');
        end
    end
end
