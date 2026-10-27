classdef tgamedwc < matlab.unittest.TestCase
%TGAMEDWC  Regression tests for the disturbance-multiplicative channel Dwc in
%   the COUPLED Markov-jump game backend MVST_GAME_GAMMA_CT_MJLS (game counterpart of TDWCMODES / TDWCSYN).  A latent boundedness
%   defect accepted NaN-limping ode45 runs as "bounded" once Dwc was present,
%   collapsing the gamma-bisection to the floor (~2e-3); the fix adds an
%   isfinite guard to the feasibility test, mirroring MVST_GAME_GAMMA_CT.  This
%   suite pins:
%     (1) BACKWARD COMPAT: the diffusion-free coupled critical level is
%         unchanged (bit-stable on the noise-free path);
%     (2) Dwc is CARRIED and SANE: with the disturbance-multiplicative channel
%         the coupled level is finite, well above the floor, and NOT below the
%         noise-free level (a diffusion channel degrades attenuation, so it can
%         only raise the critical disturbance level);
%     (3) SENSITIVITY: perturbing Dwc moves the coupled level.
%
%   See also MVST_GAME_GAMMA_CT_MJLS, MVST_GAME_GAMMA_CT, TDWCSYN, TDWCMODES.

    methods (TestClassSetup)
        function addPaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root); addpath(fullfile(root, 'lib'));
        end
    end

    properties (Constant)
        A1  = [0 1;-2 -1];   A2 = [0 1;-3 -1.5];
        Bu  = [0;1];         Bw = [1;0];
        Pig = [-1 1; 2 -2];
    end

    methods (Test)
        % (1) BACKWARD COMPAT: diffusion-free coupled level unchanged.
        function coupledNoChanUnchanged(tc)
            g = mvst_game_gamma_ct_mjls({tc.A1,tc.A2}, tc.Bu, tc.Bw, ...
                                        eye(2), 1, tc.Pig, 30);
            tc.verifyEqual(g, 1.3343, 'AbsTol', 5e-3, ...
                'Noise-free coupled MJLS game level must be unchanged.');
        end

        % (2) Dwc CARRIED and SANE: finite, above the floor, >= noise-free level.
        function coupledDwcSane(tc)
            Ci = 0.2*eye(2); Dwc = {0.1*[1;0]};
            g0 = mvst_game_gamma_ct_mjls({tc.A1,tc.A2}, tc.Bu, tc.Bw, ...
                                         eye(2), 1, tc.Pig, 30);
            g  = mvst_game_gamma_ct_mjls({tc.A1,tc.A2}, tc.Bu, tc.Bw, ...
                                         eye(2), 1, tc.Pig, 30, {{Ci},{Ci}}, Dwc);
            tc.verifyGreaterThan(g, 0.1, ...
                'Dwc coupled level must not collapse to the bisection floor.');
            tc.verifyGreaterThanOrEqual(g, g0 - 5e-3, ...
                'A diffusion channel can only raise the critical level.');
            tc.verifyLessThan(g, 10, 'Coupled level must stay finite/sane.');
        end

        % (3) SENSITIVITY: perturbing Dwc moves the coupled level.
        function coupledDwcSensitivity(tc)
            Ci = 0.2*eye(2);
            g1 = mvst_game_gamma_ct_mjls({tc.A1,tc.A2}, tc.Bu, tc.Bw, ...
                     eye(2), 1, tc.Pig, 30, {{Ci},{Ci}}, {0.1*[1;0]});
            g3 = mvst_game_gamma_ct_mjls({tc.A1,tc.A2}, tc.Bu, tc.Bw, ...
                     eye(2), 1, tc.Pig, 30, {{Ci},{Ci}}, {0.3*[1;0]});
            tc.verifyGreaterThan(abs(g3-g1), 1e-2, ...
                'Perturbing Dwc must change the coupled game level.');
        end
    end
end
