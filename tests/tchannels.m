classdef tchannels < matlab.unittest.TestCase
%TCHANNELS  Unit tests for the channel port and the critical-parameter
%   wrapper (auto-discovered by RUN_ALL_TESTS; solver-free).
%
%   Covers:
%     Channels the Brownian-multiplicative channels (Ac state-mult, Dwc
%              disturbance-mult) ported into the exact L2-gain and game backends
%              (mvst_l2gain_ct / _ct_mjls / _sd_hybrid and mvst_game_gamma_ct /
%              _ct_mjls), cross-checked against the independent Poissonian-jump
%              reference backend (with no jumps) and the collapse of a one-mode
%              MJLS to the ct case;
%     Critical mvst_stab_critical: bisection on the mean-square stability margin
%              over a named scalar parameter (dwell-time and a setter-rebuilt
%              noise scale), plus the no-bracket error.
%
%   Backward compatibility (diffusion-free inputs reproduce the pre-port
%   numbers) is guarded by the hard-coded expected values in TSYSTEM / TPERF /
%   TSTABILITY and by the golden regression, not here.
%
%   See also MVST_L2GAIN, MVST_GAME_GAMMA, MVST_STAB_CRITICAL, RUN_ALL_TESTS.

    methods (TestClassSetup)
        function addPaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root); addpath(fullfile(root, 'lib'));
        end
    end

    methods (Test)
        % ================================================================= %
        %  continuous-time diffusion L2-gain fold
        % ================================================================= %
        function ctDiffusionMatchesMjlsOneMode(tc)
            % INDEPENDENT cross-check: the ct diffusion fold is compared against
            % the one-mode MJLS coupled-BRE backend (a separate implementation,
            % NOT the Poisson backend the ct path routes through), so the check
            % is not a tautology.
            A = [0 1; -2 -1]; Bw = [0; 1]; C = [1 0]; Dw = 0.5;
            Cc = {0.15 * eye(2)}; Dwc = {[0.1; 0]};
            g    = mvst_l2gain_ct(A, Bw, C, Dw, Cc, Dwc);
            gref = mvst_l2gain_ct_mjls({A}, Bw, {C}, Dw, 0, 60, {Cc}, Dwc);
            tc.verifyEqual(g, gref, 'RelTol', 3e-3, ...
                'ct diffusion gain must equal the one-mode MJLS backend.');
            tc.verifyTrue(isfinite(g) && g > 0);
        end

        function ctDispatcherRoutesDiffusion(tc)
            % A single-mode ct struct carrying .Ac/.Dwc must fold the noise;
            % cross-checked against the independent one-mode MJLS backend.
            A = [0 1; -2 -1]; Bw = [0; 1]; C = [1 0]; Dw = 0.5;
            Cc = 0.15 * eye(2); Dwc = [0.1; 0];
            sys = struct('A', A, 'Bw', Bw, 'C', C, 'Dw', Dw, 'Ac', Cc, 'Dwc', Dwc);
            r = mvst_l2gain(sys);
            gref = mvst_l2gain_ct_mjls({A}, Bw, {C}, Dw, 0, 60, {{Cc}}, {Dwc});
            tc.verifyEqual(r.method, 'ct');
            tc.verifyEqual(r.gamma, gref, 'RelTol', 3e-3);
        end

        % ================================================================= %
        %  one-mode MJLS diffusion vs the no-jump Poisson backend
        %  (grounds the MJLS backend used as the ct reference in an independent
        %   third implementation)
        % ================================================================= %
        function mjlsDiffusionMatchesPoisson(tc)
            A = [0 1; -2 -1]; Bw = [0; 1]; C = [1 0]; Dw = 0.5;
            Cc = 0.15 * eye(2); Dwc = [0.1; 0];
            gm = mvst_l2gain_ct_mjls({A}, Bw, {C}, Dw, 0, 60, {{Cc}}, {Dwc});
            gp = mvst_l2gain_ct_poisson(A, Bw, C, Dw, {Cc}, {Dwc}, {}, [], {});
            tc.verifyEqual(gm, gp, 'RelTol', 3e-3, ...
                'one-mode MJLS diffusion gain must equal the no-jump Poisson gain.');
        end

        % ================================================================= %
        %  Dispatcher routing of .Ac/.Dwc for mjls and sd structs
        % ================================================================= %
        function mjlsDispatcherRoutesDiffusion(tc)
            % A multi-mode struct carrying .Ac (n x n x q x M) and .Dwc must be
            % routed to the mjls backend with the per-mode diffusion folded;
            % verified against the direct backend call.
            A1 = [-1 1; 0 -2]; A2 = [-2 0; 1 -1]; Bw = [0; 1]; C = [1 0]; Dw = 0.3;
            Pig = [-1 1; 2 -2];
            Ac4 = cat(4, 0.15*eye(2), 0.1*eye(2));   % n x n x 1 x 2 (one channel)
            Dwc = [0.05; 0];
            sys = struct('A0', cat(3,A1,A2), 'Bw', Bw, 'C', cat(3,C,C), 'Dw', Dw, ...
                         'Pig', Pig, 'Ac', Ac4, 'Dwc', Dwc);
            r = mvst_l2gain(sys);
            Chan = {{Ac4(:,:,1,1)}, {Ac4(:,:,1,2)}};
            gref = mvst_l2gain_ct_mjls({A1,A2}, Bw, {C,C}, Dw, Pig, [], Chan, {Dwc});
            tc.verifyEqual(r.method, 'mjls');
            tc.verifyEqual(r.gamma, gref, 'RelTol', 3e-3);
        end

        function sdDispatcherRoutesDiffusion(tc)
            % A single-mode struct with .Tbar and .Ac/.Dwc must be routed to the
            % sd backend with the flow diffusion folded; verified against the
            % direct backend call.
            A = [0 1; -1 -1]; Bw = [0; 1]; C = [1 0]; Tbar = 0.3;
            Ac = 0.2*eye(2); Dwc = [0.1; 0];
            sys = struct('A', A, 'Bw', Bw, 'C', C, 'Tbar', Tbar, 'Ac', Ac, 'Dwc', Dwc);
            r = mvst_l2gain(sys);
            gref = mvst_l2gain_sd_hybrid(A, Bw, C, Tbar, [], {Ac}, {Dwc});
            tc.verifyEqual(r.method, 'sd');
            tc.verifyEqual(r.gamma, gref, 'RelTol', 3e-3);
        end

        % ================================================================= %
        %  sampled-data hybrid diffusion raises the mean-square gain
        % ================================================================= %
        function sdDiffusionRaisesGain(tc)
            A = [0 1; -1 -1]; Bw = [0; 1]; C = [1 0]; Tbar = 0.3;
            g0 = mvst_l2gain_sd_hybrid(A, Bw, C, Tbar);
            g1 = mvst_l2gain_sd_hybrid(A, Bw, C, Tbar, [], {0.2 * eye(2)}, {});
            tc.verifyTrue(isfinite(g1) && g1 > 0);
            tc.verifyGreaterThanOrEqual(g1, g0 - 1e-3, ...
                'flow diffusion must not lower the hybrid mean-square gain.');
        end

        % ================================================================= %
        %  game (HJI) diffusion fold, independent cross-checks
        % ================================================================= %
        function gameCtDiffusionMatchesMjlsOneMode(tc)
            % INDEPENDENT cross-check: the ct game DRE (mvst_game_dre_ct) fold is
            % compared against the one-mode MJLS game backend, which uses its own
            % coupled-Riccati rhs (a separate implementation), not the same path.
            % The channel is sized so a wrong fold shifts the level well beyond
            % the tolerance (the two correct backends agree to ~1e-4): a 1.5x fold
            % perturbation was verified to move the ct level by ~4% > RelTol.
            A = [-0.5 1; -1 -0.5]; Bu = [0; 1]; Bw = [0; 1];
            Q = eye(2); Ru = 1; Cc = {0.6 * eye(2)};
            g    = mvst_game_gamma_ct(A, Bu, Bw, Q, Ru, 30, Cc);
            gref = mvst_game_gamma_ct_mjls({A}, Bu, Bw, Q, Ru, 0, 30, {Cc});
            tc.verifyEqual(g, gref, 'RelTol', 5e-3, ...
                'ct game diffusion level must equal the one-mode MJLS game level.');
            tc.verifyTrue(isfinite(g) && g > 0);
        end

        function gameMjlsDiffusionMatchesPoisson(tc)
            % Grounds the MJLS game backend in the independent no-jump Poisson game.
            A = [-0.5 1; -1 -0.5]; Bu = [0; 1]; Bw = [0; 1];
            Q = eye(2); Ru = 1; Cc = 0.6 * eye(2);
            gm = mvst_game_gamma_ct_mjls({A}, Bu, Bw, Q, Ru, 0, 30, {{Cc}});
            gp = mvst_game_gamma_ct_poisson(A, Bu, Bw, {Cc}, {}, [], Q, Ru, 30);
            tc.verifyEqual(gm, gp, 'RelTol', 5e-3, ...
                'one-mode MJLS game diffusion level must equal the no-jump Poisson game.');
        end

        function gameDwcRaisesLevel(tc)
            % Disturbance-multiplicative noise cannot lower the critical level.
            A = [-0.5 1; -1 -0.5]; Bu = [0; 1]; Bw = [0; 1];
            Q = eye(2); Ru = 1; Cc = {0.2 * eye(2)};
            g0 = mvst_game_gamma_ct(A, Bu, Bw, Q, Ru, 30, Cc, {zeros(2,1)});
            g1 = mvst_game_gamma_ct(A, Bu, Bw, Q, Ru, 30, Cc, {[0.1; 0]});
            tc.verifyTrue(isfinite(g1));
            tc.verifyGreaterThanOrEqual(g1, g0 - 1e-2);
        end

        % ================================================================= %
        %  Critical parameter on the stability margin
        % ================================================================= %
        function critFindsImpulsiveDwellTime(tc)
            % Unstable flow contracted by a reset at each dwell: stable for small
            % dwell-time, unstable for large.  The critical dwell-time zeroes the
            % monodromy log-radius.
            sys = struct('A', [0.3 1; 0 0.3], 'J', 0.5 * eye(2));
            o   = struct('class', 'impulsive');
            res = mvst_stab_critical(sys, 'T', [0.01 5], o);
            tc.verifyEqual(res.margin_at, 0, 'AbsTol', 1e-4, ...
                'margin at the critical dwell-time must be ~0.');
            % stability flips across the critical value.
            below = mvst_stability(sys, struct('class','impulsive','T',res.crit*0.9));
            above = mvst_stability(sys, struct('class','impulsive','T',res.crit*1.1));
            tc.verifyTrue(below.stable);
            tc.verifyFalse(above.stable);
        end

        function critNoBracketErrors(tc)
            sys = struct('A', [0.3 1; 0 0.3], 'J', 0.5 * eye(2));
            o   = struct('class', 'impulsive');
            % both endpoints small -> both stable -> no crossing bracketed.
            tc.verifyError(@() mvst_stab_critical(sys, 'T', [0.01 0.02], o), ...
                'mvst:stab_critical:nobracket');
        end

        function critSetterSweepsNoiseScale(tc)
            % Sweep a state-multiplicative noise scale via a builder: a Hurwitz
            % drift is destabilized in mean square as the noise grows, so the
            % lifted spectral abscissa crosses zero at a critical scale.
            A = [0 1; -1 -0.5]; N = [0 0; 1 0];
            build = @(~, c) mvst.System.ct(A, zeros(2,0), 'Ac', c * N);
            o = struct('setter', build);
            res = mvst_stab_critical([], 'noise', [0 5], o);
            tc.verifyEqual(res.margin_at, 0, 'AbsTol', 1e-4);
            tc.verifyGreaterThan(res.crit, 0);
            % check the flip
            lo = build([], res.crit * 0.9); hi = build([], res.crit * 1.1);
            tc.verifyTrue(lo.stability().stable);
            tc.verifyFalse(hi.stability().stable);
        end
    end
end
