classdef tsystem < matlab.unittest.TestCase
%TSYSTEM  Unit tests for the +mvst package.
%   Run from the matlab/tests directory:
%       results = runtests('tsystem');
%
%   Covers: every static factory (good and bad inputs), the FROMSTRUCT
%   round trip (thin mode and thin system reproduced field-for-field), the
%   RICH bridge against MVST_TORICH, and the L2-gain method sugar against a
%   direct MVST_L2GAIN call on the equivalent struct.

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root);                          % expose the +mvst package
            addpath(fullfile(root, 'lib'));         % expose the engine layer
        end
    end

    methods (Test)
        % ---------------------------------------------------------------- %
        %  Factories: construction and derived properties                   %
        % ---------------------------------------------------------------- %
        function ctConstructs(tc)
            S = mvst.System.ct([-1 0.5; 0 -2], [1; 1], 'Q', eye(2), 'R', 1);
            tc.verifyEqual(S.n, 2);
            tc.verifyEqual(S.m, 1);
            tc.verifyEqual(S.M, 1);
            tc.verifyEqual(S.tag, 'ct');
        end

        function ctRejectsNonSquare(tc)
            tc.verifyError(@() mvst.System.ct([1 2 3; 4 5 6], [1; 1]), ...
                'mvst:System:notSquare');
        end

        function ctRejectsInputDim(tc)
            tc.verifyError(@() mvst.System.ct([-1 0; 0 -1], [1; 1; 1]), ...
                'mvst:System:dim');
        end

        function mjlsConstructs(tc)
            Pi = [-1 1; 2 -2];
            S = mvst.System.mjls({[-1 0; 0 -2], [-3 1; 0 -1]}, ...
                                 {[1; 0], [0; 1]}, Pi, 'Q', eye(2), 'R', 1);
            tc.verifyEqual(S.M, 2);
            tc.verifyEqual(S.tag, 'mjls');
        end

        function mjlsRejectsBadGeneratorSize(tc)
            tc.verifyError(@() mvst.System.mjls({[-1 0; 0 -2], [-3 1; 0 -1]}, ...
                {[1; 0], [0; 1]}, [-1 1 0; 0 -1 1]), 'mvst:System:generator');
        end

        function mjlsRejectsNonMetzler(tc)
            % negative off-diagonal entry violates the Metzler property
            tc.verifyError(@() mvst.System.mjls({[-1 0; 0 -2], [-3 1; 0 -1]}, ...
                {[1; 0], [0; 1]}, [-1 -1; 2 -2]), 'mvst:System:generator');
        end

        function mjlsRejectsNonzeroRowsum(tc)
            tc.verifyError(@() mvst.System.mjls({[-1 0; 0 -2], [-3 1; 0 -1]}, ...
                {[1; 0], [0; 1]}, [-1 2; 2 -2]), 'mvst:System:generator');
        end

        function mjlsRejectsListLength(tc)
            tc.verifyError(@() mvst.System.mjls({[-1 0; 0 -2]}, ...
                {[1; 0], [0; 1]}, [-1 1; 2 -2]), 'mvst:System:dim');
        end

        function sampledConstructs(tc)
            S = mvst.System.sampled([-1 1; 0 -1], [0; 1], eye(2), 0.2);
            tc.verifyEqual(S.tag, 'sampled');
            tc.verifyEqual(S.modes{1}.Tbar, 0.2);
        end

        function sampledRejectsNonpositivePeriod(tc)
            tc.verifyError(@() mvst.System.sampled([-1 0; 0 -1], [0; 1], eye(2), -1), ...
                'MATLAB:validators:mustBePositive');
        end

        function poissonConstructs(tc)
            S = mvst.System.poisson([-1 0; 0 -1], [1; 0], 0.5 * eye(2), 2);
            tc.verifyEqual(S.tag, 'poisson');
            tc.verifyEqual(S.modes{1}.lamf, 2);
        end

        function poissonRejectsRateMismatch(tc)
            J = cat(3, 0.5 * eye(2), 0.7 * eye(2));      % two jump channels
            tc.verifyError(@() mvst.System.poisson([-1 0; 0 -1], [1; 0], J, 2), ...
                'mvst:System:dim');                       % only one rate given
        end

        function switchedConstructs(tc)
            S = mvst.System.switched({[-1 0; 0 -2], [-2 1; 0 -1]}, ...
                                     {[1; 0], [0; 1]}, 'Td', 0.5);
            tc.verifyEqual(S.M, 2);
            tc.verifyEqual(S.tag, 'switched');
        end

        function impulsiveConstructs(tc)
            S = mvst.System.impulsive([0 1; -1 0], [0; 1], 0.5 * eye(2));
            tc.verifyEqual(S.tag, 'impulsive');
        end

        % ---------------------------------------------------------------- %
        %  Impulsive is the full-channel master factory                  %
        % ---------------------------------------------------------------- %
        function impulsiveFullChannelsRoundtrip(tc)
            % Every channel of the general impulsive class, u and w carried.
            A = [0 1; -1 0]; B = [0; 1]; J = 0.5 * eye(2);
            S = mvst.System.impulsive(A, B, J, ...
                'Q', eye(2), 'R', 1, 'Ac', 0.1 * eye(2), 'Bc', [0; 0.2], ...
                'Dwc', [0.1; 0], 'Bw', [0; 1], 'C', [1 0], 'Dw', 0.3, ...
                'Jf', 1.1 * eye(2), 'lamf', 0.5, 'BJf', [0; 0.4], ...
                'BJ', [0; 0.7], 'Bd', [0; 0.2], 'Bdf', [0.1; 0], ...
                'CJ', [1 0], 'DJ', 0.2);
            tc.verifyEqual(S.tag, 'impulsive');
            % thin() round-trips losslessly through fromstruct.
            md = S.thin();
            tc.verifyTrue(isequal(mvst.System.fromstruct(md).thin(), md), ...
                'impulsive thin() must round-trip field-for-field.');
            % rich() carries every extended channel through.
            R = S.rich();
            for f = {'Ac','Bc','Dwc','Bw','C','Dw','J','Jf','lamf','BJf', ...
                     'BJ','Bd','Bdf','CJ','DJ'}
                tc.verifyTrue(isfield(R, f{1}), ...
                    sprintf('rich() must carry the .%s channel.', f{1}));
            end
        end

        function impulsiveRejectsBadRate(tc)
            % lamf count must match the Jf channel count.
            Jf = cat(3, 1.1 * eye(2), 0.9 * eye(2));      % two jump channels
            tc.verifyError(@() mvst.System.impulsive([0 1;-1 0], [0;1], eye(2), ...
                'Jf', Jf, 'lamf', 0.5), 'mvst:System:dim');
        end

        function impulsiveRejectsBadDwcRows(tc)
            % Dwc must have n rows to match the state.
            tc.verifyError(@() mvst.System.impulsive([0 1;-1 0], [0;1], eye(2), ...
                'Ac', 0.1*eye(2), 'Bw', [0;1], 'Dwc', [0.1 0 0]'), ...
                'mvst:System:dim');
        end

        function impulsiveRejectsBadCJcols(tc)
            % CJ must have n columns to match the state.
            tc.verifyError(@() mvst.System.impulsive([0 1;-1 0], [0;1], eye(2), ...
                'CJ', [1 0 0]), 'mvst:System:dim');
        end

        % ---------------------------------------------------------------- %
        %  fromstruct round trip                                            %
        % ---------------------------------------------------------------- %
        function fromstructModeRoundtrip(tc)
            md = mvst_mode([-1 0.5; 0 -2], [1; 1], eye(2), 1, ...
                           'S', [0.1; 0.1], 'X0', 2 * eye(2));
            S  = mvst.System.fromstruct(md);
            tc.verifyTrue(isequal(S.thin(), md), ...
                'thin() must reproduce the source mode field-for-field.');
        end

        function fromstructSystemRoundtrip(tc)
            m1  = mvst_mode([-1 0; 0 -2], [1; 0], eye(2), 1);
            m2  = mvst_mode([-3 1; 0 -1], [0; 1], eye(2), 1);
            sys = mvst_system({m1, m2}, 'Pi', [-1 1; 2 -2], ...
                              'pi0', [0.5 0.5], 'ZT', eye(2));
            S   = mvst.System.fromstruct(sys);
            tc.verifyTrue(isequal(S.thin(), sys), ...
                'thin() must reproduce the source system field-for-field.');
        end

        function fromstructSampledRoundtrip(tc)
            % A thin mode carrying the sampled-data extension fields.
            md = mvst_mode([-1 1; 0 -1], [0; 1], eye(2), 1, 'J', eye(2));
            md.Tbar = 0.3;
            S  = mvst.System.fromstruct(md);
            tc.verifyTrue(isequal(S.thin(), md));
            tc.verifyEqual(S.tag, 'sampled');
        end

        % ---------------------------------------------------------------- %
        %  rich bridge equals mvst_torich(thin)                             %
        % ---------------------------------------------------------------- %
        function richEqualsTorichSingle(tc)
            % Single-mode RICH is now the LOSSLESS flat-mode bridge (finding #1):
            % it equals MVST_TORICH applied to the FLAT thin mode, NOT to the
            % system wrapper (which injects the trivial 1x1 generator .Pi/.Pig).
            S = mvst.System.ct([-1 0.5; 0 -2], [1; 1], 'Q', eye(2), 'R', 1);
            tc.verifyTrue(isequal(S.rich(), mvst_torich(S.thin().modes{1})));
        end

        function richEqualsTorichMulti(tc)
            S = mvst.System.mjls({[-1 0; 0 -2], [-3 1; 0 -1]}, ...
                                 {[1; 0], [0; 1]}, [-1 1; 2 -2]);
            tc.verifyTrue(isequal(S.rich(), mvst_torich(S.thin())));
        end

        % ---------------------------------------------------------------- %
        %  method sugar equals the engine on the equivalent struct          %
        % ---------------------------------------------------------------- %
        function l2gainSugarMatchesEngine(tc)
            % small Hurwitz 2x2 continuous-time system with disturbance / output
            % (the exact backend expects matching output / disturbance dims)
            Sstruct = struct('A', [-2 1; 0 -3], 'Bw', [1; 0.5], ...
                             'C', [1 0], 'Dw', 0);
            expected = mvst_l2gain(Sstruct);           % run the actual engine
            S   = mvst.System.fromstruct(Sstruct);
            got = S.l2gain();                          % method sugar
            tc.verifyEqual(got.gamma, expected.gamma, 'RelTol', 1e-6);
            tc.verifyEqual(got.method, expected.method);
        end

        % ---------------------------------------------------------------- %
        %  finding #2: 'Ac' may be supplied without 'Bc' (padded to zeros)  %
        % ---------------------------------------------------------------- %
        function ctAcOnlyConstructsAndRuns(tc)
            % The documented state-noise-only example must construct AND run.
            S = mvst.System.ct([0 1; -2 -1], zeros(2, 0), 'Ac', 0.15 * eye(2));
            tc.verifyEqual(S.n, 2);
            tc.verifyEqual(S.tag, 'ct');
            tc.verifyEqual(size(S.modes{1}.Ac, 3), 1);   % one state-noise channel
            tc.verifyEqual(size(S.modes{1}.Bc, 3), 1);   % Bc padded to match
            r = S.stability();
            tc.verifyEqual(r.class, 'ct');
            tc.verifyTrue(islogical(r.stable));
        end

        function modeBcOnlyPadsAc(tc)
            % symmetric case: 'Bc' supplied without 'Ac' pads Ac to zeros.
            md = mvst_mode([-1 0; 0 -1], [1; 0], eye(2), 1, 'Bc', [0.1; 0]);
            tc.verifyEqual(md.Nc, 1);
            tc.verifyEqual(size(md.Ac, 3), 1);
            tc.verifyEqual(md.Ac, zeros(2, 2, 1));
        end

        % ---------------------------------------------------------------- %
        %  finding #4: DT transition P -> tag 'dtmjls' (not 'mjls')         %
        % ---------------------------------------------------------------- %
        function fromstructDTMjlsTag(tc)
            m1  = mvst_mode([1.1 0.3; 0 1.05], [0; 1], eye(2), 1);
            m2  = mvst_mode([1.05 -0.4; 0.2 1.1], [0; 1], eye(2), 1);
            sys = mvst_system({m1, m2}, 'P', [0.7 0.3; 0.4 0.6]);
            S   = mvst.System.fromstruct(sys);
            tc.verifyEqual(S.tag, 'dtmjls');
            tc.verifyEqual(mvst_stability(S).class, 'dtmjls');
        end

        % ---------------------------------------------------------------- %
        %  finding #1: single-mode bridge is LOSSLESS (correct backend)     %
        % ---------------------------------------------------------------- %
        function poissonL2gainBridge(tc)
            % CRITICAL repro: the Poisson factory object must route to the
            % Poisson backend (gamma = 0.571289), not silently to 'ct' (=1).
            A = -eye(2); B = [1; 0];
            Sp = mvst.System.poisson(A, B, 0.5 * eye(2), 2, 'Bw', [1; 0], 'C', [1 0]);
            got = Sp.l2gain();
            tc.verifyEqual(got.method, 'poisson');
            tc.verifyEqual(got.gamma, 0.571289, 'AbsTol', 1e-5);
            % same number via the direct engine on the FLAT thin mode struct
            eng = mvst_l2gain(Sp.thin().modes{1});
            tc.verifyEqual(eng.method, 'poisson');
            tc.verifyEqual(got.gamma, eng.gamma, 'RelTol', 1e-9);
        end

        function levyL2gainBridge(tc)
            % Diffusion + Poisson jump (levy tag) must route to the poisson
            % backend with the diffusion channels carried, both paths equal.
            S = mvst.System.poisson([-1 0; 0 -1], [1; 0], 0.5 * eye(2), 1, ...
                                    'Ac', 0.15 * eye(2), 'Bw', [1; 0], 'C', [1 0]);
            tc.verifyEqual(size(S.modes{1}.Ac, 3), 1);   % diffusion channel carried
            got = S.l2gain();
            eng = mvst_l2gain(S.thin().modes{1});
            tc.verifyEqual(got.method, 'poisson');       % routes through the poisson backend
            tc.verifyEqual(got.gamma, eng.gamma, 'RelTol', 1e-9);
        end

        function sampledL2gainBridge(tc)
            % Sampled-data factory must route to the 'sd' backend, both equal.
            S = mvst.System.sampled([-2 1; 0 -3], [0; 1], eye(2), 0.3, ...
                                    'Bw', [1; 0.5], 'C', [1 0]);
            got = S.l2gain();
            eng = mvst_l2gain(S.thin().modes{1});
            tc.verifyEqual(got.method, 'sd');
            tc.verifyEqual(got.gamma, eng.gamma, 'RelTol', 1e-9);
        end

        function impulsiveL2gainBridge(tc)
            % Impulsive factory: l2gain has no impulsive class, both paths agree
            % (the reset is ignored consistently: the bridge drops nothing else).
            S = mvst.System.impulsive([-2 1; 0 -3], [0; 1], 0.5 * eye(2), ...
                                      'Bw', [1; 0.5], 'C', [1 0]);
            got = S.l2gain();
            eng = mvst_l2gain(S.thin().modes{1});
            tc.verifyEqual(got.gamma, eng.gamma, 'RelTol', 1e-9);
            tc.verifyEqual(got.method, eng.method);
        end

        function perfBridge(tc)
            % perf('l2') on a Poisson object equals the engine on the flat mode.
            S = mvst.System.poisson(-eye(2), [1; 0], 0.5 * eye(2), 2, ...
                                    'Bw', [1; 0], 'C', [1 0]);
            got = S.perf('l2');
            eng = mvst_perf(S.thin().modes{1}, 'l2');
            tc.verifyEqual(got.gamma, eng.gamma, 'RelTol', 1e-9);
        end

        function gamegammaBridge(tc)
            % gamegamma on a Poisson object equals the engine on the flat mode
            % (with the control-input / weight aliases the game backend reads).
            S = mvst.System.poisson(-eye(2), [0; 1], 0.5 * eye(2), 2, ...
                                    'Bw', [1; 0], 'Q', eye(2), 'R', 1);
            got = S.gamegamma();
            md  = S.thin().modes{1};
            md.Bu = md.B; md.Ru = md.R;               % names the game backend reads
            eng = mvst_game_gamma(md);
            tc.verifyEqual(got.method, 'poisson');
            tc.verifyEqual(got.gammastar, eng.gammastar, 'RelTol', 1e-8);
        end

        % ---------------------------------------------------------------- %
        %  finding #1: multi-mode per-mode I/O is refused, not dropped      %
        % ---------------------------------------------------------------- %
        function permodeIORejected(tc)
            m1 = mvst_mode([-1 0; 0 -2], [1; 0], eye(2), 1); m1.C = [1 0]; m1.Bw = [1; 0];
            m2 = mvst_mode([-3 1; 0 -1], [0; 1], eye(2), 1); m2.C = [0 1]; m2.Bw = [1; 0];
            sys = mvst_system({m1, m2}, 'Pi', [-1 1; 2 -2]);
            S   = mvst.System.fromstruct(sys);
            tc.verifyError(@() S.l2gain(), 'mvst:System:permodeIO');
        end

        % ---------------------------------------------------------------- %
        %  finding #5: multi-mode-with-reset is refused by the pos bridge   %
        % ---------------------------------------------------------------- %
        function posSwitchedResetUnsupported(tc)
            S = mvst.System.switched({[-1 0; 3 -1], [-2 1; 0 -1]}, ...
                                     {[1; 0], [0; 1]}, 'J', eye(2));
            tc.verifyError(@() S.pos('stability'), 'mvst:System:posUnsupported');
        end
    end
end
