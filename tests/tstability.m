classdef tstability < matlab.unittest.TestCase
%TSTABILITY  Unit tests for the generic stability / stabilization dispatchers.
%   Run from the matlab/tests directory:
%       results = runtests('tstability');
%
%   For every system class the test rebuilds the SAME small instance used by
%   the corresponding example script and asserts that MVST_STABILITY reproduces
%   the example's inline exact number (spectral abscissa or monodromy radius),
%   recomputed independently here from MVST_LIFTSYM to a 1e-8 tolerance.  The
%   stabilization tests (ct + mjls + sampled) certify that the closed loop
%   returned by MVST_STABILIZE is declared stable by its own MVST_STABILITY check.

    properties (Constant)
        AbsTolNum = 1e-8;
    end

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root);                          % +mvst package
            addpath(fullfile(root, 'lib'));         % engine layer
            addpath(fullfile(root, 'examples'));    % ex_data
        end
    end

    methods (Test)
        % ================================================================ %
        %  Perron-Frobenius irreducibility (opts.pf)                        %
        %  prop:stab:PFirr + finite Burnside test prop:stab:PFburnside      %
        % ================================================================ %
        function pfReducibleDiagonal(tc)
            % A single diagonal generator leaves span{e1} invariant: reducible,
            % and its unital algebra is only the diagonal matrices (dim n<n^2).
            r = mvst_stability(mvst.System.ct(diag([-1 -2]),zeros(2,0)), ...
                               struct('pf','irreducible'));
            tc.verifyFalse(r.irreducible);
            tc.verifyLessThan(r.pf.algdim, r.pf.full);
        end
        function pfIrreducibleFlow(tc)
            % Rotation + distinct diagonal have no common invariant subspace:
            % irreducible, algebra = all of C^{n x n} (dim n^2).
            r = mvst_stability(mvst.System.ct([0 1;-1 0],zeros(2,0),'Ac',diag([1 2])), ...
                               struct('pf','irreducible'));
            tc.verifyTrue(r.irreducible);
            tc.verifyEqual(r.pf.algdim, r.pf.full);
        end
        function pfIrreducibleJump(tc)
            % Swap + diagonal jump family: irreducible (prop:stab:PFirr(a)).
            r = mvst_stability(mvst.System.ct([0 1;1 0],zeros(2,0),'Ac',diag([1 2])), ...
                               struct('pf','irreducible'));
            tc.verifyTrue(r.irreducible);
        end
        function pfReducibleSharedSubspace3(tc)
            % Two 3x3 generators sharing invariant span{e1,e2}: reducible.
            M1=[1 1 0;0 1 0;0 0 2]; M2=[0 1 0;1 0 0;0 0 1];
            r = mvst_stability(mvst.System.ct(M1,zeros(3,0),'Ac',M2), ...
                               struct('pf','irreducible'));
            tc.verifyFalse(r.irreducible);
            tc.verifyEqual(r.irreducible, r.pf.algdim==r.pf.full);
        end
        function pfIrreducibleCyclic3(tc)
            % 3x3 cyclic shift + distinct diagonal: irreducible (algdim=9).
            r = mvst_stability(mvst.System.ct([0 1 0;0 0 1;1 0 0],zeros(3,0), ...
                               'Ac',diag([1 2 3])), struct('pf','irreducible'));
            tc.verifyTrue(r.irreducible);
            tc.verifyEqual(r.pf.algdim, 9);
        end

        % ================================================================ %
        %  Analysis: each class reproduces its example's exact number       %
        % ================================================================ %
        function ctAbscissa(tc)
            % Folded Levy operating point (ex_stab_levy): pure diffusion lift.
            A0 = [0 1; -2 -0.5]; A1 = [0 1; 1 0]; c = 0.30;
            expected = max(real(eig(mvst_liftsym(A0, {sqrt(c) * A1}))));
            S = mvst.System.ct(A0, [0; 0], 'Ac', sqrt(c) * A1, 'Bc', [0; 0]);
            r = mvst_stability(S);
            tc.verifyEqual(r.class, 'ct');
            tc.verifyEqual(r.margin, expected, 'AbsTol', tc.AbsTolNum);
            tc.verifyEqual(r.stable, expected < 0);
        end

        function dtLogRadius(tc)
            % Discrete-time stochastic system (ex_stabz_dt open loop).
            A0 = [1.1 0.3; 0 1.05]; A1 = 0.2 * eye(2);
            radOpen = max(abs(eig(mvst_liftsym(zeros(2), {A0, A1}))));
            sys = struct('A', A0, 'Ac', A1);
            r = mvst_stability(sys, struct('class', 'dt'));
            tc.verifyEqual(r.margin, log(radOpen), 'AbsTol', tc.AbsTolNum);
            tc.verifyEqual(r.stable, radOpen < 1);
        end

        function mjlsCoupledAbscissa(tc)
            % MJLS with diffusion (ex00_stability panel (a)).
            A1 = [0 1; -2 -1]; A2 = [0 1; -3 -1.5];
            Cd = 0.15 * eye(2); Pig = [-1 1; 2 -2];
            nr = 3; Lc = {mvst_liftsym(A1, {Cd}), mvst_liftsym(A2, {Cd})};
            Lmj = zeros(2 * nr);
            for i = 1:2
                for j = 1:2
                    blk = Pig(j, i) * eye(nr);
                    if i == j; blk = blk + Lc{i}; end
                    Lmj((i-1)*nr+1:i*nr, (j-1)*nr+1:j*nr) = blk;
                end
            end
            expected = max(real(eig(Lmj)));
            S = mvst.System.mjls({A1, A2}, {[0; 0], [0; 0]}, Pig, 'Ac', Cd, 'Bc', [0; 0]);
            r = mvst_stability(S);
            tc.verifyEqual(r.class, 'mjls');
            tc.verifyEqual(r.margin, expected, 'AbsTol', tc.AbsTolNum);
        end

        function dtmjlsCoupledLogRadius(tc)
            % DT-MJLS (ex_stabz_dtmjls open loop): coupled lifted transition.
            n = 2; M = 2;
            A  = cat(3, [1.1 0.3; 0 1.05], [1.05 -0.4; 0.2 1.1]);
            Bc = cat(3, 0.15 * eye(2), 0.15 * eye(2));
            P  = [0.7 0.3; 0.4 0.6];
            Us = mvst_svecbasis(n); d = size(Us, 2); Tj = cell(1, M); G = zeros(M * d);
            for j = 1:M
                Tj{j} = Us' * (kron(A(:,:,j), A(:,:,j)) + kron(Bc(:,:,j), Bc(:,:,j))) * Us;
            end
            for i = 1:M
                for j = 1:M
                    G((i-1)*d+(1:d), (j-1)*d+(1:d)) = P(j, i) * Tj{j};
                end
            end
            expected = max(abs(eig(G)));
            m1 = mvst_mode(A(:,:,1), [0; 0], zeros(n), 1, 'Ac', Bc(:,:,1), 'Bc', [0; 0]);
            m2 = mvst_mode(A(:,:,2), [0; 0], zeros(n), 1, 'Ac', Bc(:,:,2), 'Bc', [0; 0]);
            sys = mvst_system({m1, m2}, 'P', P);
            r = mvst_stability(sys);
            tc.verifyEqual(r.class, 'dtmjls');
            tc.verifyEqual(r.margin, log(expected), 'AbsTol', tc.AbsTolNum);
        end

        function poissonAbscissa(tc)
            % Pure Poissonian-jump generator L0 + lam (Jlift - I).
            A = [0 1; -2 -1]; J1 = 1.15 * eye(2); lam = 0.5; nr = 3;
            expected = max(real(eig(mvst_liftsym(A, {}) + ...
                lam * (mvst_liftsym(zeros(2), {J1}) - eye(nr)))));
            S = mvst.System.poisson(A, zeros(2, 0), J1, lam);
            r = mvst_stability(S);
            tc.verifyEqual(r.class, 'poisson');
            tc.verifyEqual(r.margin, expected, 'AbsTol', tc.AbsTolNum);
        end

        function levyAbscissa(tc)
            % GBM + Poissonian jumps (ex00_stability panel (b)): diffusion + jump.
            A = [0 1; -2 -1]; C1 = 0.15 * eye(2); J1 = 1.15 * eye(2);
            lam = 0.5; nr = 3;
            L0 = mvst_liftsym(A, {C1});
            LJ = mvst_liftsym(zeros(2), {J1}) - eye(nr);
            expected = max(real(eig(L0 + lam * LJ)));
            S = mvst.System.poisson(A, [0; 0], J1, lam, 'Ac', C1, 'Bc', [0; 0]);  % diffusion + jump
            r = mvst_stability(S);
            tc.verifyEqual(r.class, 'levy');
            tc.verifyEqual(r.margin, expected, 'AbsTol', tc.AbsTolNum);
        end

        function sampledMonodromy(tc)
            % Impulsive monodromy MJ*expm(Tbar*Mfo) (ex_stabz_imp open loop).
            A0 = [0.8 1; 0 0.5]; A1 = 0.2 * eye(2); J = 0.7 * eye(2); Tbar = 0.6;
            Mfo = mvst_liftsym(A0, {A1}); MJ = mvst_liftsym(zeros(2), {J});
            expected = max(abs(eig(MJ * expm(Tbar * Mfo))));
            S = mvst.System.sampled(A0, [0; 0], J, Tbar, 'Ac', A1, 'Bc', [0; 0]);
            r = mvst_stability(S);
            tc.verifyEqual(r.class, 'sampled');
            tc.verifyEqual(r.margin, log(expected), 'AbsTol', tc.AbsTolNum);
            tc.verifyEqual(r.T, Tbar);
        end

        function switchedMonodromyScalarAndGrid(tc)
            % Switched monodromy over one cycle (ex00_stability panel (d)).
            Asw1 = [-1 4; 0 -1]; Asw2 = [-1 0; 4 -1];
            L1 = mvst_liftsym(Asw1, {}); L2 = mvst_liftsym(Asw2, {});
            Td = 1.0;
            expected = max(abs(eig(expm(Td * L2) * expm(Td * L1))));
            S = mvst.System.switched({Asw1, Asw2}, {zeros(2, 0), zeros(2, 0)}, 'Td', Td);
            r = mvst_stability(S, struct('T', Td));
            tc.verifyEqual(r.class, 'switched');
            tc.verifyEqual(r.margin, log(expected), 'AbsTol', tc.AbsTolNum);

            % dwell grid -> stability-vs-dwell curve
            grid = [0.5 1.0 2.0];
            rg = mvst_stability(S, struct('T', grid));
            tc.verifyTrue(isfield(rg, 'curve'));
            tc.verifyEqual(rg.curve.T, grid);
            for k = 1:numel(grid)
                rad = max(abs(eig(expm(grid(k) * L2) * expm(grid(k) * L1))));
                tc.verifyEqual(rg.curve.radius(k), rad, 'AbsTol', tc.AbsTolNum);
            end
        end

        % ================================================================ %
        %  Robustness of the class resolution       %
        % ================================================================ %
        function bareStructAnalyzesCt(tc)
            % finding 3a: a bare struct('A',...) WITHOUT an .Ac field must
            % analyze as 'ct' (mode_channels must not fabricate an empty channel).
            A = [0 1; -2 -1];
            r = mvst_stability(struct('A', A));
            tc.verifyEqual(r.class, 'ct');
            tc.verifyEqual(r.margin, max(real(eig(mvst_liftsym(A, {})))), ...
                           'AbsTol', tc.AbsTolNum);
        end

        function switchedResetDropped(tc)
            % finding 3b: a single-mode system carrying a common reset would
            % have the reset silently dropped and be analyzed as 'ct'; it must
            % ERROR instead, asking for the dwell-time / an explicit class.
            S = mvst.System.switched({[-1 0; 0 -2]}, {[1; 0]}, 'J', eye(2));
            tc.verifyError(@() mvst_stability(S), 'mvst:stability:resetDropped');
        end

        % ================================================================ %
        %  Stabilization: closed loop certifies through its own check        %
        % ================================================================ %
        function ctStabilizeCertifies(tc)
            % Open-loop unstable CT diffusion (ex_data D.ct family).
            S = mvst.System.ct([0 1; -1 0.2], [0; 1], 'Q', eye(2), 'R', 1, ...
                               'Ac', 0.3 * [0 0; 1 0], 'Bc', [0; 0.2]);
            r = S.stabilize();
            tc.verifyEqual(r.class, 'ct');
            tc.verifyTrue(r.check.stable, 'CT closed loop must certify stable.');
            tc.verifyLessThan(r.check.margin, 0);
        end

        function mjlsStabilizeCertifies(tc)
            D = ex_data();
            r = mvst_stabilize(D.mjls);
            tc.verifyEqual(r.class, 'mjls');
            tc.verifyTrue(iscell(r.K) && numel(r.K) == 2);
            tc.verifyTrue(r.check.stable, 'MJLS closed loop must certify stable.');
        end

        function sampledStabilizeCertifies(tc)
            D = ex_data();
            r = mvst_stabilize(D.sd);
            tc.verifyEqual(r.class, 'sampled');
            tc.verifyTrue(r.check.stable, 'Sampled-data closed loop must certify stable.');
            tc.verifyLessThan(r.check.margin, 0);
        end
    end
end
