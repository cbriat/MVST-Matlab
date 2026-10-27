classdef tsdrange < matlab.unittest.TestCase
%TSDRANGE  Tests for the range (aperiodic) sampled-data refresh-gain synthesis
%   backend MVST_STABZ_SD_RANGE (Stabilization chapter, Proposition
%   prop:stabz:sdrange scheduled design, Corollary cor:stabz:sdrange:const
%   constant design, and the dilated/Finsler SLACK constant-gain design with a
%   period-dependent certificate), on the sampled-data master (double integrator
%   with the actuated diffusion channel, reused from EX_STABZ_SD).  It verifies:
%     (1) FEASIBILITY SANITY: at h=0 (the degenerate window = periodic case) the
%         scheduled and constant programs coincide and are BOTH feasible; at a
%         representative nonzero window the scheduled program stays feasible.
%     (2) SCHEDULED >= CONSTANT: the scheduled SOS program tolerates a strictly
%         wider window than the constant one -- feasible at a representative
%         window where the constant (frozen-certificate) program is infeasible,
%         and its critical half-width is >= the constant critical half-width.
%     (3) MONODROMY CERTIFICATION: the recovered scheduled refresh gain drives
%         the exact one-period monodromy spectral radius below 1 across the
%         representative window (independent lifted-space check).
%     (4) GOLDEN COMPARE: the recovered representative scheduled gain K_J(Tc)
%         reproduces the frozen golden value in data/ex_stabz_sdrange.tex.
%     (5-8) SLACK MODE: at h=0 the dilated design is feasible and at least as
%         good as constant (lossless Finsler at a single period); at a nonzero
%         window it stays feasible where constant is not (wider tolerance); its
%         recovered CONSTANT gain drives the exact monodromy radius below 1; and
%         it reproduces the golden slack gain StabzSdrKslack.
%
%   See also MVST_STABZ_SD_RANGE, EX_STABZ_SD, EX_STABZ_SDRANGE, RUN_ALL_TESTS.

    properties (Constant)
        A  = [0 1; 2 0.1];
        B  = [0; 1];
        Ac = 0.2*eye(2);
        B1 = 0.1*[0; 1];
        Tc = 0.25;
        deg = 4;
        tol = 1e-7;
    end

    methods (TestClassSetup)
        function addPaths(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);                 % matlab/
            addpath(root); addpath(fullfile(root, 'lib'));
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS (sossolve) not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods
        function r = solve(tc, h, mode)
            plant = struct('A',tc.A,'Bu',tc.B,'Ac',tc.Ac,'Bc',tc.B1);
            r = mvst_stabz_sd_range(plant, struct('Tmin',tc.Tc-h,'Tmax',tc.Tc+h, ...
                                                  'mode',mode,'deg',tc.deg));
        end
    end

    methods (Test)
        % (1) at h=0 scheduled and constant coincide and are both feasible.
        function periodicCoincide(tc)
            rs = tc.solve(0,'scheduled');
            rc = tc.solve(0,'constant');
            tc.verifyGreaterThan(rs.slack, tc.tol, 'scheduled must be feasible at h=0.');
            tc.verifyGreaterThan(rc.slack, tc.tol, 'constant must be feasible at h=0.');
            tc.verifyEqual(rs.slack, rc.slack, 'AbsTol', 1e-4, ...
                'at h=0 the two programs are the same LMI, so their slacks coincide.');
        end

        % (2) scheduled tolerates a wider window than constant.
        function scheduledWiderThanConstant(tc)
            hrep = 0.2*tc.Tc;                       % representative nonzero window
            rs = tc.solve(hrep,'scheduled');
            rc = tc.solve(hrep,'constant');
            tc.verifyGreaterThan(rs.slack, tc.tol, ...
                'scheduled must stay feasible at the representative window.');
            tc.verifyLessThanOrEqual(rc.slack, tc.tol, ...
                'constant (frozen certificate) must be infeasible at a nonzero window.');
        end

        % (3) the recovered scheduled gain certifies by the exact monodromy.
        function scheduledCertifiesMonodromy(tc)
            hrep = 0.2*tc.Tc; Tmin = tc.Tc-hrep; Tmax = tc.Tc+hrep;
            rs = tc.solve(hrep,'scheduled');
            rad = tsdrange.windowRadius(tc.A,tc.B,tc.Ac,tc.B1, rs.KJfun, Tmin, Tmax, false);
            tc.verifyLessThan(rad, 1, ...
                'recovered scheduled gain must drive the monodromy radius below 1.');
        end

        % (4) golden compare of the representative scheduled gain K_J(Tc).
        function goldenGain(tc)
            here = fileparts(mfilename('fullpath'));
            gold = fullfile(here,'golden','data','ex_stabz_sdrange.tex');
            tc.assumeTrue(isfile(gold), 'golden file present');
            KsTc_gold = tsdrange.readMacroRow(gold, 'StabzSdrKsTc');
            rs = tc.solve(0.2*tc.Tc,'scheduled');
            KsTc = rs.KJfun(tc.Tc);
            tc.verifyEqual(KsTc(:).', KsTc_gold(:).', 'AbsTol', 5e-3, ...
                'recovered scheduled K_J(Tc) must reproduce the golden value.');
        end

        % (5) SLACK at h=0: feasible and at least as good as constant (Finsler
        %     dilation is lossless at a single period, so it cannot be worse).
        function slackAtLeastAsGoodAtZero(tc)
            rx = tc.solve(0,'slack');
            rc = tc.solve(0,'constant');
            tc.verifyGreaterThan(rx.slack, tc.tol, 'slack must be feasible at h=0.');
            tc.verifyGreaterThanOrEqual(rx.slack, rc.slack - 1e-6, ...
                'slack (dilated) must be at least as good as constant at h=0.');
        end

        % (6) SLACK tolerates a wider window than constant: feasible at a
        %     representative nonzero window where the constant program is not.
        function slackWiderThanConstant(tc)
            hrep = 0.2*tc.Tc;
            rx = tc.solve(hrep,'slack');
            rc = tc.solve(hrep,'constant');
            tc.verifyGreaterThan(rx.slack, tc.tol, ...
                'slack must stay feasible at the representative window.');
            tc.verifyLessThanOrEqual(rc.slack, tc.tol, ...
                'constant (frozen certificate) must be infeasible at a nonzero window.');
        end

        % (7) the recovered CONSTANT slack gain certifies by the exact monodromy
        %     across the representative window (independent lifted-space check).
        function slackCertifiesMonodromy(tc)
            hrep = 0.2*tc.Tc; Tmin = tc.Tc-hrep; Tmax = tc.Tc+hrep;
            rx = tc.solve(hrep,'slack');
            rad = tsdrange.windowRadius(tc.A,tc.B,tc.Ac,tc.B1, rx.KJfun, Tmin, Tmax, true);
            tc.verifyLessThan(rad, 1, ...
                'recovered slack constant gain must drive the monodromy radius below 1.');
        end

        % (8) golden compare of the representative slack constant gain K_J.
        function slackGoldenGain(tc)
            here = fileparts(mfilename('fullpath'));
            gold = fullfile(here,'golden','data','ex_stabz_sdrange.tex');
            tc.assumeTrue(isfile(gold), 'golden file present');
            K_gold = tsdrange.readMacroRow(gold, 'StabzSdrKslack');
            rx = tc.solve(0.2*tc.Tc,'slack');
            K = rx.KJfun(tc.Tc);
            tc.verifyEqual(K(:).', K_gold(:).', 'AbsTol', 5e-3, ...
                'recovered slack K_J must reproduce the golden value.');
        end
    end

    methods (Static)
        function rad = windowRadius(A,B,Ac,B1, KJfun, Tmin, Tmax, isConst)
            n = size(A,1); m = size(B,2); N = n+m;
            Atil  = [A B; zeros(m,N)];  Actil = [Ac B1; zeros(m,N)];
            J0 = [eye(n) zeros(n,m); zeros(m,N)];  BJ = [zeros(n,m); eye(m)];
            Mf = mvst_liftsym(Atil, {Actil});  Ud = mvst_svecbasis(N);
            Ts = linspace(Tmin, Tmax, 25); rad = 0;
            for T = Ts
                if isConst; K = KJfun(Tmin); else; K = KJfun(T); end
                G = J0 + BJ*K;  MJ = Ud' * kron(G,G) * Ud;
                rad = max(rad, max(abs(eig(MJ*expm(T*Mf)))));
            end
        end

        function v = readMacroRow(file, name)
            % Parse "\def\NAME{\begin{bmatrix} a & b & c \end{bmatrix}}" or a
            % scalar "\def\NAME{x}" into a numeric row vector.
            txt = fileread(file);
            % NOTE: MATLAB regexp '.' matches newlines, so restrict the capture
            % to a single line with [^\r\n] to avoid spanning adjacent macros.
            tok = regexp(txt, ['\\def\\' name '\{([^\r\n]*)\}'], 'tokens', 'once');
            assert(~isempty(tok), 'macro %s not found in %s', name, file);
            body = tok{1};
            body = regexprep(body, '\\begin\{bmatrix\}|\\end\{bmatrix\}', '');
            parts = regexp(body, '&', 'split');
            v = cellfun(@(s) str2double(strtrim(s)), parts);
        end
    end
end
