classdef System
%SYSTEM  Unified value-type wrapper for a matrix-valued dynamical system.
%   S = MVST.SYSTEM.CT(A,B,...) and the sibling static factories build an
%   immutable object that carries ONE model in the thin OC/Riccati convention
%   (MVST_MODE / MVST_SYSTEM) and bridges, on demand, to the rich SOS-layer
%   convention (MVST_TORICH).  Every analysis and synthesis routine of the
%   companion library is exposed as a method that delegates to the existing
%   engine through THIN or RICH, so the legacy solvers keep working unchanged.
%
%   The object is built through a STATIC FACTORY (the constructor is private),
%   one per system class of the book:
%       MVST.SYSTEM.CT         continuous-time diffusion
%       MVST.SYSTEM.MJLS       Markov-jump linear system (CT generator)
%       MVST.SYSTEM.SAMPLED    sampled-data (zero-order hold) system
%       MVST.SYSTEM.POISSON    Poissonian-jump diffusion
%       MVST.SYSTEM.SWITCHED   dwell-time switched system
%       MVST.SYSTEM.IMPULSIVE  impulsive system with a scheduled reset
%       MVST.SYSTEM.FROMSTRUCT wrap any legacy thin/rich struct
%
%   PROPERTIES (immutable)
%     modes  1 x M cell of mode structs in the thin convention (MVST_MODE).
%     Pi     M x M continuous-time Markov generator coupling the modes.
%     pi0    1 x M initial mode distribution.
%     ZT     terminal weight for the finite-horizon solver.
%     P      M x M discrete-time transition (optional coupling extension).
%     J      n x n common switch / impulse reset (optional extension).
%   DERIVED (dependent): n (state dimension), m (input dimension),
%     M (number of modes).  The read-only TAG names the resolved class.
%
%   BRIDGE METHODS
%     thin()  the legacy MVST_SYSTEM / MVST_MODE struct the Riccati layer
%             consumes (byte-compatible: MVST_TORICH(obj.thin()) round-trips).
%     rich()  the SOS-layer struct, obtained by MVST_TORICH on the thin form.
%
%   ANALYSIS / SYNTHESIS SUGAR (one-line delegations)
%     l2gain, hinfnorm, passivity, kyp, dissipativity, steer,
%     gamegamma, ocfinite, ocinfinite, ocaverage, stability, stabilize, dissyn,
%     pos, perf.  Every OPTS argument accepts a plain struct or an MVST.OPTIONS.
%
%   Example (build a continuous-time system and price a finite horizon):
%       S = mvst.System.ct([-1 0.5; 0 -2], [1;1], 'Q', eye(2), 'R', 1);
%       r = S.ocfinite(3);
%       r.cost
%
%   See also MVST.OPTIONS, MVST.SOLVERCHECK, MVST_MODE, MVST_SYSTEM,
%   MVST_TORICH, MVST_L2GAIN, MVST_SOLVE_FINITE.

    properties (SetAccess = immutable)
        modes            % 1 x M cell of thin mode structs
        Pi               % M x M CT generator
        pi0              % 1 x M initial mode distribution
        ZT               % terminal weight
        P                % M x M DT transition (optional)
        J                % common reset (optional)
    end

    properties (SetAccess = immutable)
        tag              % resolved class tag (char): 'ct','mjls',...,'custom'
    end

    properties (Access = private)
        IsSystem (1,1) logical = true   % true -> thin() is a MVST_SYSTEM struct
    end

    properties (Dependent)
        n                % state dimension
        m                % input dimension
        M                % number of modes
    end

    % ==================================================================== %
    %  Construction                                                         %
    % ==================================================================== %
    methods (Access = private)
        function obj = System(modes, Pi, pi0, ZT, P, J, tag, isSystem)
            %SYSTEM  Private constructor; use the static factories instead.
            obj.modes    = modes;
            obj.Pi       = Pi;
            obj.pi0      = pi0;
            obj.ZT       = ZT;
            obj.P        = P;
            obj.J        = J;
            obj.tag      = tag;
            obj.IsSystem = isSystem;
        end
    end

    methods (Static)
        function obj = ct(A, B, opts)
            %CT  Continuous-time controlled diffusion system.
            %   S = MVST.SYSTEM.CT(A,B) builds the single-mode system
            %       dx = (A x + B u) dt + sum_l (Ac_l x + Bc_l u) dW_l,
            %   with quadratic cost blocks Q, R, S and the optional
            %   disturbance / performance I/O blocks used by the L2-gain and
            %   game engines.  Name/value options:
            %       'Q','R','S','W','X0'  cost / covariance blocks (MVST_MODE)
            %       'Ac','Bc'             state / input noise channels
            %       'Bw','C','Dw'         disturbance in, output, feedthrough
            %       'ZT'                  terminal weight (finite horizon)
            %
            %   Example:
            %       S = mvst.System.ct([-1 0.5; 0 -2], [1;1], 'Q', eye(2));
            %       S.n, S.tag
            arguments
                A  double
                B  double
                opts.Q  double = zeros(size(A,1))
                opts.R  double = eye(size(B,2))
                opts.S  double = []
                opts.W  double = []
                opts.X0 double = []
                opts.Ac double = []
                opts.Bc double = []
                opts.Bw double = []
                opts.C  double = []
                opts.Dw double = []
                opts.ZT double = []
            end
            n = mvst.System.checkSquare(A, 'A');
            mvst.System.checkRows(B, n, 'B');
            md  = mvst.System.buildMode(A, B, opts);
            ZT  = mvst.System.defaultZT(opts.ZT, n);
            obj = mvst.System({md}, zeros(1), 1, ZT, [], [], 'ct', true);
        end

        function obj = mjls(Alist, Blist, Pi, opts)
            %MJLS  Markov-jump linear system with a continuous-time generator.
            %   S = MVST.SYSTEM.MJLS(ALIST,BLIST,PI) builds an M-mode system
            %   whose modes switch under the M x M CT Markov generator PI
            %   (Metzler, zero row sums).  ALIST and BLIST are 1 x M cells of
            %   the per-mode drift / input matrices.  Name/value options
            %   ('Q','R','S','W','X0','Ac','Bc') are SHARED across modes; for
            %   per-mode data build the modes yourself and wrap them with
            %   MVST.SYSTEM.FROMSTRUCT(MVST_SYSTEM(...)).  'pi0' sets the
            %   initial distribution, 'ZT' the terminal weight.
            %
            %   Example:
            %       Pi = [-1 1; 2 -2];
            %       S  = mvst.System.mjls({[-1 0;0 -2],[-3 1;0 -1]}, ...
            %                             {[1;0],[0;1]}, Pi, 'Q', eye(2));
            arguments
                Alist cell
                Blist cell
                Pi    double
                opts.Q   double = []
                opts.R   double = []
                opts.S   double = []
                opts.W   double = []
                opts.X0  double = []
                opts.Ac  double = []
                opts.Bc  double = []
                opts.pi0 double = []
                opts.ZT  double = []
            end
            M = numel(Alist);
            if numel(Blist) ~= M
                error('mvst:System:dim', ...
                    'ALIST and BLIST must have the same length (%d vs %d).', ...
                    M, numel(Blist));
            end
            mvst.System.checkGenerator(Pi, M);
            n = mvst.System.checkSquare(Alist{1}, 'A');
            modes = cell(1, M);
            for i = 1:M
                mvst.System.checkSquare(Alist{i}, sprintf('A{%d}', i));
                mvst.System.checkRows(Blist{i}, n, sprintf('B{%d}', i));
                modes{i} = mvst.System.buildMode(Alist{i}, Blist{i}, opts);
            end
            pi0 = mvst.System.defaultPi0(opts.pi0, M);
            ZT  = mvst.System.defaultZT(opts.ZT, n);
            obj = mvst.System(modes, Pi, pi0, ZT, [], [], 'mjls', true);
        end

        function obj = sampled(A, B, J, T, opts)
            %SAMPLED  Sampled-data (zero-order-hold) system.
            %   S = MVST.SYSTEM.SAMPLED(A,B,J,T) builds a single-mode CT flow
            %   with sampling period T (>0) and sample-and-hold reset map J.
            %   The period is carried as the .Tbar field the L2-gain / steering
            %   engines read.  Same name/value options as MVST.SYSTEM.CT.
            %
            %   Example:
            %       S = mvst.System.sampled([-1 1; 0 -1], [0;1], eye(2), 0.2);
            arguments
                A  double
                B  double
                J  double
                T  (1,1) double {mustBePositive}
                opts.Q  double = zeros(size(A,1))
                opts.R  double = eye(size(B,2))
                opts.S  double = []
                opts.W  double = []
                opts.X0 double = []
                opts.Ac double = []
                opts.Bc double = []
                opts.Bw double = []
                opts.C  double = []
                opts.Dw double = []
                opts.ZT double = []
            end
            n = mvst.System.checkSquare(A, 'A');
            mvst.System.checkRows(B, n, 'B');
            mvst.System.checkSquare(J, 'J');
            md = mvst.System.buildMode(A, B, opts);
            md.J    = J;
            md.Tbar = T;
            ZT  = mvst.System.defaultZT(opts.ZT, n);
            obj = mvst.System({md}, zeros(1), 1, ZT, [], [], 'sampled', true);
        end

        function obj = poisson(A, B, J, lambda, opts)
            %POISSON  Poissonian-jump diffusion system.
            %   S = MVST.SYSTEM.POISSON(A,B,J,LAMBDA) builds a single-mode CT
            %   diffusion with spontaneous (Poisson) jump maps J at rates
            %   LAMBDA.  J is n x n x pf (or n x n for one channel) and LAMBDA
            %   is 1 x pf.  These populate the .Jf / .lamf fields read by the
            %   Poisson backends.  Same name/value options as MVST.SYSTEM.CT.
            %
            %   Example:
            %       S = mvst.System.poisson([-1 0;0 -1], [1;0], ...
            %                               0.5*eye(2), 2);
            arguments
                A      double
                B      double
                J      double
                lambda double
                opts.Q  double = zeros(size(A,1))
                opts.R  double = eye(size(B,2))
                opts.S  double = []
                opts.W  double = []
                opts.X0 double = []
                opts.Ac double = []
                opts.Bc double = []
                opts.Bw double = []
                opts.C  double = []
                opts.Dw double = []
                opts.ZT double = []
            end
            n  = mvst.System.checkSquare(A, 'A');
            mvst.System.checkRows(B, n, 'B');
            pf = size(J, 3);
            if size(J,1) ~= n || size(J,2) ~= n
                error('mvst:System:dim', 'Jump map J must be %d x %d [x pf].', n, n);
            end
            if numel(lambda) ~= pf
                error('mvst:System:dim', ...
                    'LAMBDA must have one rate per jump channel (%d).', pf);
            end
            md = mvst.System.buildMode(A, B, opts);
            md.Jf   = J;
            md.lamf = reshape(lambda, 1, pf);
            ZT  = mvst.System.defaultZT(opts.ZT, n);
            obj = mvst.System({md}, zeros(1), 1, ZT, [], [], 'poisson', true);
        end

        function obj = switched(Alist, Blist, opts)
            %SWITCHED  Dwell-time switched system.
            %   S = MVST.SYSTEM.SWITCHED(ALIST,BLIST) builds an M-mode system
            %   with NO stochastic coupling (the generator is zero); the modes
            %   are governed by a dwell-time schedule.  Name/value options:
            %   the shared cost/noise blocks as in MVST.SYSTEM.MJLS, plus
            %       'Td' dwell-time carried for the dwell backends,
            %       'J'  common reset applied at a switch,
            %       'pi0','ZT'.
            %
            %   Example:
            %       S = mvst.System.switched({[-1 0;0 -2],[-2 1;0 -1]}, ...
            %                                {[1;0],[0;1]}, 'Td', 0.5);
            arguments
                Alist cell
                Blist cell
                opts.Q   double = []
                opts.R   double = []
                opts.S   double = []
                opts.W   double = []
                opts.X0  double = []
                opts.Ac  double = []
                opts.Bc  double = []
                opts.Td  double = []
                opts.J   double = []
                opts.pi0 double = []
                opts.ZT  double = []
            end
            M = numel(Alist);
            if numel(Blist) ~= M
                error('mvst:System:dim', ...
                    'ALIST and BLIST must have the same length (%d vs %d).', ...
                    M, numel(Blist));
            end
            n = mvst.System.checkSquare(Alist{1}, 'A');
            modes = cell(1, M);
            for i = 1:M
                mvst.System.checkSquare(Alist{i}, sprintf('A{%d}', i));
                mvst.System.checkRows(Blist{i}, n, sprintf('B{%d}', i));
                modes{i} = mvst.System.buildMode(Alist{i}, Blist{i}, opts);
            end
            if ~isempty(opts.Td); modes{1}.Td = opts.Td; end   % carried for dwell backends
            pi0 = mvst.System.defaultPi0(opts.pi0, M);
            ZT  = mvst.System.defaultZT(opts.ZT, n);
            Jc  = [];
            if ~isempty(opts.J); mvst.System.checkSquare(opts.J, 'J'); Jc = opts.J; end
            obj = mvst.System(modes, zeros(M), pi0, ZT, [], Jc, 'switched', true);
        end

        function obj = impulsive(A, B, J, opts)
            %IMPULSIVE  Impulsive system: the book's most general single-mode class.
            %   S = MVST.SYSTEM.IMPULSIVE(A,B,J) builds a single-mode CT flow
            %   with a scheduled (dwell-time) jump / reset map J attached as the
            %   .J field the SOS layer reads.  Beyond MVST.SYSTEM.CT this factory
            %   exposes the FULL channel set of the general impulsive class of
            %   Section 2.7 (flow drift + Brownian diffusion + spontaneous
            %   Poissonian jumps + scheduled resets), each channel carrying both
            %   the control u and the exogenous disturbance w:
            %       'Q','R','S','W','X0'  cost / covariance blocks (MVST_MODE)
            %       'Ac','Bc'             flow state / input multiplicative noise
            %       'Dwc'                 flow disturbance-multiplicative noise
            %       'Bw','C','Dw'         flow disturbance in, output, feedthrough
            %       'Jf','lamf','BJf'     spontaneous (Poisson) jump maps, rates,
            %                             and jump-control selectors
            %       'BJ'                  scheduled-jump control selector
            %       'Bd','Bdf'            scheduled- / spontaneous-jump disturbance
            %       'CJ','DJ'             jump performance output and feedthrough
            %       'ZT'                  terminal weight (finite horizon)
            %   Every channel is dimension-validated (mvst:System:dim); the
            %   rich() / fromstruct() round trip carries them all losslessly.
            %
            %   Example:
            %       S = mvst.System.impulsive([0 1; -1 0], [0;1], 0.5*eye(2), ...
            %              'Ac', 0.1*eye(2), 'Bw', [0;1], 'C', [1 0], ...
            %              'Jf', 1.1*eye(2), 'lamf', 0.5, 'Bd', [0;0.2]);
            arguments
                A    double
                B    double
                J    double
                opts.Q    double = zeros(size(A,1))
                opts.R    double = eye(size(B,2))
                opts.S    double = []
                opts.W    double = []
                opts.X0   double = []
                opts.Ac   double = []
                opts.Bc   double = []
                opts.Dwc  double = []
                opts.BJ   double = []
                opts.Bw   double = []
                opts.C    double = []
                opts.Dw   double = []
                opts.Jf   double = []
                opts.lamf double = []
                opts.BJf  double = []
                opts.Bd   double = []
                opts.Bdf  double = []
                opts.CJ   double = []
                opts.DJ   double = []
                opts.ZT   double = []
            end
            n = mvst.System.checkSquare(A, 'A');
            m = size(B, 2);
            mvst.System.checkRows(B, n, 'B');
            if size(J,1) ~= n || size(J,2) ~= n
                error('mvst:System:dim', 'Reset map J must be %d x %d [x p].', n, n);
            end
            md = mvst.System.buildMode(A, B, opts);
            md.J = J;
            % ---- validate and attach the full impulsive channel set ----
            mvst.System.checkImpulsiveChannels(opts, n, m, size(J,3));
            extra = {'Dwc','BJ','Jf','lamf','BJf','Bd','Bdf','CJ','DJ'};
            for f = extra
                if isfield(opts, f{1}) && ~isempty(opts.(f{1}))
                    md.(f{1}) = opts.(f{1});
                end
            end
            ZT  = mvst.System.defaultZT(opts.ZT, n);
            obj = mvst.System({md}, zeros(1), 1, ZT, [], [], 'impulsive', true);
        end

        function obj = fromstruct(sys)
            %FROMSTRUCT  Wrap a legacy thin or rich struct as a SYSTEM.
            %   S = MVST.SYSTEM.FROMSTRUCT(SYS) accepts any of
            %     * a thin mode  (MVST_MODE output, has .A, no .modes),
            %     * a thin system(MVST_SYSTEM output, has .modes),
            %     * a rich struct(SOS-layer, has .A0),
            %   and returns the equivalent SYSTEM.  For a thin input the round
            %   trip is exact: MVST.SYSTEM.FROMSTRUCT(S).THIN() reproduces S.
            %
            %   Example:
            %       md = mvst_mode([-1 0;0 -2], [1;1], eye(2), 1);
            %       S  = mvst.System.fromstruct(md);
            %       isequal(S.thin(), md)   % true
            arguments
                sys struct
            end
            if isfield(sys, 'modes')                         % thin system
                modes = sys.modes;
                Msz   = numel(modes);
                Pi    = mvst_getf(sys, 'Pi',  zeros(Msz));
                pi0   = mvst_getf(sys, 'pi0', mvst.System.defaultPi0([], Msz));
                nn    = mvst.System.dimOf(modes{1});
                ZT    = mvst_getf(sys, 'ZT',  zeros(nn));
                P     = mvst_getf(sys, 'P',   []);
                J     = mvst_getf(sys, 'J',   []);
                tag   = mvst.System.inferClass(modes, Pi, P, J);
                obj   = mvst.System(modes, Pi, pi0, ZT, P, J, tag, true);

            elseif isfield(sys, 'A0')                        % rich struct
                [modes, Pi, pi0, ZT, P, J, isSys] = mvst.System.richToModel(sys);
                tag = mvst.System.inferClass(modes, Pi, P, J);
                obj = mvst.System(modes, Pi, pi0, ZT, P, J, tag, isSys);

            elseif isfield(sys, 'A')                         % thin bare mode
                modes = {sys};
                tag   = mvst.System.inferClass(modes, [], [], []);
                obj   = mvst.System(modes, zeros(1), 1, [], [], [], tag, false);

            else
                error('mvst:System:fromstruct', ...
                    'SYS must carry .A (thin mode), .modes (thin system) or .A0 (rich).');
            end
        end
    end

    % ==================================================================== %
    %  Dependent getters                                                    %
    % ==================================================================== %
    methods
        function v = get.n(obj); v = mvst.System.dimOf(obj.modes{1}); end

        function v = get.m(obj)
            m1 = obj.modes{1};
            if     isfield(m1, 'm');  v = m1.m;
            elseif isfield(m1, 'B');  v = size(m1.B, 2);
            elseif isfield(m1, 'B0'); v = size(m1.B0, 2);
            else;                     v = 0;
            end
        end

        function v = get.M(obj); v = numel(obj.modes); end
    end

    % ==================================================================== %
    %  Bridge methods                                                       %
    % ==================================================================== %
    methods
        function s = thin(obj)
            %THIN  Legacy MVST_SYSTEM / MVST_MODE struct for the Riccati layer.
            %   S = THIN(OBJ) returns the byte-compatible thin struct the
            %   optimal-control / steering / estimation engines consume.  A
            %   single bare mode round-trips to itself; a system round-trips to
            %   the MVST_SYSTEM struct.  MVST_TORICH(THIN(OBJ)) is idempotent.
            if obj.IsSystem
                args = {'Pi', obj.Pi, 'pi0', obj.pi0, 'ZT', obj.ZT};
                if ~isempty(obj.P); args = [args, {'P', obj.P}]; end
                if ~isempty(obj.J); args = [args, {'J', obj.J}]; end
                s = mvst_system(obj.modes, args{:});
            else
                s = obj.modes{1};
            end
        end

        function r = rich(obj)
            %RICH  SOS-layer struct (.A0 naming) the SOS / L2-gain / game engines consume.
            %   R = RICH(OBJ) bridges the object to the rich convention through
            %   MVST_TORICH.  The bridge is LOSSLESS in both cases:
            %
            %     * SINGLE mode (M==1): R is MVST_TORICH applied to the FLAT thin
            %       mode itself (MODE_TO_RICH), so EVERY field is carried through
            %       (A->A0, B->B0, plus .Ac/.Bc, the jump/reset .J/.BJ, the
            %       Poisson data .Jf/.lamf/.BJf, and all dissipativity I/O blocks
            %       .Bw/.C/.Dw/.Dwc/.Bd/.Bdf/.CJ/.DJ).  This is the field set the
            %       exact L2-gain / game / performance backends read, so a Poisson
            %       / Levy / sampled / impulsive factory object routes to the
            %       correct backend rather than silently to the 'ct' one.
            %
            %     * MULTI mode (M>1): R is MVST_TORICH applied to THIN() (the
            %       mode-indexed SYSTEM_TO_RICH: A0 n x n x M, stacked Ac/B0/Bc,
            %       the generator .Pi/.Pig, the DT transition .P and the common
            %       reset .J), with the mode-SHARED dissipativity I/O and Poisson
            %       jump data re-attached from the first mode.  Modes carrying
            %       DIFFERING per-mode .Bw/.C/.Dw raise 'mvst:System:permodeIO'
            %       (per-mode disturbance / output data is not supported).
            if obj.M == 1
                r = mvst_torich(obj.modes{1});      % flat, lossless single mode
            else
                r = mvst_torich(obj.thin());        % mode-indexed system
                r = obj.reattachShared(r);          % + shared I/O / Poisson data
            end
        end
    end

    % ==================================================================== %
    %  Analysis / synthesis sugar (delegations to the engine layer)         %
    % ==================================================================== %
    methods
        function res = l2gain(obj, opts)
            %L2GAIN  (Mean-square) L2-gain, via MVST_L2GAIN.
            %   RES = L2GAIN(OBJ) computes the L2-gain of the system, routing
            %   to the matching exact backend.  Requires the disturbance input
            %   .Bw and output .C (supply them through the factory).  OPTS is
            %   forwarded to MVST_L2GAIN.
            %
            %   Example:
            %       S = mvst.System.ct([-1 0.5; 0 -2], zeros(2,0), ...
            %                          'Bw', [1;0], 'C', [1 0]);
            %       S.l2gain().gamma
            arguments
                obj  (1,1) mvst.System
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_l2gain(obj.routeStruct(opts), opts);
        end

        function res = steer(obj, Xf, T, opts)
            %STEER  Finite-time optimal covariance steering, via MVST_STEER.
            %   RES = STEER(OBJ,XF,T) steers the state covariance to XF over
            %   [0,T].  OPTS is forwarded to MVST_STEER (e.g. .X0, .M0, .tol).
            %
            %   Example:
            %       S = mvst.System.ct([-1 0;0 -1], [1;1], 'X0', eye(2));
            %       r = S.steer(0.5*eye(2), 2);
            arguments
                obj  (1,1) mvst.System
                Xf   double
                T    (1,1) double {mustBePositive}
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_steer(obj.thin(), Xf, T, opts);
        end


        function res = gamegamma(obj, opts)
            %GAMEGAMMA  H-infinity game level gamma*, via MVST_GAME_GAMMA.
            %   RES = GAMEGAMMA(OBJ) computes the critical disturbance level of
            %   the state-feedback zero-sum game.  Uses the control input .Bu
            %   (= B), control weight .Ru (= R), state weight .Q, and the
            %   disturbance input .Bw (supply .Bw through the factory).  OPTS is
            %   forwarded to MVST_GAME_GAMMA.
            %
            %   Example:
            %       S = mvst.System.ct([-1 1; 0 -2], [0;1], ...
            %                          'Q', eye(2), 'R', 1, 'Bw', [1;0]);
            %       S.gamegamma().gammastar
            arguments
                obj  (1,1) mvst.System
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_game_gamma(obj.solverStruct(), opts);
        end

        function res = ocfinite(obj, T, opts)
            %OCFINITE  Finite-horizon linear-cost control, via MVST_SOLVE_FINITE.
            %   RES = OCFINITE(OBJ,T) solves the finite-horizon I-GRE on [0,T].
            %
            %   Example:
            %       S = mvst.System.ct([-1 0;0 -2], [1;1], 'Q', eye(2), 'R', 1);
            %       S.ocfinite(3).cost
            arguments
                obj  (1,1) mvst.System
                T    (1,1) double {mustBePositive}
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_solve_finite(obj.thin(), T, opts);
        end

        function res = ocinfinite(obj, opts)
            %OCINFINITE  Infinite-horizon linear-cost control, via MVST_SOLVE_INFINITE.
            %   RES = OCINFINITE(OBJ) solves the algebraic (coupled) I-GRE.
            %
            %   Example:
            %       S = mvst.System.ct([-1 0;0 -2], [1;1], 'Q', eye(2), 'R', 1);
            %       S.ocinfinite().cost
            arguments
                obj  (1,1) mvst.System
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_solve_infinite(obj.thin(), opts);
        end

        function res = ocaverage(obj, opts)
            %OCAVERAGE  Long-term average-cost control, via MVST_SOLVE_AVERAGE.
            %   RES = OCAVERAGE(OBJ) returns the optimal average cost rate .rho
            %   (weighted by the stationary distribution) on top of the
            %   infinite-horizon solution.
            %
            %   Example:
            %       S = mvst.System.ct([-1 0;0 -2], [1;1], 'W', eye(2));
            %       S.ocaverage().rho
            arguments
                obj  (1,1) mvst.System
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_solve_average(obj.thin(), opts);
        end

        function res = stability(obj, opts)
            %STABILITY  Exact mean-square stability analysis, via MVST_STABILITY.
            %   RES = STABILITY(OBJ) reads the stability of the (autonomous)
            %   system from the exact spectrum of its lifted second-moment
            %   operator, routing to the criterion of the resolved class.  OPTS
            %   is forwarded to MVST_STABILITY (e.g. .T dwell grid, .class, .tol).
            %
            %   Example:
            %       S = mvst.System.ct([0 1;-2 -1], zeros(2,0), 'Ac', 0.15*eye(2));
            %       S.stability().stable
            arguments
                obj  (1,1) mvst.System
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_stability(obj, opts);
        end

        function res = stabilize(obj, opts)
            %STABILIZE  Exact state-feedback stabilization, via MVST_STABILIZE.
            %   RES = STABILIZE(OBJ) designs a mean-square stabilizing gain by
            %   the closed-form Riccati route of the resolved class and certifies
            %   the closed loop through MVST_STABILITY.  OPTS is forwarded to
            %   MVST_STABILIZE.  Discrete / dwell-time-constrained synthesis is
            %   handled by MVST_STABZ_SOS instead.
            %
            %   Example:
            %       S = mvst.System.ct([0 1;-1 0.2], [0;1], 'Q', eye(2), 'R', 1);
            %       S.stabilize().check.stable
            arguments
                obj  (1,1) mvst.System
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_stabilize(obj, opts);
        end

        function res = dissyn(obj, supply, dwell, mode, opts)
            %DISSYN  Dissipative state-feedback synthesis, via MVST_DISSYN_SOS.
            %   RES = DISSYN(OBJ) synthesizes a state-feedback gain rendering the
            %   closed loop dissipative with respect to SUPPLY over the DWELL-time
            %   regime, in MODE ('gamma' optimizes the supply level, 'feas' only
            %   certifies feasibility).  SUPPLY defaults to the L2-gain supply,
            %   DWELL to the arbitrary (fixed-sequence / CT LTI) regime.  All three
            %   arguments are forwarded to MVST_DISSYN_SOS.
            %
            %   Example:
            %       S = mvst.System.ct([0 1;-2 -1], [0;1], ...
            %                          'Bw', eye(2), 'C', [1 0;0 0]);
            %       S.dissyn().gammastar
            arguments
                obj    (1,1) mvst.System
                supply struct = struct('type', 'l2gain')
                dwell  struct = struct('type', 'arbitrary')
                mode   char   = 'gamma'
                opts   = struct()       % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_dissyn_sos(obj.solverStruct(), supply, dwell, mode, opts);
        end

        function res = pos(obj, what, opts)
            %POS  Positive-systems copositive-LP analysis / synthesis, via the
            %   MVST_POS_* family.  RES = POS(OBJ,WHAT) routes to
            %       'stability'  MVST_POS_STABILITY  (Metzler / dwell-time LP),
            %       'gain'       MVST_POS_GAIN       (induced L1 / Linf gain),
            %       'stabilize'  MVST_POS_STABILIZE  (positivity-preserving gain),
            %       'supply'     MVST_POS_SUPPLY     (general linear supply),
            %   passing OPTS through unchanged.  An unknown WHAT raises
            %   'mvst:System:pos'.
            %
            %   Unlike the Riccati / SOS engines, the positive-systems family
            %   reads the drift and I/O data at the TOP level of a FLAT struct
            %   (SYS.A a matrix for a single mode, a 1xN cell of drifts for a
            %   multi-mode switched / Markov-jump system), not through the
            %   .modes wrapper THIN() returns.  POS bridges the object to that
            %   contract through POSSTRUCT before delegating.
            %
            %   Example:
            %       S = mvst.System.ct([-1 0; 3 -1], zeros(2,0));
            %       S.pos('stability', struct('class','ct')).stable
            arguments
                obj  (1,1) mvst.System
                what char
                opts = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            s = obj.posStruct();
            switch lower(what)
                case 'stability'; res = mvst_pos_stability(s, opts);
                case 'gain';      res = mvst_pos_gain(s, opts);
                case 'stabilize'; res = mvst_pos_stabilize(s, opts);
                case 'supply';    res = mvst_pos_supply(s, opts);
                case 'oc';        res = mvst_pos_oc(s, opts);
                otherwise
                    error('mvst:System:pos', ...
                        'Unknown WHAT ''%s''; use stability|gain|stabilize|supply|oc.', ...
                        what);
            end
        end

        function res = perf(obj, measure, opts)
            %PERF  Further-performance-measure analysis, via MVST_PERF.
            %   RES = PERF(OBJ,MEASURE) computes the named performance MEASURE
            %   ('l2','tracegain','e2p','peak','h2', ...), routing to the matching
            %   exact engine of the resolved class.  OPTS is forwarded to
            %   MVST_PERF.
            %
            %   Example:
            %       S = mvst.System.ct([-1 0.5; 0 -2], zeros(2,0), ...
            %                          'Bw', [1;0], 'C', [1 0]);
            %       S.perf('h2').gamma
            arguments
                obj     (1,1) mvst.System
                measure char
                opts    = struct()      % struct or mvst.Options
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_perf(obj.solverStruct(), measure, opts);
        end

        function res = hinfnorm(obj, opts)
            %HINFNORM  H-infinity norm / L2-gain with certificate, via MVST_HINFNORM.
            %   RES = HINFNORM(OBJ) returns the uniform struct of MVST_HINFNORM
            %   (.gamma, .certificate, .class, .method, .diagnostics) for the
            %   resolved class.  For a discrete-time system pass OPTS.class = 'dt'
            %   (single mode, drift = one-step map); a system with a DT transition
            %   matrix (tag 'dtmjls') is routed to the DT-MJLS backend.  OPTS may be
            %   a struct or an MVST.OPTIONS.
            %
            %   Example:
            %       S = mvst.System.ct([-1 0.5; 0 -2], zeros(2,0), ...
            %                          'Bw', [1;0], 'C', [1 0]);
            %       r = S.hinfnorm();  r.gamma
            arguments
                obj  (1,1) mvst.System
                opts = struct()
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_hinfnorm(obj.routeStruct(opts), opts);
        end

        function res = passivity(obj, opts)
            %PASSIVITY  Passivity test with ISP/OSP indices, via MVST_PASSIVITY.
            %   RES = PASSIVITY(OBJ) returns struct(passive, nu_isp, rho_osp,
            %   certificate, method, class, diagnostics); the CT LTI default is the
            %   exact Hamiltonian test, other classes route to the SOS engines.
            %   The index convention is the book's (supply He(w'z) - nu|w|^2).
            %   Needs a square port (size(C,1) == size(Bw,2)).
            %
            %   Example:
            %       S = mvst.System.ct([0 1;-1 -1], zeros(2,0), 'Bw', [0;1], ...
            %                          'C', [1 0], 'Dw', 0.5);
            %       r = S.passivity();  r.passive, r.nu_isp
            arguments
                obj  (1,1) mvst.System
                opts = struct()
            end
            opts = mvst.System.normOpts(opts);
            res = mvst_passivity(obj.routeStruct(opts), opts);
        end

        function res = kyp(obj, opts)
            %KYP  Operator-domain KYP engine, via MVST_KYP (deterministic classes).
            %   RES = KYP(OBJ, OPTS) forwards the single-mode deterministic LTI
            %   system to MVST_KYP.  OPTS.domain is 'ct' (default), 'dt' or
            %   'hybrid' (needs the sampled / impulsive tag; OPTS.Tp defaults to
            %   the sampling period, OPTS.J to the reset); OPTS.supply, .gamma, .M
            %   are those of MVST_KYP.  Stochastic or multi-mode classes (noise
            %   channels, MJLS, switched, Poisson) are refused with
            %   'mvst:System:kyp' since the KYP lemma of the book is stated for
            %   the deterministic operator.
            %
            %   Example:
            %       S = mvst.System.ct([-1 0.5; 0 -2], zeros(2,0), ...
            %                          'Bw', [1;0], 'C', [1 0]);
            %       r = S.kyp(struct('supply','boundedreal','gamma',2));  r.feas
            arguments
                obj  (1,1) mvst.System
                opts = struct()
            end
            opts = mvst.System.normOpts(opts);
            ok = {'ct','sampled','impulsive'};
            if obj.M ~= 1 || ~any(strcmp(obj.tag, ok))
                error('mvst:System:kyp', ...
                    ['KYP needs a single-mode deterministic system (tags ct, sampled, ' ...
                     'impulsive); this system has tag ''%s'' and %d mode(s).'], obj.tag, obj.M);
            end
            m1 = obj.modes{1};
            for f = {'Ac','Bc','Dwc','Jf','lamf'}
                if isfield(m1, f{1}) && ~isempty(m1.(f{1})) && any(m1.(f{1})(:) ~= 0)
                    error('mvst:System:kyp', ...
                        'KYP is deterministic: the system carries a nonzero .%s channel.', f{1});
                end
            end
            dom = lower(mvst_getf(opts, 'domain', 'ct'));
            if strcmp(dom, 'hybrid')
                if ~any(strcmp(obj.tag, {'sampled','impulsive'}))
                    error('mvst:System:kyp', 'domain ''hybrid'' needs a sampled or impulsive system.');
                end
                if ~isfield(opts, 'Tp') && isfield(m1, 'Tbar'); opts.Tp = m1.Tbar; end
                if ~isfield(opts, 'J') && ~isempty(obj.J);       opts.J  = obj.J;   end
            end
            sys = struct('A', m1.A, 'Bw', mvst_getf(m1, 'Bw', zeros(size(m1.A,1), 0)), ...
                         'C', mvst_getf(m1, 'C', zeros(0, size(m1.A,1))));
            sys.Dw = mvst_getf(m1, 'Dw', zeros(size(sys.C,1), size(sys.Bw,2)));
            res = mvst_kyp(sys, opts);
        end

        function res = dissipativity(obj, supply, dwell, mode, opts)
            %DISSIPATIVITY  SOS dissipativity analysis for a named or (Q,S,R) supply.
            %   RES = DISSIPATIVITY(OBJ, SUPPLY) forwards to MVST_DISS_SOS (single
            %   mode) or MVST_DISS_SOS_MODES (MJLS / switched).  SUPPLY is either
            %   the string 'l2gain' | 'passivity', a struct with .type (and
            %   .gamma, .nu, .rho), or the QSR CONTRACT
            %       struct('type','qsr','Q',Q,'S',S,'R',R)
            %   meaning the quadratic supply [z;w]'*[Q S;S' R]*[z;w] with
            %   z = C x + Dw w (book sign convention: l2gain is Q=-I, S=0,
            %   R=gamma^2 I; passivity is Q=-rho I, S=I, R=-nu I).  The shapes
            %   Q (p x p symmetric), S (p x m), R (m x m symmetric) are validated
            %   here, error 'mvst:System:supply'.  DWELL (struct) defaults to the
            %   class regime (arbitrary; ctmjls for tag mjls), MODE to 'feas'.
            %
            %   Example:
            %       S = mvst.System.ct([-1 0.5; 0 -2], zeros(2,0), ...
            %                          'Bw', [1;0], 'C', [1 0]);
            %       r = S.dissipativity(struct('type','qsr','Q',-1,'S',0,'R',9));
            %       r.feas
            arguments
                obj    (1,1) mvst.System
                supply = 'l2gain'
                dwell  = []
                mode   char = 'feas'
                opts   = struct()
            end
            opts = mvst.System.normOpts(opts);
            r = obj.solverStruct();
            C = mvst_getf(r, 'C', []);
            Bw = mvst_getf(r, 'Bw', []);
            if isempty(C) || isempty(Bw)
                error('mvst:System:supply', 'The system needs a disturbance input .Bw and output .C.');
            end
            p = size(C,1); m = size(Bw,2);
            if ischar(supply) || isstring(supply)
                supply = struct('type', char(lower(supply)));
            end
            if ~isstruct(supply) || ~isfield(supply, 'type')
                error('mvst:System:supply', 'SUPPLY must be ''l2gain''|''passivity'' or a struct with .type.');
            end
            switch lower(supply.type)
                case 'qsr'
                    mvst_qsr_supply(supply, p, m);   % single-sourced validation (ids, tolerances)
                    if any(isfield(supply, {'Qk','Sk','Rk'}))
                        mvst_qsr_supply(supply, p, m, 'jump');
                    end
                case {'l2gain','passivity'}
                otherwise
                    error('mvst:System:supply', 'Unknown supply type ''%s''.', supply.type);
            end
            if isempty(dwell)
                switch obj.tag
                    case 'mjls';     dwell = struct('type', 'ctmjls');
                    case 'switched'; dwell = struct('type', 'switched', 'Tbar', mvst_getf(obj.modes{1}, 'Td', 1));
                    otherwise;       dwell = struct('type', 'arbitrary');
                end
            end
            if obj.M > 1
                res = mvst_diss_sos_modes(r, supply, dwell, mode, opts);
            else
                res = mvst_diss_sos(r, supply, dwell, mode, opts);
            end
        end
    end

    % ==================================================================== %
    %  Private helpers used by the sugar methods                            %
    % ==================================================================== %
    methods (Access = private)
        function r = routeStruct(obj, opts)
            %ROUTESTRUCT  solverStruct() with the discrete-time markers for the
            %   hinfnorm / passivity routers: a DT-MJLS system (tag 'dtmjls')
            %   gets .domain = 'dt' and its transition matrix .Ptr; OPTS.class
            %   'dt' / 'dtmjls' marks a single-mode DT system.
            r = obj.solverStruct();
            isdt = strcmp(obj.tag, 'dtmjls') || (isfield(opts, 'class') && ...
                   any(strcmpi(opts.class, {'dt','dtmjls'})));
            if isdt
                r.domain = 'dt';
                if ~isempty(obj.P); r.Ptr = obj.P; end
            end
        end

        function r = solverStruct(obj)
            %SOLVERSTRUCT  Rich struct with the I/O and game data re-attached.
            %   The L2-gain / game / performance backends read .Bu/.Bw/.Q/.Ru
            %   (and the class-specific .Jf/.lamf/.Ac/.Dwc/.Tbar/.Td) at the top
            %   level.  RICH() is now LOSSLESS (it carries the disturbance /
            %   output / Poisson jump data for both the single- and multi-mode
            %   cases), so this helper only adds the cost data (.Q/.Ru) and the
            %   control-input alias (.Bu = .B0) the game backend expects, plus the
            %   dwell / sampling period (.Tbar/.Td) that MVST_TORICH does not map.
            r  = obj.rich();
            m1 = obj.modes{1};
            if isfield(r, 'B0') && ~isfield(r, 'Bu'); r.Bu = r.B0; end
            if isfield(m1, 'R') && ~isfield(r, 'Ru'); r.Ru = m1.R;  end
            if isfield(m1, 'Q') && ~isfield(r, 'Q');  r.Q  = m1.Q;  end
            passthru = {'Bw', 'C', 'Dw', 'Dwc', 'Tbar', 'Td'};
            for f = passthru
                if isfield(m1, f{1}) && ~isfield(r, f{1}); r.(f{1}) = m1.(f{1}); end
            end
        end

        function r = reattachShared(obj, r)
            %REATTACHSHARED  Re-attach the mode-shared I/O / Poisson jump fields
            %   that SYSTEM_TO_RICH drops from a multi-mode rich struct.
            %   MVST_TORICH stacks only the drift / noise / coupling data of a
            %   multi-mode system; the disturbance / output blocks and the
            %   spontaneous-jump data live on the modes and are dropped.  This
            %   helper copies them from the first mode, first ASSERTING that the
            %   per-mode I/O is identical across modes (per-mode disturbance /
            %   output data is unsupported: silent dropping is forbidden).
            ms = obj.modes;
            m1 = ms{1};
            iocheck = {'Bw', 'C', 'Dw', 'Dwc', 'Bd', 'Bdf', 'CJ', 'DJ'};
            for f = iocheck
                a = mvst_getf(m1, f{1}, []);
                for i = 2:numel(ms)
                    b = mvst_getf(ms{i}, f{1}, []);
                    if ~isequal(a, b)
                        error('mvst:System:permodeIO', ...
                            ['Modes carry differing per-mode .%s; the rich ', ...
                             'bridge only supports mode-shared disturbance / ', ...
                             'output data.'], f{1});
                    end
                end
            end
            shared = {'Bw', 'C', 'Dw', 'Dwc', 'Bd', 'Bdf', 'CJ', 'DJ', ...
                      'Jf', 'lamf', 'BJf'};
            for f = shared
                if isfield(m1, f{1}) && ~isempty(m1.(f{1})) && ~isfield(r, f{1})
                    r.(f{1}) = m1.(f{1});
                end
            end
        end

        function s = posStruct(obj)
            %POSSTRUCT  Flat struct the MVST_POS_* copositive-LP family reads.
            %   The positive-systems engines read the drift and I/O blocks at
            %   the TOP level (SYS.A a matrix for a single mode, a 1xN cell of
            %   per-mode drifts for a switched / Markov-jump system), NOT through
            %   the .modes wrapper THIN() produces.  This helper conforms to
            %   that contract: for a single-mode system it returns the (already
            %   flat) thin mode struct, whose .A / .J / .Tbar / .Bw / .C / .Dw
            %   sit at the top level and whose empty noise .Bc falls through to
            %   the flow .B in the engines (MVST_GETF treats empty as absent);
            %   for a multi-mode system it gathers the per-mode drifts into the
            %   cell arrangement MVST_POS_STABILITY expects and re-attaches the
            %   generator (MJLS coupled program) and any common reset.
            if obj.M == 1
                s = obj.modes{1};
            else
                % A multi-mode system with a COMMON RESET has no home in the
                % positive-systems family: MVST_POS_STABILITY has no
                % switched/Markov-jump-with-reset class, and emitting both a cell
                % .A and a .J would misroute to the single-mode 'impulsive' LP
                % (which then fails on the cell .A).  Refuse explicitly instead of
                % silently dropping the reset.
                if ~isempty(obj.J)
                    error('mvst:System:posUnsupported', ...
                        ['The positive-systems layer has no multi-mode ', ...
                         '(switched / Markov-jump) with-reset class; a ', ...
                         'multi-mode system carrying a common reset J cannot ', ...
                         'be routed to the MVST_POS_* family.  Drop the reset, ', ...
                         'or use a single-mode impulsive system for the ', ...
                         'reset-based dwell-time program.']);
                end
                s = struct();
                s.A = cellfun(@(md) md.A, obj.modes, 'UniformOutput', false);
                if ~isempty(obj.Pi) && any(abs(obj.Pi(:)) > 0)
                    s.Pi = obj.Pi;                 % coupled MJLS copositive LP
                end
            end
        end
    end

    % ==================================================================== %
    %  Static validation / construction utilities                           %
    % ==================================================================== %
    methods (Static, Access = private)
        function opts = normOpts(opts)
            %NORMOPTS  Accept an MVST.OPTIONS wherever a struct is expected.
            %   An MVST.OPTIONS object becomes the plain struct (.solver, .reltol,
            %   .verbosity) through MVST_OPTS; a struct is returned unchanged.
            if isa(opts, 'mvst.Options')
                opts = mvst_opts(opts, false);
            elseif isempty(opts)
                opts = struct();
            elseif ~isstruct(opts)
                error('mvst:System:opts', 'OPTS must be a struct or an mvst.Options object.');
            end
        end

        function n = checkSquare(A, name)
            %CHECKSQUARE  Assert A is square; return its dimension.
            if size(A,1) ~= size(A,2)
                error('mvst:System:notSquare', ...
                    '%s must be square; got %d x %d.', name, size(A,1), size(A,2));
            end
            n = size(A,1);
        end

        function checkRows(B, n, name)
            %CHECKROWS  Assert B has N rows (matching the state dimension).
            if ~isempty(B) && size(B,1) ~= n
                error('mvst:System:dim', ...
                    '%s must have %d rows to match the state dimension.', name, n);
            end
        end

        function checkImpulsiveChannels(opts, n, m, p)
            %CHECKIMPULSIVECHANNELS  Validate the general impulsive channel set.
            %   Each channel is checked only when supplied (empty = absent).  NW
            %   (disturbance count) is read from .Bw when present; column checks
            %   that depend on it are skipped when .Bw is absent.  Errors use the
            %   mvst:System:dim identifier.
            dimerr = @(msg, varargin) error('mvst:System:dim', msg, varargin{:});
            has = @(f) isfield(opts, f) && ~isempty(opts.(f));
            nw = [];
            if has('Bw'); nw = size(opts.Bw, 2); end
            q = 0; if has('Ac'); q = size(opts.Ac, 3); end
            % spontaneous (Poisson) jump maps and rates
            if has('Jf')
                if size(opts.Jf,1) ~= n || size(opts.Jf,2) ~= n
                    dimerr('Jf must be %d x %d [x pf].', n, n);
                end
            end
            pf = 0; if has('Jf'); pf = size(opts.Jf,3); end
            if has('lamf')
                if has('Jf') && numel(opts.lamf) ~= pf
                    dimerr('lamf must have one rate per Jf channel (%d).', pf);
                end
                if ~has('Jf'); pf = numel(opts.lamf); end
            end
            if has('BJf')
                if size(opts.BJf,1) ~= n
                    dimerr('BJf must have %d rows to match the state.', n);
                end
                if size(opts.BJf,2) ~= m
                    dimerr('BJf must have %d columns to match the input.', m);
                end
                if pf > 0 && size(opts.BJf,3) ~= pf
                    dimerr('BJf must have %d jump channels to match Jf/lamf.', pf);
                end
            end
            % scheduled-jump control selector
            if has('BJ')
                if size(opts.BJ,1) ~= n
                    dimerr('BJ must have %d rows to match the state.', n);
                end
                if size(opts.BJ,2) ~= m
                    dimerr('BJ must have %d columns to match the input.', m);
                end
                if size(opts.BJ,3) ~= 1 && size(opts.BJ,3) ~= p
                    dimerr('BJ must have 1 or %d jump channels to match J.', p);
                end
            end
            % flow disturbance-multiplicative noise (paired with Ac diffusion)
            if has('Dwc')
                if size(opts.Dwc,1) ~= n
                    dimerr('Dwc must have %d rows to match the state.', n);
                end
                if ~isempty(nw) && size(opts.Dwc,2) ~= nw
                    dimerr('Dwc must have %d columns to match the disturbance.', nw);
                end
                if q > 0 && size(opts.Dwc,3) ~= q
                    dimerr('Dwc must have %d channels to match Ac.', q);
                end
            end
            % jump disturbance channels (scheduled Bd, spontaneous Bdf)
            if has('Bd')
                if size(opts.Bd,1) ~= n
                    dimerr('Bd must have %d rows to match the state.', n);
                end
                if ~isempty(nw) && size(opts.Bd,2) ~= nw
                    dimerr('Bd must have %d columns to match the disturbance.', nw);
                end
            end
            if has('Bdf')
                if size(opts.Bdf,1) ~= n
                    dimerr('Bdf must have %d rows to match the state.', n);
                end
                if ~isempty(nw) && size(opts.Bdf,2) ~= nw
                    dimerr('Bdf must have %d columns to match the disturbance.', nw);
                end
            end
            % jump performance output / feedthrough
            if has('CJ') && size(opts.CJ,2) ~= n
                dimerr('CJ must have %d columns to match the state.', n);
            end
            if has('CJ') && has('DJ') && size(opts.DJ,1) ~= size(opts.CJ,1)
                dimerr('DJ must have %d rows to match CJ.', size(opts.CJ,1));
            end
        end

        function checkGenerator(Pi, M)
            %CHECKGENERATOR  Assert PI is an M x M CT Markov generator.
            if ~isequal(size(Pi), [M M])
                error('mvst:System:generator', 'PI must be %d x %d.', M, M);
            end
            off = Pi - diag(diag(Pi));
            if any(off(:) < -1e-9)
                error('mvst:System:generator', ...
                    'PI must be Metzler (nonnegative off-diagonal entries).');
            end
            if any(abs(sum(Pi, 2)) > 1e-8)
                error('mvst:System:generator', ...
                    'PI rows must sum to zero (continuous-time generator).');
            end
        end

        function md = buildMode(A, B, opts)
            %BUILDMODE  MVST_MODE call plus optional I/O passthrough attach.
            args = {};
            for f = {'Ac', 'Bc', 'S', 'W', 'X0'}
                if isfield(opts, f{1}) && ~isempty(opts.(f{1}))
                    args = [args, {f{1}, opts.(f{1})}]; %#ok<AGROW>
                end
            end
            Q = mvst.System.pickOr(opts, 'Q', zeros(size(A,1)));
            R = mvst.System.pickOr(opts, 'R', eye(size(B,2)));
            md = mvst_mode(A, B, Q, R, args{:});
            for f = {'Bw', 'C', 'Dw'}
                if isfield(opts, f{1}) && ~isempty(opts.(f{1}))
                    md.(f{1}) = opts.(f{1});
                end
            end
        end

        function v = pickOr(opts, f, d)
            if isfield(opts, f) && ~isempty(opts.(f)); v = opts.(f); else; v = d; end
        end

        function ZT = defaultZT(ZT, n)
            if isempty(ZT); ZT = zeros(n); end
        end

        function pi0 = defaultPi0(pi0, M)
            if isempty(pi0); pi0 = [1, zeros(1, M-1)]; end
        end

        function nn = dimOf(md)
            %DIMOF  State dimension of a thin or rich mode struct.
            if     isfield(md, 'n');  nn = md.n;
            elseif isfield(md, 'A');  nn = size(md.A, 1);
            elseif isfield(md, 'A0'); nn = size(md.A0, 1);
            else;                     nn = 0;
            end
        end

        function tag = inferClass(modes, Pi, P, J) %#ok<INUSD>
            %INFERCLASS  Best-effort class tag from the assembled model.
            m1 = modes{1};
            M  = numel(modes);
            has = @(f) isfield(m1, f) && ~isempty(m1.(f));
            hasJf   = has('Jf') || has('lamf');
            hasAc   = has('Ac') && size(m1.Ac, 3) > 0;
            hasTbar = has('Tbar');
            hasReset= has('J');
            PiNZ    = ~isempty(Pi) && any(abs(Pi(:)) > 0);
            if     hasJf && hasAc;        tag = 'levy';
            elseif hasJf;                 tag = 'poisson';
            elseif hasTbar;               tag = 'sampled';
            elseif M > 1 && ~isempty(P);  tag = 'dtmjls';   % DT transition -> discrete MJLS
            elseif M > 1 && PiNZ;         tag = 'mjls';
            elseif M > 1;                 tag = 'switched';
            elseif hasReset;              tag = 'impulsive';
            else;                         tag = 'ct';
            end
        end

        function [modes, Pi, pi0, ZT, P, J, isSys] = richToModel(sys)
            %RICHTOMODEL  Permissive rich (.A0) struct -> thin modes + coupling.
            n   = size(sys.A0, 1);
            Msz = size(sys.A0, 3);
            modes = cell(1, Msz);
            for i = 1:Msz
                modes{i} = mvst.System.richMode(sys, i);
            end
            P  = mvst_getf(sys, 'P', []);
            J  = mvst_getf(sys, 'J', []);
            ZT = mvst_getf(sys, 'ZT', zeros(n));
            if Msz == 1
                Pi = zeros(1); pi0 = 1; isSys = false;
            else
                Pi  = mvst_getf(sys, 'Pi', mvst_getf(sys, 'Pig', zeros(Msz)));
                pi0 = [1, zeros(1, Msz-1)];
                isSys = true;
            end
        end

        function md = richMode(sys, i)
            %RICHMODE  Extract mode I of a rich struct as a thin mode struct.
            n = size(sys.A0, 1);
            A = sys.A0(:, :, min(i, size(sys.A0, 3)));
            if isfield(sys, 'B0') && ~isempty(sys.B0)
                B = sys.B0(:, :, min(i, size(sys.B0, 3)));
            else
                B = zeros(n, 0);
            end
            m = size(B, 2);
            args = {};
            if isfield(sys, 'Ac') && ~isempty(sys.Ac)
                args = [args, {'Ac', mvst.System.slice4(sys.Ac, i)}];
            end
            if isfield(sys, 'Bc') && ~isempty(sys.Bc)
                args = [args, {'Bc', mvst.System.slice4(sys.Bc, i)}];
            end
            md = mvst_mode(A, B, zeros(n), eye(m), args{:});
            % copy jump / reset / dissipativity passthrough (mode-shared)
            for f = {'J','BJ','Jf','lamf','BJf','Bw','Dwc','C','Dw','Bdf','Bd','CJ','DJ','Tbar','Td'}
                if isfield(sys, f{1}) && ~isempty(sys.(f{1}))
                    md.(f{1}) = sys.(f{1});
                end
            end
        end

        function X = slice4(X, i)
            %SLICE4  Mode slice of a possibly 4-D noise array (n x n x q x M).
            if ndims(X) >= 4
                X = X(:, :, :, min(i, size(X, 4)));
            end
        end
    end
end
