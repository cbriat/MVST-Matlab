classdef Options
%OPTIONS  Solver options for the MVST package.
%   O = MVST.OPTIONS() returns the default options.  O = MVST.OPTIONS(Name,
%   Value,...) sets any of the fields below with validation.
%
%   PROPERTIES
%     solver     'sedumi' (default) | 'sdpt3' | 'mosek' (SDP solvers) | 'nosdp'
%                (SOSNOSDP, no SDP solver; SOS engines only, not the legacy YALMIP engines).
%     reltol     relative tolerance passed to the solver (default 1e-8).
%     verbosity  solver verbosity level, 0 = silent (default 0).
%
%   THREADING.  The .solver choice is threaded through every SDP / SOS engine
%   in the package.  Pass an MVST.OPTIONS (or a plain struct with a .solver
%   field) as the trailing OPTS argument of an engine and it selects the SDP
%   backend for that solve:
%
%       mvst_stab_sos(sys, dwell, o)          mvst_diss_sos(sys,sup,dwell,mode,o)
%       mvst_stabz_sos(sys, dwell, o)         mvst_diss_sos_modes(...,mode,o)
%       mvst_stab_sos_modes(sys, dwell, o)    mvst_game_sos(sys,game,dwell,mode,o)
%       mvst_oc_sos(sys, cost, dwell, o)      mvst_dof_sos(sys, dwell, o)
%       mvst_oc_sos_modes(sys, cost, dwell, o)
%
%   Each engine normalizes OPTS through MVST_OPTS, which validates that the
%   requested backend is installed and errors ('mvst:opts:solverMissing')
%   otherwise, with NO silent fallback.  SOSTOOLS receives the choice through
%   sossolve's option field .solver; YALMIP users through sdpsettings('solver').
%   The class fixes the release-facing surface and validates inputs.
%
%   Example:
%       o = mvst.Options('solver', 'mosek', 'verbosity', 1);
%       o.solver
%       res = mvst_stab_sos(sys, dwell, o);   % solved with MOSEK if installed
%
%   See also MVST.SOLVERCHECK, MVST.SYSTEM, MVST_OPTS.

    properties
        solver    (1,:) char {mustBeMember(solver, {'sedumi','sdpt3','mosek','nosdp'})} = 'sedumi'
        reltol    (1,1) double {mustBePositive} = 1e-8
        verbosity (1,1) double {mustBeNonnegative, mustBeInteger} = 0
    end

    methods
        function obj = Options(opts)
            %OPTIONS  Construct a validated options object.
            %   O = MVST.OPTIONS('solver','sdpt3','reltol',1e-6).
            arguments
                opts.solver    (1,:) char   = 'sedumi'
                opts.reltol    (1,1) double = 1e-8
                opts.verbosity (1,1) double = 0
            end
            obj.solver    = opts.solver;
            obj.reltol    = opts.reltol;
            obj.verbosity = opts.verbosity;
        end
    end

    methods (Static)
        function o = default()
            %DEFAULT  The default options object.
            %   O = MVST.OPTIONS.DEFAULT() is equivalent to MVST.OPTIONS().
            o = mvst.Options();
        end
    end
end
