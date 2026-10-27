function out = mvst_pos_lp(lp, opts)
%MVST_POS_LP  Internal LINPROG wrapper for the positive-systems LP family.
%   OUT = MVST_POS_LP(LP) solves the linear program described by the struct LP
%   with MATLAB's LINPROG and returns a uniform result.  It is the single
%   solver entry point shared by MVST_POS_STABILITY, MVST_POS_GAIN,
%   MVST_POS_STABILIZE and MVST_POS_SUPPLY, so that solver selection, option
%   pass-through and the Optimization-Toolbox availability check live in one
%   place.
%
%   LP FIELDS (LINPROG naming; empty or absent fields default sensibly):
%     .f     objective vector (default: all zeros, i.e. a pure feasibility LP)
%     .Aub   inequality matrix  Aub x <= bub          (default [])
%     .bub   inequality vector                        (default [])
%     .Aeq   equality matrix    Aeq x =  beq          (default [])
%     .beq   equality vector                          (default [])
%     .lb    lower bounds                             (default [])
%     .ub    upper bounds                             (default [])
%     .nvar  number of variables (used when .f is omitted)
%
%   OPTS FIELDS (all optional):
%     .solver     solver name; only 'linprog' is supported (default 'linprog').
%     .lpalgo     LINPROG algorithm (default 'dual-simplex'; the simplex
%                 vertex solution reproduces the closed-form LP optima of the
%                 monograph exactly, unlike an interior point).
%     .lpopts     an OPTIMOPTIONS('linprog', ...) object passed through
%                 verbatim; overrides .lpalgo when given.
%     .display    LINPROG 'Display' value (default 'off').
%
%   OUT FIELDS:
%     .x        optimizer (LINPROG x), [] if infeasible;
%     .fval     optimal objective value, [] if infeasible;
%     .feas     logical feasibility (LINPROG exitflag == 1);
%     .exitflag raw LINPROG exit flag;
%     .lp       the resolved LP struct actually solved (for RES.lp echo).
%
%   Errors with id 'mvst:pos:noOptimToolbox' when LINPROG is unavailable.
%
%   This is an internal helper of the positive-systems solver family; it
%   implements no mathematics of its own.  See also MVST_POS_STABILITY,
%   MVST_POS_GAIN, MVST_POS_STABILIZE, MVST_POS_SUPPLY.

    if nargin < 2 || isempty(opts); opts = struct(); end

    solver = mvst_getf(opts, 'solver', 'linprog');
    if ~strcmpi(solver, 'linprog')
        error('mvst:pos:solver', ...
              'MVST_POS_LP supports only the ''linprog'' solver (got ''%s'').', ...
              solver);
    end
    if exist('linprog', 'file') ~= 2
        error('mvst:pos:noOptimToolbox', ...
              ['The positive-systems LP solvers require MATLAB''s ', ...
               'Optimization Toolbox (LINPROG), which is not available.']);
    end

    Aub = mvst_getf(lp, 'Aub', []);
    bub = mvst_getf(lp, 'bub', []);
    Aeq = mvst_getf(lp, 'Aeq', []);
    beq = mvst_getf(lp, 'beq', []);
    lb  = mvst_getf(lp, 'lb',  []);
    ub  = mvst_getf(lp, 'ub',  []);
    f   = mvst_getf(lp, 'f',   []);
    if isempty(f)
        nvar = mvst_getf(lp, 'nvar', size(Aub, 2));
        f = zeros(nvar, 1);
    end

    lpopts = mvst_getf(opts, 'lpopts', []);
    if isempty(lpopts)
        algo = mvst_getf(opts, 'lpalgo', 'dual-simplex');
        disp = mvst_getf(opts, 'display', 'off');
        lpopts = optimoptions('linprog', 'Algorithm', algo, 'Display', disp);
    end

    [x, fval, exitflag] = linprog(f, Aub, bub, Aeq, beq, lb, ub, lpopts);

    out = struct();
    out.exitflag = exitflag;
    out.feas     = (exitflag == 1);
    if out.feas
        out.x    = x;
        out.fval = fval;
    else
        out.x    = [];
        out.fval = [];
    end
    out.lp = struct('f', f, 'Aub', Aub, 'bub', bub, ...
                    'Aeq', Aeq, 'beq', beq, 'lb', lb, 'ub', ub);
end
