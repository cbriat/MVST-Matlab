function [prog, info] = mvst_sossolve(prog, opts)
%MVST_SOSSOLVE  Solve a SOSTOOLS program with the backend selected by OPTS.SOLVER.
%   [PROG, INFO] = MVST_SOSSOLVE(PROG, OPTS) replaces every direct SOSSOLVE call
%   in the engines.  OPTS.SOLVER selects the backend:
%
%     'sedumi' | 'sdpt3' | 'mosek'  ->  SOSTOOLS sossolve(PROG, OPTS) (classical
%                                       SDP route; OPTS passed unchanged).
%     'nosdp'                       ->  sossolve_nosdp(PROG) (SOSNOSDP, no SDP
%                                       solver).  Options forwarded: every field of
%                                       OPTS except .solver and .params, plus the
%                                       fields of OPTS.nosdp (struct, wins on
%                                       clashes); reltol -> tol and verbosity -> verbose (mvst.Options
%                                       fields); OPTS.params.fid == 0 gives
%                                       verbose = 0 unless set explicitly.
%
%   INFO always carries the SeDuMi-style fields the engines test, i.e. the clean
%   solve test (INFO.pinf == 0 && INFO.dinf == 0 && INFO.numerr < 2).  For the
%   'nosdp' route they are derived from SOSNOSDP's INFO.status (the original
%   status is kept in INFO.status and in INFO.nosdp_status):
%
%     status                                  pinf dinf numerr feasratio
%     'optimal', 'feasible'                    0    0    0       1
%     'feasible_max_iter', 'feasible_<other>'  0    0    1       1   (certified, reduced accuracy)
%     'infeasible', 'infeasible_equalities'    0    1    0      -1
%     'uncertified_*', 'max_iter', 'stalled'   0    0    2       0   (not trusted)
%     anything else                            0    0    2       0
%
%   (Infeasible SOS programs may report pinf = 1 or dinf = 1 depending on the backend; the
%   nosdp mapping uses dinf = 1, as SeDuMi does for a primal-infeasible program.)
%   When no solution was produced PROG.solinfo.x is set to NaN so that SOSGETSOL
%   returns NaN rather than erroring, as for the classical solvers.
%
%   See also SOSSOLVE, SOSSOLVE_NOSDP, MVST_OPTS.
if nargin < 2 || isempty(opts), opts = struct(); end
if isa(opts, 'mvst.Options')
    opts = struct('solver', opts.solver, 'reltol', opts.reltol, 'verbosity', opts.verbosity);
end
if isfield(opts, 'solver') && strcmpi(char(opts.solver), 'nosdp')
    no = struct();
    fn = fieldnames(opts);
    for k = 1:numel(fn)
        if ~any(strcmp(fn{k}, {'solver', 'params', 'nosdp'}))
            no.(fn{k}) = opts.(fn{k});
        end
    end
    if isfield(no, 'reltol'),    if ~isfield(no, 'tol'),     no.tol = no.reltol;         end, no = rmfield(no, 'reltol');    end
    if isfield(no, 'verbosity'), if ~isfield(no, 'verbose'), no.verbose = no.verbosity > 0; end, no = rmfield(no, 'verbosity'); end
    if isfield(opts, 'params') && isstruct(opts.params) && isfield(opts.params, 'fid') ...
            && opts.params.fid == 0 && ~isfield(no, 'verbose')
        no.verbose = 0;
    end
    if isfield(opts, 'nosdp') && isstruct(opts.nosdp)
        fn = fieldnames(opts.nosdp);
        for k = 1:numel(fn), no.(fn{k}) = opts.nosdp.(fn{k}); end
    end
    [prog, info] = sossolve_nosdp(prog, no);
    st = info.status;
    info.nosdp_status = st;
    if any(strcmp(st, {'optimal', 'feasible'}))
        pinf = 0; dinf = 0; numerr = 0; fr = 1;
    elseif strncmp(st, 'feasible', 8)
        pinf = 0; dinf = 0; numerr = 1; fr = 1;
    elseif any(strcmp(st, {'infeasible', 'infeasible_equalities'}))
        pinf = 0; dinf = 1; numerr = 0; fr = -1;
    else
        pinf = 0; dinf = 0; numerr = 2; fr = 0;
    end
    info.pinf = pinf; info.dinf = dinf; info.numerr = numerr; info.feasratio = fr;
    if ~isfield(info, 'iter'), info.iter = info.iters; end
    if isempty(prog.solinfo.x)
        nU = prog.var.idx{end} - 1;
        prog.solinfo.x = nan(nU, 1); prog.solinfo.RRx = nan(nU, 1);
    end
    prog.solinfo.info = info;
else
    % classical backends: reltol -> the backend's own tolerance parameter (unless set explicitly)
    if isfield(opts, 'reltol') && isfield(opts, 'solver')
        if ~isfield(opts, 'params') || ~isstruct(opts.params), opts.params = struct(); end
        switch lower(char(opts.solver))
            case 'sedumi'; tolf = 'eps';
            case 'sdpt3';  tolf = 'gaptol';
            case 'mosek';  tolf = 'MSK_DPAR_INTPNT_CO_TOL_REL_GAP';
            otherwise;     tolf = '';
        end
        if ~isempty(tolf) && ~isfield(opts.params, tolf), opts.params.(tolf) = opts.reltol; end
    end
    [prog, info] = sossolve(prog, opts);
end
end
