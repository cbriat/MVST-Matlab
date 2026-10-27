function res = mvst_pos_gain_sos(sys, opts)
%MVST_POS_GAIN_SOS  Timer-dependent copositive-SOS induced-gain programs for
%   positive impulsive / switched systems (minimum and range dwell-time), the
%   timer-dependent branch of MVST_POS_GAIN.  Implements the book's
%   displayed sliced conditions LITERALLY (forward timer tau = t - t_k of
%   eq:pos:timer):
%
%     Linf x ellinf  (peak gain, cor:pos:Linfslice, eq:pos:LinfLP) -- FORWARD
%       TUBE  xi(tau) >= 0:
%         dxi/dtau      >= A(tau) xi + Ec 1            on (0,T)
%         xi(0+)        >= J xi(theta) + Ed 1          jump
%         Cc xi(tau) + Fc 1 <= gamma 1                 on (0,T)
%         Cd xi(theta) + Fd 1 <= gamma 1               jump
%       minimum dwell-time adds the frozen peak tail (eq:pos:pktail):
%         A(Tmin) xi(Tmin) + Ec 1 <= 0.
%
%     L1 x ell1  (trace gain, cor:pos:L1slice, eq:pos:L1LP) -- BACKWARD ADJOINT
%       lambda(tau) >= 0:
%         dlambda/dtau + A(tau)' lambda + Cc' 1 <= 0   on (0,T)
%         Ec' lambda + Fc' 1 <= gamma 1                on (0,T)
%         J' lambda(0+) - lambda(theta) + Cd' 1 <= 0   jump
%         Ed' lambda(0+) + Fd' 1 <= gamma 1            jump
%       minimum dwell-time adds the frozen trace tail (eq:pos:trtail):
%         A(Tmin)' lambda(Tmin) + Cc' 1 <= 0.
%
%   The two programs are GENUINELY different (backward adjoint vs forward tube),
%   not transposes; the book states them separately.  NO differentiability pin
%   is imposed (the book uses the algebraic frozen tail only); range dwell-time
%   imposes the jump rows for all theta in [Tmin,Tmax] (a robust family) through
%   the interval Positivstellensatz multiplier (theta-Tmin)(Tmax-theta).
%
%   CLASSES:
%     impulsive : single Metzler A(tau), nonnegative J,Ec,Ed,Cc,Cd,Fc,Fd.
%     switched  : per-mode Metzler A_i(tau) (cells), state continuous at a
%                 switch (identity map, no disturbance injection), so the book
%                 jump row specializes to  xi_new(0) >= xi_old(theta)  (L1:
%                 lambda_new... see below) for every admissible transition
%                 old -> new; certificate is one tube per mode.
%
%   SYS FIELDS (impulsive): .A (n x n, or a function handle @(tau)->n x n for a
%     timer-varying Metzler flow), .Ec/.E (n x mc), .Cc/.C (qc x n),
%     .Fc/.F (qc x mc; default 0), .J (n x n), .Ed (n x md; default 0),
%     .Cd (qd x n; default Cc), .Fd (qd x md; default 0).
%   SYS FIELDS (switched): .A a 1xN cell of Metzler drifts (or handles),
%     per-mode .Ec,.Cc,.Fc cells (or shared matrices broadcast to all modes),
%     optional .trans (K x 2 admissible ordered transitions [old new];
%     default all old ~= new).
%
%   OPTS FIELDS:
%     .norm    'L1' | 'Linf' (required).
%     .class   'impulsive' | 'switched'.
%     .dwell   'min' | 'range' (required here).
%     .Tmin,.Tmax  dwell bounds ('min' uses Tmin; 'range' uses both).
%     .deg     timer-polynomial degree of the certificate (default 4).
%     .eps     strict positivity floor of the certificate (default 1e-4).
%     .solver  SDP backend for SOSTOOLS (default 'sedumi'), via MVST_OPTS.
%
%   RES FIELDS:
%     .gamma   the certified induced gain (min gamma), [] if infeasible;
%     .weights the certificate: a function handle @(tau)->vector (lambda or xi),
%              or a 1xN cell of such handles for 'switched';
%     .feas    logical feasibility (finite gamma returned);
%     .norm, .class, .dwell, .info (raw SOS solution / degrees).
%
%   See also MVST_POS_GAIN, MVST_SOSPOSMATRVAR, MVST_SOSMATRVAR, MVST_POS_LP.

    if nargin < 2; opts = struct(); end
    oo   = mvst_opts(opts);
    nrm  = lower(mvst_getf(opts, 'norm',  ''));
    cls  = lower(mvst_getf(opts, 'class', 'impulsive'));
    dwl  = lower(mvst_getf(opts, 'dwell', ''));
    deg  = mvst_getf(opts, 'deg', 4);
    eps0 = mvst_getf(opts, 'eps', 1e-4);
    if ~any(strcmp(nrm, {'l1','linf'}))
        error('mvst:pos:norm', 'OPTS.norm must be ''L1'' or ''Linf''.');
    end
    if ~any(strcmp(dwl, {'min','range','constant'}))
        error('mvst:pos:dwell', ...
              'MVST_POS_GAIN_SOS handles dwell ''constant'', ''min'' or ''range''.');
    end

    Tmin = mvst_getf(opts, 'Tmin', mvst_getf(sys, 'Tmin', []));
    Tmax = mvst_getf(opts, 'Tmax', mvst_getf(sys, 'Tmax', []));
    if strcmp(dwl,'constant')
        % exact dwell-time Tbar: timer certificate on [0,Tbar], single jump at
        % theta = Tbar, NO frozen tail (book eq:pos:traceper).
        Tbar = mvst_getf(opts,'Tbar', mvst_getf(sys,'Tbar', Tmin));
        if isempty(Tbar); error('mvst:pos:Tbar','dwell ''constant'' needs Tbar.'); end
        T = Tbar; Tmin = Tbar; Tmax = Tbar;
    elseif strcmp(dwl,'min')
        if isempty(Tmin); error('mvst:pos:Tmin','dwell ''min'' needs Tmin.'); end
        T = Tmin;
    else
        if isempty(Tmin) || isempty(Tmax)
            error('mvst:pos:Trange','dwell ''range'' needs Tmin and Tmax.');
        end
        T = Tmax;
    end

    res = struct('gamma', [], 'weights', [], 'feas', false, ...
                 'norm', nrm, 'class', cls, 'dwell', dwl, 'info', []);

    pvar th;
    gv   = dpvar('gam');
    prog = sosprogram(th);
    prog = sosdecvar(prog, gv);
    ctx  = struct('th', th, 'gv', gv, 'deg', deg, 'eps', eps0, ...
                  'T', T, 'Tmin', Tmin, 'Tmax', Tmax, 'dwl', dwl, 'nrm', nrm);

    switch cls
        case 'impulsive'
            [prog, wfun] = build_impulsive(prog, sys, ctx);
            weights = wfun;
        case 'switched'
            [prog, wfun] = build_switched(prog, sys, ctx);
            weights = wfun;
        otherwise
            error('mvst:pos:class', ...
                  'MVST_POS_GAIN_SOS class must be ''impulsive'' or ''switched''.');
    end

    prog = sossetobj(prog, gv);
    so = mvst_sosopts(oo);
    prog = mvst_sossolve(prog, so);

    gval = double(sosgetsol(prog, gv));
    res.info = struct('deg', deg, 'T', T, 'solver', oo.solver);
    if isfinite(gval)
        res.feas    = true;
        res.gamma   = gval;
        res.weights = weights;
    end
end

% ======================================================================= %
%  Impulsive builder                                                       %
% ======================================================================= %
function [prog, wfun] = build_impulsive(prog, sys, ctx)
    th = ctx.th;
    A  = mvst_getf(sys, 'A', []);
    Ec = mvst_getf(sys, 'Ec', mvst_getf(sys, 'E', []));
    Cc = mvst_getf(sys, 'Cc', mvst_getf(sys, 'C', []));
    if isempty(A) || isempty(Ec) || isempty(Cc)
        error('mvst:pos:gaindata', 'SYS needs .A, .Ec/.E and .Cc/.C.');
    end
    n  = matdim(A); mc = size(Ec,2); qc = size(Cc,1);
    Fc = mvst_getf(sys, 'Fc', mvst_getf(sys, 'F', zeros(qc, mc)));
    J  = mvst_getf(sys, 'J', zeros(n));
    Ed = mvst_getf(sys, 'Ed', zeros(n, 0));
    Cd = mvst_getf(sys, 'Cd', Cc);
    md = size(Ed,2); qd = size(Cd,1);
    Fd = mvst_getf(sys, 'Fd', zeros(qd, md));

    Ath  = evalM(A,  ctx.T*th);      % A(tau) at tau = T*th (normalized timer)
    AT   = evalM(A,  ctx.Tmin);      % frozen flow at Tmin (min-dwell tail)

    mons = monomials(th, 0:ctx.deg);
    [prog, w] = mvst_sosmatrvar(prog, mons, [n 1]);   % lambda or xi
    w0 = subs(w, th, 0);

    prog = gain_rows(prog, ctx, w, w0, Ath, AT, Ec, Cc, Fc, J, Ed, Cd, Fd, n);
    wfun = @(tau) double(subs(w, th, tau/ctx.T));
end

% ======================================================================= %
%  Switched builder (state continuous at a switch: identity map, no Ed)    %
% ======================================================================= %
function [prog, wfun] = build_switched(prog, sys, ctx)
    th = ctx.th;
    Acell = mvst_getf(sys, 'A', []);
    if ~iscell(Acell); error('mvst:pos:switched','switched needs a cell .A.'); end
    N = numel(Acell); n = matdim(Acell{1});
    Ec = ascellN(mvst_getf(sys,'Ec',mvst_getf(sys,'E',[])), N);
    Cc = ascellN(mvst_getf(sys,'Cc',mvst_getf(sys,'C',[])), N);
    Fc = ascellN(mvst_getf(sys,'Fc',mvst_getf(sys,'F',[])), N);
    trans = mvst_getf(sys, 'trans', all_trans(N));

    mons = monomials(th, 0:ctx.deg);
    w = cell(1,N); w0 = cell(1,N);
    for i = 1:N
        [prog, w{i}] = mvst_sosmatrvar(prog, mons, [n 1]);
        w0{i} = subs(w{i}, th, 0);
    end
    % Per-mode flow / output / positivity / frozen tail (no per-mode jump self).
    for i = 1:N
        Eci = Ec{i}; Cci = Cc{i};
        qc = size(Cci,1); mc = size(Eci,2);
        if isempty(Fc{i}); Fci = zeros(qc,mc); else; Fci = Fc{i}; end
        Ath = evalM(Acell{i}, ctx.T*th);  AT = evalM(Acell{i}, ctx.Tmin);
        prog = flow_output_rows(prog, ctx, w{i}, Ath, AT, Eci, Cci, Fci, n);
    end
    % Switch coupling old -> new: book jump specialized to identity map / no Ed.
    for k = 1:size(trans,1)
        io = trans(k,1); jn = trans(k,2);          % old -> new
        prog = switch_rows(prog, ctx, w{io}, w0{jn}, n);
    end
    wfun = cellfun(@(p) @(tau) double(subs(p, th, tau/ctx.T)), w, 'UniformOutput', false);
end

% ======================================================================= %
%  Row builders                                                            %
% ======================================================================= %
function prog = gain_rows(prog, ctx, w, w0, Ath, AT, Ec, Cc, Fc, J, Ed, Cd, Fd, n)
% Full impulsive gain rows (flow + output + positivity + jump + tail).
    md = size(Ed,2); qd = size(Cd,1);
    gv = ctx.gv;
    prog = flow_output_rows(prog, ctx, w, Ath, AT, Ec, Cc, Fc, n);
    % ---- jump rows (impulsive: at theta) -----------------------------------
    if strcmp(ctx.nrm, 'linf')
        wj = jumpval(ctx, w);                          % xi(theta)
        ej  = J*wj + Ec_agg(Ed,md) - w0;               % J xi(theta) + Ed1 - xi(0) <= 0
        eoj = Cd*wj + Ec_agg(Fd,md) - gv*ones(qd,1);   % Cd xi(theta) + Fd1 <= gamma1
        prog = jump_ineq(prog, ctx, -ej,  n);
        prog = jump_ineq(prog, ctx, -eoj, qd);
    else
        wj = jumpval(ctx, w);                          % lambda(theta)
        ej  = J.'*w0 - wj + Cc_agg(Cd);                % J' lam(0) - lam(theta) + Cd'1 <= 0
        eoj = Ed.'*w0 + Cc_agg(Fd) - gv*ones(md,1);    % Ed' lam(0) + Fd'1 <= gamma1
        prog = jump_ineq(prog, ctx, -ej,  n);
        prog = jump_ineq(prog, ctx, -eoj, md);
    end
end

function prog = flow_output_rows(prog, ctx, w, Ath, AT, Ec, Cc, Fc, n)
% Flow, output-flow, positivity on (0,T), and the minimum-dwell frozen tail.
    th = ctx.th; T = ctx.T; gv = ctx.gv;
    mc = size(Ec,2); qc = size(Cc,1);
    wdot = diff(w, th)/T;        % d/dtau on normalized timer th = tau/T
    if strcmp(ctx.nrm, 'linf')
        ef  = Ath*w + Ec_agg(Ec,mc) - wdot;            % A xi + Ec1 - dxi <= 0
        eof = Cc*w + Ec_agg(Fc,mc) - gv*ones(qc,1);    % Cc xi + Fc1 <= gamma1
    else
        ef  = wdot + Ath.'*w + Cc_agg(Cc);             % dlam + A'lam + Cc'1 <= 0
        eof = Ec.'*w + Cc_agg(Fc) - gv*ones(mc,1);     % Ec'lam + Fc'1 <= gamma1
    end
    prog = interval_ineq(prog, ctx, -ef,  n);          % flow
    prog = interval_ineq(prog, ctx, -eof, numel_rows(eof));  % output flow
    prog = interval_ineq(prog, ctx, w - ctx.eps*ones(n,1), n);  % positivity
    % ---- frozen tail (minimum dwell only) ----------------------------------
    if strcmp(ctx.dwl, 'min')
        wT = subs(w, th, 1);        % tau = Tmin (= T) at th = 1
        if strcmp(ctx.nrm, 'linf')
            et = AT*wT + Ec_agg(Ec,mc);                % A(Tmin) xi(Tmin) + Ec1 <= 0
        else
            et = AT.'*wT + Cc_agg(Cc);                 % A(Tmin)'lam(Tmin) + Cc'1 <= 0
        end
        prog = point_ineq(prog, -et, n);
    end
end

function prog = switch_rows(prog, ctx, wold, w0new, n)
% Switched jump (book row, identity map, no disturbance):
%   Linf:  xi_new(0) >= xi_old(theta)  ->  xi_old(theta) - xi_new(0) <= 0
%   L1  :  lambda_old(theta) - lambda_new(0) <= 0  (adjoint continuity)
    wj = jumpval(ctx, wold);
    e  = wj - w0new;
    prog = jump_ineq(prog, ctx, -e, n);
end

% ======================================================================= %
%  Positivstellensatz imposition                                          %
% ======================================================================= %
function prog = interval_ineq(prog, ctx, expr, m)
% expr(tau) >= 0 (each of the m scalar rows) on tau in [0,T], via the interval
% multiplier g = tau (T - tau) and a fresh SOS multiplier per row.
    th = ctx.th; g = th*(1 - th);          % normalized timer interval [0,1]
    for i = 1:m
        [prog, S] = mvst_sosposmatrvar(prog, 1, max(0,ctx.deg-2), th);
        prog = sosineq(prog, expr(i) - g*S);
    end
end

function prog = jump_ineq(prog, ctx, expr, m)
% For 'min' the jump is a single endpoint (theta = Tmin): a plain scalar >= 0.
% For 'range' the jump holds for all theta in [Tmin,Tmax]: interval multiplier
% gm = (tau - Tmin)(Tmax - tau).
    if ~strcmp(ctx.dwl, 'range')
        prog = point_ineq(prog, expr, m);           % constant/min: single endpoint
    else
        th = ctx.th; gm = (th - ctx.Tmin/ctx.T)*(1 - th);   % theta in [Tmin/T,1]
        for i = 1:m
            [prog, S] = mvst_sosposmatrvar(prog, 1, max(0,ctx.deg-2), th);
            prog = sosineq(prog, expr(i) - gm*S);
        end
    end
end

function prog = point_ineq(prog, expr, m)
% expr >= 0 at a fixed timer value (constant / decision-variable-affine).
    for i = 1:m; prog = sosineq(prog, expr(i)); end
end

function v = jumpval(ctx, w)
% Pre-jump certificate value at the dwell instant theta: for 'min' the fixed
% endpoint w(Tmin); for 'range' the polynomial w(tau) (imposed over [Tmin,Tmax]).
    if strcmp(ctx.dwl, 'range'); v = w; else; v = subs(w, ctx.th, 1); end
end

% ======================================================================= %
%  Small utilities                                                         %
% ======================================================================= %
function y = Ec_agg(M, m)                 % M * ones(m,1) (aggregate disturbance)
    if isempty(M) || m == 0; y = zeros(size(M,1),1); else; y = M*ones(m,1); end
end
function y = Cc_agg(M)                     % M' * ones(rows(M),1) (adjoint aggregate)
    if isempty(M); y = 0; else; y = M.'*ones(size(M,1),1); end
end
function d = matdim(A); if iscell(A); d = size(A{1},1); else; d = size(evalM(A,0),1); end; end
function M = evalM(A, tau)                 % constant matrix or handle @(tau)
    if isa(A,'function_handle'); M = A(tau); else; M = A; end
end
function c = ascellN(X, N)
    if iscell(X); c = X; else; c = repmat({X}, 1, N); end
end
function tr = all_trans(N)
    tr = zeros(N*(N-1),2); k = 0;
    for i = 1:N; for j = 1:N; if i~=j; k=k+1; tr(k,:)=[i j]; end; end; end
end
function r = numel_rows(e); r = size(e,1); end
