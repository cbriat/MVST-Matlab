function res = mvst_pos_stabilize_sos(sys, opts)
%MVST_POS_STABILIZE_SOS  Timer-dependent copositive-SOS state-feedback synthesis
%   for positive impulsive systems under constant / minimum / range dwell-time,
%   the timer-dependent branch of MVST_POS_STABILIZE.  It designs a
%   flow gain Kc(tau) and a jump gain Kd (or Kd(theta)) rendering the closed
%   loop BOTH positive and exponentially stable, via the copositive change of
%   variables of Lemma lem:pos:cov (X = diag(lambda), Uc = Kc X, Ud = Kd X), so
%   the design conditions are linear in (lambda, Uc, Ud) at each timer value.
%
%   Two objectives:
%     * STABILIZATION (default): the timer-dependent conditions of
%       Thm th:pos:stabzfix(b) (constant) and Prop. prop:pos:stabzdwell
%       (minimum eq:pos:stabzMDT / range eq:pos:stabzRDT).  The routine
%       maximizes the copositive decrease slack and certifies feasibility.
%     * GAIN-OPTIMAL (when OPTS.norm = 'Linf' and disturbance/output data are
%       present): the SAME positivity change of variables together with the
%       closed-loop peak-gain rows of Cor. cor:pos:Linfslice (forward tube
%       xi = lambda), minimizing the induced Linf x ellinf gain gamma.
%
%   ROWS (exactly as displayed; forward timer tau = t - t_k of eq:pos:timer):
%     positivity  A(tau) X + Bc Uc + alpha I >= 0        (entrywise, on the
%                 dwell interval),  J X(theta) + Bd Ud >= 0  (entrywise);
%     decrease    [-dX/dtau + A X + Bc Uc] 1 <= -eps 1   (stabilization)
%                 or  [-dX/dtau + (A X + Bc Uc) 1] + Ec 1 <= 0  (gain, the tube
%                 flow of the closed loop);
%     jump        [J X(theta) + Bd Ud - X(0)] 1 <= -eps 1  (stabilization)
%                 or  <= -Ed 1  (gain);
%     output      (Cc X + Dc Uc) 1 + Fc 1 <= gamma 1,  (Cd X(theta) + Dd Ud) 1
%                 + Fd 1 <= gamma 1   (gain only);
%     tail        minimum dwell-time adds the frozen tail
%                 [A(Tmin) X(Tmin) + Bc Uc(Tmin)] 1 <= (-eps 1 | -Ec 1)
%                 (eq:pos:szm:tail).
%   No differentiability pin (book form).  Range dwell-time imposes the jump
%   rows for all theta in [Tmin,Tmax] via the interval multiplier, with a
%   timer-dependent jump gain Ud(theta).
%
%   OPTIONAL script-side regularization (default OFF; NOT a book condition):
%   OPTS.Kmax caps the recovered flow gain magnitude by the copositive rows
%   Kmax X - Uc >= 0 and Uc + Kmax X >= 0 (the NAHS C5/C6 rows); documented as
%   a numerical regularizer only.
%
%   SYS FIELDS: .A (Metzler; n x n or a handle @(tau)), .Bc/.B (n x mc),
%     .J (n x n), .Bd (n x md; default Bc); gain mode also reads
%     .Ec/.E,.Cc/.C,.Fc/.F,.Ed,.Cd,.Fd,.Dc,.Dd.
%   OPTS FIELDS: .dwell 'constant'|'min'|'range'; .Tbar/.Tmin/.Tmax; .norm
%     ('Linf' selects gain-optimal); .deg (timer degree, default 2); .eps
%     (positivity floor, default 1e-4); .margin (decrease slack, stabilization);
%     .Kmax (optional gain cap, default []); .solver.
%
%   RES FIELDS: .feas; .gamma (gain mode); .slack (stabilization); .Kc @(tau),
%     .Kd (md x n, or @(theta) for range); .lambda @(tau); .Uc,.Ud (raw
%     polynomials); .class,.dwell,.info.
%
%   See also MVST_POS_STABILIZE, MVST_POS_GAIN_SOS, MVST_SOSMATRVAR,
%   MVST_SOSPOSMATRVAR.

    if nargin < 2; opts = struct(); end
    oo   = mvst_opts(opts);
    dwl  = lower(mvst_getf(opts,'dwell',''));
    deg  = mvst_getf(opts,'deg', 2);
    eps0 = mvst_getf(opts,'eps', 1e-4);
    mrg  = mvst_getf(opts,'margin', 1e-3);
    Kmax = mvst_getf(opts,'Kmax', []);
    if ~any(strcmp(dwl,{'constant','min','range'}))
        error('mvst:pos:dwell','dwell must be ''constant'',''min'' or ''range''.');
    end

    A  = mvst_getf(sys,'A',[]);
    Bc = mvst_getf(sys,'Bc', mvst_getf(sys,'B',[]));
    if isempty(A) || isempty(Bc); error('mvst:pos:syn','SYS needs .A and .Bc/.B.'); end
    n  = matdim(A); mc = size(Bc,2);
    J  = mvst_getf(sys,'J', zeros(n));
    Bd = mvst_getf(sys,'Bd', Bc); md = size(Bd,2);

    % gain-optimal mode iff Linf norm requested AND disturbance/output present
    Ec = mvst_getf(sys,'Ec', mvst_getf(sys,'E',[]));
    Cc = mvst_getf(sys,'Cc', mvst_getf(sys,'C',[]));
    gainmode = strcmpi(mvst_getf(opts,'norm',''),'Linf') && ~isempty(Ec) && ~isempty(Cc);

    Tmin = mvst_getf(opts,'Tmin', mvst_getf(sys,'Tmin',[]));
    Tmax = mvst_getf(opts,'Tmax', mvst_getf(sys,'Tmax',[]));
    switch dwl
        case 'constant'
            T = mvst_getf(opts,'Tbar', mvst_getf(sys,'Tbar', Tmin));
            if isempty(T); error('mvst:pos:Tbar','dwell ''constant'' needs Tbar.'); end
        case 'min'
            if isempty(Tmin); error('mvst:pos:Tmin','dwell ''min'' needs Tmin.'); end
            T = Tmin;
        case 'range'
            if isempty(Tmin)||isempty(Tmax); error('mvst:pos:Trange','range needs Tmin,Tmax.'); end
            T = Tmax;
    end

    if gainmode
        [Fc,Ed,Cd,Fd,Dc,Dd,qc,qd,~] = gain_data(sys, Ec, Cc, n, mc, md);
    end

    pvar th;
    prog = sosprogram(th);
    mons  = monomials(th, 0:deg);
    [prog, lam] = mvst_sosmatrvar(prog, mons, [n 1]);
    [prog, Uc ] = mvst_sosmatrvar(prog, mons, [mc n]);
    if strcmp(dwl,'range')
        [prog, Ud] = mvst_sosmatrvar(prog, mons, [md n]);   % Ud(theta)
    else
        [prog, Ud] = mvst_sosmatrvar(prog, monomials(th,0), [md n]);  % constant Ud
    end
    [prog, al] = mvst_sosposmatrvar(prog, 1, 0, th);        % alpha >= 0
    if gainmode
        gv = dpvar('gam'); prog = sosdecvar(prog, gv);
    else
        sl = dpvar('sl');  prog = sosdecvar(prog, sl);       % decrease slack
        prog = sosineq(prog, 1 - sl);                        % slack <= 1: scaling a feasible (lambda,U,sl) by c>1 stays feasible (monotonicity), so the maximized slack is unbounded without this cap
    end

    onen = ones(n,1); onemc = ones(mc,1); onemd = ones(md,1);
    g  = th*(1-th);
    lamdot = diff(lam,th)/T;
    lam0 = subs(lam,th,0); lamT = subs(lam,th,1); UcT = subs(Uc,th,1);
    Ath  = evalM(A, T*th);  AT = evalM(A, Tmin);
    AX = @(L,M) M .* (onen*L.');                 % (M diag(L))_{ij} = M_ij L_j

    % ---- positivity (entrywise) on the dwell interval --------------------- %
    Mmetz = AX(lam,Ath) + Bc*Uc + al*eye(n);
    prog  = intv_mat(prog, Mmetz, g, th, deg);
    % ---- flow / decrease on the dwell interval ---------------------------- %
    AXBU  = AX(lam,Ath) + Bc*Uc;
    if gainmode
        flow = -lamdot + AXBU*onen + Ec*onemc;                    % <= 0
        prog = intv_vec(prog, -flow, g, th, deg, n);
        outf = Cc*lam + Dc*(Uc*onen) + Fc*onemc - gv*ones(qc,1);  % <= 0
        prog = intv_vec(prog, -outf, g, th, deg, qc);
    else
        flow = -lamdot + AXBU*onen + sl*onen;                     % <= 0 (slack)
        prog = intv_vec(prog, -flow, g, th, deg, n);
    end
    % ---- positivity of the certificate ------------------------------------ %
    prog = intv_vec(prog, lam - eps0*onen, g, th, deg, n);
    % ---- optional gain cap (script-side regularizer, default off) --------- %
    if ~isempty(Kmax)
        Xrow = onemc*lam.';                        % broadcast lambda over mc rows
        prog = intv_mat(prog, Kmax*Xrow - Uc, g, th, deg);
        prog = intv_mat(prog, Uc + Kmax*Xrow, g, th, deg);
    end

    % ---- jump rows (theta = T for constant/min; [Tmin,Tmax] for range) ---- %
    if strcmp(dwl,'range')
        Mj  = AX(lam,J) + Bd*Ud;                  % J X(theta) + Bd Ud(theta)
        prog = jmp_mat(prog, Mj, th, Tmin, T, deg);
        if gainmode
            jt = Mj*onen - lam0 + Ed*onemd;       % <= 0
            prog = jmp_vec(prog, -jt, th, Tmin, T, deg, n);
            oj = Cd*lam + Dd*(Ud*onen) + Fd*onemd - gv*ones(qd,1);
            prog = jmp_vec(prog, -oj, th, Tmin, T, deg, qd);
        else
            jt = Mj*onen - lam0 + sl*onen;        % <= 0 (slack)
            prog = jmp_vec(prog, -jt, th, Tmin, T, deg, n);
        end
    else
        Mj  = AX(lamT,J) + Bd*Ud;                 % J X(T) + Bd Ud
        for i=1:n; for j=1:n; prog = siq(prog, Mj(i,j)); end; end
        if gainmode
            jt = Mj*onen - lam0 + Ed*onemd;  for i=1:n; prog=siq(prog,-jt(i)); end
            oj = Cd*lamT + Dd*(Ud*onen) + Fd*onemd - gv*ones(qd,1);
            for i=1:qd; prog=siq(prog,-oj(i)); end
        else
            jt = Mj*onen - lam0 + sl*onen;   for i=1:n; prog=siq(prog,-jt(i)); end
        end
    end

    % ---- minimum-dwell frozen tail (eq:pos:szm:tail) ---------------------- %
    if strcmp(dwl,'min')
        tail = AX(lamT,AT)*onen + Bc*(UcT*onen);
        if gainmode; tail = tail + Ec*onemc; rhs = 0;
        else;        tail = tail + sl*onen;  rhs = 0; end %#ok<NASGU>
        for i=1:n; prog = siq(prog, -tail(i)); end
    end

    % ---- objective -------------------------------------------------------- %
    if gainmode; prog = sossetobj(prog, gv); else; prog = sossetobj(prog, -sl); end
    so = mvst_sosopts(oo);
    [prog, sinfo] = mvst_sossolve(prog, so);
    % Clean-solve test (mvst_kyp convention): numerr < 2 on every backend. On nosdp, numerr == 2
    % (uncertified_*, stalled, max_iter) still returns the design but labeled uncertified.
    isno = strcmpi(char(oo.solver),'nosdp');
    solved = (sinfo.pinf == 0 && sinfo.dinf == 0);
    certified = (sinfo.numerr < 2);
    if solved && ~certified
        if ~isno; solved = false; end
    end

    lamv = sosgetsol(prog, lam); Ucv = sosgetsol(prog, Uc); Udv = sosgetsol(prog, Ud);
    res = struct('feas',false,'gamma',[],'slack',[],'Kc',[],'Kd',[], ...
                 'lambda',[],'Uc',Ucv,'Ud',Udv,'class','impulsive','dwell',dwl,'info',[]);
    res.info = struct('deg',deg,'T',T,'solver',oo.solver,'certified',certified);
    if solved && ~certified
        warning('mvst:uncertified', 'mvst_pos_stabilize_sos: the nosdp optimum is uncertified (status %s); the design is returned but info.certified = false.', sinfo.nosdp_status);
    end
    res.lambda = @(tau) double(subs(lamv, th, tau/T));
    res.Kc = @(tau) gainrow(Ucv, lamv, th, tau/T, mc);
    if strcmp(dwl,'range'); res.Kd = @(theta) gainrow(Udv, lamv, th, theta/T, md);
    else;                   res.Kd = gainrow(Udv, lamv, th, 1, md); end
    if gainmode
        gval = double(sosgetsol(prog, gv));
        if solved && isfinite(gval); res.feas = true; res.gamma = gval; end
    else
        sval = double(sosgetsol(prog, sl));
        res.slack = sval; res.feas = solved && sval > 1e-7; %#ok<*PROP>
    end
end

% ======================================================================= %
function [Fc,Ed,Cd,Fd,Dc,Dd,qc,qd,mc] = gain_data(sys, Ec, Cc, n, mc, md)
    qc = size(Cc,1); mcx = size(Ec,2);
    Fc = mvst_getf(sys,'Fc', mvst_getf(sys,'F', zeros(qc,mcx)));
    Ed = mvst_getf(sys,'Ed', zeros(n,md));
    Cd = mvst_getf(sys,'Cd', Cc); qd = size(Cd,1);
    Fd = mvst_getf(sys,'Fd', zeros(qd, size(Ed,2)));
    Dc = mvst_getf(sys,'Dc', zeros(qc, mc));
    Dd = mvst_getf(sys,'Dd', zeros(qd, md));
    mc = mcx;
end

function K = gainrow(U, lam, th, thv, m)
% Kc(tau) = U(tau) * diag(lambda(tau))^{-1}  (rational, stored as the U,X pair;
% this evaluator returns the numeric gain at a timer value).
    Uv = double(subs(U, th, thv));  Uv = reshape(Uv, m, []);
    lv = double(subs(lam, th, thv)); lv = lv(:).';
    K = Uv ./ (ones(m,1)*lv);
end

function prog = intv_mat(prog, M, g, th, deg, varargin)
    [p,q] = size(M);
    for i=1:p; for j=1:q
        [prog,S] = mvst_sosposmatrvar(prog, 1, max(0,deg), th);
        prog = sosineq(prog, M(i,j) - g*S);
    end; end
end
function prog = intv_vec(prog, v, g, th, deg, m)
    for i=1:m
        [prog,S] = mvst_sosposmatrvar(prog, 1, max(0,deg), th);
        prog = sosineq(prog, v(i) - g*S);
    end
end
function prog = jmp_mat(prog, M, th, Tmin, T, deg)
    gm = (th - Tmin/T)*(1 - th); [p,q] = size(M);
    for i=1:p; for j=1:q
        [prog,S] = mvst_sosposmatrvar(prog, 1, max(0,deg), th);
        prog = sosineq(prog, M(i,j) - gm*S);
    end; end
end
function prog = jmp_vec(prog, v, th, Tmin, T, deg, m)
    gm = (th - Tmin/T)*(1 - th);
    for i=1:m
        [prog,S] = mvst_sosposmatrvar(prog, 1, max(0,deg), th);
        prog = sosineq(prog, v(i) - gm*S);
    end
end
function prog = siq(prog, e)
% Guarded scalar SOS inequality: SOSTOOLS 4.00 builds a DEGENERATE SDP block
% for a symbolically-zero argument (e.g. a structurally-zero output column),
% so skip a row that is PROVABLY the constant 0 (a vacuous 0 >= 0); impose any
% row carrying a polynomial or decision variable.
    if isZeroConst(e); return; end
    prog = sosineq(prog, e);
end
function tf = isZeroConst(e)
    try d = double(e); tf = all(d(:)==0); catch; tf = false; end
end
function d = matdim(A); if isa(A,'function_handle'); d=size(A(0),1); else; d=size(A,1); end; end
function M = evalM(A, tau); if isa(A,'function_handle'); M=A(tau); else; M=A; end; end
