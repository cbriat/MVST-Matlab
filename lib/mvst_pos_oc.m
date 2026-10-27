function res = mvst_pos_oc(sys, opts)
%MVST_POS_OC  Linear-cost optimal control of linear positive systems by the
%   copositive costate programs of Chapter "Linear Positive Systems"
%   (Section sec:pos:oc).  A single dispatcher over the value/costate results:
%   the value of a positive system under a LINEAR cost is itself linear,
%   V(x) = p'x, and the optimal costate p solves a copositive program.
%
%   RES = MVST_POS_OC(SYS, OPTS) selects a result with OPTS.mode.
%
%   The cost (weights entrywise nonnegative, q,q_d,q_T >= 0) is
%     J_T = \int_{t0}^{t0+T}(q'x + r_c'u_c)dt + \sum_k(q_d'x + r_d'u_d) + q_T'x(T)
%   (eq:pos:occost / eq:pos:occostu); r_c,r_d may be sign-indefinite.
%
%   MODES (OPTS.mode):
%     'finite'  Finite-horizon, no control penalty (Prop. prop:pos:ocfinite):
%               optimal cost J* = p(t0)'x1_0 with the adjoint costate recursion
%                 -dp/dt = A(t)'p + q,   p(t_k^-) = J'p(t_k^+) + q_d,  p(T)=q_T
%               (eq:pos:ocLP with the input terms Xi = Xi_d = 0).  Jump times in
%               OPTS.tk (default none = pure continuous-time).
%     'inf'     Infinite-horizon LTI, no penalty (Prop. prop:pos:ocinf): the
%               bounded costate is p = -(A')^{-1} q (A Metzler-Hurwitz), giving
%               J* = p'x1_0.
%     'dwell'   Periodic dwell-time T_p (Prop. prop:pos:ocdwell(a)): the unique
%               T_p-periodic costate p(tau) on [0,T_p] with the flow recursion
%               and the jump p(0^+) = J'p(T_p) + q_d; J* = p(0)'x1_0.
%     'bangbang' Constrained inputs u_c in U_c, u_d in U_d (Prop.
%               prop:pos:ocbangbang, corrected revpurple forms): the adjoint
%               costate (eq:pos:ocadjoint), switching functions
%                 sigma_c(t) = r_c(t) + B_c'p(t),  sigma_d(k) = r_d + B_d'p(t_k^+)
%               (eq:pos:ocswitch), the VALUE IDENTITY eq:pos:ocvalueid and the
%               optimal cost eq:pos:ocvalueopt
%                 J* = p(t0)'x0 + \int min_{u in U_c} sigma_c'u dt
%                             + \sum_k min_{u in U_d} sigma_d'u,
%               with the box bang-bang rule eq:pos:ocbangrule.
%     'stateconstr' State-scaled amplitude constraint |u_c| <= E_c x
%               (Prop. prop:pos:ocstateconstr, corrected revpurple forms): the
%               value is linear V(t,x) = ptil(t)'x with the BELLMAN COSTATE
%               eq:pos:ocexcostate
%                 -dptil/dt = A'ptil + q - E_c'|sigma_c(ptil)|,  ptil(T)=q_T,
%               sigma_c(ptil) = r_c + B_c'ptil, and the linear saturating
%               feedback rows eq:pos:ocKrow
%                 (K_c)_i = sign(sigma_{c,i}(ptil)) (E_c)_i,
%               which satisfies |K_c| <= E_c (eq:pos:ocKbound).
%
%   SYS FIELDS: .A (n x n Metzler; handle @(t) allowed for 'finite'/'bangbang'),
%     .J (reset; default none), .Bc/.B (n x mc), .Bd (n x md; default Bc),
%     .Ec (n x mc, 'stateconstr' amplitude scaling).
%   OPTS FIELDS: .mode; cost .q (n x 1 or @(t)), .qd (n x 1), .qT (n x 1),
%     .rc (mc x 1 or @(t)), .rd (md x 1); .x0/.X10 (n x 1) initial state;
%     .T (horizon) / .Tp (period); .tk (jump times, row); input boxes
%     .ubc,.ubd (upper bounds, lower bound 0) for 'bangbang'; .ngrid (time grid,
%     default 400); .solver.
%
%   RES FIELDS: .Jopt (optimal cost p(t0)'x0), .p (@(t)->costate) or .ptil,
%     .p0 (costate at t0), .sigmac/.sigmad (@(t)/vector switching functions,
%     bangbang), .Kc (@(t)-> saturating gain rows, stateconstr), .feas, .mode.
%
%   See also MVST_POS_GAIN, MVST_POS_STABILIZE, MVST_POS_STABILITY, MVST_POS.

    arguments
        sys  (1,1) struct
        opts (1,1) struct = struct()
    end

    mode = lower(mvst_getf(opts,'mode',''));
    if isempty(mode); error('mvst:pos:ocmode','OPTS.mode is required.'); end
    x0 = mvst_getf(opts,'x0', mvst_getf(opts,'X10', []));

    res = struct('Jopt',[],'p',[],'p0',[],'ptil',[],'sigmac',[],'sigmad',[], ...
                 'Kc',[],'feas',false,'mode',mode);

    switch mode
        case 'inf'
            A = getf(sys,'A'); q = numvec(mvst_getf(opts,'q',[]), size(A,1));
            if max(real(eig(A))) >= 0
                error('mvst:pos:ocinf','''inf'' needs a Metzler-Hurwitz A.');
            end
            p = -(A.')\q;                       % bounded costate (eq:pos:ocinfLP)
            res.p0 = p; res.feas = all(p >= -1e-9);
            res.p = @(t) p;
            if ~isempty(x0); res.Jopt = p.'*x0; end

        case 'finite'
            [p0, pfun] = costate_backward(sys, opts);
            res.p0 = p0; res.p = pfun; res.feas = all(p0 >= -1e-9);
            if ~isempty(x0); res.Jopt = p0.'*x0; end

        case 'dwell'
            [p0, pfun] = costate_periodic(sys, opts);
            res.p0 = p0; res.p = pfun; res.feas = all(p0 >= -1e-9);
            if ~isempty(x0); res.Jopt = p0.'*x0; end

        case 'bangbang'
            res = oc_bangbang(sys, opts, res, x0);

        case 'stateconstr'
            res = oc_stateconstr(sys, opts, res, x0);

        otherwise
            error('mvst:pos:ocmode', ...
                  'OPTS.mode must be finite|inf|dwell|bangbang|stateconstr.');
    end
end

% ======================================================================= %
%  Costate backward recursion (no control penalty): -dp = A'p + q         %
% ======================================================================= %
function [p0, pfun] = costate_backward(sys, opts)
    A  = getf(sys,'A'); n = matdim(A);
    q  = mvst_getf(opts,'q', zeros(n,1));
    qd = numvec(mvst_getf(opts,'qd', zeros(n,1)), n);
    qT = numvec(mvst_getf(opts,'qT', zeros(n,1)), n);
    J  = mvst_getf(sys,'J', []);
    T  = mvst_getf(opts,'T', []);
    if isempty(T); error('mvst:pos:ocT','''finite'' needs a horizon OPTS.T.'); end
    tk = sort(mvst_getf(opts,'tk', []));
    tk = tk(tk > 0 & tk < T);
    ng = mvst_getf(opts,'ngrid', 400);
    tgrid = linspace(0, T, ng);
    events = [tk(:).', T];
    p = qT; tsamp = T; psamp = qT;               % backward accumulation
    tprev = T;
    for e = numel(events):-1:1
        tr = events(e);
        if e == 1; tl = 0; else; tl = events(e-1); end
        % integrate backward on (tl, tr): sample within the interval
        seg = flip(tgrid(tgrid >= tl - 1e-12 & tgrid <= tr + 1e-12));
        for tt = seg
            dt = tr - tt;
            pt = costate_flow(A, q, p, dt, tr);   % p(tt) from p(tr)
            tsamp(end+1) = tt; psamp(:,end+1) = pt; %#ok<AGROW>
        end
        % move to the left endpoint and apply the jump p(tl^-) = J'p(tl^+)+qd
        p = costate_flow(A, q, p, tr - tl, tr);
        if e > 1 && ~isempty(J)
            p = J.'*p + qd;
        end
    end
    p0 = p;
    [tsu, iu] = unique(tsamp); psu = psamp(:, iu);
    pfun = @(t) interp1(tsu, psu.', min(max(t,0),T), 'linear').';
end

function pt = costate_flow(A, q, pr, dt, tr)
% Solve -dp/dt = A'p + q backward over a length-dt interval ending at p(tr)=pr.
% For constant A this is the exact matrix-exponential map; a handle A(t) is
% frozen at the interval's right endpoint tr (piecewise-constant on the grid).
     At = evalM(A, tr); n = size(At,1);
    E  = expm(At.'*dt);
    pt = E*pr + (At.')\(E - eye(n))*numvec(q, n);
end

% ======================================================================= %
%  Periodic costate (Prop. prop:pos:ocdwell(a))                           %
% ======================================================================= %
function [p0, pfun] = costate_periodic(sys, opts)
    A = getf(sys,'A'); n = matdim(A);
    q = numvec(mvst_getf(opts,'q', zeros(n,1)), n);
    qd = numvec(mvst_getf(opts,'qd', zeros(n,1)), n);
    J = mvst_getf(sys,'J', eye(n));
    Tp = mvst_getf(opts,'Tp', mvst_getf(opts,'Tbar', mvst_getf(opts,'T',[])));
    if isempty(Tp); error('mvst:pos:ocTp','''dwell'' needs a period OPTS.Tp.'); end
    % One-period backward map: p(0^+) = M p(Tp) + c, with p(0^+) periodic ==
    % p(Tp) after applying the jump: p(Tp) = J' p(0^+) + q_d, and the flow map
    % p(0^+) = e^{A'Tp} p(Tp) + (A')^{-1}(e^{A'Tp}-I) q.  Solve the linear
    % periodicity system for the periodic costate.
    E = expm(A.'*Tp); c = (A.')\(E - eye(n))*q;
    % p(0^+) = E*(J'*p(0^+)+q_d) + c  =>  (I - E J') p(0^+) = E q_d + c
    p0plus = (eye(n) - E*J.') \ (E*qd + c);
    p0 = p0plus;
    ng = mvst_getf(opts,'ngrid', 200); tgrid = linspace(0,Tp,ng);
    P = zeros(n, ng);
    for i = 1:ng
        dt = Tp - tgrid(i);
        pTp = J.'*p0plus + qd;                   % p(Tp) pre-jump value
        P(:,i) = expm(A.'*dt)*pTp + (A.')\(expm(A.'*dt)-eye(n))*q;
    end
    pfun = @(t) interp1(tgrid, P.', min(max(t,0),Tp), 'linear').';
end

% ======================================================================= %
%  Bang-bang optimal control (corrected revpurple forms)                  %
% ======================================================================= %
function res = oc_bangbang(sys, opts, res, x0)
    A  = getf(sys,'A'); n = matdim(A);
    Bc = getf(sys,'Bc','B'); mc = size(Bc,2);
    T  = mvst_getf(opts,'T',[]);
    if isempty(T); error('mvst:pos:ocT','''bangbang'' needs OPTS.T.'); end
    [p0, pfun] = costate_backward(sys, opts);
    res.p0 = p0; res.p = pfun;
    rc = mvst_getf(opts,'rc', zeros(mc,1));
    ubc = numvec(mvst_getf(opts,'ubc', ones(mc,1)), mc);   % box [0, ubc]
    sigc = @(t) numvec(evalv(rc,t),mc) + evalM(Bc,t).'*pfun(t);   % eq:pos:ocswitch
    res.sigmac = sigc;
    % Optimal cost eq:pos:ocvalueopt: J* = p(t0)'x0 + \int min_{u in Uc} sig'u.
    ng = mvst_getf(opts,'ngrid', 400); tg = linspace(0,T,ng);
    integ = zeros(1,ng);
    for i = 1:ng
        s = sigc(tg(i));
        integ(i) = sum(min(0, s).*ubc);          % vertex min over the box
    end
    Jstar = trapz(tg, integ);
    if ~isempty(x0); res.Jopt = p0.'*x0 + Jstar; end
    res.feas = true;
    % Box bang-bang gain rule eq:pos:ocbangrule as a function of time.
    res.Kc = @(t) bangrule(sigc(t), ubc);
end

function u = bangrule(s, ub)
    u = zeros(size(s)); u(s < 0) = ub(s < 0);    % u_i = ub if sigma_i<0, else 0
end

% ======================================================================= %
%  State-constrained OC |u_c| <= E_c x (corrected revpurple forms)        %
% ======================================================================= %
function res = oc_stateconstr(sys, opts, res, x0)
    A  = getf(sys,'A'); n = matdim(A);
    Bc = getf(sys,'Bc','B'); mc = size(Bc,2);
    Ec = getf(sys,'Ec');
    T  = mvst_getf(opts,'T',[]);
    if isempty(T); error('mvst:pos:ocT','''stateconstr'' needs OPTS.T.'); end
    q  = numvec(mvst_getf(opts,'q', zeros(n,1)), n);
    qT = numvec(mvst_getf(opts,'qT', zeros(n,1)), n);
    rc = numvec(mvst_getf(opts,'rc', zeros(mc,1)), mc);
    % Bellman costate eq:pos:ocexcostate: -dptil = A'ptil + q - Ec'|sigma_c|,
    % sigma_c = rc + Bc'ptil.  Backward nonlinear ODE, ptil(T)=qT.
    % integrate backward: let s = T - t, dptil/ds = A'ptil + q - Ec'|sigma|
    fwd = @(s,pt) A.'*pt + q - Ec.'*abs(rc + Bc.'*pt);
    ng = mvst_getf(opts,'ngrid', 400); sg = linspace(0,T,ng);
    [~, P] = ode45(fwd, sg, qT);
    ptil0 = P(end,:).';                          % s=T corresponds to t=0
    res.ptil = ptil0; res.p0 = ptil0;
    res.p = @(t) interp1(T - sg, flipud(P), min(max(t,0),T), 'linear').';
    % Linear saturating feedback rows eq:pos:ocKrow: (Kc)_i = sign(sigma_i) (Ec)_i
    res.Kc = @(t) diag(sign(rc + Bc.'*res.p(t)))*Ec;
    res.feas = true;
    if ~isempty(x0); res.Jopt = ptil0.'*x0; end
end

% ======================================================================= %
%  Small utilities                                                        %
% ======================================================================= %
function v = getf(s, f, alt)
    v = mvst_getf(s, f, []);
    if isempty(v) && nargin>2; v = mvst_getf(s, alt, []); end
    if isempty(v); error('mvst:pos:ocdata','SYS/OPTS missing field ''%s''.', f); end
end
function d = matdim(A); if isa(A,'function_handle'); d=size(A(0),1); else; d=size(A,1); end; end
function M = evalM(A, t); if isa(A,'function_handle'); M=A(t); else; M=A; end; end
function v = evalv(x, t); if isa(x,'function_handle'); v=x(t); else; v=x; end; end
function y = numvec(x, n)
    if isa(x,'function_handle'); y = x(0); else; y = x; end
    if isempty(y); y = zeros(n,1); end; y = y(:);
end
