function res = mvst_kyp(sys, opts)
%MVST_KYP  Operator-domain Kalman-Yakubovich-Popov engine for matrix-valued
%   dissipativity (Chapter "Dissipativity", Section sec:diss:kyp): the general
%   KYP lemma for an arbitrary-sign supply matrix M, in continuous time, discrete
%   time, and the hybrid (periodic impulsive) setting.
%
%   Default engine: SOSTOOLS (primal SDPs with the witness as primal variable,
%   solved through MVST_SOSSOLVE: SeDuMi, SDPT3, MOSEK or the SDP-free 'nosdp').
%   OPTS.engine='yalmip' selects the legacy YALMIP route (kept one release; it
%   needs YALMIP and a classical SDP solver, error mvst:opts:yalmipMissing
%   otherwise).  The regularity Gramians use base MATLAB (sylvester, Kronecker
%   solve), no Control System Toolbox.
%
%   RES = MVST_KYP(SYS, OPTS).  OPTS.domain in {'ct'(default),'dt','hybrid'};
%   OPTS.engine in {'sostools'(default),'yalmip'}.
%
%   SYSTEM.  The LTI (matrix-valued covariance) system with state channel A,
%   disturbance channel Bw and, for the named supplies, output z=C x + Dw w:
%       ct :  d/dt X1 = F(Z),   F(Z)=A X1+X1 A'+Bw X2'+X2 Bw',   X1=E Z E'
%       dt :  X1+ = D(Z),       D(Z)=[A Bw] Z [A Bw]'
%   on the composite second moment Z=[X1 X2; X2' X3] (state X1 n x n, cross X2
%   n x m, disturbance X3 m x m), E=[I_n 0].  The flow adjoint (eq:diss:kyp:split)
%   is F*(P)=[A'P+P A, P Bw; Bw'P, 0]; the discrete one is D*(P)=[A Bw]'P[A Bw].
%
%   SUPPLY (OPTS.supply): 'boundedreal' (L2-gain gamma; OPTS.gamma, or the engine
%   bisects the least feasible gamma) with M=[-C'C, -C'Dw; -Dw'C, gamma^2 I-Dw'Dw];
%   'passivity' with M=[0 C'; C Dw+Dw']; or a general matrix
%   OPTS.M (n+m square, arbitrary sign) with OPTS.Mk the jump supply (hybrid).
%
%   CT KYP (th:diss:kypCT + lem:diss:kyp:dc).  Under F_1 exp. stable and the
%   noise-channel reachability (R) (some Ybar>0 in ker F; for this deterministic
%   lift Phi(0,.)=0, and (R) holds iff (A,Bw) is controllable, see CT_REGULARITY), the
%   equivalent conditions are:
%     (a) LMI feasibility:  F*(P)+M >= 0  for some P=P';
%     (b) zero-frequency inequality:  <M,Y> >= 0  for all Y in ker F cap PSD.
%   The engine solves both, reports which certifies, and returns the
%   frequency-domain WITNESS Y (in the section, with <M,Y> < 0) when infeasible.
%   If (R) fails the theorem does not apply and the engine REFUSES to conclude
%   (res.regular=false, res.feas=[]): the reachability hypothesis is required.
%
%   DT KYP (th:diss:kypDT + lem:diss:kypDT:dc).  Under D_1 Schur and (Rd) (the DT
%   DC Gramian (I-D_1)^{-1}Phi_D(0,I)>0): (a) D*(P)-E'P E+M >= 0; (b) <M,Y> >= 0
%   for Y in {D(Y)=E Y E'} cap PSD.
%
%   HYBRID KYP (th:diss:kyp:hyb:nec/loss, prop:diss:kyp:hyb:collapse).  Periodic
%   impulsive data (flow A over [0,Tp], jump J at Tp).  The engine covers:
%     - the monodromy Schur check (M = J e^{A Tp} on the covariance operator);
%     - (R_hyb) via the one-period reachability Gramian (def:diss:kyp:hyb:reach);
%     - the COLLAPSED case (M,P constant): the whole characterization collapses to
%       the monodromy DT-KYP  Mono*(P0)-E'P0 E+Mper >= 0 (prop collapse), solved by
%       the DT path with the monodromy operator in place of D_1;
%     - the rank-one hybrid-harmonic sweep, PERMITTED ONLY in the collapsed or
%       jump-free cases (cor:diss:kyp:hyb:sweeplossless); OUTSIDE them the sweep is
%       NON-EXHAUSTIVE (prop:diss:kyp:hyb:nonexh) and the engine will not use it.
%
%   RES FIELDS: .domain, .supply, .feas (LMI feasible; [] if regularity fails),
%     .P (certificate), .fdi (zero/unit-frequency inequality holds), .certifies
%     ('both'|'neither'|'lmi'|'fdi'), .witness (section point Y with <M,Y><0 when
%     infeasible), .regular ((R)/(Rd)/(R_hyb) holds), .gamma (bounded-real least
%     feasible gain), .M, .monoSchur (hybrid), .info (.engine; .dual = SDP dual
%     of the LMI program of the last solve when the SOSTOOLS backend is a
%     classical SDP solver, empty for 'nosdp' and the YALMIP route; .val =
%     optimal <M,Y> of the zero-frequency program).
%
%   See also MVST_L2GAIN, MVST_PASSIVITY_CT, MVST_DISS_SOS.

    if nargin < 2; opts = struct(); end
    dom = lower(mvst_getf(opts,'domain','ct'));
    engine = lower(mvst_getf(opts,'engine','sostools'));
    if ~any(strcmp(engine,{'sostools','yalmip'}))
        error('mvst:kyp:engine','OPTS.engine must be ''sostools'' (default) or ''yalmip''.');
    end
    if strcmp(engine,'yalmip')
        oo = mvst_opts(opts, true, 'yalmip');
        yml = sdpsettings('solver',oo.solver,'verbose',0);
    else
        oo = mvst_opts(opts, true, 'sos');
        yml = [];
    end
    ctx = struct('engine',engine,'yml',yml,'sos',mvst_sosopts(oo),'dual',[],'val',[]);

    A  = mvst_getf(sys,'A',[]);  if isempty(A); error('mvst:kyp:A','SYS needs .A.'); end
    n  = size(A,1);
    Bw = mvst_getf(sys,'Bw', mvst_getf(sys,'B', zeros(n,0)));
    m  = size(Bw,2);
    C  = mvst_getf(sys,'C', zeros(0,n));
    Dw = mvst_getf(sys,'Dw', zeros(size(C,1), m));

    res = struct('domain',dom,'supply',mvst_getf(opts,'supply','M'),'feas',[], ...
                 'P',[],'fdi',[],'certifies','','witness',[],'regular',[], ...
                 'gamma',[],'M',[],'monoSchur',[],'info',[]);

    switch dom
        case 'ct';     [res,ctx] = kyp_ct(A,Bw,C,Dw,n,m,opts,ctx,res);
        case 'dt';     [res,ctx] = kyp_dt(A,Bw,C,Dw,n,m,opts,ctx,res);
        case 'hybrid'; [res,ctx] = kyp_hybrid(A,Bw,C,Dw,n,m,opts,ctx,res);
        otherwise; error('mvst:kyp:domain','OPTS.domain must be ct|dt|hybrid.');
    end
    if isempty(res.info); res.info = struct(); end
    res.info.engine = engine; res.info.dual = ctx.dual; res.info.val = ctx.val;
end

% ======================================================================= %
%  Supply matrices                                                        %
% ======================================================================= %
function M = supply_matrix(A,Bw,C,Dw,n,m,opts,gam)
    sup = lower(mvst_getf(opts,'supply','M'));
    switch sup
        case 'boundedreal'
            if nargin < 8 || isempty(gam); gam = mvst_getf(opts,'gamma',1); end
            M = [ -C.'*C,        -C.'*Dw;
                  -Dw.'*C,        gam^2*eye(m) - Dw.'*Dw ];
        case 'passivity'
            M = [ zeros(n),       C.';
                  C,              Dw + Dw.' ];
        otherwise
            M = mvst_getf(opts,'M',[]);
            if isempty(M); error('mvst:kyp:M','General supply needs OPTS.M.'); end
    end
    M = (M+M.')/2;
end

% ======================================================================= %
%  Continuous-time KYP (th:diss:kypCT)                                     %
% ======================================================================= %
function [res,ctx] = kyp_ct(A,Bw,C,Dw,n,m,opts,ctx,res)
    % Bounded-real gain bisection: least gamma with the LMI feasible.
    if strcmpi(mvst_getf(opts,'supply',''),'boundedreal') && isempty(mvst_getf(opts,'gamma',[]))
        res.gamma = brl_bisect(@(g) lmi_feasible('ct',A,Bw,supply_matrix(A,Bw,C,Dw,n,m,opts,g),n,m,ctx), ...
                               norm(Dw)+1e-6);
        opts.gamma = res.gamma;
    end
    M = supply_matrix(A,Bw,C,Dw,n,m,opts); res.M = M;
    res.regular = ct_regularity(A,Bw,n,m);
    [res.feas, res.P, ctx] = lmi_feasible_ctx('ct',A,Bw,M,n,m,ctx);
    [res.fdi, res.witness, ctx] = fdi_solve('ct',A,Bw,M,n,m,ctx);
    res.certifies = verdict(res);
    if ~res.regular
        res.feas = []; res.certifies = 'refused';     % (R) fails: theorem N/A
    end
end

function ok = ct_regularity(A,Bw,n,m)
    % (R): there is Ybar > 0 with F(Ybar) = 0.  For this deterministic lift
    % Phi(X2,X3) = Bw X2' + X2 Bw', so Phi(0,.) = 0 and the DC-Gramian sufficient
    % condition of the book never applies; (R) is produced by the CROSS block.
    % With A Hurwitz, (R) holds iff (A,Bw) is controllable, i.e. iff the
    % controllability Gramian Wc (A Wc + Wc A' + Bw Bw' = 0) is positive definite:
    %   if Wc > 0, take Ybar1 = Wc, Ybar2 = Bw/2 (then A Ybar1 + Ybar1 A'
    %   + Bw Ybar2' + Ybar2 Bw' = -Bw Bw' + Bw Bw' = 0) and
    %   Ybar3 = I_m + Ybar2' Wc^{-1} Ybar2 (Schur complement I_m > 0), so Ybar > 0;
    %   conversely, F(Ybar) = 0 gives Ybar1 = int_0^inf e^{At}(Bw Ybar2' + Ybar2 Bw')
    %   e^{A't} dt, and a left eigenvector v of A with v'Bw = 0 would give
    %   v'Ybar1 v = 0, contradicting Ybar1 > 0.
    if m == 0; ok = false; return; end
    if max(real(eig(A))) >= 0; ok = false; return; end
    Wc = sylvester(A, A.', -Bw*Bw.');                          % A Wc + Wc A' + Bw Bw' = 0
    ok = min(real(eig((Wc+Wc.')/2))) > 1e-9;
end

% ======================================================================= %
%  LMI and zero-frequency (section) kernels, two engines                   %
%  ctx.engine = 'sostools' (default, MVST_SOSSOLVE; the witness is the     %
%  PRIMAL PSD variable of the zero-frequency program) | 'yalmip' (legacy). %
% ======================================================================= %
function ok = lmi_feasible(dom,A,Bw,M,n,m,ctx)
    ok = lmi_feasible_ctx(dom,A,Bw,M,n,m,ctx);
end

function [feas, P, ctx] = lmi_feasible_ctx(dom,A,Bw,M,n,m,ctx)
% LMI:  ct F*(P)+M >= 0,   dt D*(P)-E'PE+M >= 0,   P = P'.
    M = (M+M.')/2;
    if strcmp(ctx.engine,'yalmip')
        P = sdpvar(n,n,'symmetric');
        switch dom
            case 'ct'
                Cn = [ A.'*P + P*A,  P*Bw; Bw.'*P, zeros(m) ] + M;
            case 'dt'
                AB = [A Bw]; E = [eye(n) zeros(n,m)];
                Cn = AB.'*P*AB - E.'*P*E + M;
        end
        d = optimize(Cn >= 0, [], ctx.yml);
        feas = (d.problem == 0); if feas; P = value(P); else; P = []; end
        return;
    end
    pvar th;
    prog = sosprogram(th);
    [prog,Pv] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
    switch dom
        case 'ct'
            Cn = [A.'*Pv + Pv*A, Pv*Bw; Bw.'*Pv, zeros(m)] + M;
        case 'dt'
            AB = [A Bw]; E = [eye(n) zeros(n,m)];
            Cn = AB.'*Pv*AB - E.'*Pv*E + M;
    end
    prog = sosmatrixineq(prog, Cn);
    [prog,i1] = mvst_sossolve(prog, ctx.sos);
    % Clean-solve test: numerr == 0 on the classical backends (the YALMIP route counts
    % numerr == 1 as a numerical problem, i.e. infeasible, which keeps the bisected gamma
    % on the feasible side of the boundary); the SOSNOSDP map reports reduced-accuracy
    % certified feasibility as numerr == 1, which is accepted there.
    nmax = 1 + strcmp(ctx.sos.solver,'nosdp');
    feas = (i1.pinf == 0 && i1.dinf == 0 && i1.numerr < nmax);
    P = []; if feas; P = double(sosgetsol(prog,Pv)); end
    if isfield(prog.solinfo,'y') && ~isempty(prog.solinfo.y) && isfield(prog.solinfo,'var') ...
            && isfield(prog.solinfo.var,'dual')
        ctx.dual = struct('y',prog.solinfo.y, 'gram', {prog.solinfo.var.dual});
    else
        ctx.dual = [];
    end
end

function [ok, Y, ctx] = fdi_solve(dom,A,Bw,M,n,m,ctx)
% min <M,Z> s.t. Z >= 0, F(Z)=0 (ct) / D(Z)=E Z E' (dt), tr Z = 1;  ok iff optimum >= -tol.
    M = (M+M.')/2; AB = [A Bw]; E = [eye(n) zeros(n,m)]; ix = 1:n; iu = n+1:n+m;
    Y = [];
    if strcmp(ctx.engine,'yalmip')
        Z = sdpvar(n+m,n+m,'symmetric');
        switch dom
            case 'ct'
                X1 = Z(ix,ix); X2 = Z(ix,iu);
                Cn = [Z >= 0, A*X1 + X1*A.' + Bw*X2.' + X2*Bw.' == 0, trace(Z) == 1];
            case 'dt'
                Cn = [Z >= 0, AB*Z*AB.' == E*Z*E.', trace(Z) == 1];
        end
        d = optimize(Cn, trace(M*Z), ctx.yml);
        val = value(trace(M*Z));
        ok = (d.problem ~= 0) || (val >= -1e-7);       % infeasible section => vacuously ok
        if d.problem == 0 && val < -1e-7; Y = value(Z); end
        ctx.val = val;
        return;
    end
    pvar th;
    prog = sosprogram(th);
    [prog,Z] = sospolymatrixvar(prog, monomials(th,0), [n+m n+m], 'symmetric');
    prog = sosmatrixineq(prog, Z);
    switch dom
        case 'ct'
            X1 = Z(ix,ix); X2 = Z(ix,iu);
            R = A*X1 + X1*A.' + Bw*X2.' + X2*Bw.';
        case 'dt'
            R = AB*Z*AB.' - E*Z*E.';
    end
    for i = 1:n
        for j = i:n
            prog = soseq(prog, R(i,j));
        end
    end
    prog = soseq(prog, trace(Z) - 1);
    prog = sossetobj(prog, trace(M*Z));
    [prog,i2] = mvst_sossolve(prog, ctx.sos);
    solved = (i2.pinf == 0 && i2.dinf == 0 && i2.numerr < 2);
    val = NaN;
    if solved
        Zv = double(sosgetsol(prog,Z)); val = trace(M*Zv);
    end
    ok = ~solved || (val >= -1e-7);                    % infeasible section => vacuously ok
    if solved && val < -1e-7; Y = Zv; end
    ctx.val = val;
end


% ======================================================================= %
%  Discrete-time KYP (th:diss:kypDT)                                       %
% ======================================================================= %
function [res,ctx] = kyp_dt(A,Bw,C,Dw,n,m,opts,ctx,res)
    if strcmpi(mvst_getf(opts,'supply',''),'boundedreal') && isempty(mvst_getf(opts,'gamma',[]))
        res.gamma = brl_bisect(@(g) lmi_feasible('dt',A,Bw,supply_matrix(A,Bw,C,Dw,n,m,opts,g),n,m,ctx), ...
                               norm(Dw)+1e-6);
        opts.gamma = res.gamma;
    end
    M = supply_matrix(A,Bw,C,Dw,n,m,opts); res.M = M;
    res.regular = dt_regularity(A,Bw,n,m);
    [res.feas, res.P, ctx] = lmi_feasible_ctx('dt',A,Bw,M,n,m,ctx);
    [res.fdi, res.witness, ctx] = fdi_solve('dt',A,Bw,M,n,m,ctx);
    res.certifies = verdict(res);
    if ~res.regular; res.feas = []; res.certifies = 'refused'; end
end

function ok = dt_regularity(A,Bw,n,m)
    % (Rd): DT DC Gramian (I-D_1)^{-1}Phi_D(0,I) > 0, = DT Lyapunov solution (Kronecker solve, base MATLAB).
    if m == 0; ok = false; return; end
    if max(abs(eig(A))) >= 1; ok = false; return; end
    Wc = reshape((eye(n*n) - kron(conj(A), A)) \ reshape(Bw*Bw.', [], 1), n, n);  % A Wc A' - Wc + Bw Bw' = 0
    ok = min(real(eig((Wc+Wc.')/2))) > 1e-9;
end

% ======================================================================= %
%  Hybrid KYP (collapse case + monodromy + sweep guard)                    %
% ======================================================================= %
function [res,ctx] = kyp_hybrid(A,Bw,C,Dw,n,m,opts,ctx,res)
    Tp = mvst_getf(opts,'Tp',[]);
    J  = mvst_getf(opts,'J', eye(n));
    if isempty(Tp); error('mvst:kyp:Tp','hybrid needs OPTS.Tp.'); end
    % Monodromy super-operator on the covariance state: Mono(X)=Jt X Jt', Jt=J e^{A Tp}.
    Jt = J*expm(A*Tp);
    res.monoSchur = max(abs(eig(Jt))) < 1;              % Schur on the covariance operator
    res.regular = dt_regularity(Jt, ct_reach_input(A,Bw,Tp), n, size(ct_reach_input(A,Bw,Tp),2));
    collapsed = mvst_getf(opts,'collapsed', true);      % M,P constant over the period
    if collapsed
        % prop:diss:kyp:hyb:collapse: the characterization collapses to the
        % monodromy DT-KYP with the per-period aggregated supply Mper.
        Mper = supply_matrix(Jt, ct_reach_input(A,Bw,Tp), C, Dw, n, ...
                             size(ct_reach_input(A,Bw,Tp),2), opts);
        res.M = Mper;
        BwT = ct_reach_input(A,Bw,Tp); mT = size(BwT,2);
        [res.feas, res.P, ctx] = lmi_feasible_ctx('dt', Jt, BwT, Mper, n, mT, ctx);
        [res.fdi, res.witness, ctx] = fdi_solve('dt', Jt, BwT, Mper, n, mT, ctx);
        res.certifies = verdict(res);
        if ~res.monoSchur || ~res.regular
            res.feas = []; res.certifies = 'refused';
        end
    else
        % Genuine non-collapsed case: solve the periodic differential-difference
        % LMI (th:diss:kyp:hyb:nec / eq:diss:kyp:hyb:lmi) with a timer-varying
        % storage P(tau) via the SOS-timer route.  The rank-one hybrid-harmonic
        % SWEEP is non-exhaustive here (prop:diss:kyp:hyb:nonexh), but the
        % periodic LMI itself is the correct object and IS solved.
        res = kyp_hybrid_nc(A,Bw,C,Dw,n,m,Tp,J,opts,res);
        if ~res.monoSchur || ~res.regular
            res.feas = []; res.certifies = 'refused';
        end
    end
end

% ----------------------------------------------------------------------- %
% Non-collapsed periodic differential-difference LMI (SOS-timer).         %
%   C_M(tau,P(tau)) = [dotP + A'P + PA, P Bw; Bw'P, 0] + M   >= 0 on [0,Tp]
%   D_M(P(0+),P(Tp-)) = [J F]' P(0+) [J F] - Eb' P(Tp-) Eb + M_k >= 0      %
% (eq:IMP:LM, eq:diss:kyp:hyb:lmi; deterministic rank-one reading).       %
% ----------------------------------------------------------------------- %
function res = kyp_hybrid_nc(A,Bw,C,Dw,n,m,Tp,J,opts,res)
    F  = mvst_getf(opts,'F', zeros(n,0));   mu = size(F,2);
    Mk = mvst_getf(opts,'Mk', zeros(n+mu)); Mk = (Mk+Mk.')/2;
    deg = mvst_getf(opts,'deg',4);
    br  = strcmpi(mvst_getf(opts,'supply',''),'boundedreal');
    if br && isempty(mvst_getf(opts,'gamma',[]))
        res.gamma = brl_bisect(@(g) nc_feasible(A,Bw,C,Dw,n,m,Tp,J,F,mu,Mk, ...
                        supply_matrix(A,Bw,C,Dw,n,m,opts,g),deg,opts), norm(Dw)+1e-6);
        opts.gamma = res.gamma;
    end
    M = supply_matrix(A,Bw,C,Dw,n,m,opts); res.M = M;
    [feas, P, sl] = nc_solve(A,Bw,n,m,Tp,J,F,mu,Mk,M,deg,opts);
    res.feas = feas; res.P = P; res.fdi = [];
    res.info = struct('slack',sl,'deg',deg,'noncollapsed',true);
    res.certifies = verdict(res);
end

function ok = nc_feasible(A,Bw,C,Dw,n,m,Tp,J,F,mu,Mk,M,deg,opts) %#ok<INUSD>
    ok = nc_solve(A,Bw,n,m,Tp,J,F,mu,Mk,M,deg,opts);
end

function [feas, Pv, slv] = nc_solve(A,Bw,n,m,Tp,J,F,mu,Mk,M,deg,opts)
    pvar th; gt = th*(1-th);
    sl = dpvar('sl'); prog = sosprogram(th); prog = sosdecvar(prog,sl);
    [prog,P] = sospolymatrixvar(prog,monomials(th,0:deg),[n n],'symmetric');
    Pd = diff(P,th)/Tp;
    CM = [Pd+A.'*P+P*A, P*Bw; Bw.'*P, zeros(m)] + M;          % C_M(tau,P)
    prog = matSOS(prog, CM - sl*eye(n+m), gt, th, deg, n+m);  % >= sl I on [0,Tp]
    P0 = subs(P,th,0); PT = subs(P,th,1);
    Eb = [eye(n) zeros(n,mu)]; JF = [J F];
    DM = JF.'*P0*JF - Eb.'*PT*Eb + Mk;                        % D_M(P(0+),P(Tp-))
    prog = sosmatrixineq(prog, DM - sl*eye(n+mu));            % >= sl I
    prog = sossetobj(prog, -sl);                              % maximize slack
    oo = mvst_opts(opts, true, 'sos'); o = mvst_sosopts(oo); prog = mvst_sossolve(prog, o);
    slv = double(sosgetsol(prog,sl)); feas = slv > -1e-6;
    if feas; Pv = sosgetsol(prog,P); else; Pv = []; end
end

function prog = matSOS(prog,Mx,gmul,th,deg,sz)
% Mx(tau) matrix-SOS (PSD) on [0,1]:  Mx = S0 + gmul*S1, S0,S1 matrix-SOS.
    [prog,S] = sospolymatrixvar(prog,monomials(th,0:max(0,deg-2)),[sz sz],'symmetric');
    prog = sosmatrixineq(prog, Mx - gmul*S); prog = sosmatrixineq(prog, S);
end

function BwT = ct_reach_input(A,Bw,Tp)
    % One-period disturbance-input operator: the reachability map over [0,Tp].
    n = size(A,1); BwT = integral(@(t) expm(A*t)*Bw, 0, Tp, 'ArrayValued', true);
    if isempty(Bw); BwT = zeros(n,0); end
end

% ======================================================================= %
%  Bounded-real bisection / verdict                                        %
% ======================================================================= %
function g = brl_bisect(feasfun, glo)
    ghi = max(glo, 1);
    while ~feasfun(ghi) && ghi < 1e8; ghi = ghi*2; end
    for it = 1:50
        gm = 0.5*(glo+ghi);
        if feasfun(gm); ghi = gm; else; glo = gm; end
    end
    g = ghi;
end

function v = verdict(res)
    a = ~isempty(res.feas) && res.feas; b = ~isempty(res.fdi) && res.fdi;
    if a && b; v='both'; elseif a; v='lmi'; elseif b; v='fdi'; else; v='neither'; end
end
