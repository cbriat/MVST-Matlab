function res = mvst_stab_delay(sys, opts)
%MVST_STAB_DELAY  Delay-independent stability of matrix-valued (and positive)
%   time-delay systems (Chapter "Stability", Section sec:stab:delay, and Chapter
%   "Linear Positive Systems", prop:pos:classncs).  Certifies stability for ALL
%   delays h_i >= 0 at once (delay-independent), by the linear copositive
%   Lyapunov-Krasovskii LMI and, equivalently, the zero-delay operator spectrum.
%
%   RES = MVST_STAB_DELAY(SYS, OPTS).  OPTS.domain in {'ct'(default),'dt','positive'}.
%
%   CONTINUOUS TIME (th:stab:delay, eq:stab:delaysys).  Flow
%   F0(X)=A0 X+X A0'+sum_l Al X Al', delayed channels Hi(X)=Bi X Bi'
%   (eq:stab:Hops).  Equivalent conditions (a)<=>(b)<=>(c)<=>(d):
%     (c) LMI  eq:stab:delayLMI:  P,Qi in S^n_{>0},
%           F0*(P)+sum_i Qi < 0,   Hi*(P)-Qi <= 0,
%         F0*(P)=A0'P+PA0+sum_l Al'P Al,  Hi*(P)=Bi'P Bi;
%     (d) spectrum on S^n of the zero-delay operator F0+sum_i Hi in C_{<0}.
%   SYS: .A0/.A (drift), .Al (cell, drift channels), .Ad (cell of delayed Bi).
%
%   DISCRETE TIME (th:stab:delayDT, eq:stab:delaysysDT).  J0(X)=sum_l Jl X Jl',
%   delayed Ti(X)=Hi X Hi'.  LMI eq:stab:delayLMIDT: J0*(P)-P+sum_i Qi < 0,
%   Ti*(P)-Qi <= 0; spectrum of J0+sum_i Ti inside the unit disk.
%   SYS: .Jl (cell, drift channels), .Hd (cell of delayed Hi).
%
%   POSITIVE (prop:pos:classncs, eq:pos:delaysys).  Vector positive delay system
%   dx=A0 x+A1 x(t-h), A0 Metzler, A1>=0.  Delay-independently stable IFF the
%   zero-delay matrix A0+A1 is Hurwitz, certified by the linear copositive
%   Lyapunov-Krasovskii functional V=p'x+int_{t-h}^t p'A1 x(s)ds with
%   (A0+A1)'p<0, p>0 (the copositive LP of MVST_POS_STABILITY on A0+A1).
%   SYS: .A0 (Metzler), .A1 (nonnegative).
%
%   OPTS.engine selects the LMI backend of the CT/DT certificates: 'sostools'
%   (default; degree-zero matrix SOS variables via MVST_SOSSOLVE, so every
%   backend of OPTS.solver including 'nosdp' works) or 'yalmip' (legacy route,
%   kept for one release for cross-checks; needs YALMIP, no 'nosdp').
%
%   RES: .stable, .margin (spectral abscissa CT / log-radius DT / abscissa of
%   A0+A1 positive), .spectrum, .P/.Q (LMI certificate), .p (positive LKF
%   vector), .domain, .kind, .info.
%
%   See also MVST_STABILITY, MVST_POS_STABILITY, MVST_LIFTSYM.

    if nargin < 2; opts = struct(); end
    domain = lower(mvst_getf(opts,'domain','ct'));
    oo = mvst_opts(opts,false); solver = oo.solver; eps0 = mvst_getf(opts,'eps',1e-6);
    engine = lower(mvst_getf(opts,'engine','sostools'));
    if ~any(strcmp(engine,{'sostools','yalmip'}))
        error('mvst:delay:engine','OPTS.engine must be ''sostools'' (default) or ''yalmip''.');
    end
    res = struct('stable',false,'margin',[],'spectrum',[],'P',[],'Q',[], ...
                 'p',[],'domain',domain,'kind','','info',[]);

    switch domain
        case 'positive'
            A0 = mvst_getf(sys,'A0', mvst_getf(sys,'A',[]));
            A1 = mvst_getf(sys,'A1', []);
            if isempty(A0)||isempty(A1); error('mvst:delay:pos','positive needs .A0,.A1.'); end
            M = A0 + A1;                                   % zero-delay Metzler matrix
            pr = mvst_pos_stability(struct('A',M), struct('class','ct'));
            res.stable = pr.stable; res.margin = max(real(eig(M)));
            res.spectrum = eig(M); res.p = mvst_getf(pr,'p',[]);
            res.kind = 'positive:copositive-LK (A0+A1)''p<0';
            res.info = pr;

        case 'dt'
            Jl = tocell(mvst_getf(sys,'Jl', mvst_getf(sys,'J',{})));
            Hd = tocell(mvst_getf(sys,'Hd', mvst_getf(sys,'Ad',{})));
            n  = size(Jl{1},1);
            G  = mvst_liftsym(zeros(n), [Jl, Hd]);         % zero-delay operator J0+sum Ti
            res.spectrum = eig(G); res.margin = log(max(abs(res.spectrum)));
            [res.stable,res.P,res.Q] = delay_lmi_dt(Jl,Hd,n,eps0,oo,engine);
            res.kind = 'dt:delay-independent-LMI (+log-spectral-radius)';

        otherwise    % 'ct'
            A0 = mvst_getf(sys,'A0', mvst_getf(sys,'A',[]));
            Al = tocell(mvst_getf(sys,'Al',{}));
            Ad = tocell(mvst_getf(sys,'Ad', mvst_getf(sys,'Hd',{})));
            n  = size(A0,1);
            G  = mvst_liftsym(A0, [Al, Ad]);               % zero-delay operator F0+sum Hi
            res.spectrum = eig(G); res.margin = max(real(res.spectrum));
            [res.stable,res.P,res.Q] = delay_lmi_ct(A0,Al,Ad,n,eps0,oo,engine);
            res.kind = 'ct:delay-independent-LMI (+spectral-abscissa)';
    end
end

% ----------------------------------------------------------------------- %
function [ok,Pv,Qv] = delay_lmi_ct(A0,Al,Ad,n,eps0,solver,engine)
    if strcmp(engine,'sostools'); [ok,Pv,Qv] = delay_sos(A0,Al,Ad,n,eps0,solver,'ct'); return; end
    solver = solver.solver; mvst_opts(struct('solver',solver), true, 'yalmip'); yalmip('clear'); yml = sdpsettings('solver',solver,'verbose',0);
    P = sdpvar(n,n,'symmetric'); N = numel(Ad);
    F0s = A0.'*P + P*A0; for l=1:numel(Al); F0s = F0s + Al{l}.'*P*Al{l}; end
    Q = cell(1,N); sumQ = zeros(n); C = [P >= eps0*eye(n)];
    for i=1:N
        Q{i} = sdpvar(n,n,'symmetric'); C = [C, Q{i} >= eps0*eye(n)];
        C = [C, Ad{i}.'*P*Ad{i} - Q{i} <= 0];             % Hi*(P)-Qi <= 0
        sumQ = sumQ + Q{i};
    end
    C = [C, F0s + sumQ <= -eps0*eye(n)];                  % F0*(P)+sum Qi < 0
    d = optimize(C, [], yml); ok = (d.problem == 0);
    Pv = []; Qv = {}; if ok; Pv = value(P); Qv = cellfun(@value,Q,'uni',0); end
end

function [ok,Pv,Qv] = delay_lmi_dt(Jl,Hd,n,eps0,solver,engine)
    if strcmp(engine,'sostools'); [ok,Pv,Qv] = delay_sos(Jl,[],Hd,n,eps0,solver,'dt'); return; end
    solver = solver.solver; mvst_opts(struct('solver',solver), true, 'yalmip'); yalmip('clear'); yml = sdpsettings('solver',solver,'verbose',0);
    P = sdpvar(n,n,'symmetric'); N = numel(Hd);
    J0s = zeros(n); for l=1:numel(Jl); J0s = J0s + Jl{l}.'*P*Jl{l}; end
    Q = cell(1,N); sumQ = zeros(n); C = [P >= eps0*eye(n)];
    for i=1:N
        Q{i} = sdpvar(n,n,'symmetric'); C = [C, Q{i} >= eps0*eye(n)];
        C = [C, Hd{i}.'*P*Hd{i} - Q{i} <= 0];             % Ti*(P)-Qi <= 0
        sumQ = sumQ + Q{i};
    end
    C = [C, J0s - P + sumQ <= -eps0*eye(n)];              % J0*(P)-P+sum Qi < 0
    d = optimize(C, [], yml); ok = (d.problem == 0);
    Pv = []; Qv = {}; if ok; Pv = value(P); Qv = cellfun(@value,Q,'uni',0); end
end

function [ok,Pv,Qv] = delay_sos(X0,Al,Ad,n,eps0,solver,dom)
% SOSTOOLS (degree-zero matrix SOS) form of eq:stab:delayLMI / eq:stab:delayLMIDT.
% The conditions are a cone (homogeneous in (P,Qi)), so the scale is normalized:
%   P >= I, Qi >= I, Hi*(P) <= Qi, and
%   ct: F0*(P)+sum Qi <= -eps I      dt: J0*(P)-P+sum Qi <= -eps I
% (eps relative to that scale).  The returned certificate is re-verified a
% posteriori against these inequalities, because SeDuMi can report a
% tiny-scale spurious solution (numerr = 1) on a truly infeasible instance.
% X0 = A0 (ct) or the cell Jl (dt); Al = cell of drift channels (ct only).
    mvst_opts(solver, true, 'sos');
    pvar th; prog = sosprogram(th); N = numel(Ad); I = eye(n);
    [prog,P] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
    prog = psd_ineq(prog, P - I);
    Q = cell(1,N); sumQ = zeros(n);
    for i = 1:N
        [prog,Q{i}] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric');
        prog = psd_ineq(prog, Q{i} - I);
        prog = psd_ineq(prog, Q{i} - Ad{i}.'*P*Ad{i});     % Hi*(P) - Qi <= 0
        sumQ = sumQ + Q{i};
    end
    if strcmp(dom,'ct')
        F0s = X0.'*P + P*X0; for l = 1:numel(Al); F0s = F0s + Al{l}.'*P*Al{l}; end
        prog = psd_ineq(prog, -(F0s + sumQ) - eps0*I);      % F0*(P)+sum Qi < 0
    else
        J0s = zeros(n); for l = 1:numel(X0); J0s = J0s + X0{l}.'*P*X0{l}; end
        prog = psd_ineq(prog, -(J0s - P + sumQ) - eps0*I);  % J0*(P)-P+sum Qi < 0
    end
    prog = psd_flush(prog);
    o = mvst_sosopts(solver);
    [prog,info] = mvst_sossolve(prog, o);
    ok = (info.pinf == 0 && info.dinf == 0 && info.numerr < 2);
    Pv = []; Qv = {};
    if ok
        Pv = double(sosgetsol(prog,P)); Qv = cell(1,N);
        for i = 1:N; Qv{i} = double(sosgetsol(prog,Q{i})); end
        % a posteriori verification of the rows (tolerance relative to the scale)
        tol = 1e-7*max(1, norm(Pv,'fro'));
        Sq = zeros(n); chk = min(eig((Pv+Pv.')/2)) >= 1 - tol;
        for i = 1:N
            Sq = Sq + Qv{i};
            R = Qv{i} - Ad{i}.'*Pv*Ad{i};
            chk = chk && min(eig((R+R.')/2)) >= -tol && min(eig((Qv{i}+Qv{i}.')/2)) >= 1 - tol;
        end
        if strcmp(dom,'ct')
            F = X0.'*Pv + Pv*X0; for l = 1:numel(Al); F = F + Al{l}.'*Pv*Al{l}; end
        else
            F = -Pv; for l = 1:numel(X0); F = F + X0{l}.'*Pv*X0{l}; end
        end
        F = F + Sq;
        chk = chk && max(eig((F+F.')/2)) <= -eps0 + tol;
        if ~chk; ok = false; Pv = []; Qv = {}; end
    end
end

function prog = psd_ineq(prog, EXPR)
% EXPR >= 0 as a matrix SOS constraint.  1x1 constraints are collected in
% prog.mvst_pending and flushed as ONE diagonal block (see PSD_FLUSH): programs
% programs whose Gram blocks are ALL 1x1 crash SeDuMi's pretransfo (still the case;
% the former SOSNOSDP crash on mixed 1x1/larger blocks, nsd_gram>apply_A, is fixed in
% the vendored sosnosdp, so this merge is now needed for SeDuMi only).
    if size(EXPR,1) == 1
        if ~isfield(prog,'mvst_pending'); prog.mvst_pending = {}; end
        prog.mvst_pending{end+1} = EXPR;
    else
        prog = sosmatrixineq(prog, EXPR);
    end
end

function prog = psd_flush(prog)
    if ~isfield(prog,'mvst_pending') || isempty(prog.mvst_pending); return; end
    L = prog.mvst_pending; k = numel(L); if k == 1; k = 2; end
    D = 0 * L{1} * zeros(k);
    for i = 1:numel(L); Ei = zeros(k); Ei(i,i) = 1; D = D + L{i}*Ei; end
    if numel(L) == 1; D = D + 1*[0 0;0 1]; end
    prog = rmfield(prog,'mvst_pending');
    prog = sosmatrixineq(prog, D);
end

function c = tocell(X)
    if iscell(X); c=X; elseif isempty(X); c={}; else; p=size(X,3); c=cell(1,p);
        for k=1:p; c{k}=X(:,:,k); end; end
end
