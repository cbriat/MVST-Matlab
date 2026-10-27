function res = mvst_perf_dwell(sys, measure, opts)
%MVST_PERF_DWELL  Dwell-time and deterministic-peak backends for MVST_PERF
%   (Chapter "Further Performance Measures").  All timer-dependent tubes use
%   the house SOS-timer machinery (SOSTOOLS); the decay-rate certificates use a
%   two-parameter (lambda,nu) bisection whose inner LMIs are constant-degree
%   SOSTOOLS programs (OPTS.engine = 'sostools', default; every OPTS.solver
%   backend including 'nosdp') or, as a legacy cross-check route kept for one
%   release, YALMIP (OPTS.engine = 'yalmip').
%
%   RES = MVST_PERF_DWELL(SYS, MEASURE, OPTS).  MEASURE and OPTS.dwell:
%     'peak' / 'e2p' + OPTS.dwell in {min,range}: the supersolution reachability
%       tube Xi(tau)>=0 of prop:perf:peakgaindwell / th:perf:e2p:MDT/RDT,
%         dot Xi >= A Xi + Xi A' + Bw Bw'   (flow, on [0,T]),
%         Xi(0+) >= J Xi(theta) J'          (jump; theta=Tmin, or the whole
%                                            family [Tmin,Tmax] for range),
%         C Xi C' <= s I                    (output),  minimize s,
%       with the frozen tail  A Xi(Tmin)+Xi(Tmin)A'+Bw Bw' <= 0  for minimum
%       dwell-time.  peak = s (= intensity peak-to-peak, prop:perf:ptpint:dwell,
%       same tube); e2p = sqrt(s).
%     'h2' + OPTS.dwell in {min,range}: the observability-side budget of
%       th:perf:h2:MDT/RDT, linear in (Qbar(tau),eta):
%         dot Qbar + A' Qbar + Qbar A + C'C <= 0     (flow),
%         J' Qbar(0+) J - Qbar(theta) + Cd'Cd <= 0   (jump),
%         int_0^T tr(Qbar Bw Bw') dtau <= eta*T,  tr(Qbar(Tmin) Bw Bw') <= eta,
%       plus the frozen tail for minimum dwell-time; H2 = sqrt(eta).
%     'p2p': the DETERMINISTIC peak-to-peak (star-norm) of th:perf:ptp:seq /
%       th:perf:ptp:dwell, the (lambda,nu)-parameterized invariant-set LMIs
%         [dot P+A'P+P A+lambda P, P Bw; Bw'P, -lambda I] <= 0,
%         J'P(0+)J <= (1-nu) P(theta),   C'C <= gamma^2 P,   (Cd'Cd <= gamma^2 P),
%       minimized over (P,gamma^2) for fixed (lambda,nu) and bisected over
%       (lambda,nu).  With OPTS.timer=true the storage is a TIMER-VARYING
%       P(tau) (SOS-timer) for the constant/MDT/RDT families
%       (th:perf:ptp:dwell (a)/(b)/(c)); it warm-starts (lambda,nu) from the
%       constant-P solve and is never worse than that slice (dominance).
%       Companion lower bound gammaLB = ||C A^{-1} Bw||_2 (prop:perf:ptplower)
%       for the CT LTI slice.
%     'peakfull': the amplitude-bounded peak, full-block dual certificate of
%       prop:perf:peakgainfull, same (lambda,nu) LMIs with the LINEAR output
%       C'C <= gamma P (gamma, not gamma^2); minimize gamma.
%
%   OPTS: .dwell, .Tmin, .Tmax (or .Tbar for constant), .deg (SOS timer degree,
%   default 4), .J (reset, default I), .Cd (jump output), .lamgrid/.nugrid
%   (bisection grids), .solver, .engine ('sostools' | 'yalmip', p2p/peakfull LMIs).
%   RES: .gamma, .measure, .class, .dwell,
%   .gammaLB (p2p/peakfull), .lambda/.nu (chosen), .info.
%
%   See also MVST_PERF, MVST_STAB_SOS, MVST_KYP.

    if nargin < 3; opts = struct(); end
    oo = mvst_opts(opts, false); solver = oo.solver;
    A  = mvst_getf(sys,'A', mvst_getf(sys,'A0',[])); A = A(:,:,1);
    n  = size(A,1);
    Bw = mvst_getf(sys,'Bw',[]);  C = mvst_getf(sys,'C',[]); C = C(:,:,1);
    if isempty(Bw) || isempty(C); error('mvst:perf:io','need .Bw and .C.'); end
    J  = mvst_getf(sys,'J', mvst_getf(opts,'J', eye(n)));
    Cd = mvst_getf(sys,'Cd', mvst_getf(opts,'Cd', zeros(0,n)));
    dwl = lower(mvst_getf(opts,'dwell',''));
    Tmin = mvst_getf(opts,'Tmin', mvst_getf(sys,'Tmin', mvst_getf(opts,'Tbar', mvst_getf(sys,'Tbar',[]))));
    Tmax = mvst_getf(opts,'Tmax', mvst_getf(sys,'Tmax', Tmin));
    deg  = mvst_getf(opts,'deg', 4);

    engine = lower(mvst_getf(opts,'engine','sostools'));
    if ~any(strcmp(engine,{'sostools','yalmip'}))
        error('mvst:perf:engine','OPTS.engine must be ''sostools'' (default) or ''yalmip''.');
    end
    res = struct('gamma',[],'measure',measure,'class','dwell','dwell',dwl, ...
                 'gammaLB',[],'lambda',[],'nu',[],'info',[]);

    % dependency guard (SOS routes: SOSTOOLS, any solver; YALMIP routes: no 'nosdp')
    switch lower(measure)
        case {'peak','e2p','intp2p','h2'}
            mvst_opts(oo, true, 'sos');
        case {'p2p','peakfull'}
            if strcmp(engine,'yalmip'); mvst_opts(oo, true, 'yalmip'); else; mvst_opts(oo, true, 'sos'); end
            if strcmpi(measure,'p2p') && mvst_getf(opts,'timer',false); mvst_opts(oo, true, 'sos'); end
    end

    switch lower(measure)
        case {'peak','e2p','intp2p'}
            s = peak_tube(A,Bw,C,J,dwl,Tmin,Tmax,deg,oo);
            if strcmpi(measure,'e2p'); res.gamma = sqrt(max(s,0)); else; res.gamma = s; end
        case 'h2'
            eta = h2_budget(A,Bw,C,Cd,J,dwl,Tmin,Tmax,deg,oo);
            res.gamma = sqrt(max(eta,0));
        case 'p2p'
            if mvst_getf(opts,'timer',false)
                % Timer-varying certificate P(tau) (th:perf:ptp:dwell (a)/(b)/(c)):
                % SOS-timer route for the constant/MDT/RDT families.
                [g,lam,nu] = ptp_timer_bisect(A,Bw,C,Cd,J,dwl,Tmin,Tmax,deg,oo,opts);
            else
                [g,lam,nu] = ptp_bisect(A,Bw,C,Cd,J,dwl,Tmin,Tmax,oo,opts,false);
            end
            res.gamma = g; res.lambda = lam; res.nu = nu;
            if max(real(eig(A)))<0; res.gammaLB = norm(C*(A\Bw)); end
        case 'peakfull'
            [g,lam,nu] = ptp_bisect(A,Bw,C,Cd,J,dwl,Tmin,Tmax,oo,opts,true);
            res.gamma = g; res.lambda = lam; res.nu = nu;
            if max(real(eig(A)))<0; res.gammaLB = norm(C*(A\Bw)); end
        otherwise
            error('mvst:perf:measure','mvst_perf_dwell: unknown measure ''%s''.',measure);
    end
end

% ======================================================================= %
%  Peak / e2p / intensity-p2p reachability tube (SOS-timer)                %
% ======================================================================= %
function s = peak_tube(A,Bw,C,J,dwl,Tmin,Tmax,deg,solver)
    n=size(A,1); q=size(C,1);
    if strcmp(dwl,'range'); T=Tmax; else; T=Tmin; end
    pvar th; g=dpvar('g'); prog=sosprogram(th); prog=sosdecvar(prog,g);
    [prog,Xi]=sospolymatrixvar(prog,monomials(th,0:deg),[n n],'symmetric');
    gt=th*(1-th); Xid=diff(Xi,th)/T;
    % flow supersolution:  Xidot - (A Xi+Xi A'+Bw Bw') >= 0  on [0,1]
    prog=mat_ge(prog, Xid-(A*Xi+Xi*A.'+Bw*Bw.'), gt, th, deg, n);
    prog=mat_ge(prog, Xi, gt, th, deg, n);                         % Xi >= 0
    prog=mat_ge(prog, g*eye(q)-C*Xi*C.', gt, th, deg, q);         % C Xi C' <= g I
    Xi0=subs(Xi,th,0); XiT=subs(Xi,th,1);
    if strcmp(dwl,'range')
        gm=(th-Tmin/T)*(1-th);
        prog=mat_ge(prog, Xi0-J*Xi*J.', gm, th, deg, n);          % jump for all theta
        prog=mat_ge(prog, g*eye(q)-C*Xi*C.', gm, th, deg, q);     % output at theta
    else
        prog=sosmatrixineq(prog, Xi0-J*XiT*J.');                  % jump at Tmin
        prog=sosmatrixineq(prog, -(A*XiT+XiT*A.'+Bw*Bw.'));       % frozen tail
        prog=sosmatrixineq(prog, g*eye(q)-C*XiT*C.');
    end
    prog=sossetobj(prog,g); o=mvst_sosopts(solver);
    prog=mvst_sossolve(prog, o); s=double(sosgetsol(prog,g));
end

% ======================================================================= %
%  H2 dwell budget (SOS-timer observability side, linear in (Qbar,eta))    %
% ======================================================================= %
function eta = h2_budget(A,Bw,C,Cd,J,dwl,Tmin,Tmax,deg,solver)
    n=size(A,1);
    if strcmp(dwl,'range'); T=Tmax; else; T=Tmin; end
    pvar th; et=dpvar('eta'); prog=sosprogram(th); prog=sosdecvar(prog,et);
    [prog,Q]=sospolymatrixvar(prog,monomials(th,0:deg),[n n],'symmetric');
    gt=th*(1-th); Qd=diff(Q,th)/T;
    prog=mat_ge(prog, -(Qd+A.'*Q+Q*A+C.'*C), gt, th, deg, n);     % flow DLMI <=0
    prog=mat_ge(prog, Q, gt, th, deg, n);                          % Q >= 0
    Q0=subs(Q,th,0); QT=subs(Q,th,1);
    CdCd = zeros(n); if ~isempty(Cd); CdCd=Cd.'*Cd; end
    if strcmp(dwl,'range')
        gm=(th-Tmin/T)*(1-th);
        prog=mat_ge(prog, -(J.'*Q0*J-Q+CdCd), gm, th, deg, n);     % jump for all theta
    else
        prog=sosmatrixineq(prog, -(J.'*Q0*J-QT+CdCd));            % jump at Tmin
        prog=sosmatrixineq(prog, -(A.'*QT+QT*A+C.'*C));           % frozen tail
    end
    % budget: int_0^T tr(Q Bw Bw') dtau <= eta T ; tr(Q(Tmin) Bw Bw') <= eta.
    % The integrand tr(Q(th) Bw Bw') is affine in the certificate coefficients;
    % integrate over tau=T*th by trapezoidal quadrature (each node a linear
    % dpvar expression), so the budget is an LP row in (Q-coeffs, eta).
    BB=Bw*Bw.'; Nq=40; ths=linspace(0,1,Nq); wq=ones(1,Nq); wq(1)=0.5; wq(end)=0.5; wq=wq/(Nq-1);
    intval = 0;
    for i=1:Nq; intval = intval + wq(i)*trace(subs(Q,th,ths(i))*BB); end   % ~ int_0^1
    prog=sosineq(prog, et*T - T*intval);                          % int_0^T <= eta T
    prog=sosineq(prog, et - trace(subs(Q,th,1)*BB));             % tr(Q(Tmin)BwBw')<=eta
    prog=sossetobj(prog,et); o=mvst_sosopts(solver);
    prog=mvst_sossolve(prog, o); eta=double(sosgetsol(prog,et));
end

% ======================================================================= %
%  Deterministic p2p (star-norm) / full-block peak: (lambda,nu) bisection  %
% ======================================================================= %
function [gbest,lambest,nubest] = ptp_bisect(A,Bw,C,Cd,J,dwl,Tmin,Tmax,solver,opts,fullblock)
    ctx = struct('engine',mvst_getf(opts,'engine','sostools'),'solver',solver,'yml',[]);
    if strcmp(ctx.engine,'yalmip'); ctx.yml = sdpsettings('solver',solver.solver,'verbose',0); end
    yml = ctx;
    jumpfree = isequal(J,eye(size(J,1))) && isempty(dwl);
    lamg = mvst_getf(opts,'lamgrid', logspace(-1.5,1.5,13));       % flow decay-rate grid
    if jumpfree; nug = 0; else; nug = mvst_getf(opts,'nugrid', linspace(0,1,6)); end
    gbest=inf; lambest=[]; nubest=[];
    for lam=lamg
        for nu=nug
            g = ptp_leastgamma(A,Bw,C,Cd,J,lam,nu,yml,fullblock);
            if g < gbest; gbest=g; lambest=lam; nubest=nu; end
        end
    end
    % local refine on lambda (golden-section) around the best grid point.
    if isfinite(gbest) && ~isempty(lambest)
        lo=lambest/2; hi=lambest*2;
        for it=1:12
            m1=lo+(hi-lo)/3; m2=hi-(hi-lo)/3;
            g1=ptp_leastgamma(A,Bw,C,Cd,J,m1,nubest,yml,fullblock);
            g2=ptp_leastgamma(A,Bw,C,Cd,J,m2,nubest,yml,fullblock);
            if g1<g2; hi=m2; if g1<gbest; gbest=g1; lambest=m1; end
            else;     lo=m1; if g2<gbest; gbest=g2; lambest=m2; end; end
        end
    end
end

function g = ptp_leastgamma(A,Bw,C,Cd,J,lam,nu,yml,fullblock)
    % Least gamma for fixed (lambda,nu).  The output condition C'C <= gamma^2 P
    % (p2p) / gamma P (full-block) is BILINEAR in (gamma,P) -- a GEVP.
    % 'yalmip' engine: bisect gamma, each step a genuine LMI in P (legacy).
    % 'sostools' engine: with mu = 1/gamma^2 (p2p) / 1/gamma (full-block) the
    %   output rows become P >= mu C'C (and P >= mu Cd'Cd), LINEAR in (P,mu), so
    %   the least gamma is ONE SDP: maximize mu (same optimum as the bisection).
    if strcmp(yml.engine,'sostools')
        g = ptp_gevp_sos(A,Bw,C,Cd,J,lam,nu,yml.solver,fullblock); return;
    end
    feas = @(gm) ptp_feas(A,Bw,C,Cd,J,lam,nu,gm,yml,fullblock);
    glo=0; ghi=1; while ~feas(ghi) && ghi<1e6; ghi=ghi*2; end
    if ~feas(ghi); g=inf; return; end
    for it=1:36; gm=0.5*(glo+ghi); if feas(gm); ghi=gm; else; glo=gm; end; end
    g=ghi;
end

function g = ptp_gevp_sos(A,Bw,C,Cd,J,lam,nu,solver,fullblock)
    n=size(A,1); m=size(Bw,2);
    % output weights normalized to unit scale (mu is then O(1)); gamma is rescaled back.
    sc = max([norm(C,'fro'), norm(Cd,'fro')]);
    if sc == 0; g = 0; return; end
    CC = (C.'*C)/sc^2; CdCd = (Cd.'*Cd)/sc^2;
    % disturbance scale: with Bw = b*Bh, ||Bh||_F = 1, the substitution P = Ph/b^2, mu = muh/b^2
    % (congruence diag(I,1/b) on the flow LMI) maps the problem EXACTLY to Bh with gamma = b*gamma_h,
    % so mu stays O(1) for any scale of Bw.
    b = norm(Bw,'fro'); if b == 0; g = 0; return; end
    Bw = Bw/b;
    cap = 1e4;                                                         % mu <= cap; raised while active
    for itcap = 1:6
        pvar th; prog=sosprogram(th); mu=dpvar('mu'); prog=sosdecvar(prog,mu);
        [prog,P]=sospolymatrixvar(prog,monomials(th,0),[n n],'symmetric');
        prog=sosmatrixineq(prog, P-1e-6*eye(n));
        Flow=[A.'*P+P*A+lam*P, P*Bw; Bw.'*P, -lam*eye(m)];
        prog=sosmatrixineq(prog, -Flow);                                   % flow invariant-set LMI
        if ~isequal(J,eye(n))                                              % jump
            prog=sosmatrixineq(prog, (1-nu)*P-J.'*P*J);
        end
        prog=sosmatrixineq(prog, P-mu*CC);                                 % P >= mu C'C
        if ~isempty(Cd); prog=sosmatrixineq(prog, P-mu*CdCd); end
        prog=sosmatrixineq(prog, cap-mu);                                  % mu bounded
        prog=sossetobj(prog,-mu);
        o=mvst_sosopts(solver);
        [prog,info]=mvst_sossolve(prog,o);
        ok = (info.pinf==0 && info.dinf==0 && info.numerr<2);
        if ~ok; g = inf; return; end
        mv = double(sosgetsol(prog,mu));
        if ~isfinite(mv) || mv <= 1e-12; g = inf; return; end
        if mv < 0.99*cap; break; end                                       % cap inactive: done
        if itcap == 6
            warning('mvst:perf:mucap','gamma bound limited by the mu cap (mu=%g).',cap);
        end
        cap = cap*1e3;
    end
    % mu = sc^2/gamma^2 (p2p) / sc^2/gamma (full-block)
    if fullblock; g = sc^2/mv; else; g = sc/sqrt(mv); end
    % gamma^2 scales as b^2 (p2p), gamma as b^2 (full-block, both output rows share the rescaled P)
    if fullblock; g = b^2*g; else; g = b*g; end
end

function ok = ptp_feas(A,Bw,C,Cd,J,lam,nu,gam,yml,fullblock)
    % Feasibility of the invariant-set LMIs in P at a fixed decay rate and gain.
    yml = yml.yml;
    n=size(A,1); m=size(Bw,2);
    P=sdpvar(n,n,'symmetric');
    Cn=[P>=1e-6*eye(n), ...
        [A.'*P+P*A+lam*P, P*Bw; Bw.'*P, -lam*eye(m)] <= 0];   % flow invariant-set LMI
    if ~isequal(J,eye(n)); Cn=[Cn, J.'*P*J <= (1-nu)*P]; end   % jump J'P(0)J <= (1-nu)P
    outc = gam^2; if fullblock; outc = gam; end                % gamma^2 (p2p) / gamma (full)
    Cn=[Cn, C.'*C <= outc*P];
    if ~isempty(Cd); Cn=[Cn, Cd.'*Cd <= outc*P]; end
    d=optimize(Cn, [], yml); ok=(d.problem==0);
end

% ======================================================================= %
%  Timer-varying deterministic p2p (star-norm): SOS-timer P(tau)           %
%  (th:perf:ptp:dwell (a) constant, (b) minimum, (c) range dwell-time)     %
% ======================================================================= %
function [gbest,lambest,nubest] = ptp_timer_bisect(A,Bw,C,Cd,J,dwl,Tmin,Tmax,deg,solver,opts)
% Warm-start (lambda,nu) from the fast constant-P solve, then refine the gain
% with the timer-varying certificate at that (and neighboring) decay rates.
% The timer certificate is never worse than the constant-P slice (dominance).
    [gc,lam0,nu0] = ptp_bisect(A,Bw,C,Cd,J,dwl,Tmin,Tmax,solver,opts,false);
    if ~isfinite(gc) || isempty(lam0)
        gbest=gc; lambest=lam0; nubest=nu0; return;
    end
    F  = mvst_getf(opts,'F', zeros(size(A,1),0));
    lamset = mvst_getf(opts,'lamtimer', lam0);              % warm-started decay rate
    gbest = gc; lambest = lam0; nubest = nu0;       % constant slice is the fallback
    % The constant-slice nu0 can be a non-unique optimum (backend-dependent warm start): probe the
    % timer certificate on the whole nu grid as well, and keep the best (lambda,nu).
    nuset = unique([nu0, mvst_getf(opts,'nugrid', linspace(0,1,6))]);
    if isequal(J,eye(size(J,1))) && isempty(dwl); nuset = nu0; end
    for l = lamset
        for nu = nuset
            g = ptp_timer_leastgamma(A,Bw,C,Cd,J,F,dwl,Tmin,Tmax,l,nu,deg,solver);
            if g < gbest; gbest=g; lambest=l; nubest=nu; end
        end
    end
end

function g = ptp_timer_leastgamma(A,Bw,C,Cd,J,F,dwl,Tmin,Tmax,lam,nu,deg,solver)
    feas = @(gm) ptp_timer_feas(A,Bw,C,Cd,J,F,dwl,Tmin,Tmax,lam,nu,gm,deg,solver);
    glo=0; ghi=1; while ~feas(ghi) && ghi<1e6; ghi=ghi*2; end
    if ~feas(ghi); g=inf; return; end
    for it=1:14; gm=0.5*(glo+ghi); if feas(gm); ghi=gm; else; glo=gm; end; end
    g=ghi;
end

function ok = ptp_timer_feas(A,Bw,C,Cd,J,F,dwl,Tmin,Tmax,lam,nu,gam,deg,solver)
% SOS-timer feasibility of the (lambda,nu) invariant-set LMIs with a
% timer-varying storage P(tau) (eq:perf:ptpflow/ptpjump/ptpout, rank-one).
    n=size(A,1); m=size(Bw,2); mu=size(F,2);
    if strcmp(dwl,'range'); T=Tmax; else; T=Tmin; end
    pvar th; gt=th*(1-th); sl=dpvar('sl');
    prog=sosprogram(th); prog=sosdecvar(prog,sl);
    prog=sosineq(prog,1-sl);                                       % slack <= 1: scaling a feasible (P,sl) by c>1 stays feasible (monotonicity), so the maximized slack is unbounded without this cap
    [prog,P]=sospolymatrixvar(prog,monomials(th,0:deg),[n n],'symmetric');
    Pd=diff(P,th)/T;
    Flow=[Pd+A.'*P+P*A+lam*P, P*Bw; Bw.'*P, -lam*eye(m)];         % eq:perf:ptpflow
    prog=mat_ge(prog, -Flow - sl*eye(n+m), gt, th, deg, n+m);      % Flow <= -sl I
    prog=mat_ge(prog, P - 1e-6*eye(n), gt, th, deg, n);            % P > 0
    prog=mat_ge(prog, gam^2*P - C.'*C - sl*eye(n), gt, th, deg, n);% eq:perf:ptpout (mu=1)
    P0=subs(P,th,0); PT=subs(P,th,1);
    JF=[J F]; Eb=[eye(n) zeros(n,mu)]; Ebot=[zeros(mu,n) eye(mu)];
    if strcmp(dwl,'range')
        gm=(th-Tmin/T)*(1-th);
        Jm=JF.'*P0*JF - (1-nu)*(Eb.'*P*Eb) - nu*(Ebot.'*Ebot);    % all theta in [Tmin,Tmax]
        prog=mat_ge(prog, -Jm - sl*eye(n+mu), gm, th, deg, n+mu);
    else
        Jm=JF.'*P0*JF - (1-nu)*(Eb.'*PT*Eb) - nu*(Ebot.'*Ebot);   % theta = Tmin / Tp
        prog=sosmatrixineq(prog, -Jm - sl*eye(n+mu));
        if strcmp(dwl,'min')                                       % constant-certificate tail
            Tail=[A.'*PT+PT*A+lam*PT, PT*Bw; Bw.'*PT, -lam*eye(m)];
            prog=sosmatrixineq(prog, -Tail - sl*eye(n+m));
        end
    end
    if ~isempty(Cd); prog=sosmatrixineq(prog, gam^2*PT - Cd.'*Cd); end
    prog=sossetobj(prog,-sl); o=mvst_sosopts(solver);
    [prog,si]=mvst_sossolve(prog, o);
    okSolve = (si.pinf==0 && si.dinf==0 && si.numerr<2);             % failed solve = infeasible probe
    ok = okSolve && double(sosgetsol(prog,sl)) > -1e-6;
end

% ======================================================================= %
%  SOS / polynomial helpers                                                %
% ======================================================================= %
function prog = mat_ge(prog, M, gmul, th, deg, sz)
% M(th) >= 0 on the interval where gmul>=0:  M = gmul*S + SOS-matrix, S SOS-matrix.
    [prog,S]=sospolymatrixvar(prog, monomials(th,0:2*deg), [sz sz], 'symmetric');
    prog=sosmatrixineq(prog, M - gmul*S);
    prog=sosmatrixineq(prog, S);
end
