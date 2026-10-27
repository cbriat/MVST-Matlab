function res = mvst_oc_dwell(md, jump, opts)
%MVST_OC_DWELL  Causal minimum / range dwell-time optimal control of a
%   single-cone linear impulsive system (Chapter "Linear Cost", Section
%   sec:OC:dwell): the timer-dependent forward-Riccati / jump-GRE programs
%   Theorem th:OC:MDT (minimum dwell-time) and Theorem th:OC:RDT (range
%   dwell-time), with the average-cost bounds th:IMP:AC:dwell:MDT /
%   th:IMP:AC:dwell:RDT as a secondary mode.
%
%   RES = MVST_OC_DWELL(MD, JUMP, OPTS) where MD is a mode struct (MVST_MODE:
%   flow data A,B,Q,R and GBM channels Ac,Bc), JUMP carries the reset/jump-cost
%   data, and OPTS selects the dwell family and mode.
%
%   METHOD.  The dwell certificate P_s over the timer window is the fixed point
%       P_s(0^+) = FlowBack( JumpRiccati(P_s(0^+)), Delta ),
%   where FlowBack integrates the flow differential I-GRE (MVST_DRE_CT) backward
%   over the dwell window and JumpRiccati is the discrete jump I-GRE that
%   eliminates the jump input by its Schur complement.  This is the flow-control
%   generalization of the sampled-data periodic solver MVST_IGRE_PERIODIC_SD.
%   The certified cost bound is J_T(X*,X10) <= <P_s(0), X10> for every admissible
%   dwell sequence (Thm th:OC:MDT / th:OC:RDT): an UPPER bound realized by the
%   recovered policy.
%     * 'min'  : window [0,Tmin]; the flow gain is frozen at K(Tmin) beyond, and
%                the end-of-dwell LMI (iii) of th:OC:MDT,
%                   G_inf(P_s(Tmin)) = mvst_ric_ct(P_s(Tmin), md_inf) <= 0,
%                is verified (the frozen flow dissipates for every Delta>=Tmin).
%     * 'range': window [0,Tmax]; the jump gain K_d(Delta) is scheduled on the
%                just-observed dwell (causal), and the jump I-GRE is required
%                across the family Delta in [Tmin,Tmax].
%
%   AVERAGE MODE (OPTS.mode='average'): returns the long-term average-cost rate
%   upper bound (eq:IMP:AC:dwell:MDT:rho / eq:IMP:AC:dwell:RDT:rho) using the
%   timer-occupation weighting, given the forcing covariances W(tau), W_inf,
%   W_jp and the jump rate nu (or dwell distribution).
%
%   MD FIELDS (MVST_MODE): .A,.B,.Q,.R,.Ac,.Bc,.S,.n,.m.  For 'min' the flow
%   data are time-invariant beyond Tmin (eq:MDT:LTIassumption); a time-varying
%   flow may be supplied as MD a function handle @(tau)->mode (frozen at Tmin).
%   JUMP FIELDS: .J (reset map, n x n), .Bd (jump input selector, n x md;
%   default 0 = no jump control), .Qd (pre-jump state cost, n x n; default 0),
%   .Rd (jump input cost, md x md; default I).
%   OPTS FIELDS: .dwell 'min'|'range' (required); .mode 'cost'(default)|'average';
%   .Tmin, .Tmax; .X10 (initial covariance, for the numeric cost bound);
%   .tol/.maxit (fixed point, defaults 1e-9/200); .ngrid (timer grid, 200);
%   average-mode .W (@(tau)->W or matrix), .Winf, .Wjp, .nu (jump rate),
%   .mu (dwell pdf @(Delta) for 'range', default uniform on [Tmin,Tmax]).
%
%   RES FIELDS: .Jbound (<P_s(0),X10>, the certified cost upper bound; [] if
%   X10 absent), .P0 (=P_s(0^+)), .Pfun (@(tau)->P_s(tau)), .Kc (@(tau)->flow
%   gain), .Kd (jump gain: matrix for 'min', @(Delta) for 'range'), .feas
%   (fixed point converged AND, for 'min', the end-of-dwell LMI holds),
%   .lmiMDT (max eigenvalue of the (iii) residual; <=tol certifies), .rho
%   (average-cost bound, average mode), .iters,.residual,.dwell,.mode.
%
%   See also MVST_IGRE_PERIODIC_SD, MVST_DRE_CT, MVST_CARE_CT, MVST_RIC_CT.

    if nargin < 3; opts = struct(); end
    dwl  = lower(mvst_getf(opts,'dwell',''));
    mode = lower(mvst_getf(opts,'mode','cost'));
    if ~any(strcmp(dwl,{'min','range'}))
        error('mvst:ocdwell:dwell','OPTS.dwell must be ''min'' or ''range''.');
    end
    tol   = mvst_getf(opts,'tol',1e-9);
    maxit = mvst_getf(opts,'maxit',200);
    ng    = mvst_getf(opts,'ngrid',200);
    Tmin  = mvst_getf(opts,'Tmin',[]);
    Tmax  = mvst_getf(opts,'Tmax',[]);
    % Dfp = jump-closure horizon (the SHORTEST admissible dwell: the jump I-GRE
    % must close at the shortest dwell so the jump Schur contribution
    % <JumpRiccati(P0)-P_s(theta),X> stays <=0 for every theta>=Tmin, keeping the
    % bound sound); Dwin = certificate window ([0,Tmin] for MDT, extended to
    % [0,Tmax] for RDT so the causal theta-scheduled gain is available).
    if strcmp(dwl,'min')
        if isempty(Tmin); error('mvst:ocdwell:Tmin','dwell ''min'' needs Tmin.'); end
        Dfp = Tmin; Dwin = Tmin;
    else
        if isempty(Tmin)||isempty(Tmax); error('mvst:ocdwell:Tr','range needs Tmin,Tmax.'); end
        Dfp = Tmin; Dwin = Tmax;
    end

    md0 = evalmd(md, Dfp);                         % frozen mode (const data)
    n   = md0.n;
    J   = mvst_getf(jump,'J', eye(n));
    Bd  = mvst_getf(jump,'Bd', zeros(n,0)); mdj = size(Bd,2);
    Qd  = mvst_getf(jump,'Qd', zeros(n));
    Rd  = mvst_getf(jump,'Rd', eye(max(mdj,1)));

    % ---- fixed point  P0 = FlowBack( JumpRiccati(P0), Delta ) ------------- %
    P0 = eye(n); delta = inf; it = 0;
    while it < maxit && delta > tol
        it = it + 1;
        Pm = jump_riccati(P0, J, Bd, Qd, Rd);      % pre-jump P_s(Delta^-)
        [~,~,P0new] = mvst_dre_ct(md0, Dfp, Pm);   % jump closure at shortest dwell
        P0new = (P0new+P0new.')/2;
        delta = norm(P0new - P0,'fro'); P0 = P0new;
    end
    if delta > max(tol,1e-6)
        % The closure is a boundary-value problem in P(0^+): the forward
        % Riccati eq:forwardRiccati is ANTI-stable about P_inf, so it is solved
        % here as the well-conditioned backward fixed point
        % P(0^+) = FlowBack(JumpRiccati(P(0^+)), Delta) (FlowBack and the
        % forward flow are exact inverses of dP/dtau = -G(P)).  Non-convergence
        % is reported as infeasible rather than returning a bogus certificate.
        warning('mvst:ocdwell:converge', ...
            ['Dwell I-GRE fixed-point residual %.2e exceeds tol %.1e ', ...
             '(returning feas=false).'], delta, max(tol,1e-6));
    end
    Pm = jump_riccati(P0, J, Bd, Qd, Rd);
    [~, Kd] = jump_riccati(P0, J, Bd, Qd, Rd);

    % Flow certificate over the window (dense, accurately interpolated) and the
    % EXACT optimal flow gain read from it.  The gain is recomputed from the
    % certificate at each tau, NOT interpolated: interpolating K(tau) directly
    % makes the applied policy slightly suboptimal, which (over many cycles)
    % pushes the realized cost above the certified bound; the exact gain keeps
    % the flow term <C_Z,X*>=0 along the optimal trajectory, so the bound holds.
    [tt, Pt] = build_cert(md, md0, Dfp, Dwin, Pm);
    Pfun  = @(tau) interp3d(tt, Pt, min(max(tau,0),Dwin));   % pchip on P_s
    Kcfun = @(tau) exact_gain(Pfun(min(max(tau,0),Dwin)), md0);

    res = struct('Jbound',[],'P0',P0,'Pfun',Pfun,'Kc',Kcfun,'Kd',[], ...
                 'feas',false,'lmiMDT',[],'lmiRDT',[],'rho',[],'iters',it, ...
                 'residual',delta,'dwell',dwl,'mode',mode);

    conv = delta <= max(tol,1e-6);
    switch dwl
        case 'min'
            % end-of-dwell LMI (iii): frozen Riccati residual <= 0
            Gtail = mvst_ric_ct(Pm, md0);
            res.lmiMDT = max(real(eig((Gtail+Gtail.')/2)));
            res.feas = conv && res.lmiMDT <= 1e-6;
            res.Kd = Kd;                            % single jump gain (dwell=Tmin)
        case 'range'
            % causal jump gain scheduled on the observed dwell Delta_{k-1}
            res.Kd = @(D) jump_gain_at(P0, md0, J, Bd, Qd, Rd, min(max(D,Tmin),Tmax));
            % theta-FAMILY jump-GRE check (th:OC:RDT (ii), verified not assumed):
            % the jump at dwell theta contributes <D_Z(P0,P_s(theta)),X*> to the
            % key identity; with the optimal jump gain this equals the Schur
            % complement of D_Z, namely JumpRiccati(P0)-P_s(theta) = Pm-P_s(theta).
            % The bound J<=<P_s(0),X10> holds only if this is <=0 for EVERY
            % admissible theta, i.e. P_s(theta) >= Pm on [Tmin,Tmax].  The forward
            % Riccati is not monotone in general (a dipping Z(tau) can push
            % P_s(theta) below Pm at an interior theta), so this is checked on a
            % grid and folded into feasibility -- never assumed.
            thg = linspace(Tmin, Tmax, max(ng,25));
            viol = -inf;
            for th = thg
                Sc = Pm - Pfun(th);                 % Schur complement of D_Z(P0,P_s(theta))
                viol = max(viol, max(real(eig((Sc+Sc.')/2))));
            end
            res.lmiRDT = viol;                      % worst family-condition residual
            res.feas = conv && (viol <= 1e-6);
    end

    X10 = mvst_getf(opts,'X10', mvst_getf(opts,'x0',[]));
    if ~isempty(X10) && res.feas
        res.Jbound = trace(P0.'*X10);               % <P_s(0), X10>
    end

    % ---- average-cost bound (secondary mode) ----------------------------- %
    if strcmp(mode,'average') && res.feas
        res.rho = avg_rate(res, md, opts, Tmin, Tmax, dwl, ng, P0, Pm);
    end
end

% ======================================================================= %
%  Discrete jump I-GRE (post-jump P^+ -> pre-jump P^-), Schur elimination  %
% ======================================================================= %
function [Pm, Kd] = jump_riccati(Pp, J, Bd, Qd, Rd)
    Pp = (Pp+Pp.')/2;
    if size(Bd,2) == 0
        Pm = J.'*Pp*J + Qd; Pm = (Pm+Pm.')/2; Kd = zeros(0,size(J,1)); return;
    end
    Th = Rd + Bd.'*Pp*Bd; Th = (Th+Th.')/2;
    Sd = Bd.'*Pp*J;
    Kd = -(Th\Sd);
    Pm = J.'*Pp*J + Qd - Sd.'*(Th\Sd); Pm = (Pm+Pm.')/2;
end

function Kd = jump_gain_at(P0, md0, J, Bd, Qd, Rd, D) %#ok<INUSD>
% Causal RDT jump gain: the jump-optimal gain eliminating the jump input by
% its Schur complement uses the POST-jump certificate P_s(0^+)=P0 (the jump
% Bellman gain u_d = -(Rd+Bd'P0 Bd)^{-1} Bd'P0 J x), so it is the same for
% every observed dwell D in [Tmin,Tmax]; the argument D is accepted for the
% causal interface.
    [~, Kd] = jump_riccati(P0, J, Bd, Qd, Rd);
end

% ======================================================================= %
%  Average-cost rate bounds (eq:IMP:AC:dwell:MDT:rho / RDT:rho)            %
% ======================================================================= %
function rho = avg_rate(res, md, opts, Tmin, Tmax, dwl, ng, P0, Pm)
    W    = mvst_getf(opts,'W', []);
    Winf = mvst_getf(opts,'Winf', []);
    Wjp  = mvst_getf(opts,'Wjp', zeros(size(P0)));
    nu   = mvst_getf(opts,'nu', []);
    Pfun = res.Pfun; ip = @(P,W) trace(P.'*W);
    if strcmp(dwl,'min')
        if isempty(nu); nu = 1/Tmin; end            % default constant-dwell rate
        tg = linspace(0,Tmin,ng); vals = zeros(1,ng);
        for i=1:ng; vals(i) = ip(Pfun(tg(i)), evalW(W,tg(i))); end
        flowInt = trapz(tg, vals);
        rho = nu*flowInt + (1-nu*Tmin)*ip(Pm, def(Winf, evalW(W,Tmin))) ...
              + nu*ip(P0, Wjp);
    else
        if isempty(nu); nu = 1/((Tmin+Tmax)/2); end
        mu = mvst_getf(opts,'mu', @(D) 1/(Tmax-Tmin));   % uniform default
        Dg = linspace(Tmin,Tmax,ng); outer = zeros(1,ng);
        for j=1:ng
            D = Dg(j); tg = linspace(0,D,ng); vals = zeros(1,ng);
            for i=1:ng; vals(i) = ip(Pfun(tg(i)), evalW(W,tg(i))); end
            outer(j) = trapz(tg,vals)*mu(D);
        end
        rho = nu*trapz(Dg,outer) + nu*ip(P0, Wjp);
    end
end

% ======================================================================= %
%  Small utilities                                                        %
% ======================================================================= %
function md = evalmd(md, tau)
    if isa(md,'function_handle'); md = md(tau); end
end
function [tt, Pt] = build_cert(md_orig, md0, Dfp, Dwin, Pm)
% Certificate P_s(tau) on [0,Dwin].  On [0,Dfp] it is the backward flow I-GRE
% from the pre-jump value Pm at tau=Dfp (MVST_DRE_CT), using the mode frozen at
% Dfp (const flow data on the binding window).  For range dwell-time (Dwin>Dfp)
% it is continued past Dfp by integrating the forward Riccati dP/dtau = -G(P)
% from Pm over [Dfp,Dwin] using the ACTUAL (possibly timer-varying) flow data
% md_orig(tau) -- so a genuinely non-monotone P_s (e.g. a dipping Z(tau)) is
% represented and the theta-family jump check sees it.  The forward flow is
% anti-stable about P_inf, but the window is short and bounded.
    [tt, Pt] = mvst_dre_ct(md0, Dfp, Pm, 1e-10);          % [0,Dfp]
    if Dwin > Dfp + 1e-12
        n = size(Pm,1);
        odef = @(tau,pv) reshape(-mvst_ric_ct(reshape(pv,n,n), evalmd(md_orig,tau)), [], 1);
        tg = linspace(Dfp, Dwin, 60);
        [ts, Pv] = ode45(odef, tg, Pm(:), odeset('RelTol',1e-9,'AbsTol',1e-11));
        for k = 2:numel(ts)
            Pk = reshape(Pv(k,:), n, n); Pt(:,:,end+1) = (Pk+Pk.')/2; %#ok<AGROW>
        end
        tt = [tt(:); ts(2:end)];
    end
end
function M = interp3d(tt, Pt, tau)
    n = size(Pt,1); M = zeros(n);
    for i=1:n; for j=1:n
        M(i,j) = interp1(tt, squeeze(Pt(i,j,:)), tau, 'pchip');
    end; end
    M = (M+M.')/2;
end
function K = exact_gain(P, md)
    [~, K] = mvst_ric_ct((P+P.')/2, md);
end
function W = evalW(W, tau)
    if isa(W,'function_handle'); W = W(tau); end
    if isempty(W); W = 0; end
end
function y = def(x, d); if isempty(x); y = d; else; y = x; end; end
