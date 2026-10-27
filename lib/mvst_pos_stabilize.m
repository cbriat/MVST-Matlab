function res = mvst_pos_stabilize(sys, opts)
%MVST_POS_STABILIZE  Positivity-preserving state-feedback stabilization by
%   copositive linear programming (Chapter "Linear Positive Systems" of the
%   monograph).  Designs a state feedback that renders the closed loop BOTH
%   positive and exponentially stable, via the copositive change of variables
%   X = diag(lambda), U = K X of Lemma lem:pos:cov, so the design conditions
%   are linear programs in (lambda, U).
%
%   RES = MVST_POS_STABILIZE(SYS, OPTS) selects a class with OPTS.class.
%
%   SUPPORTED CLASSES (OPTS.class):
%     'ct'        continuous-time LTI positive system  x' = A x + Bc u.
%                 Implements the flow conditions of Thm th:pos:stabzfix
%                 (the static, arbitrary-dwell-time positive-system synthesis;
%                 vector level of Aitrami:07): find lambda>0, Uc such that
%                   (A X + Bc Uc)_{ij} >= 0  for i ~= j   (closed loop Metzler),
%                   (A X + Bc Uc) 1 <= -eps 1             (copositive decrease),
%                 with X = diag(lambda).  The gain is Kc = Uc X^{-1} and the
%                 closed loop A + Bc Kc is Metzler-Hurwitz.
%     'impulsive' impulsive LTI positive system with reset, constant
%                 dwell-time Tbar.  Implements Thm th:pos:stabzfix(b),
%                 eq:pos:stabzconst, with a constant (timer-independent)
%                 diagonal certificate X (Xdot = 0):
%                   (A X + Bc Uc)_{ij} >= 0  (i ~= j),
%                   (J X + Bd Ud)      >= 0  (entrywise),
%                   (A X + Bc Uc) 1        <= -eps 1,
%                   (J X + Bd Ud - X) 1    <= -eps 1,
%                 giving flow gain Kc = Uc X^{-1} and jump gain Kd = Ud X^{-1},
%                 closed loop positive and stable under constant dwell-time.
%
%   The timer-dependent minimum / range dwell-time families
%   (Prop. prop:pos:stabzdwell) need a piecewise-C1 certificate X(tau) over
%   the dwell window and are not implemented here (their finite-dimensional
%   reduction requires a timer discretization outside the book's statements);
%   the constant-dwell-time class above is the exact finite-dimensional case.
%
%   SYS FIELDS:
%     .A            Metzler-controllable drift            (n x n)
%     .Bc or .B     flow input selector                   (n x mc)
%     .J            reset map (class 'impulsive')         (n x n)
%     .Bd           jump input selector (class 'impulsive'; default Bc)
%                                                          (n x md)
%     .Tbar         constant dwell-time (class 'impulsive'; or OPTS.Tbar)
%
%   OPTS FIELDS:
%     .class    'ct' | 'impulsive' (required; inferred: SYS.J -> 'impulsive',
%              else 'ct').
%     .eps     copositive decrease margin (default 1e-3).
%     .pmin    lower bound on the diagonal weight lambda (default 1; sets c1).
%     .pmax    upper bound on lambda (default [] = unbounded; sets c2).
%     .Tbar    constant dwell-time (class 'impulsive').
%     .certify logical; re-certify the closed loop with MVST_POS_STABILITY
%              (default true).
%     plus any field consumed by MVST_POS_LP (.solver, .lpalgo, .lpopts).
%
%   RES FIELDS:
%     .K        the stabilizing flow gain Kc (mc x n); [] if infeasible.
%     .Kd       the stabilizing jump gain Kd (class 'impulsive'); else [].
%     .lambda   the diagonal certificate weight lambda (n x 1); [] if infeas.
%     .lp       the LP data actually solved (LINPROG naming; see MVST_POS_LP).
%     .feas     logical LP feasibility.
%     .cert     the MVST_POS_STABILITY result certifying the closed loop
%               (present when OPTS.certify=true and the design is feasible).
%     .class    the resolved class tag.
%
%   Example (unstable drift made Metzler-Hurwitz by full-actuation feedback):
%     sys = struct('A',[1 2; 0 1],'Bc',eye(2));
%     r = mvst_pos_stabilize(sys, struct('class','ct'));
%     % r.feas == true; A + Bc*r.K is Metzler and Hurwitz; r.cert.stable == true.
%   (A single shared input column Bc=[1;1] cannot make this A Metzler AND
%   Hurwitz, so full or adequate actuation is required for feasibility.)
%
%   See also MVST_POS_STABILITY, MVST_POS_GAIN, MVST_POS_SUPPLY, MVST_POS_LP.

    arguments
        sys  (1,1) struct
        opts (1,1) struct = struct()
    end

    epsm = mvst_getf(opts, 'eps',  1e-3);
    pmin = mvst_getf(opts, 'pmin', 1);
    pmax = mvst_getf(opts, 'pmax', []);

    A = mvst_getf(sys, 'A', []);
    if isempty(A); error('mvst:pos:noA', 'SYS must carry a drift .A.'); end
    Bc = mvst_getf(sys, 'Bc', mvst_getf(sys, 'B', []));
    n  = size(A,1);
    if isempty(Bc)
        % L1 output-injection synthesis (opts.norm='L1') needs no state-feedback
        % channel; every other path requires .Bc.
        if strcmpi(mvst_getf(opts,'norm',''),'L1'); Bc = zeros(n,0);
        else; error('mvst:pos:noB', 'SYS must carry a flow input .Bc.'); end
    end
    mc = size(Bc,2);

    cls = mvst_getf(opts, 'class', '');
    if isempty(cls)
        if ~isempty(mvst_getf(sys,'J',[])); cls = 'impulsive'; else; cls = 'ct'; end
    end
    cls = lower(cls);

    % ---- genuinely timer-dependent dwell families -> SOS synthesis -------- %
    dwl = lower(mvst_getf(opts, 'dwell', 'arbitrary'));
    if strcmp(cls,'impulsive') && any(strcmp(dwl,{'constant','min','range'}))
        r = mvst_pos_stabilize_sos(sys, opts);
        res = struct('K', r.Kc, 'Kd', r.Kd, 'lambda', r.lambda, 'lp', [], ...
                     'feas', r.feas, 'cert', [], 'class', cls, 'dwell', dwl, ...
                     'gamma', r.gamma, 'slack', r.slack, 'info', r.info, ...
                     'Uc', r.Uc, 'Ud', r.Ud);
        return;
    end

    % ---- static gain-optimal synthesis (opts.norm), continuous-time slice -- %
    % Linf  : L-infinity-optimal STATE feedback (primal tube LP, prop:pos:ebihara
    %         dual / th:pos:stabzfix with the gain minimized).
    % L1    : L1-optimal OUTPUT INJECTION, computed by the design duality
    %         prop:pos:Linfdesign -- L1-OI on (A,Cinj,Ew,Cz) equals L-infinity-
    %         optimal STATE feedback on the transposed data (A',Cinj',Cz',Ew'),
    %         with the injection gain recovered as L = K'.
    nrm = lower(mvst_getf(opts,'norm',''));
    if strcmp(cls,'ct') && any(strcmp(nrm,{'linf','l1'}))
        res = pos_gainopt_static(sys, A, Bc, nrm, opts);
        return;
    end

    res = struct('K', [], 'Kd', [], 'lambda', [], 'lp', [], ...
                 'feas', false, 'cert', [], 'class', cls);

    switch cls
        case 'ct'
            % variables x = [lambda(n); uc(mc*n)]
            nv = n + mc*n;
            [Aub, bub] = flow_rows(A, Bc, n, mc, nv, epsm);
            lb = [pmin*ones(n,1); -inf(mc*n,1)];
            ub = ubvec(pmax, n, mc*n);
            lp = struct('nvar', nv, 'Aub', Aub, 'bub', bub, 'lb', lb, 'ub', ub);
            o  = mvst_pos_lp(lp, opts);
            res.lp = o.lp;
            if o.feas
                lam = o.x(1:n);
                Uc  = reshape(o.x(n+(1:mc*n)), mc, n);
                res.feas   = true;
                res.lambda = lam;
                res.K      = Uc / diag(lam);
            end

        case 'impulsive'
            J = mvst_getf(sys, 'J', []);
            if isempty(J); error('mvst:pos:impulsive', ...
                    'Class ''impulsive'' requires the reset SYS.J.'); end
            Bd = mvst_getf(sys, 'Bd', Bc);
            md = size(Bd,2);
            Tbar = mvst_getf(opts, 'Tbar', mvst_getf(sys, 'Tbar', []));
            if isempty(Tbar); error('mvst:pos:impulsive', ...
                    'Class ''impulsive'' requires the dwell-time Tbar.'); end
            % variables x = [lambda(n); uc(mc*n); ud(md*n)]
            nv = n + mc*n + md*n;
            [Af, bf] = flow_rows(A, Bc, n, mc, nv, epsm);
            [Aj, bj] = jump_rows(J, Bd, n, mc, md, nv, epsm);
            Aub = [Af; Aj]; bub = [bf; bj];
            lb = [pmin*ones(n,1); -inf(mc*n+md*n,1)];
            ub = ubvec(pmax, n, mc*n+md*n);
            lp = struct('nvar', nv, 'Aub', Aub, 'bub', bub, 'lb', lb, 'ub', ub);
            o  = mvst_pos_lp(lp, opts);
            res.lp = o.lp;
            if o.feas
                lam = o.x(1:n);
                Uc  = reshape(o.x(n+(1:mc*n)), mc, n);
                Ud  = reshape(o.x(n+mc*n+(1:md*n)), md, n);
                res.feas   = true;
                res.lambda = lam;
                res.K      = Uc / diag(lam);
                res.Kd     = Ud / diag(lam);
            end

        otherwise
            error('mvst:pos:class', ...
                  'OPTS.class must be ''ct'' or ''impulsive'' (got ''%s'').', cls);
    end

    % ---- certify the closed loop -------------------------------------- %
    if res.feas && mvst_getf(opts, 'certify', true)
        switch cls
            case 'ct'
                Acl = A + Bc*res.K;
                res.cert = mvst_pos_stability(struct('A', Acl), ...
                                              struct('class', 'ct'));
            case 'impulsive'
                Acl = A + Bc*res.K;
                Jcl = mvst_getf(sys,'J',[]) + mvst_getf(sys,'Bd',Bc)*res.Kd;
                res.cert = mvst_pos_stability( ...
                    struct('A', Acl, 'J', Jcl, 'Tbar', ...
                           mvst_getf(opts,'Tbar',mvst_getf(sys,'Tbar',[]))), ...
                    struct('class', 'impulsive'));
        end
    end
end

% ===================================================================== %
%  Row builders (linear in [lambda; uc; ...])                            %
% ===================================================================== %
function [Aub, bub] = flow_rows(A, Bc, n, mc, nv, epsm)
% Off-diagonal Metzler rows  -(A X + Bc Uc)_{ij} <= 0  (i ~= j)
% and copositive decrease rows  [(A X + Bc Uc) 1]_i <= -eps.
    % Metzler: n*(n-1) rows
    Am = zeros(n*(n-1), nv); r = 0;
    for i = 1:n
        for j = 1:n
            if i == j; continue; end
            r = r + 1;
            row = zeros(1, nv);
            row(j) = A(i,j);                       % coeff of lambda_j in (AX)_ij
            for s = 1:mc
                row(n + (j-1)*mc + s) = Bc(i,s);   % coeff of uc(s,j)
            end
            Am(r,:) = -row;                        % >= 0  ->  -row <= 0
        end
    end
    bm = zeros(n*(n-1), 1);
    % Decrease: n rows, row i = sum_j (AX+BcUc)_ij
    Ad = zeros(n, nv);
    for i = 1:n
        row = zeros(1, nv);
        for j = 1:n
            row(j) = row(j) + A(i,j);
            for s = 1:mc
                row(n + (j-1)*mc + s) = row(n + (j-1)*mc + s) + Bc(i,s);
            end
        end
        Ad(i,:) = row;
    end
    bd = -epsm*ones(n,1);
    Aub = [Am; Ad]; bub = [bm; bd];
end

function [Aub, bub] = jump_rows(J, Bd, n, mc, md, nv, epsm)
% All-entry nonnegativity  -(J X + Bd Ud)_{ij} <= 0
% and jump decrease  [(J X + Bd Ud - X) 1]_i <= -eps.
    off = n + mc*n;                                % ud starts here
    % Nonnegativity: n*n rows
    An = zeros(n*n, nv); r = 0;
    for i = 1:n
        for j = 1:n
            r = r + 1;
            row = zeros(1, nv);
            row(j) = J(i,j);                       % coeff of lambda_j in (JX)_ij
            for s = 1:md
                row(off + (j-1)*md + s) = Bd(i,s); % coeff of ud(s,j)
            end
            An(r,:) = -row;
        end
    end
    bn = zeros(n*n, 1);
    % Decrease: row i = sum_j (JX+BdUd)_ij - lambda_i
    Ad = zeros(n, nv);
    for i = 1:n
        row = zeros(1, nv);
        for j = 1:n
            row(j) = row(j) + J(i,j);
            for s = 1:md
                row(off + (j-1)*md + s) = row(off + (j-1)*md + s) + Bd(i,s);
            end
        end
        row(i) = row(i) - 1;                       % - lambda_i
        Ad(i,:) = row;
    end
    bd = -epsm*ones(n,1);
    Aub = [An; Ad]; bub = [bn; bd];
end

function ub = ubvec(pmax, n, ntail)
    if isempty(pmax)
        ub = [];
    else
        ub = [pmax*ones(n,1); inf(ntail,1)];
    end
end

% ======================================================================= %
%  Static gain-optimal synthesis: Linf state feedback and,  %
%  by the design duality prop:pos:Linfdesign, L1 output injection.         %
% ======================================================================= %
function res = pos_gainopt_static(sys, A, Bc, nrm, opts)
    n = size(A,1);
    Ec = mvst_getf(sys,'Ec', mvst_getf(sys,'Bw', mvst_getf(sys,'E', zeros(n,0))));
    C  = mvst_getf(sys,'Cc', mvst_getf(sys,'C',  zeros(0,n)));
    Fc = mvst_getf(sys,'Fc', mvst_getf(sys,'F',  zeros(size(C,1),size(Ec,2))));
    res = struct('K',[],'L',[],'lambda',[],'gamma',[],'feas',false, ...
                 'class','ct','norm',nrm,'Acl',[]);
    if strcmp(nrm,'linf')
        % L-infinity-optimal state feedback on (A,Bc,Ec,C).
        [g,K,lam,ok] = linf_sf(A,Bc,Ec,C,Fc,opts);
        res.feas=ok; res.gamma=g; res.K=K; res.lambda=lam; res.Acl=A+Bc*K;
    else
        % L1-optimal output injection on (A,Cinj,Ew,Cz):
        %   Cinj = measured output (.Cinj / .C), Ew = disturbance (.Bw/.Ec/.E),
        %   Cz = performance output (.Cz / .C).  By prop:pos:Linfdesign this is
        %   L-infinity-optimal state feedback on (A',Cinj',Cz',Ew'); L = K'.
        Cinj = mvst_getf(sys,'Cinj', C);
        Cz   = mvst_getf(sys,'Cz', C);
        Ew   = Ec;
        Fz   = mvst_getf(sys,'Fz', zeros(size(Cz,1), size(Ew,2)));
        [g,K,lam,ok] = linf_sf(A.', Cinj.', Cz.', Ew.', Fz.', opts);
        res.feas=ok; res.gamma=g; res.L=K.'; res.lambda=lam;
        if ok; res.Acl = A + res.L*Cinj; end
    end
end

function [g,K,lam,ok] = linf_sf(A,Bc,Ec,C,Fc,opts)
% L-infinity-optimal state feedback via the primal reachability-tube LP:
%   min gamma  s.t.  lambda>0, Uc, alpha>=0,
%     (A diag(lambda)+Bc Uc)_{ij} >= 0  (i~=j)   [closed loop Metzler]
%     (A diag(lambda)+Bc Uc) 1 + Ec 1 <= 0        [tube supersolution]
%     C lambda + Fc 1 <= gamma                     [peak output bound]
%   K = Uc diag(lambda)^{-1}  (lem:pos:cov).
    n=size(A,1); mc=size(Bc,2); qc=size(C,1); mE=size(Ec,2);
    nv=n+mc*n+1; gi=nv; Aub=[]; bub=[];
    for i=1:n; for j=1:n; if i==j; continue; end
        row=zeros(1,nv); row(j)=A(i,j);
        for s=1:mc; row(n+(j-1)*mc+s)=Bc(i,s); end
        Aub=[Aub;-row]; bub=[bub;0]; %#ok<AGROW>
    end; end
    for i=1:n
        row=zeros(1,nv);
        for j=1:n; row(j)=row(j)+A(i,j);
            for s=1:mc; row(n+(j-1)*mc+s)=row(n+(j-1)*mc+s)+Bc(i,s); end
        end
        Aub=[Aub;row]; bub=[bub;-sum(Ec(i,:))]; %#ok<AGROW>
    end
    for i=1:qc
        row=zeros(1,nv); for j=1:n; row(j)=C(i,j); end; row(gi)=-1;
        Aub=[Aub;row]; bub=[bub;-sum(Fc(i,:))]; %#ok<AGROW>
    end
    lb=[1e-6*ones(n,1); -inf(mc*n,1); 0]; f=[zeros(nv-1,1);1];
    o=mvst_pos_lp(struct('f',f,'Aub',Aub,'bub',bub,'lb',lb), opts);
    ok=o.feas; g=[]; K=[]; lam=[];
    if ok; x=o.x; lam=x(1:n); Uc=reshape(x(n+(1:mc*n)),mc,n); K=Uc/diag(lam); g=x(gi); end
end
