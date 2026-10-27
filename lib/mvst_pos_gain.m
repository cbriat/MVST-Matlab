function res = mvst_pos_gain(sys, opts)
%MVST_POS_GAIN  Class-dispatched induced-gain copositive programs for linear
%   positive systems (Chapter "Linear Positive Systems").  Computes the induced
%   L1 x ell1 (trace) gain or the induced Linf x ellinf (peak-to-peak) gain of a
%   stable positive system by the sliced copositive programs of the chapter,
%   routing to the matching per-class backend and returning a uniform result.
%
%   RES = MVST_POS_GAIN(SYS, OPTS).  OPTS.norm selects the gain ('L1'|'Linf');
%   OPTS.class the system class; OPTS.dwell the dwell-time family.
%
%   TWO DUAL PROGRAMS (they are genuinely different, not transposes):
%     'L1'   trace gain, BACKWARD-ADJOINT weight lambda (uses A',J',C',E'), the
%            sliced program of Cor. cor:pos:L1slice (eq:pos:L1LP).
%     'Linf' peak gain, FORWARD-TUBE xi >= 0 (uses A,J,E,C directly), the sliced
%            program of Cor. cor:pos:Linfslice (eq:pos:LinfLP).
%
%   CLASSES (OPTS.class; inferred when omitted) and dwell-time families:
%     'ct'        continuous-time  x' = A x + E w, z = C x + F w  (A
%                 Metzler-Hurwitz; E,C,F >= 0).  Static: closed-form vertex LP,
%                 the induced gain equals the DC gain C(-A)^{-1}E (+F).
%     'dt'        discrete-time  x+ = J x + Ed w, z = Cd x + Fd w  (J >= 0
%                 Schur).  Peak: (J-I) xi + Ed1 <= 0, Cd xi + Fd1 <= gamma1.
%     'impulsive' flow (A,Ec,Cc,Fc) + reset (J,Ed,Cd,Fd), positive impulsive
%                 system eq:pos:perfsys.  OPTS.dwell:
%                   'arbitrary' (default): static common certificate LP (the
%                       arbitrary-dwell-time copositive gain, YALMIP baseline),
%                   'constant'|'min'|'range': the genuinely timer-dependent
%                       certificate lambda(tau)/xi(tau) via the SOSTOOLS
%                       timer-polynomial route (MVST_POS_GAIN_SOS), the sliced
%                       programs eq:pos:L1LP / eq:pos:LinfLP.
%     'switched'  x' = A_sigma x (modes A_i Metzler, identity state map, timer
%                 reset), per-mode (Ec_i,Cc_i,Fc_i).  'arbitrary': common
%                 certificate LP; 'min'|'range': per-mode timer certificate,
%                 switch coupling xi_new(0) >= xi_old(theta) (SOS).
%     'mjls'      Markov jump  x' = A_r x, generator Pi = (pi_ij).  Mean-reading
%                 mode-indexed coupled copositive rows (eq:pos:mjlscop): the
%                 flow row of mode i carries + sum_j pi_ij p_j.
%     'poisson'   x' = A x + Poissonian jumps of rate lambda_r and map J_r; the
%                 flow rows carry the rate term sum_r lambda_r (J_r - I) (Linf,
%                 primal) / sum_r lambda_r (J_r' - I) (L1, adjoint).
%
%   SYS FIELDS (thin or rich naming accepted):
%     .A/.A0   drift (n x n; a 1xN cell for 'switched'/'mjls'; a function handle
%              @(tau) for a timer-varying Metzler flow in the SOS dwell modes)
%     .E/.Bw/.Ec  disturbance input (n x mc; per-mode cell for switched/mjls)
%     .C/.Cc      performance output (p x n; per-mode cell)
%     .F/.Dw/.Fc  feedthrough (p x mc; default 0)
%     .J,.Ed,.Cd,.Fd  reset / discrete data ('dt','impulsive')
%     .Pi/.Pig    generator ('mjls');  .Jf (cell) + .lamf ('poisson')
%     .Tbar ('constant'), .Tmin,.Tmax ('min','range'), .trans ('switched')
%
%   OPTS FIELDS:
%     .norm 'L1'|'Linf' (required); .class; .dwell; .Tbar/.Tmin/.Tmax;
%     .deg (SOS timer degree, default 4); .eps (SOS positivity floor);
%     .solver (SOS backend); plus any field consumed by MVST_POS_LP.
%
%   RES FIELDS:
%     .gamma   the induced gain (min gamma); [] if infeasible.
%     .weights the copositive certificate (lambda for 'L1', xi for 'Linf'): an
%              n-vector (LP branches), a 1xN cell (mjls/switched), or a function
%              handle @(tau) (SOS dwell branches).
%     .feas    logical feasibility.
%     .norm,.class,.dwell,.method  the resolved tags.
%     .lp      LP data (LINPROG branches);  .info (SOS branches).
%     .gammaDC the DC gain C(-A)^{-1}E (+F) for 'ct' (a closed-form cross-check).
%
%   Example (compartmental L1 gain, Ex. ex:pos:gain):
%     r = mvst_pos_gain(struct('A',[-0.9 0;0.5 -0.6],'E',[1;0],'C',[1 1]), ...
%                       struct('norm','L1'));   % r.gamma == 2.0370...
%
%   See also MVST_POS_GAIN_SOS, MVST_POS_STABILITY, MVST_POS_STABILIZE,
%   MVST_POS_LP, MVST_POS.

    arguments
        sys  (1,1) struct
        opts (1,1) struct = struct()
    end

    nrm = lower(mvst_getf(opts, 'norm', ''));
    if ~any(strcmp(nrm, {'l1','linf'}))
        error('mvst:pos:norm', 'OPTS.norm must be ''L1'' or ''Linf''.');
    end
    cls = resolve_class(sys, opts);
    dwl = lower(mvst_getf(opts, 'dwell', 'arbitrary'));

    res = struct('gamma', [], 'weights', [], 'feas', false, 'norm', nrm, ...
                 'class', cls, 'dwell', dwl, 'method', '', 'lp', [], ...
                 'info', [], 'gammaDC', []);

    % ---- genuinely timer-dependent dwell families -> SOS engine ----------- %
    if any(strcmp(cls, {'impulsive','switched'})) && ...
       any(strcmp(dwl, {'constant','min','range'}))
        o = opts; o.norm = nrm; o.class = cls; o.dwell = dwl;
        r = mvst_pos_gain_sos(sys, o);
        res.gamma = r.gamma; res.weights = r.weights; res.feas = r.feas;
        res.info  = r.info;  res.method = ['sos:' cls ':' dwl];
        return;
    end

    % ---- static / arbitrary LINPROG branches ------------------------------ %
    switch cls
        case 'ct';        res = branch_ct(sys, nrm, opts, res);
        case 'dt';        res = branch_dt(sys, nrm, opts, res);
        case 'impulsive'; res = branch_imp_arb(sys, nrm, opts, res);
        case 'switched';  res = branch_switched_arb(sys, nrm, opts, res);
        case 'mjls';      res = branch_mjls(sys, nrm, opts, res);
        case 'poisson';   res = branch_poisson(sys, nrm, opts, res);
        otherwise
            error('mvst:pos:class', 'Unknown positive-gain class ''%s''.', cls);
    end
    res.method = ['lp:' cls];
end

% ======================================================================= %
%  Class inference                                                         %
% ======================================================================= %
function cls = resolve_class(sys, opts)
    cls = lower(mvst_getf(opts, 'class', ''));
    if ~isempty(cls); return; end
    if ~isempty(mvst_getf(sys, 'Pi', mvst_getf(sys, 'Pig', [])))
        cls = 'mjls';
    elseif ~isempty(mvst_getf(sys, 'Jf', [])) || ~isempty(mvst_getf(sys, 'lamf', []))
        cls = 'poisson';
    elseif iscell(mvst_getf(sys, 'A', mvst_getf(sys, 'A0', [])))
        cls = 'switched';
    elseif ~isempty(mvst_getf(sys, 'J', []))
        cls = 'impulsive';
    else
        cls = 'ct';
    end
end

% ======================================================================= %
%  ct : continuous-time (closed-form vertex LP, DC-gain cross-check)       %
% ======================================================================= %
function res = branch_ct(sys, nrm, opts, res)
    A = getA1(sys); E = getE(sys); C = getC(sys);
    n = size(A,1); m = size(E,2); p = size(C,1);
    F = getF(sys, p, m);
    if strcmp(nrm,'l1')
        [Aub,bub,lb] = lp_l1_flow(A, E, C, F, n, m, p);
        o = mvst_pos_lp(struct('f',objf(n),'Aub',Aub,'bub',bub,'lb',lb), opts);
        res = finish_single(res, o, n, ...
              @() max(E.'*((-A.')\(C.'*ones(p,1))) + F.'*ones(p,1)));
    else
        [Aub,bub,lb] = lp_linf_flow(A, E, C, F, n, m, p);
        o = mvst_pos_lp(struct('f',objf(n),'Aub',Aub,'bub',bub,'lb',lb), opts);
        res = finish_single(res, o, n, ...
              @() max(C*((-A)\(E*ones(m,1))) + F*ones(m,1)));
    end
end

% ======================================================================= %
%  dt : discrete-time pure jump                                            %
% ======================================================================= %
function res = branch_dt(sys, nrm, opts, res)
    J  = mvst_getf(sys,'J', mvst_getf(sys,'A', mvst_getf(sys,'A0',[])));
    if isempty(J); error('mvst:pos:dt','Class ''dt'' needs a transition .J (or .A).'); end
    n  = size(J,1);
    Ed = mvst_getf(sys,'Ed', mvst_getf(sys,'E', mvst_getf(sys,'Ec',[])));
    if isempty(Ed); error('mvst:pos:dt','Class ''dt'' needs a disturbance .Ed.'); end
    md = size(Ed,2);
    Cd = mvst_getf(sys,'Cd', mvst_getf(sys,'C', mvst_getf(sys,'Cc',[])));
    if isempty(Cd); error('mvst:pos:dt','Class ''dt'' needs an output .Cd.'); end
    qd = size(Cd,1);
    Fd = mvst_getf(sys,'Fd', zeros(qd,md));
    if isscalar(Fd) && (qd~=1 || md~=1); Fd = Fd*ones(qd,md); end
    if strcmp(nrm,'l1')
        % (J'-I) lam + Cd'1 <= 0 ; Ed' lam + Fd'1 <= gamma1
        Aub = [ (J.'-eye(n)), zeros(n,1) ; Ed.', -ones(md,1) ];
        bub = [ -Cd.'*ones(qd,1) ; -Fd.'*ones(qd,1) ];
    else
        % (J-I) xi + Ed1 <= 0 ; Cd xi + Fd1 <= gamma1
        Aub = [ (J-eye(n)), zeros(n,1) ; Cd, -ones(qd,1) ];
        bub = [ -Ed*ones(md,1) ; -Fd*ones(md,1) ];
    end
    lb = [zeros(n,1); -inf];
    o  = mvst_pos_lp(struct('f',objf(n),'Aub',Aub,'bub',bub,'lb',lb), opts);
    res = finish_single(res, o, n, []);
end

% ======================================================================= %
%  impulsive, arbitrary dwell-time (static common certificate)            %
% ======================================================================= %
function res = branch_imp_arb(sys, nrm, opts, res)
    A = getA1(sys); Ec = getE(sys); Cc = getC(sys);
    n = size(A,1); mc = size(Ec,2); qc = size(Cc,1);
    Fc = getF(sys, qc, mc);
    J  = mvst_getf(sys,'J', zeros(n));
    Ed = mvst_getf(sys,'Ed', zeros(n,0)); md = size(Ed,2);
    Cd = mvst_getf(sys,'Cd', Cc);          qd = size(Cd,1);
    Fd = mvst_getf(sys,'Fd', zeros(qd,md));
    if strcmp(nrm,'l1')
        Aub = [ A.',            zeros(n,1) ;      % A'lam + Cc'1 <= 0
                (J.'-eye(n)),   zeros(n,1) ;      % (J'-I)lam + Cd'1 <= 0
                Ec.',           -ones(mc,1) ;     % Ec'lam + Fc'1 <= g1
                Ed.',           -ones(md,1) ];    % Ed'lam + Fd'1 <= g1
        bub = [ -Cc.'*ones(qc,1) ; -Cd.'*ones(qd,1) ; ...
                -Fc.'*ones(qc,1) ; -Fd.'*ones(qd,1) ];
    else
        Aub = [ A,             zeros(n,1) ;       % A xi + Ec1 <= 0
                (J-eye(n)),    zeros(n,1) ;       % (J-I)xi + Ed1 <= 0
                Cc,            -ones(qc,1) ;      % Cc xi + Fc1 <= g1
                Cd,            -ones(qd,1) ];     % Cd xi + Fd1 <= g1
        bub = [ -Ec*ones(mc,1) ; -Ed*ones(md,1) ; ...
                -Fc*ones(mc,1) ; -Fd*ones(md,1) ];
    end
    lb = [zeros(n,1); -inf];
    o  = mvst_pos_lp(struct('f',objf(n),'Aub',Aub,'bub',bub,'lb',lb), opts);
    res = finish_single(res, o, n, []);
end

% ======================================================================= %
%  switched, arbitrary switching (common certificate over the modes)      %
% ======================================================================= %
function res = branch_switched_arb(sys, nrm, opts, res)
    Ac = ascell(getAcell(sys)); N = numel(Ac); n = size(Ac{1},1);
    Ec = ascellN(getE(sys), N); Cc = ascellN(getC(sys), N);
    Fc = ascellN(mvst_getf(sys,'F',mvst_getf(sys,'Fc',[])), N);
    Aub = []; bub = [];
    for i = 1:N
        [Ai,bi] = mode_gain_rows(Ac{i}, Ec{i}, Cc{i}, Fc{i}, nrm, n);
        Aub = [Aub; Ai]; bub = [bub; bi]; %#ok<AGROW>
    end
    lb = [zeros(n,1); -inf];
    o  = mvst_pos_lp(struct('f',objf(n),'Aub',Aub,'bub',bub,'lb',lb), opts);
    res = finish_single(res, o, n, []);
end

% ======================================================================= %
%  mjls : mean-reading coupled copositive rows (+ sum_j pi_ij p_j)        %
% ======================================================================= %
function res = branch_mjls(sys, nrm, opts, res)
    Ac = ascell(getAcell(sys)); N = numel(Ac); n = size(Ac{1},1);
    Pi = mvst_getf(sys,'Pi', mvst_getf(sys,'Pig', []));
    if isempty(Pi); error('mvst:pos:mjls','Class ''mjls'' needs generator .Pi.'); end
    Ec = ascellN(getE(sys), N); Cc = ascellN(getC(sys), N);
    Fc = ascellN(mvst_getf(sys,'F',mvst_getf(sys,'Fc',[])), N);
    nv = N*n + 1; gcol = nv;             % [p_1;...;p_N; gamma]
    Aub = []; bub = [];
    for i = 1:N
        Eci = Ec{i}; Cci = Cc{i}; mc = size(Eci,2); qc = size(Cci,1);
        if isempty(Fc{i}); Fci = zeros(qc,mc); else; Fci = Fc{i}; end
        % flow row block (n rows): [.. A_i(') at i, pi_ij I at j ..] p + out'1 <= 0
        rf = zeros(n, nv);
        if strcmp(nrm,'l1'); rf(:, blk(i,n)) = Ac{i}.'; else; rf(:, blk(i,n)) = Ac{i}; end
        % Generator coupling: L1 (backward adjoint) uses pi_ij; Linf (forward
        % mean tube) uses pi_ji -- mass flows INTO mode i from j at rate pi_ji.
        for j = 1:N
            if strcmp(nrm,'l1'); cpl = Pi(i,j); else; cpl = Pi(j,i); end
            rf(:, blk(j,n)) = rf(:, blk(j,n)) + cpl*eye(n);
        end
        if strcmp(nrm,'l1'); bf = -Cci.'*ones(qc,1); else; bf = -Eci*ones(mc,1); end
        % output row block
        if strcmp(nrm,'l1')
            ro = zeros(mc, nv); ro(:, blk(i,n)) = Eci.'; ro(:,gcol) = -ones(mc,1);
            bo = -Fci.'*ones(qc,1);
        else
            ro = zeros(qc, nv); ro(:, blk(i,n)) = Cci;   ro(:,gcol) = -ones(qc,1);
            bo = -Fci*ones(mc,1);
        end
        Aub = [Aub; rf; ro]; bub = [bub; bf; bo]; %#ok<AGROW>
    end
    lb = [zeros(N*n,1); -inf];
    o  = mvst_pos_lp(struct('f',[zeros(N*n,1);1],'Aub',Aub,'bub',bub,'lb',lb), opts);
    res.lp = o.lp;
    if o.feas
        res.feas = true; res.gamma = o.x(end);
        res.weights = splitp(o.x(1:N*n), N, n);
    end
end

% ======================================================================= %
%  poisson : rate terms sum_r lambda_r (J_r - I) in the flow rows         %
% ======================================================================= %
function res = branch_poisson(sys, nrm, opts, res)
    A = getA1(sys); Ec = getE(sys); Cc = getC(sys);
    n = size(A,1); mc = size(Ec,2); qc = size(Cc,1);
    Fc = getF(sys, qc, mc);
    Jc  = tocell(mvst_getf(sys,'Jf', mvst_getf(sys,'Jr', {})));
    lam = mvst_getf(sys,'lamf', mvst_getf(sys,'lam', []));
    R = A;                                       % flow generator incl. rate terms
    for r = 1:numel(Jc); R = R + lam(r)*(Jc{r} - eye(n)); end
    if strcmp(nrm,'l1')
        Aub = [ R.', zeros(n,1) ; Ec.', -ones(mc,1) ];
        bub = [ -Cc.'*ones(qc,1) ; -Fc.'*ones(qc,1) ];
    else
        Aub = [ R, zeros(n,1) ; Cc, -ones(qc,1) ];
        bub = [ -Ec*ones(mc,1) ; -Fc*ones(mc,1) ];
    end
    lb = [zeros(n,1); -inf];
    o  = mvst_pos_lp(struct('f',objf(n),'Aub',Aub,'bub',bub,'lb',lb), opts);
    res = finish_single(res, o, n, []);
end

% ======================================================================= %
%  Shared row builders / utilities                                        %
% ======================================================================= %
function [Aub,bub,lb] = lp_l1_flow(A, E, C, F, n, m, p)
    Aub = [ A.', zeros(n,1) ; E.', -ones(m,1) ];
    bub = [ -C.'*ones(p,1) ; -F.'*ones(p,1) ];
    lb  = [ zeros(n,1) ; -inf ];
end
function [Aub,bub,lb] = lp_linf_flow(A, E, C, F, n, m, p)
    Aub = [ A, zeros(n,1) ; C, -ones(p,1) ];
    bub = [ -E*ones(m,1) ; -F*ones(m,1) ];
    lb  = [ zeros(n,1) ; -inf ];
end
function [Ai,bi] = mode_gain_rows(A, Ec, Cc, Fc, nrm, n)
    mc = size(Ec,2); qc = size(Cc,1);
    if isempty(Fc); Fc = zeros(qc,mc); end
    if strcmp(nrm,'l1')
        Ai = [ A.', zeros(n,1) ; Ec.', -ones(mc,1) ];
        bi = [ -Cc.'*ones(qc,1) ; -Fc.'*ones(qc,1) ];
    else
        Ai = [ A, zeros(n,1) ; Cc, -ones(qc,1) ];
        bi = [ -Ec*ones(mc,1) ; -Fc*ones(mc,1) ];
    end
end

function res = finish_single(res, o, n, dcfun)
    res.lp = o.lp;
    if o.feas
        res.feas = true; res.gamma = o.x(end); res.weights = o.x(1:n);
        if ~isempty(dcfun); res.gammaDC = dcfun(); end
    end
end

function f = objf(n); f = [zeros(n,1); 1]; end
function idx = blk(i,n); idx = (i-1)*n + (1:n); end
function A = getA1(sys)
    A = mvst_getf(sys,'A', mvst_getf(sys,'A0',[]));
    if iscell(A); A = A{1}; end
    if isempty(A); error('mvst:pos:noA','SYS must carry a drift .A.'); end
end
function A = getAcell(sys); A = mvst_getf(sys,'A', mvst_getf(sys,'A0',[])); end
function E = getE(sys)
    E = mvst_getf(sys,'E', mvst_getf(sys,'Ec', mvst_getf(sys,'Bw',[])));
    if isempty(E); error('mvst:pos:noE','SYS must carry a disturbance .E/.Ec.'); end
end
function C = getC(sys)
    C = mvst_getf(sys,'C', mvst_getf(sys,'Cc',[]));
    if isempty(C); error('mvst:pos:noC','SYS must carry an output .C/.Cc.'); end
end
function F = getF(sys, p, m)
    F = mvst_getf(sys,'F', mvst_getf(sys,'Fc', mvst_getf(sys,'Dw', zeros(p,m))));
    if isscalar(F) && (p~=1 || m~=1); F = F*ones(p,m); end
end
function P = splitp(x, N, n); P = cell(1,N); for i=1:N; P{i}=x(blk(i,n)); end; end
function c = ascell(A); if iscell(A); c=A; else; c={A}; end; end
function c = ascellN(X, N); if iscell(X); c=X; else; c=repmat({X},1,N); end; end
function c = tocell(X); if iscell(X); c=X; elseif isempty(X); c={}; else; c={X}; end; end
