function res = mvst_steer(sys, Xf, T, opts)
%MVST_STEER  Class-dispatched finite-time optimal covariance steering.
%   RES = MVST_STEER(SYS,XF,T) steers the state (second-moment) covariance of
%   the controlled system SYS from its initial covariance to the target XF over
%   the horizon [0,T], by the forward-backward (FB) algorithm: for a symmetric
%   terminal multiplier M it runs the BACKWARD I-GRE from P(T)=M to obtain the
%   time-varying gains, propagates the closed-loop covariance FORWARD, and finds
%   the M for which Sigma(T;M)=XF by the damped Newton core MVST_STEER_NEWTON.
%   RES = MVST_STEER(SYS,XF,T,OPTS) passes options (see below).
%
%   This is the class-dispatched wrapper around the shared Newton core
%   MVST_STEER_NEWTON.  It absorbs the per-class terminal-gap closures and
%   forward covariance propagators that previously lived inline in the steering
%   example scripts, and retires the dead single-class MVST_STEER_CT.
%
%   CLASS SELECTION.  OPTS.class chooses the system class; when omitted it is
%   inferred from SYS.  Supported classes and the SYS layout each expects:
%
%     'ct'       single-mode continuous-time diffusion.  SYS is an MVST_MODE
%                or a single-mode MVST_SYSTEM.  Backward: MVST_DRE_CT.
%     'mjls'     Markov-jump diffusion.  SYS is a multi-mode MVST_SYSTEM
%                (.modes,.Pi,.pi0).  Backward: the coupled MVST_SOLVE_FINITE
%                (or MVST_SOLVE_FINITE_TV when the modes are time-varying).
%                Forward: the coupled per-mode moment equation, aggregate
%                Sigma = sum_i X_i started from pi0(i)*X0.
%     'switched' dwell-time switched diffusion.  SYS = struct('modes',{modes},
%                'Tdwell',Td).  Backward: the mode-dependent I-GRE integrated
%                through the schedule (MVST_RIC_CT / MVST_RIC_CT_TV).
%     'sd'       sampled-data (zero-order-hold) control.  SYS = struct('A','B',
%                'Ac','Bc','Q','R','Tbar').  Backward: the impulsive I-GRE on
%                the augmented state z=[x;u]; forward: augmented covariance with
%                sample-and-hold resets.  A time-varying flow is supplied by
%                making SYS.A a function handle @(t).
%     'poisson'  GBM + Poissonian-jump diffusion.  SYS = struct('A','B','Cc',
%                'Dc','Jc','lam','Q','R') with Cc,Dc,Jc cell arrays of the
%                diffusion / jump channels and lam the rate vector.  Backward:
%                MVST_DRE_CT_POISSON.
%
%   Time-varying variants ('mjls','switched','sd') are detected automatically
%   when the relevant drift is a function handle, or forced by OPTS.tv.
%
%   OPTS fields (all optional):
%     .class  one of the strings above (default: inferred from SYS).
%     .tv     force the time-varying path (default: auto-detected).
%     .X0     initial state covariance.  Defaults to the mode's .X0 for 'ct'
%             and 'mjls'; REQUIRED for 'switched', 'sd', 'poisson'.
%     .M0     initial Newton multiplier (default zeros(n)).  For hard targets
%             the isotropic default M0=0 can leave the damped Newton stranded at
%             a large terminal gap (the trial step drives the backward Riccati
%             into finite escape, then the line search stalls); passing an M0 in
%             the solution basin restores quadratic convergence.  The
%             sampled-data and Poissonian-jump steering examples set it for this
%             reason.
%     .tol    Newton gap tolerance (default 1e-7, via MVST_STEER_NEWTON).
%     .maxit  Newton iteration cap (default 60, via MVST_STEER_NEWTON).
%
%   RES fields (mirroring the retired MVST_STEER_CT):
%     .M      optimal terminal multiplier,
%     .t      time grid of the returned covariance trajectory,
%     .Sigma  n x n x nt covariance trajectory,
%     .K      optimal gain object (shape depends on class: an m x n x nt array
%             for 'ct'/'switched'/'poisson', a 1 x M cell of such arrays for
%             'mjls', a 1 x N cell of the sample gains for 'sd'),
%     .res    final terminal gap norm,
%     .iters  Newton iterations taken,
%     .hist   gap-norm history,
%     .class  the resolved class string.
%
%   See also MVST_STEER_NEWTON, MVST_STEER_ELLIPSE, MVST_STEER_COMPONENTS.
if nargin < 4; opts = struct(); end
if ~isstruct(sys) || ~(isfield(sys,'A') || isfield(sys,'modes'))
    error('mvst:steer:input', ...
        'SYS must be a mode/system struct carrying .A or .modes.');
end
if ~(isscalar(T) && isreal(T) && T > 0)
    error('mvst:steer:horizon', 'T must be a positive scalar.');
end
if size(Xf,1) ~= size(Xf,2)
    error('mvst:steer:target', 'XF (target covariance) must be square.');
end
cls = lower(mvst_getf(opts,'class', infer_class(sys)));

switch cls
    case 'ct'
        if isfield(sys,'modes'); md = sys.modes{1}; else; md = sys; end
        n = md.n; X0 = mvst_getf(opts,'X0', md.X0);
        term = @(M) term_ct(M, md, T, X0, n);
    case 'mjls'
        n = sys.n; X0 = mvst_getf(opts,'X0', sys.modes{1}.X0);
        tv = mvst_getf(opts,'tv', isa(sys.modes{1}.A,'function_handle'));
        if tv; term = @(M) term_mjls_tv(M, sys, T, X0, n);
        else;  term = @(M) term_mjls(M, sys, T, X0, n); end
    case 'switched'
        modes = sys.modes; Td = mvst_getf(sys,'Tdwell', mvst_getf(sys,'Td',[]));
        n = modes{1}.n; X0 = req_X0(opts, cls);
        tv = mvst_getf(opts,'tv', isa(modes{1}.A,'function_handle'));
        if tv; term = @(M) term_sw_tv(M, modes, Td, T, X0, n);
        else;  term = @(M) term_sw(M, modes, Td, T, X0, n); end
    case 'sd'
        n = size(mvst_evalc(sys.A,0),1); m = size(sys.B,2);
        X0 = req_X0(opts, cls);
        tv = mvst_getf(opts,'tv', isa(sys.A,'function_handle'));
        if tv; term = @(M) term_sd_tv(M, sys, T, X0, n, m);
        else;  term = @(M) term_sd(M, sys, T, X0, n, m); end
    case 'poisson'
        n = size(sys.A,1); X0 = req_X0(opts, cls);
        term = @(M) term_poisson(M, sys, T, X0, n);
    otherwise
        error('mvst:steer:class','Unknown steering class ''%s''.', cls);
end

gapfn = @(mv) mvst_vech( term(mvst_unvech(mv,n)) - Xf );
M0 = mvst_getf(opts,'M0', zeros(n));
nopts = struct();
if isfield(opts,'tol');   nopts.tol   = opts.tol;   end
if isfield(opts,'maxit'); nopts.maxit = opts.maxit; end
[M, hist, it] = mvst_steer_newton(gapfn, M0, n, nopts);

[~, tg, Sig, K] = term(M);
res.M = M; res.t = tg; res.Sigma = Sig; res.K = K;
res.res = hist(end); res.iters = it; res.hist = hist; res.class = cls;
end

% ======================================================================= %
%  Class inference / small utilities                                       %
% ======================================================================= %
function cls = infer_class(sys)
if isfield(sys,'Jc') || isfield(sys,'lam')
    cls = 'poisson';
elseif isfield(sys,'Tbar')
    cls = 'sd';
elseif isfield(sys,'Tdwell') || isfield(sys,'Td')
    cls = 'switched';
elseif isfield(sys,'modes') && numel(sys.modes) > 1
    cls = 'mjls';
else
    cls = 'ct';
end
end

function X0 = req_X0(opts, cls)
if ~isfield(opts,'X0') || isempty(opts.X0)
    error('mvst:steer:X0','OPTS.X0 (initial covariance) is required for class ''%s''.', cls);
end
X0 = opts.X0;
end

function S = symm(v, n)
S = reshape(v, n, n);  S = (S + S.')/2;
end

% ======================================================================= %
%  CT: single-mode continuous-time diffusion                               %
% ======================================================================= %
function [SigT, tg, Sig, K] = term_ct(M, md, T, X0, n)
rt = 1e-8; m = md.m;
[tb, ~, ~, Kt] = mvst_dre_ct(md, T, M, rt);
Kmat = reshape(Kt, m*n, numel(tb)).';
Kfun = @(t) reshape(interp1(tb, Kmat, min(max(t,tb(1)),tb(end)), 'pchip'), m, n);
tg = linspace(0, T, 200);
[~, Sv] = ode45(@(t,sv) ct_fwd(sv, Kfun(t), md, n), tg, X0(:), ...
                odeset('RelTol',rt,'AbsTol',rt*1e-2));
Sig = zeros(n,n,numel(tg));
for k = 1:numel(tg); Sig(:,:,k) = symm(Sv(k,:).', n); end
SigT = Sig(:,:,end);
if nargout > 3
    K = zeros(m,n,numel(tg));
    for k = 1:numel(tg); K(:,:,k) = Kfun(tg(k)); end
end
end

function ds = ct_fwd(sv, K, md, n)
S = symm(sv, n);
Acl = md.A + md.B*K; dS = Acl*S + S*Acl.';
for l = 1:md.Nc
    Lc = md.Ac(:,:,l) + md.Bc(:,:,l)*K; dS = dS + Lc*S*Lc.';
end
ds = dS(:);
end

% ======================================================================= %
%  MJLS: Markov-jump diffusion (constant)                                  %
% ======================================================================= %
function [SigT, tg, Sig, K] = term_mjls(M, sys, T, X0, n)
Mn = numel(sys.modes); m = sys.m; sys.ZT = M;
rf = mvst_solve_finite(sys, T);                 % coupled K_i(t)
tg = linspace(0, T, 160);
y0 = zeros(Mn*n*n,1);
for i = 1:Mn; Xi0 = sys.pi0(i)*X0; y0((i-1)*n*n+1:i*n*n) = Xi0(:); end
[~, Y] = ode45(@(t,y) mjls_rhs(t,y,sys,rf,n,m,Mn), tg, y0, ...
               odeset('RelTol',1e-7,'AbsTol',1e-9, ...
                      'Events',@(tt,yy) mvst_odecap(yy,1e8)));
if size(Y,1) < numel(tg)
    Y = [Y; repmat(Y(end,:), numel(tg) - size(Y,1), 1)];
end
Sig = zeros(n,n,numel(tg));
for k = 1:numel(tg)
    Sk = zeros(n);
    for i = 1:Mn; Xi = reshape(Y(k,(i-1)*n*n+1:i*n*n),n,n); Sk = Sk + Xi; end
    Sig(:,:,k) = (Sk+Sk.')/2;
end
SigT = Sig(:,:,end);
if nargout > 3; K = mjls_gains(rf, tg, m, n, Mn); end
end

function dy = mjls_rhs(t, y, sys, rf, n, m, Mn)
Pi = sys.Pi; X = cell(1,Mn);
for i = 1:Mn; Xi = reshape(y((i-1)*n*n+1:i*n*n),n,n); X{i} = (Xi+Xi.')/2; end
dy = zeros(Mn*n*n,1);
tt = min(max(t,rf.t(1)),rf.t(end));
for i = 1:Mn
    md = sys.modes{i};
    Ki = reshape(interp1(rf.t, reshape(rf.K{i},m*n,[]).', tt, 'pchip'), m, n);
    Acl = md.A + md.B*Ki; dX = Acl*X{i} + X{i}*Acl.';
    for l = 1:md.Nc; Lc = md.Ac(:,:,l)+md.Bc(:,:,l)*Ki; dX = dX + Lc*X{i}*Lc.'; end
    for j = 1:Mn; dX = dX + Pi(j,i)*X{j}; end     % Markov coupling
    dy((i-1)*n*n+1:i*n*n) = dX(:);
end
end

% ======================================================================= %
%  MJLS: Markov-jump diffusion (time-varying)                              %
% ======================================================================= %
function [SigT, tg, Sig, K] = term_mjls_tv(M, sys, T, X0, n)
Mn = numel(sys.modes); m = sys.m; sys.ZT = M;
rf = mvst_solve_finite_tv(sys, T);              % coupled K_i(t), time-varying
tg = linspace(0, T, 160);
y0 = zeros(Mn*n*n,1);
for i = 1:Mn; Xi0 = sys.pi0(i)*X0; y0((i-1)*n*n+1:i*n*n) = Xi0(:); end
[~, Y] = ode45(@(t,y) mjls_rhs_tv(t,y,sys,rf,n,m,Mn), tg, y0, ...
               odeset('RelTol',1e-7,'AbsTol',1e-9, ...
                      'Events',@(tt,yy) mvst_odecap(yy,1e8)));
if size(Y,1) < numel(tg)
    Y = [Y; repmat(Y(end,:), numel(tg) - size(Y,1), 1)];
end
Sig = zeros(n,n,numel(tg));
for k = 1:numel(tg)
    Sk = zeros(n);
    for i = 1:Mn; Xi = reshape(Y(k,(i-1)*n*n+1:i*n*n),n,n); Sk = Sk + Xi; end
    Sig(:,:,k) = (Sk+Sk.')/2;
end
SigT = Sig(:,:,end);
if nargout > 3; K = mjls_gains(rf, tg, m, n, Mn); end
end

function dy = mjls_rhs_tv(t, y, sys, rf, n, m, Mn)
Pi = sys.Pi; X = cell(1,Mn);
for i = 1:Mn; Xi = reshape(y((i-1)*n*n+1:i*n*n),n,n); X{i} = (Xi+Xi.')/2; end
dy = zeros(Mn*n*n,1);
tt = min(max(t,rf.t(1)),rf.t(end));
for i = 1:Mn
    md = sys.modes{i};
    Ki = reshape(interp1(rf.t, reshape(rf.K{i},m*n,[]).', tt, 'pchip'), m, n);
    A = mvst_evalc(md.A,t);                     % time-varying flow drift
    Acl = A + mvst_evalc(md.B,t)*Ki; dX = Acl*X{i} + X{i}*Acl.';
    Ac = mvst_evalc(md.Ac,t); Bc = mvst_evalc(md.Bc,t);
    for l = 1:size(Ac,3); Lc = Ac(:,:,l)+Bc(:,:,l)*Ki; dX = dX + Lc*X{i}*Lc.'; end
    for j = 1:Mn; dX = dX + Pi(j,i)*X{j}; end   % Markov coupling
    dy((i-1)*n*n+1:i*n*n) = dX(:);
end
end

function K = mjls_gains(rf, tg, m, n, Mn)
K = cell(1,Mn);
for i = 1:Mn
    Ki = zeros(m,n,numel(tg));
    Kmat = reshape(rf.K{i}, m*n, []).';
    for k = 1:numel(tg)
        tt = min(max(tg(k),rf.t(1)),rf.t(end));
        Ki(:,:,k) = reshape(interp1(rf.t, Kmat, tt, 'pchip'), m, n);
    end
    K{i} = Ki;
end
end

% ======================================================================= %
%  Switched diffusion (constant)                                           %
% ======================================================================= %
function [SigT, tg, Sig, K] = term_sw(M, modes, Td, T, X0, n)
sw = Td:Td:(T-1e-9); edges = [0 sw T];
seq = 1 + mod(0:numel(edges)-2, numel(modes));      % active mode per interval
rt = 1e-8;
oo = odeset('RelTol',rt,'AbsTol',rt*1e-2, ...
            'Events',@(tt,pv) mvst_odecap(pv, 1e8));

% backward: P continuous across switches, P(T)=M
tb = []; Pb = zeros(n,n,0); Pcur = M;
for i = numel(edges)-1:-1:1
    md = modes{seq(i)};
    [ts, Pv] = ode45(@(t,pv) reshape(-mvst_ric_ct(reshape(pv,n,n),md),[],1), ...
                     [edges(i+1) edges(i)], Pcur(:), oo);
    [ts, ix] = sort(ts); Pv = Pv(ix,:);
    for k = 1:numel(ts); Pk = reshape(Pv(k,:),n,n); Pb(:,:,end+1) = (Pk+Pk.')/2; end %#ok<AGROW>
    tb = [tb; ts];                                  %#ok<AGROW>
    Pcur = reshape(Pv(1,:),n,n); Pcur = (Pcur+Pcur.')/2;
end
[tb, ix] = sort(tb); Pb = Pb(:,:,ix);
[tb, iu] = unique(tb); Pb = Pb(:,:,iu);             % drop duplicate switch-boundary times
Pmat = reshape(Pb, n*n, []).';                      % nt x n^2 for interpolation

% forward: Sigma(0)=X0, mode-dependent closed loop using the gain from P(t)
tg = linspace(0, T, 200);
[~, Sv] = ode45(@(t,sv) sw_fwd(t,sv,modes,seq,Td,tb,Pmat,n), tg, X0(:), oo);
if size(Sv,1) < numel(tg)
    Sv = [Sv; repmat(Sv(end,:), numel(tg) - size(Sv,1), 1)];
end
Sig = zeros(n,n,numel(tg));
for k = 1:numel(tg); Sk = reshape(Sv(k,:),n,n); Sig(:,:,k) = (Sk+Sk.')/2; end
SigT = Sig(:,:,end);
if nargout > 3
    m = modes{1}.m; K = zeros(m,n,numel(tg));
    for k = 1:numel(tg)
        mi = seq(min(floor(tg(k)/Td)+1, numel(seq)));
        P = reshape(interp1(tb, Pmat, min(max(tg(k),tb(1)),tb(end)), 'pchip'), n, n);
        [~, Kk] = mvst_ric_ct((P+P.')/2, modes{mi}); K(:,:,k) = Kk;
    end
end
end

function ds = sw_fwd(t, sv, modes, seq, Td, tb, Pmat, n)
S = reshape(sv,n,n); S = (S+S.')/2;
m = seq(min(floor(t/Td)+1, numel(seq))); md = modes{m};
P = reshape(interp1(tb, Pmat, min(max(t,tb(1)),tb(end)), 'pchip'), n, n);
[~, K] = mvst_ric_ct((P+P.')/2, md);
Acl = md.A + md.B*K; dS = Acl*S + S*Acl.';
for l = 1:md.Nc; Lc = md.Ac(:,:,l)+md.Bc(:,:,l)*K; dS = dS + Lc*S*Lc.'; end
ds = dS(:);
end

% ======================================================================= %
%  Switched diffusion (time-varying)                                       %
% ======================================================================= %
function [SigT, tg, Sig, K] = term_sw_tv(M, md, Td, T, X0, n)
sw = Td:Td:(T-1e-9); edges = [0 sw T];
seq = 1 + mod(0:numel(edges)-2, numel(md));         % active mode per interval
rt = 1e-8;
oo = odeset('RelTol',rt,'AbsTol',rt*1e-2, ...
            'Events',@(tt,pv) mvst_odecap(pv, 1e8));

% backward: P continuous across switches, P(T)=M, time-varying I-GRE
tb = []; Pb = zeros(n,n,0); Pcur = M;
for i = numel(edges)-1:-1:1
    am = md{seq(i)};
    [ts, Pv] = ode45(@(t,pv) reshape(-mvst_ric_ct_tv(t,symm(pv,n),am),[],1), ...
                     [edges(i+1) edges(i)], Pcur(:), oo);
    [ts, ix] = sort(ts); Pv = Pv(ix,:);
    for k = 1:numel(ts); Pb(:,:,end+1) = symm(Pv(k,:).',n); end %#ok<AGROW>
    tb = [tb; ts];                                  %#ok<AGROW>
    Pcur = symm(Pv(1,:).',n);
end
[tb, ix] = sort(tb); Pb = Pb(:,:,ix);
[tb, iu] = unique(tb); Pb = Pb(:,:,iu);             % drop duplicate switch-boundary times
Pmat = reshape(Pb, n*n, []).';                      % nt x n^2 for interpolation

% forward: Sigma(0)=X0, mode-dependent closed loop using the gain from P(t)
tg = linspace(0, T, 200);
[~, Sv] = ode45(@(t,sv) sw_fwd_tv(t,sv,md,seq,Td,tb,Pmat,n), tg, X0(:), oo);
if size(Sv,1) < numel(tg)
    Sv = [Sv; repmat(Sv(end,:), numel(tg) - size(Sv,1), 1)];
end
Sig = zeros(n,n,numel(tg));
for k = 1:numel(tg); Sk = reshape(Sv(k,:),n,n); Sig(:,:,k) = (Sk+Sk.')/2; end
SigT = Sig(:,:,end);
if nargout > 3
    m = size(mvst_evalc(md{1}.B,0),2); K = zeros(m,n,numel(tg));
    for k = 1:numel(tg)
        mi = seq(min(floor(tg(k)/Td)+1, numel(seq)));
        P = reshape(interp1(tb, Pmat, min(max(tg(k),tb(1)),tb(end)), 'pchip'), n, n);
        [~, Kk] = mvst_ric_ct_tv(tg(k), (P+P.')/2, md{mi}); K(:,:,k) = Kk;
    end
end
end

function ds = sw_fwd_tv(t, sv, md, seq, Td, tb, Pmat, n)
S = reshape(sv,n,n); S = (S+S.')/2;
mi = seq(min(floor(t/Td)+1, numel(seq))); am = md{mi};
P = reshape(interp1(tb, Pmat, min(max(t,tb(1)),tb(end)), 'pchip'), n, n);
[~, K] = mvst_ric_ct_tv(t, (P+P.')/2, am);
A = mvst_evalc(am.A,t);                             % time-varying flow drift
Acl = A + mvst_evalc(am.B,t)*K; dS = Acl*S + S*Acl.';
Ac = mvst_evalc(am.Ac,t); Bc = mvst_evalc(am.Bc,t);
for l = 1:size(Ac,3); Lc = Ac(:,:,l)+Bc(:,:,l)*K; dS = dS + Lc*S*Lc.'; end
ds = dS(:);
end

% ======================================================================= %
%  Sampled-data (constant)                                                 %
% ======================================================================= %
function [SigT, tg, Sig, K] = term_sd(M, s, T, X0, n, m)
d = n+m; ix = 1:n; iu = n+1:d;
Ab0 = [s.A s.B; zeros(m,d)]; Abl = [s.Ac s.Bc; zeros(m,d)];
Zf  = blkdiag((s.Q+s.Q.')/2, (s.R+s.R.')/2);
tk = s.Tbar:s.Tbar:(T-1e-9); N = numel(tk); edges = [0 tk T];
rt = 1e-8;
oo = odeset('RelTol',rt,'AbsTol',rt*1e-2, ...
            'Events',@(tt,pv) mvst_odecap(pv, 1e8));   % clean stop on escape

% backward: gains K_k at each sample (P(T)=blkdiag(M,0))
P = blkdiag((M+M.')/2, zeros(m)); Kc = cell(1,N);
for i = N+1:-1:1
    [~, Pv] = ode45(@(t,pv) sd_flowback(pv,Ab0,Abl,Zf,d), [edges(i+1) edges(i)], P(:), oo);
    P = reshape(Pv(end,:),d,d); P = (P+P.')/2;
    if i >= 2
        Puu = P(iu,iu); Pux = P(iu,ix); Pxx = P(ix,ix); Pxu = P(ix,iu);
        Kc{i-1} = -(Puu \ Pux);
        sc = Pxx - Pxu*(Puu \ Pux); P = blkdiag((sc+sc.')/2, zeros(m));
    end
end

% forward: Z(0)=blkdiag(X0,0); flow + sample-and-hold resets
Z = blkdiag((X0+X0.')/2, zeros(m)); tg = []; Sig = zeros(n,n,0);
for i = 1:N+1
    ti = linspace(edges(i), edges(i+1), 30);
    [ts, Zv] = ode45(@(t,zv) sd_flowfwd(zv,Ab0,Abl,d), ti, Z(:), oo);
    for k = 1:numel(ts)
        Zk = reshape(Zv(k,:),d,d); Zk = (Zk+Zk.')/2; Sig(:,:,end+1) = Zk(ix,ix); %#ok<AGROW>
    end
    tg = [tg; ts(:)];                                         %#ok<AGROW>
    Z = reshape(Zv(end,:),d,d); Z = (Z+Z.')/2;
    if i <= N
        Jb = [eye(n) zeros(n,m); Kc{i} zeros(m,m)]; Z = Jb*Z*Jb.'; Z = (Z+Z.')/2;
    end
end
SigT = Sig(:,:,end);
if nargout > 3; K = Kc; end
end

function dp = sd_flowback(pv, Ab0, Abl, Zf, d)
P = reshape(pv,d,d); P = (P+P.')/2;
dP = Ab0.'*P + P*Ab0 + Abl.'*P*Abl + Zf;
dp = reshape(-(dP+dP.')/2, [], 1);
end

function dz = sd_flowfwd(zv, Ab0, Abl, d)
Z = reshape(zv,d,d); Z = (Z+Z.')/2;
dZ = Ab0*Z + Z*Ab0.' + Abl*Z*Abl.';
dz = reshape((dZ+dZ.')/2, [], 1);
end

% ======================================================================= %
%  Sampled-data (time-varying)                                             %
% ======================================================================= %
function [SigT, tg, Sig, K] = term_sd_tv(M, s, T, X0, n, m)
d = n+m; ix = 1:n; iu = n+1:d;
Afun = s.A;                                 % time-varying flow drift handle
Abl = [s.Ac s.Bc; zeros(m,d)];
Zf  = blkdiag((s.Q+s.Q.')/2, (s.R+s.R.')/2);
tk = s.Tbar:s.Tbar:(T-1e-9); N = numel(tk); edges = [0 tk T];
rt = 1e-8;
oo = odeset('RelTol',rt,'AbsTol',rt*1e-2, ...
            'Events',@(tt,pv) mvst_odecap(pv, 1e8));   % clean stop on escape

% backward: gains K_k at each sample (P(T)=blkdiag(M,0)), time-varying flow
P = blkdiag((M+M.')/2, zeros(m)); Kc = cell(1,N);
for i = N+1:-1:1
    [~, Pv] = ode45(@(t,pv) sd_flowback_tv(t,pv,Afun,s.B,Abl,Zf,d), [edges(i+1) edges(i)], P(:), oo);
    P = reshape(Pv(end,:),d,d); P = (P+P.')/2;
    if i >= 2
        Puu = P(iu,iu); Pux = P(iu,ix); Pxx = P(ix,ix); Pxu = P(ix,iu);
        Kc{i-1} = -(Puu \ Pux);
        sc = Pxx - Pxu*(Puu \ Pux); P = blkdiag((sc+sc.')/2, zeros(m));
    end
end

% forward: Z(0)=blkdiag(X0,0); time-varying flow + sample-and-hold resets
Z = blkdiag((X0+X0.')/2, zeros(m)); tg = []; Sig = zeros(n,n,0);
for i = 1:N+1
    ti = linspace(edges(i), edges(i+1), 30);
    [ts, Zv] = ode45(@(t,zv) sd_flowfwd_tv(t,zv,Afun,s.B,Abl,d), ti, Z(:), oo);
    for k = 1:numel(ts)
        Zk = reshape(Zv(k,:),d,d); Zk = (Zk+Zk.')/2; Sig(:,:,end+1) = Zk(ix,ix); %#ok<AGROW>
    end
    tg = [tg; ts(:)];                                         %#ok<AGROW>
    Z = reshape(Zv(end,:),d,d); Z = (Z+Z.')/2;
    if i <= N
        Jb = [eye(n) zeros(n,m); Kc{i} zeros(m,m)]; Z = Jb*Z*Jb.'; Z = (Z+Z.')/2;
    end
end
SigT = Sig(:,:,end);
if nargout > 3; K = Kc; end
end

function dp = sd_flowback_tv(t, pv, Afun, B, Abl, Zf, d)
n = size(B,1); m = size(B,2);
Ab0 = [Afun(t) B; zeros(m,n+m)];            % time-varying augmented flow
P = reshape(pv,d,d); P = (P+P.')/2;
dP = Ab0.'*P + P*Ab0 + Abl.'*P*Abl + Zf;
dp = reshape(-(dP+dP.')/2, [], 1);
end

function dz = sd_flowfwd_tv(t, zv, Afun, B, Abl, d)
n = size(B,1); m = size(B,2);
Ab0 = [Afun(t) B; zeros(m,n+m)];            % time-varying augmented flow
Z = reshape(zv,d,d); Z = (Z+Z.')/2;
dZ = Ab0*Z + Z*Ab0.' + Abl*Z*Abl.';
dz = reshape((dZ+dZ.')/2, [], 1);
end

% ======================================================================= %
%  Poissonian-jump diffusion                                               %
% ======================================================================= %
function [SigT, tg, Sig, K] = term_poisson(M, sys, T, X0, n)
A = sys.A; B = sys.B; Cc = sys.Cc; Dc = sys.Dc; Jc = sys.Jc;
lam = sys.lam; Q = sys.Q; R = sys.R; m = size(B,2);
nT = 160; tg = linspace(0, T, nT);
[tk, ~, Kt, okb] = mvst_dre_ct_poisson(A, B, Cc, Dc, Jc, lam, Q, R, T, M, nT);
if ~okb
    % backward Riccati escaped (bad trial multiplier during the Newton line
    % search): return a huge terminal moment so the step is rejected.
    SigT = 1e8*eye(n); Sig = repmat(SigT,1,1,nT);
    if nargout > 3; K = zeros(m,n,nT); end
    return
end
Kfun = @(t) poisson_interp_gain(Kt, tk, t);
oo = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
            'Events',@(tt,yy) mvst_odecap(yy,1e8));
[~, Y] = ode45(@(t,s) poisson_fwd(t,s,Kfun,A,B,Cc,Dc,Jc,lam,n), tg, X0(:), oo);
if size(Y,1) < nT
    Y = [Y; repmat(Y(end,:), nT - size(Y,1), 1)];
end
Sig = zeros(n,n,nT);
for k = 1:nT
    Sk = reshape(Y(k,:),n,n); Sig(:,:,k) = (Sk+Sk.')/2;
end
SigT = Sig(:,:,end);
if nargout > 3
    K = zeros(m,n,nT);
    for k = 1:nT; K(:,:,k) = Kfun(tg(k)); end
end
end

function K = poisson_interp_gain(Kt, tk, t)
t = min(max(t, tk(1)), tk(end));
K = zeros(size(Kt,1), size(Kt,2));
for i = 1:size(Kt,1)
    for j = 1:size(Kt,2)
        K(i,j) = interp1(tk, squeeze(Kt(i,j,:)), t);
    end
end
end

function ds = poisson_fwd(t, s, Kfun, A, B, Cc, Dc, Jc, lam, n)
S = reshape(s,n,n); S = (S+S.')/2;
K = Kfun(t); Acl = A + B*K;
dS = Acl*S + S*Acl.';
for l = 1:numel(Cc); Ccl = Cc{l} + Dc{l}*K; dS = dS + Ccl*S*Ccl.'; end
for r = 1:numel(Jc); dS = dS + lam(r)*(Jc{r}*S*Jc{r}.' - S); end
dS = (dS+dS.')/2; ds = dS(:);
end
