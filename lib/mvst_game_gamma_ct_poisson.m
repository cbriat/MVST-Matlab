function [gstar, Pinf, info] = mvst_game_gamma_ct_poisson(A, Bu, Bw, Cc, Jc, lam, Q, Ru, Tlong, BuJc, BwJc, Rdc, rtol, threl)
%MVST_GAME_GAMMA_CT_POISSON  Critical level gamma* for GBM+Poisson games.
%   Bisection on the disturbance level, feasibility decided by boundedness
%   of the backward game Riccati of MVST_GAME_DRE_CT_POISSON over [0,TLONG]
%   (finite escape <=> the maximizer wins <=> gamma < gamma*).
%
%   Optional trailing cells BuJc, BwJc, Rdc (one entry per jump map) route the
%   jump control/disturbance channels and jump-control weight into the static
%   Isaacs jump saddle of MVST_GAME_DRE_CT_POISSON; omit them for a state-only
%   reset x -> J_r x.
%
%   [GSTAR,PINF,INFO] = MVST_GAME_GAMMA_CT_POISSON(...,RDC,RTOL,THREL) switches
%   to the LONG-HORIZON criterion ( omit
%   RTOL or pass it empty to keep the former test unchanged).  A trial level
%   GAMMA is then feasible iff BOTH
%     (a) the backward game Riccati over [0,TLONG] (pass TLONG = 1000 for the
%         book examples) reaches t=0 with a finite solution and ||P|| < 1e6
%         throughout, and
%     (b) the Isaacs pivot of every jump saddle keeps the UNIFORM margin
%           lambda_min(Theta_r(P(t))) >= THREL*GAMMA^2   for all t,
%           Theta_r = GAMMA^2 I - BwJc{r}'P BwJc{r},
%         THREL = 1e-6 by default (terminal ODE event).  The flow carries no
%         disturbance-multiplicative channel, so its pivot is GAMMA^2 I.
%   The bisection stops at the relative bracket (hi-lo) <= RTOL*max(1,hi)
%   (RTOL = 1e-6 for the book examples); GSTAR = hi is the feasible end.
%
%   See also MVST_GAME_DRE_CT_POISSON, MVST_GAME_GAMMA_CT, MVST_L2GAIN_CT_POISSON.
if size(A,1) ~= size(A,2)
    error('mvst:game_gamma_ct_poisson:notSquare', 'A must be square.');
end
if size(Q,1) ~= size(A,1)
    error('mvst:game_gamma_ct_poisson:dim', 'Q must be %d x %d to match A.', size(A,1), size(A,1));
end
if size(Bu,1) ~= size(A,1) || size(Bw,1) ~= size(A,1)
    error('mvst:game_gamma_ct_poisson:dim', 'BU and BW must have %d rows to match A.', size(A,1));
end
if nargin < 9 || isempty(Tlong); Tlong = 30; end
if nargin < 10; BuJc = {}; end
if nargin < 11; BwJc = {}; end
if nargin < 12; Rdc  = {}; end
cap = 1e6;
if nargin < 13; rtol = []; end
if nargin < 14 || isempty(threl); threl = 1e-6; end
info = struct('lo',[],'hi',[],'threl',[],'thetaMin',[],'mode','legacy');
if ~isempty(rtol)
    [gstar, Pinf, info] = longhorizon(A,Bu,Bw,Cc,Jc,lam,Q,Ru,Tlong,BuJc,BwJc,Rdc,rtol,threl,cap);
    return
end

% Boundedness trials at infeasible levels drive the jump-saddle solve toward
% singularity by design; silence the expected near-boundary matrix warnings.
wstate = warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix');
cleanup = onCleanup(@() warning(wstate));

    function [ok, P0] = feasible(g)
        P0 = [];
        try
            [~, Pt] = mvst_game_dre_ct_poisson(A,Bu,Bw,Cc,Jc,lam,Q,Ru,g,[0 Tlong],[],1e-7,BuJc,BwJc,Rdc);
            nrm = zeros(1,size(Pt,3));
            for k = 1:size(Pt,3); nrm(k) = norm(Pt(:,:,k)); end
            ok = all(isfinite(nrm)) && max(nrm) < cap;
            if ok; P0 = Pt(:,:,1); end
        catch
            ok = false;
        end
    end

lo = 1e-3; hi = 1;
for kk = 1:60
    if feasible(hi); break; end
    hi = hi*2;
    if kk == 60
        warning('mvst:game_gamma_ct_poisson:nohi','No feasible gamma bracketed.');
        gstar = hi; Pinf = []; return;
    end
end
Pinf = [];
while (hi - lo) > 1e-3*max(1,hi)
    mid = 0.5*(lo+hi);
    [ok, P0] = feasible(mid);
    if ok; hi = mid; Pinf = P0; else; lo = mid; end
end
gstar = hi;
if isempty(Pinf); [~, Pinf] = feasible(hi); end
end

% ------------------------------------------------------------------------
function [gstar, Pinf, info] = longhorizon(A,Bu,Bw,Cc,Jc,lam,Q,Ru,Tlong,BuJc,BwJc,Rdc,rtol,threl,cap)
%LONGHORIZON  Bisection with the boundedness + uniform jump-Isaacs-margin test.
n = size(A,1); p = numel(Jc);
hasIO = ~isempty(BuJc) || ~isempty(BwJc);
if hasIO
    if isempty(BuJc); BuJc = repmat({zeros(n,0)},1,p); end
    if isempty(BwJc); BwJc = repmat({zeros(n,0)},1,p); end
    if isempty(Rdc);  Rdc  = repmat({[]},1,p);         end
end
BuRu = Bu*(Ru\Bu.'); BwBw = Bw*Bw.';
wstate = warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:illConditionedMatrix');
cleanup = onCleanup(@() warning(wstate));

    function thm = margin(P, g)
        thm = inf;
        if ~hasIO; return; end
        for r = 1:p
            BwJ = BwJc{r};
            if isempty(BwJ); continue; end
            Th = g^2*eye(size(BwJ,2)) - BwJ.'*P*BwJ;
            thm = min(thm, min(eig((Th+Th.')/2)));
        end
    end

    function G = rhs(pv, g)
        P = reshape(pv,n,n); P = (P+P.')/2;
        mult = zeros(n);
        for l = 1:numel(Cc); mult = mult + Cc{l}.'*P*Cc{l}; end
        pois = zeros(n);
        for r = 1:p
            Phat = P;
            if hasIO && (~isempty(BuJc{r}) || ~isempty(BwJc{r}))
                BuJ = BuJc{r}; BwJ = BwJc{r}; mwJ = size(BwJ,2);
                Rd = Rdc{r}; if isempty(Rd); Rd = eye(size(BuJ,2)); end
                if mwJ > 0
                    Th = g^2*eye(mwJ) - BwJ.'*P*BwJ;
                    % concavity guard; the margin itself is monitored by ev()
                    if min(eig((Th+Th.')/2)) <= 0; G = nan(n*n,1); return; end
                end
                Br = [BuJ, BwJ];
                Mr = Br.'*P*Br + blkdiag(Rd, -g^2*eye(mwJ));
                Phat = P - P*Br*(Mr\(Br.'*P)); Phat = (Phat+Phat.')/2;
            end
            pois = pois + lam(r)*(Jc{r}.'*Phat*Jc{r} - P);
        end
        G = A.'*P + P*A + mult + pois + Q - P*BuRu*P + (1/g^2)*P*BwBw*P;
        G = (G+G.')/2; G = G(:);
    end

    function [ok, P0, thmin] = feasible(g)
        P0 = []; thmin = NaN; eps_th = threl*g^2;
        opts = odeset('RelTol',1e-7,'AbsTol',1e-9, ...
                      'Events',@(s,pv) ev(pv,g,eps_th));
        try
            % time-invariant data: integrate in the backward time s = TLONG - t
            [ss, Y] = ode45(@(s,pv) rhs(pv,g), [0 Tlong], zeros(n*n,1), opts);
        catch
            ok = false; return
        end
        ok = (ss(end) >= Tlong - 1e-6) && all(isfinite(Y(:)));
        if ~ok; return; end
        th = zeros(numel(ss),1); nr = th;
        for k = 1:numel(ss)
            P = reshape(Y(k,:),n,n); P = (P+P.')/2;
            th(k) = margin(P,g); nr(k) = norm(P);
        end
        thmin = min(th);
        ok = (thmin >= eps_th) && (max(nr) < cap);
        if ok; P0 = reshape(Y(end,:),n,n); P0 = (P0+P0.')/2; end
    end

    function [val, term, dir] = ev(pv, g, eps_th)
        P = reshape(pv,n,n); P = (P+P.')/2;
        thm = margin(P,g); if ~isfinite(thm); thm = 1; end
        val = [cap - norm(pv); thm - eps_th];
        term = [1; 1]; dir = [-1; -1];
    end

lo = 1e-3; hi = 1;
for kk = 1:60
    if feasible(hi); break; end
    lo = hi; hi = hi*2;
    if kk == 60
        warning('mvst:game_gamma_ct_poisson:nohi','No feasible gamma bracketed.');
    end
end
if feasible(lo)
    error('mvst:game_gamma_ct_poisson:lo', 'lower bracket end %.4g is feasible', lo);
end
while (hi - lo) > rtol*max(1,hi)
    mid = 0.5*(lo+hi);
    if feasible(mid); hi = mid; else; lo = mid; end
end
gstar = hi;
[~, Pinf, thmin] = feasible(hi);
info = struct('lo',lo,'hi',hi,'threl',threl,'thetaMin',thmin,'mode','longhorizon');
end
