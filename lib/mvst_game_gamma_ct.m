function [gstar, Pinf, info] = mvst_game_gamma_ct(Afun, Bu, Bw, Q, Ru, Tlong, Cc, Dwc, rtol, threl)
%MVST_GAME_GAMMA_CT  Critical disturbance level gamma* of the H-infinity game.
%   [GSTAR,PINF] = MVST_GAME_GAMMA_CT(AFUN,BU,BW,Q,RU,TLONG) returns the
%   smallest disturbance penalty GAMMA for which the state-feedback game
%   (HJI) Riccati equation admits a bounded stabilizing solution, i.e. the
%   H-infinity-optimal level GAMMA*.  Feasibility of a trial GAMMA is tested
%   by integrating the game differential Riccati (MVST_GAME_DRE_CT) backward
%   over the long horizon [0,TLONG] from P(T)=0: for GAMMA > GAMMA* the dual
%   stays bounded and settles, while for GAMMA < GAMMA* it has a finite
%   escape (the maximizer wins) and the norm blows up.  GSTAR is found by
%   bisection; PINF is the (approximate) infinite-horizon dual at GAMMA*.
%   AFUN may be a constant matrix or a handle @(t) (its t->inf value sets
%   the asymptotic game level for time-varying drift).
%
%   [...] = MVST_GAME_GAMMA_CT(...,CC,DWC) carries the Brownian diffusion
%   channels CC (state-multiplicative) and DWC (disturbance-multiplicative) into
%   the game Riccati (see MVST_GAME_DRE_CT), returning the mean-square critical
%   level.  With CC and DWC absent or empty the deterministic level is returned
%   unchanged.
%
%   [GSTAR,PINF,INFO] = MVST_GAME_GAMMA_CT(...,CC,DWC,RTOL,THREL) switches to the
%   LONG-HORIZON criterion ( omit RTOL or
%   pass it empty to keep the former test, used by the other callers, unchanged).
%   A trial level GAMMA is then declared feasible iff BOTH
%     (a) the backward game Riccati, integrated from P(TLONG)=0 down to t=0
%         (pass TLONG = 1000 for the book examples), reaches t=0 with a finite
%         solution and ||P(t)|| < 1e6 throughout (no finite escape), and
%     (b) the Isaacs pivot keeps the UNIFORM margin
%           lambda_min(Theta(P(t))) >= THREL*GAMMA^2   for all t in [0,TLONG],
%           Theta(P) = GAMMA^2 I - sum_l DWC{l}'P DWC{l},
%         with THREL = 1e-6 by default.  The margin is monitored by a terminal
%         ODE event, so a solution that drifts onto the Theta-singular manifold
%         (P plateauing with Theta near 0) is classified infeasible.
%   The bisection then stops at the RELATIVE bracket (hi-lo) <= RTOL*max(1,hi)
%   (RTOL = 1e-6 for the book examples) and GSTAR = hi is the feasible end.
%   Without disturbance-multiplicative channels Theta = GAMMA^2 I and (b) is
%   void.  INFO returns the final bracket (lo, hi), the margin used, and the
%   smallest least eigenvalue of Theta along the feasible-end trajectory.
%
%   See also MVST_GAME_DRE_CT, MVST_GAME_GAMMA_CT_MJLS, MVST_L2GAIN_CT.
if size(Q,1) ~= size(Q,2)
    error('mvst:game_gamma_ct:notSquare', 'Q must be square.');
end
if size(Bu,1) ~= size(Q,1) || size(Bw,1) ~= size(Q,1)
    error('mvst:game_gamma_ct:dim', 'BU and BW must have %d rows to match Q.', size(Q,1));
end
if nargin < 6 || isempty(Tlong); Tlong = 30; end
if nargin < 7; Cc  = {}; end
if nargin < 8; Dwc = {}; end
if nargin < 9; rtol = []; end
if nargin < 10 || isempty(threl); threl = 1e-6; end
cap = 1e6;                                   % "blow-up" threshold on ||P||
info = struct('lo',[],'hi',[],'threl',[],'thetaMin',[],'mode','legacy');
if ~isempty(rtol)
    [gstar, Pinf, info] = longhorizon(Afun, Bu, Bw, Q, Ru, Tlong, Cc, Dwc, ...
                                      rtol, threl, cap);
    return
end

    function [ok, P0] = feasible(g)
        P0 = [];
        try
            [~, Pt] = mvst_game_dre_ct(Afun, Bu, Bw, Q, Ru, g, [0 Tlong], [], 1e-7, Cc, Dwc);
            nrm = zeros(1,size(Pt,3));
            for k = 1:size(Pt,3); nrm(k) = norm(Pt(:,:,k)); end
            ok = all(isfinite(nrm)) && max(nrm) < cap;
            if ok; P0 = Pt(:,:,1); end
        catch
            ok = false;                      % ode45 failure => finite escape
        end
    end

% bracket: grow hi until feasible, keep lo infeasible
lo = 0.05; hi = 1.0;
for it = 1:60
    [ok,~] = feasible(hi);
    if ok; break; end
    hi = hi*1.6;
    if it == 60; warning('mvst:game_gamma_ct:nohi','no feasible gamma found'); end
end
% ensure lo is infeasible
[oklo,~] = feasible(lo);
if oklo; gstar = lo; [~,Pinf] = feasible(lo*1.05); return; end
% bisection
for it = 1:60
    mid = 0.5*(lo+hi);
    [ok,~] = feasible(mid);
    if ok; hi = mid; else; lo = mid; end
    if (hi-lo) < 1e-3; break; end
end
gstar = hi;
[~, Pinf] = feasible(hi*1.02);
end

% ------------------------------------------------------------------------
function [gstar, Pinf, info] = longhorizon(Afun, Bu, Bw, Q, Ru, Tlong, Cc, Dwc, rtol, threl, cap)
%LONGHORIZON  Bisection with the boundedness + uniform Isaacs-margin test.
if ~iscell(Cc);  Cc  = tocell3(Cc);  end
if ~iscell(Dwc); Dwc = tocell3(Dwc); end
n = size(Q,1); m = size(Bw,2);
qd = max(numel(Cc), numel(Dwc));
if numel(Cc)  < qd; Cc(end+1:qd)  = {zeros(n)};   end
if numel(Dwc) < qd; Dwc(end+1:qd) = {zeros(n,m)}; end
BuRu = Bu*(Ru\Bu.');

    function [ok, P0, thmin] = feasible(g)
        P0 = []; thmin = NaN;
        eps_th = threl*g^2;
        opts = odeset('RelTol',1e-7,'AbsTol',1e-9, ...
                      'Events',@(tt,p) ev(p,g,eps_th));
        try
            [tt, Y] = ode45(@(tt,p) -rhs(tt,p,g), [Tlong 0], zeros(n*n,1), opts);
        catch
            ok = false; return                % integrator failure = escape
        end
        ok = (tt(end) <= 1e-6) && all(isfinite(Y(:)));
        if ~ok; return; end
        th = zeros(numel(tt),1); nr = th;
        for k = 1:numel(tt)
            P = reshape(Y(k,:),n,n); P = (P+P.')/2;
            th(k) = min(eig(pivot(P,g))); nr(k) = norm(P);
        end
        thmin = min(th);
        ok = (thmin >= eps_th) && (max(nr) < cap);
        if ok; P0 = reshape(Y(end,:),n,n); P0 = (P0+P0.')/2; end
    end

    function Th = pivot(P, g)
        Th = g^2*eye(m);
        for l = 1:qd; Th = Th - Dwc{l}.'*P*Dwc{l}; end
        Th = (Th+Th.')/2;
    end

    function dP = rhs(t, p, g)
        A = mvst_evalc(Afun, t);
        P = reshape(p,n,n); P = (P+P.')/2;
        mult = zeros(n); v = P*Bw;
        for l = 1:qd
            mult = mult + Cc{l}.'*P*Cc{l};
            v = v + Cc{l}.'*P*Dwc{l};
        end
        Th = pivot(P, g);
        [~,pc] = chol(Th);
        if pc ~= 0; dP = nan(n*n,1); return; end
        G = A.'*P + P*A + Q + mult - P*BuRu*P + v*(Th\v.');
        G = (G+G.')/2; dP = G(:);
    end

    function [val, term, dir] = ev(p, g, eps_th)
        P = reshape(p,n,n); P = (P+P.')/2;
        val = [cap - norm(p); min(eig(pivot(P,g))) - eps_th];
        term = [1; 1]; dir = [-1; -1];
    end

lo = 0.05; hi = 1.0;
for it = 1:60
    if feasible(hi); break; end
    lo = hi; hi = hi*1.6;
    if it == 60
        warning('mvst:game_gamma_ct:nohi','no feasible gamma found');
    end
end
if feasible(lo)
    error('mvst:game_gamma_ct:lo', 'lower bracket end %.4g is feasible', lo);
end
while (hi - lo) > rtol*max(1,hi)
    mid = 0.5*(lo+hi);
    if feasible(mid); hi = mid; else; lo = mid; end
end
gstar = hi;
[~, Pinf, thmin] = feasible(hi);
info = struct('lo',lo,'hi',hi,'threl',threl,'thetaMin',thmin,'mode','longhorizon');
end

% ------------------------------------------------------------------------
function c = tocell3(X)
%TOCELL3  Slice a 3-D array along its 3rd dimension into a row cell.
if iscell(X)
    c = X;
elseif isempty(X)
    c = {};
else
    p = size(X, 3);
    c = cell(1, p);
    for k = 1:p; c{k} = X(:,:,k); end
end
end
