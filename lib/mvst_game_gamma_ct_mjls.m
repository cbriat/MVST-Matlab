function [gstar, P0, info] = mvst_game_gamma_ct_mjls(Ac, Bu, Bw, Q, Ru, Pig, Tlong, Chan, Dwc, rtol, threl)
%MVST_GAME_GAMMA_CT_MJLS  Exact coupled critical level of the MJLS game.
%   [GSTAR,P0] = MVST_GAME_GAMMA_CT_MJLS(AC,BU,BW,Q,RU,PIG,TLONG) returns
%   the smallest disturbance level gamma for which the COUPLED
%   state-feedback game (HJI) Riccati equations
%     -Pdot_i = Ai'Pi + Pi Ai + sum_j PIG(i,j) Pj + Q
%               - Pi Bu Ru^{-1} Bu' Pi + (1/gamma^2) Pi Bw Bw' Pi
%   admit a bounded solution (backward integration over [0,TLONG] from
%   P_i(T)=0, finite escape <=> the maximizer wins).  P0 is a 1 x M cell
%   of the settled duals P_i(0) at the returned level (evaluated at the
%   last feasible bisection point).  This is the exact critical level of
%   the coupled game; it may be larger or smaller than the worst-mode
%   level max_i gamma*_i.
%
%   [...] = MVST_GAME_GAMMA_CT_MJLS(...,CHAN,DWC) adds the per-mode Brownian
%   diffusion channels CHAN (1 x M cell, mode i a 1 x q cell / n x n x q array of
%   state-multiplicative maps) and shared disturbance-multiplicative channels DWC
%   (1 x q cell): each coupled mode Riccati gains the fold
%   sum_l Chan_i{l}'P_i Chan_i{l} and, when DWC is present, the storage-dependent
%   maximizer Hessian Theta_i = gamma^2 I - sum_l Dwc{l}'P_i Dwc{l}.  With CHAN
%   and DWC absent or empty the noise-free coupled equations are integrated
%   unchanged.
%
%   [GSTAR,P0,INFO] = MVST_GAME_GAMMA_CT_MJLS(...,CHAN,DWC,RTOL,THREL) switches
%   to the LONG-HORIZON criterion ( omit
%   RTOL or pass it empty to keep the former test unchanged).  A trial level
%   GAMMA is then feasible iff BOTH
%     (a) the coupled backward integration over [0,TLONG] (pass TLONG = 1000
%         for the book examples) reaches the end with a finite solution and
%         max_i ||P_i|| < 1e6 throughout, and
%     (b) every mode keeps the UNIFORM Isaacs margin
%           lambda_min(Theta_i(P_i(t))) >= THREL*GAMMA^2   for all t,
%           Theta_i = GAMMA^2 I - sum_l DWC{l}'P_i DWC{l},
%         THREL = 1e-6 by default (terminal ODE event), so a drift onto the
%         Theta-singular manifold is classified infeasible.
%   The bisection stops at the relative bracket (hi-lo) <= RTOL*max(1,hi)
%   (RTOL = 1e-6 for the book examples); GSTAR = hi is the feasible end.
%   INFO returns the final bracket, the margin, and the smallest least
%   eigenvalue of the Theta_i along the feasible-end trajectory.
%
%   See also MVST_GAME_GAMMA_CT, MVST_L2GAIN_CT_MJLS, MVST_GAME_DRE_CT.
if ~iscell(Ac)
    error('mvst:game_gamma_ct_mjls:modes', 'AC must be a 1 x M cell of drifts.');
end
if size(Q,1) ~= size(Q,2)
    error('mvst:game_gamma_ct_mjls:notSquare', 'Q must be square.');
end
if size(Bu,1) ~= size(Q,1) || size(Bw,1) ~= size(Q,1)
    error('mvst:game_gamma_ct_mjls:dim', 'BU and BW must have %d rows to match Q.', size(Q,1));
end
if nargin < 7 || isempty(Tlong); Tlong = 30; end
M = numel(Ac); n = size(Q,1); m = size(Bw,2);
if nargin < 8; Chan = {}; end
if nargin < 9; Dwc  = {}; end
Chan = normchan(Chan, M);
if ~iscell(Dwc); Dwc = tocell3(Dwc); end
qd = numel(Dwc);
for ii = 1:M; qd = max(qd, numel(Chan{ii})); end
if qd > 0
    if numel(Dwc) < qd; Dwc(end+1:qd) = {zeros(n,m)}; end
    for ii = 1:M
        if numel(Chan{ii}) < qd; Chan{ii}(end+1:qd) = {zeros(n)}; end
    end
end
% Structurally all-zero channels stay on the deterministic exact path.
hasNoise = anynz(Dwc);
for ii = 1:M; hasNoise = hasNoise || anynz(Chan{ii}); end
BuRu = Bu*(Ru\Bu.'); BwBw = Bw*Bw.';
cap = 1e6;
if nargin < 10; rtol = []; end
if nargin < 11 || isempty(threl); threl = 1e-6; end
info = struct('lo',[],'hi',[],'threl',[],'thetaMin',[],'mode','legacy');
longMode = ~isempty(rtol);

    function [ok, Pc, thmin] = feasible(g)
        Pc = []; thmin = NaN;
        if longMode
            [ok, Pc, thmin] = feasibleLong(g); return
        end
        opts = odeset('RelTol',1e-7,'AbsTol',1e-9, ...
                      'Events',@(s,p) mvst_odecap(p,cap));
        try
            [ss, Y] = ode45(@(s,p) rhs(p,g), [0 Tlong], zeros(M*n*n,1), opts);
            % Boundedness requires reaching Tlong AND a finite trajectory: with a
            % disturbance-multiplicative channel the maximizer-Hessian Theta_i can
            % lose positivity mid-integration, and ode45 then limps to Tlong with
            % NaN entries (chol-guard returns a NaN derivative).  The ss(end) test
            % alone accepts those NaN runs as "bounded", collapsing the bisection;
            % the isfinite guard (mirroring MVST_GAME_GAMMA_CT) rejects them so the
            % coupled critical level is found correctly.  Noise-free runs are never
            % NaN, so this is bit-identical on the diffusion-free path.
            ok = (ss(end) >= Tlong - 1e-6) && all(isfinite(Y(:)));
            if ok
                Pc = cell(1,M);
                for i = 1:M
                    Pi_ = reshape(Y(end,(i-1)*n*n+1:i*n*n).',n,n);
                    Pc{i} = (Pi_+Pi_.')/2;
                end
            end
        catch
            ok = false;
        end
    end

    function [ok, Pc, thmin] = feasibleLong(g)
        % Long-horizon test: boundedness AND uniform Isaacs margin (see header).
        Pc = []; thmin = NaN; eps_th = threl*g^2;
        opts = odeset('RelTol',1e-7,'AbsTol',1e-9, ...
                      'Events',@(s,p) evLong(p,g,eps_th));
        try
            [ss, Y] = ode45(@(s,p) rhs(p,g), [0 Tlong], zeros(M*n*n,1), opts);
        catch
            ok = false; return
        end
        ok = (ss(end) >= Tlong - 1e-6) && all(isfinite(Y(:)));
        if ~ok; return; end
        th = inf(numel(ss),1); nr = zeros(numel(ss),1);
        for k = 1:numel(ss)
            for i = 1:M
                Pi_ = reshape(Y(k,(i-1)*n*n+1:i*n*n).',n,n); Pi_ = (Pi_+Pi_.')/2;
                th(k) = min(th(k), min(eig(pivot(Pi_,g))));
                nr(k) = max(nr(k), norm(Pi_));
            end
        end
        thmin = min(th);
        ok = (thmin >= eps_th) && (max(nr) < cap);
        if ok
            Pc = cell(1,M);
            for i = 1:M
                Pi_ = reshape(Y(end,(i-1)*n*n+1:i*n*n).',n,n);
                Pc{i} = (Pi_+Pi_.')/2;
            end
        end
    end

    function Th = pivot(Pi_, g)
        Th = g^2*eye(m);
        for l = 1:qd; Th = Th - Dwc{l}.'*Pi_*Dwc{l}; end
        Th = (Th+Th.')/2;
    end

    function [val, term, dir] = evLong(p, g, eps_th)
        thm = inf;
        for i = 1:M
            Pi_ = reshape(p((i-1)*n*n+1:i*n*n),n,n); Pi_ = (Pi_+Pi_.')/2;
            thm = min(thm, min(eig(pivot(Pi_,g))));
        end
        val = [cap - norm(p); thm - eps_th];
        term = [1; 1]; dir = [-1; -1];
    end

    function dp = rhs(p, g)
        P = cell(1,M);
        for i = 1:M
            Pi_ = reshape(p((i-1)*n*n+1:i*n*n),n,n); P{i} = (Pi_+Pi_.')/2;
        end
        dp = zeros(M*n*n,1);
        for i = 1:M
            cpl = zeros(n);
            for j = 1:M; cpl = cpl + Pig(i,j)*P{j}; end
            if ~hasNoise
                G = Ac{i}.'*P{i} + P{i}*Ac{i} + cpl + Q ...
                    - P{i}*BuRu*P{i} + (1/g^2)*P{i}*BwBw*P{i};
            else
                mult = zeros(n); v = P{i}*Bw; Th = g^2*eye(m);
                for l = 1:qd
                    mult = mult + Chan{i}{l}.'*P{i}*Chan{i}{l};
                    v  = v  + Chan{i}{l}.'*P{i}*Dwc{l};
                    Th = Th - Dwc{l}.'*P{i}*Dwc{l};
                end
                [~,pc] = chol((Th+Th.')/2);
                if pc ~= 0; dp = nan(M*n*n,1); return; end
                G = Ac{i}.'*P{i} + P{i}*Ac{i} + cpl + Q ...
                    - P{i}*BuRu*P{i} + mult + v*(Th\v.');
            end
            G = (G+G.')/2;
            dp((i-1)*n*n+1:i*n*n) = G(:);
        end
    end

lo = 1e-3; hi = 1; P0 = [];
for kk = 1:60
    [ok, Pc] = feasible(hi);
    if ok; P0 = Pc; break; end
    if longMode; lo = hi; end               % an infeasible hi is a valid lo
    hi = hi*2;
    if kk == 60
        warning('mvst:game_gamma_ct_mjls:nohi','No feasible gamma bracketed.');
        gstar = hi; return;
    end
end
if longMode
    if feasible(lo)
        error('mvst:game_gamma_ct_mjls:lo', 'lower bracket end %.4g is feasible', lo);
    end
    tolb = rtol;
else
    tolb = 1e-3;
end
while (hi - lo) > tolb*max(1,hi)
    mid = 0.5*(lo+hi);
    [ok, Pc] = feasible(mid);
    if ok; hi = mid; P0 = Pc; else; lo = mid; end
end
gstar = hi;
if longMode
    [~, P0, thmin] = feasible(hi);
    info = struct('lo',lo,'hi',hi,'threl',threl,'thetaMin',thmin,'mode','longhorizon');
end
end

% ------------------------------------------------------------------------
function C = normchan(Chan, M)
%NORMCHAN  Normalise per-mode diffusion channels to a 1xM cell of 1xq cells.
if isempty(Chan); C = repmat({{}}, 1, M); return; end
if ~iscell(Chan)
    error('mvst:game_gamma_ct_mjls:chan', 'CHAN must be a 1 x M cell (per mode).');
end
C = cell(1, M);
for i = 1:M
    if i <= numel(Chan); ci = Chan{i}; else; ci = {}; end
    if ~iscell(ci); ci = tocell3(ci); end
    C{i} = ci;
end
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

% ------------------------------------------------------------------------
function tf = anynz(C)
%ANYNZ  True if a channel cell carries any structurally nonzero matrix.
if ~iscell(C); tf = ~isempty(C) && nnz(C) > 0; return; end
tf = any(cellfun(@(x) ~isempty(x) && nnz(x) > 0, C));
end
