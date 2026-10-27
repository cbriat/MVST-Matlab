function [gam, info] = mvst_l2gain_sd_hybrid(A, Bw, C, Tbar, Nper, Cc, Dwc, rtol, Ntol)
%MVST_L2GAIN_SD_HYBRID  Exact hybrid L2-gain, flow disturbance / sampled output.
%   GAM = MVST_L2GAIN_SD_HYBRID(A,BW,C,TBAR,NPER) returns the hybrid
%   L2-gain of the periodic sampled-output system: continuous flow
%   xdot = A x + Bw w, sampled performance output z_k = C x(k*TBAR).  By
%   the periodic dissipativity theorem the gain is the smallest gamma for
%   which the periodic bounded-real I-D-GRE has a bounded solution, which
%   is decided by backward integration over NPER periods:
%   flow Riccati -Sdot = A'S + S A + S Bw (gamma^2 I)^{-1} Bw' S between
%   samples (disturbance credit only: the output is sampled, not
%   integrated) and the jump update S <- S + C'C at each sample; finite
%   escape <=> gamma below the hybrid gain.
%
%   GAM = MVST_L2GAIN_SD_HYBRID(...,CC,DWC) adds the flow Brownian diffusion
%   channels of the multiplicative-noise flow
%     dx = (A x + Bw w)dt + sum_l (Cc{l} x + Dwc{l} w) dW_l
%   (CC, DWC 1 x q cells; DWC padded with zeros to match CC), giving the
%   mean-square hybrid gain: the flow Riccati gains the fold sum_l Cc{l}'S Cc{l}
%   and the storage-dependent pivot Theta = gamma^2 I - sum_l Dwc{l}'S Dwc{l}.
%   With CC and DWC absent or empty the deterministic-flow equation above is
%   integrated unchanged.  (This reset-free sampled-output backend carries no
%   state reset, so a jump-disturbance channel does not apply; the Poissonian
%   jump disturbance lives in MVST_L2GAIN_CT_POISSON's Bdf channel.)
%
%   NPER AND CONVERGENCE (changed 2026-10-27).  The escape test over a FINITE
%   number of periods under-reports the gain (book example data A=[0 1;-1 -1],
%   Bw=[0;1], C=[1 0], Tbar=0.3: NPER=60 gives 2.053, NPER=200 gives 2.102, the
%   converged periodic gain is 2.108).  When NPER is omitted or empty the
%   routine is therefore ADAPTIVE: it computes the gain for NPER = 60, 120, 240,
%   ... (doubling, at most 7680) until the increment between two successive
%   values is below NTOL (relative, default 5e-4), and inside each trial the
%   recursion stops early once the storage S reaches its periodic fixed point
%   (a bounded fixed point certifies the trial gamma).  Without multiplicative
%   noise the flow Riccati is propagated in closed form by a matrix exponential
%   (conjugate-point loss of positive semidefiniteness = escape), which makes the
%   long horizons cheap; with noise the ODE integration is used.  An explicit
%   scalar NPER reproduces the former fixed-horizon behavior unchanged.
%   [GAM,INFO] = ... returns INFO.Nper (periods used), INFO.increment (last
%   relative increment), INFO.converged (increment below NTOL), INFO.adaptive.
%
%   See also MVST_L2GAIN_CT, MVST_L2GAIN_CT_POISSON, MVST_IGRE_PERIODIC_SD,
%   MVST_DISS_PT_CT.
if size(A,1) ~= size(A,2)
    error('mvst:l2gain_sd_hybrid:notSquare', 'A must be square.');
end
if size(Bw,1) ~= size(A,1)
    error('mvst:l2gain_sd_hybrid:dim', 'BW must have %d rows to match A.', size(A,1));
end
if size(C,2) ~= size(A,1)
    error('mvst:l2gain_sd_hybrid:dim', 'C must have %d columns to match A.', size(A,1));
end
if ~(isscalar(Tbar) && isreal(Tbar) && Tbar > 0)
    error('mvst:l2gain_sd_hybrid:period', 'TBAR must be a positive scalar.');
end
adaptive = (nargin < 5 || isempty(Nper));
if adaptive; Nper = 60; end
if nargin < 9 || isempty(Ntol); Ntol = 5e-4; end
if nargin < 6; Cc  = {}; end
if nargin < 7; Dwc = {}; end
% optional relative bisection tolerance RTOL (default 1e-3, added 2026-10-04);
% the finite-horizon escape test under-estimates the gain, so many periods
% NPER and a small RTOL are needed for reference values.
if nargin < 8 || isempty(rtol)
    rtol = 1e-3; if adaptive; rtol = 2e-4; end
end
n = size(A,1); m = size(Bw,2); CtC = C.'*C;
if ~iscell(Cc);  Cc  = tocell3(Cc);  end
if ~iscell(Dwc); Dwc = tocell3(Dwc); end
qd = max(numel(Cc), numel(Dwc));
if qd > 0
    if numel(Cc)  < qd; Cc(end+1:qd)  = {zeros(n)};   end
    if numel(Dwc) < qd; Dwc(end+1:qd) = {zeros(n,m)}; end
end
hasCc  = any(~cellfun(@(x) isempty(x) || nnz(x)==0, Cc));
hasDwc = any(~cellfun(@(x) isempty(x) || nnz(x)==0, Dwc));
hasNoise = hasCc || hasDwc;
NperMax = 7680;
if adaptive
    gprev = bisect_gain();
    inc = inf;
    while true
        if Nper >= NperMax; break; end
        Nper = 2*Nper;
        gam = bisect_gain();
        inc = abs(gam - gprev)/max(1, gam);
        gprev = gam;
        if inc <= Ntol; break; end
    end
    gam = gprev;
    converged = (inc <= Ntol);
    if ~converged
        warning('mvst:l2gain_sd_hybrid:nperconv', ...
            ['Adaptive NPER stopped at %d periods with relative increment %.2e > %.1e; ' ...
             'the returned gain may still under-report.'], Nper, inc, Ntol);
    end
    info = struct('Nper', Nper, 'increment', inc, 'converged', converged, 'adaptive', true);
    return;
end
gam = bisect_gain();
info = struct('Nper', Nper, 'increment', NaN, 'converged', NaN, 'adaptive', false);
return;

    function g0 = bisect_gain()
        lo = 1e-6; hi = 1;
        for kk = 1:60
            if feasible(hi); break; end
            hi = hi*2;
            if kk == 60
                warning('mvst:l2gain_sd_hybrid:nohi','No feasible gamma bracketed.');
                g0 = hi; return;
            end
        end
        while (hi - lo) > rtol*max(1,hi)
            mid = 0.5*(lo+hi);
            if feasible(mid); hi = mid; else; lo = mid; end
        end
        g0 = hi;
    end

    function ok = feasible_fast(g)
        % closed-form flow Riccati (no noise): S(tau) = Y X^{-1} with
        % [X;Y]' = [-A -Q; 0 A'][X;Y], Q = Bw Bw'/g^2, over 4 sub-steps.
        Q = (Bw*Bw.')/g^2; ns = 4;
        Hm = [-A, -Q; zeros(n), A.'];
        F = expm(Hm*(Tbar/ns)); F11 = F(1:n,1:n); F12 = F(1:n,n+1:end); F22 = F(n+1:end,n+1:end);
        S = zeros(n); ok = true;
        for k = 1:Nper
            Sold = S; S = S + CtC;
            for q = 1:ns
                X = F11 + F12*S;
                if rcond(X) < 1e-13; ok = false; return; end
                S = (F22*S)/X; S = (S+S.')/2;
                if ~all(isfinite(S(:))) || norm(S) > 1e5 || min(eig(S)) < -1e-8*(1+norm(S))
                    ok = false; return;
                end
            end
            if adaptive && norm(S - Sold) <= 1e-10*(1+norm(S)); return; end
        end
    end

    function ok = feasible(g)
        if adaptive && ~hasNoise; ok = feasible_fast(g); return; end
        Ri = (1/g^2)*eye(m);
        S = zeros(n); ok = true;
        if hasDwc                            % also stop when the pivot Theta -> 0
            ev = @(s,q) escapeNoiseSD(q, g, 1e6);
        else
            ev = @(s,q) mvst_odecap(q, 1e6);
        end
        opts = odeset('RelTol',1e-7,'AbsTol',1e-9, 'Events', ev);
        for k = 1:Nper
            Sold = S;
            S = S + CtC;                                  % sample (backward)
            [ss, Y] = ode45(@(s,q) rhs(q,Ri,g), [0 Tbar], S(:), opts);
            if ss(end) < Tbar - 1e-9 || ~all(isfinite(Y(end,:))); ok = false; return; end
            S = reshape(Y(end,:).',n,n); S = (S+S.')/2;
            if norm(S) > 1e5; ok = false; return; end
            if adaptive && norm(S - Sold) <= 1e-10*(1+norm(S)); return; end
        end
    end

    function [val, term, dir] = escapeNoiseSD(q, g, big)
        S = reshape(q,n,n); S = (S+S.')/2;
        Th = g^2*eye(m);
        for l = 1:qd; Th = Th - Dwc{l}.'*S*Dwc{l}; end
        val = [big - norm(q); min(eig((Th+Th.')/2)) - 1e-9]; term = [1;1]; dir = [-1;-1];
    end

    function dq = rhs(q, Ri, g)
        S = reshape(q,n,n); S = (S+S.')/2;
        if ~hasNoise
            G = A.'*S + S*A + S*Bw*Ri*Bw.'*S;
        else
            mult = zeros(n);
            for l = 1:qd; mult = mult + Cc{l}.'*S*Cc{l}; end
            if ~hasDwc
                G = A.'*S + S*A + mult + S*Bw*Ri*Bw.'*S;
            else
                % Theta positivity enforced by the escapeNoiseSD event.
                v = S*Bw; Th = g^2*eye(m);
                for l = 1:qd
                    v  = v  + Cc{l}.'*S*Dwc{l};
                    Th = Th - Dwc{l}.'*S*Dwc{l};
                end
                G = A.'*S + S*A + mult + v*(Th\v.');
            end
        end
        G = (G+G.')/2; dq = G(:);
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
