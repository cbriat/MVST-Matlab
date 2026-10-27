function gam = mvst_l2gain_ct_mjls(Ac, Bw, Cc, Dw, Pig, Thoriz, Chan, Dwc, rtol)
%MVST_L2GAIN_CT_MJLS  Exact coupled mean-square L2-gain of an MJLS.
%   GAM = MVST_L2GAIN_CT_MJLS(AC,BW,CC,DW,PIG,THORIZ) returns the
%   mean-square L2-gain of the Markov jump linear system with mode
%   drifts AC{i}, common disturbance input BW, mode outputs CC{i},
%   feedthrough DW, and generator PIG (PIG(i,j) >= 0 for i ~= j, zero row
%   sums).  The gain is the smallest gamma for which the COUPLED
%   bounded-real Riccati equations
%     -Sdot_i = Ai'Si + Si Ai + sum_j PIG(i,j) Sj + Ci'Ci
%               + (Si Bw + Ci'Dw)(gamma^2 I - Dw'Dw)^{-1}(Bw'Si + Dw'Ci)
%   admit a bounded solution; feasibility of a trial gamma is decided by
%   integrating the coupled equations backward from S_i(T)=0 over the
%   horizon THORIZ (default 60) and testing for finite escape.  This is
%   the exact coupled gain: it may be larger or smaller than the
%   worst-mode gain max_i gamma_i.
%
%   GAM = MVST_L2GAIN_CT_MJLS(...,CHAN,DWC) adds the per-mode Brownian
%   diffusion channels of the multiplicative-noise MJLS
%     dx = (A_i x + Bw w)dt + sum_l (Chan_i{l} x + Dwc{l} w) dW_l.
%   CHAN is a 1 x M cell (mode i carrying a 1 x q cell, or an n x n x q array,
%   of state-multiplicative channels) and DWC a 1 x q cell of shared
%   disturbance-multiplicative channels (padded with zeros to match CHAN).
%   Each mode Riccati then gains the fold sum_l Chan_i{l}'S_i Chan_i{l} and the
%   storage-dependent pivot Theta_i = gamma^2 I - Dw'Dw - sum_l Dwc{l}'S_i Dwc{l}
%   (loss of positive-definiteness marks the trial gamma infeasible).  With
%   CHAN and DWC absent or empty the noise-free equations above are integrated
%   unchanged (identical output to the six-argument call).
%
%   See also MVST_L2GAIN_CT, MVST_L2GAIN_CT_POISSON, MVST_GAME_GAMMA_CT_MJLS,
%   MVST_L2GAIN.
if ~iscell(Ac) || ~iscell(Cc)
    error('mvst:l2gain_ct_mjls:modes', 'AC and CC must be 1 x M cells.');
end
if ~isequal(size(Pig), [numel(Ac) numel(Ac)])
    error('mvst:l2gain_ct_mjls:generator', 'PIG must be %d x %d.', numel(Ac), numel(Ac));
end
if nargin < 6 || isempty(Thoriz); Thoriz = 60; end
% optional relative bisection tolerance RTOL (default 1e-3, added 2026-10-04);
% the finite-horizon escape test under-estimates the gain, so long THORIZ and
% small RTOL are needed for reference values.
if nargin < 9 || isempty(rtol); rtol = 1e-3; end
M = numel(Ac); n = size(Ac{1},1); m = size(Dw,2);
if nargin < 7; Chan = {}; end
if nargin < 8; Dwc  = {}; end
% Per-mode diffusion channels as a 1xM cell of 1xq cells; shared Dwc as 1xq.
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
dnorm = norm(Dw);
lo = dnorm + 1e-6;
hi = max(2*lo, 1.0);
for kk = 1:60
    if feasible(hi); break; end
    hi = hi*2;
    if kk == 60
        warning('mvst:l2gain_ct_mjls:nohi','No feasible gamma bracketed.');
        gam = hi; return;
    end
end
while (hi - lo) > rtol*max(1,hi)
    mid = 0.5*(lo+hi);
    if feasible(mid); hi = mid; else; lo = mid; end
end
gam = hi;

    function ok = feasible(g)
        Rg = g^2*eye(m) - Dw.'*Dw;
        [~,pch] = chol((Rg+Rg.')/2);
        if pch ~= 0; ok = false; return; end
        big = 1e6;
        if hasNoise                          % also stop when a pivot Theta_i -> 0
            ev = @(s,p) escapeNoise(p, Rg, big);
        else
            ev = @(s,p) mvst_odecap(p, big);
        end
        opts = odeset('RelTol',1e-6,'AbsTol',1e-8, 'Events', ev);
        try
            [ss,Y] = ode45(@(s,p) rhs(p,Rg), [0 Thoriz], zeros(M*n*n,1), opts);
            ok = (ss(end) >= Thoriz - 1e-6) && all(isfinite(Y(end,:)));
        catch
            ok = false;
        end
    end

    function [val, term, dir] = escapeNoise(p, Rg, big)
        % stop on norm blow-up OR loss of any mode's disturbance pivot Theta_i
        mn = inf;
        for i = 1:M
            Si = reshape(p((i-1)*n*n+1:i*n*n),n,n); Si = (Si+Si.')/2;
            Th = Rg;
            for l = 1:qd; Th = Th - Dwc{l}.'*Si*Dwc{l}; end
            mn = min(mn, min(eig((Th+Th.')/2)));
        end
        val = [big - norm(p); mn - 1e-9]; term = [1;1]; dir = [-1;-1];
    end

    function dp = rhs(p, Rg)
        S = cell(1,M);
        for i = 1:M
            Si = reshape(p((i-1)*n*n+1:i*n*n),n,n); S{i} = (Si+Si.')/2;
        end
        dp = zeros(M*n*n,1);
        for i = 1:M
            cpl = zeros(n);
            for j = 1:M; cpl = cpl + Pig(i,j)*S{j}; end
            if ~hasNoise
                v = S{i}*Bw + Cc{i}.'*Dw;
                G = Ac{i}.'*S{i} + S{i}*Ac{i} + cpl + Cc{i}.'*Cc{i} + v*(Rg\v.');
            else
                mult = zeros(n); v = S{i}*Bw + Cc{i}.'*Dw; Th = Rg;
                for l = 1:qd
                    mult = mult + Chan{i}{l}.'*S{i}*Chan{i}{l};
                    v  = v  + Chan{i}{l}.'*S{i}*Dwc{l};
                    Th = Th - Dwc{l}.'*S{i}*Dwc{l};
                end
                % Theta positivity is enforced by the escapeNoise event, which
                % stops the integration before the pivot degenerates.
                G = Ac{i}.'*S{i} + S{i}*Ac{i} + cpl + Cc{i}.'*Cc{i} ...
                    + mult + v*(Th\v.');
            end
            G = (G+G.')/2;
            dp((i-1)*n*n+1:i*n*n) = G(:);
        end
    end
end

% ------------------------------------------------------------------------
function C = normchan(Chan, M)
%NORMCHAN  Normalise per-mode diffusion channels to a 1xM cell of 1xq cells.
if isempty(Chan); C = repmat({{}}, 1, M); return; end
if ~iscell(Chan)
    error('mvst:l2gain_ct_mjls:chan', 'CHAN must be a 1 x M cell (per mode).');
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
