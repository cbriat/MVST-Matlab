function ok = mvst_diss_dwell_cert(Ac, Bw, Cc, Dw, Td, g, Chan, Dwc)
%MVST_DISS_DWELL_CERT  MDT dissipativity certificate for a switched system.
%   OK = MVST_DISS_DWELL_CERT(AC,BW,CC,DW,TD,G) attempts to certify that
%   the switched system with mode drifts AC{i}, outputs CC{i}, common
%   (BW,DW), and minimum dwell-time TD has L2-gain at most G, by
%   constructing the timer-dependent mode storages of the MDT
%   dissipativity theorem (identity switches, no jump supply):
%     (i)   S_i(tau) solves the mode-i bounded-real Riccati (equality)
%           on [0,TD], parameterized by its terminal value S_i(TD);
%     (ii)  switch condition S_j(0) <= S_i(TD) for all j ~= i;
%     (iii) end-of-dwell condition BRI_i(S_i(TD)) <= 0 (frozen tail);
%     (iv)  S_i(tau) >= 0 on [0,TD].
%   The terminal values are searched by Nelder-Mead (fminsearch, no
%   toolbox) starting from the per-mode stabilizing solutions.  OK = true
%   only when a certificate satisfying (i)-(iv) to tolerance is FOUND, so
%   a true return is a genuine upper-bound certificate; a false return
%   means only that the search failed.
%
%   OK = MVST_DISS_DWELL_CERT(...,CHAN,DWC) additionally carries the per-mode
%   Brownian diffusion channels CHAN (1 x M cell, mode i a 1 x q cell / n x n x q
%   array of state-multiplicative maps) and the shared disturbance-multiplicative
%   channels DWC (1 x q cell), folding sum_l chi_i{l}'S chi_i{l} into every mode
%   bounded-real inequality and loading the storage-dependent pivot
%   Theta_i = Rg - sum_l Dwc{l}'S Dwc{l}.  With CHAN and DWC absent or empty the
%   noise-free certificate above is constructed unchanged.
%
%   See also MVST_L2GAIN_CT_DWELL, MVST_GAME_GAMMA_CT_DWELL, MVST_L2GAIN_CT.
if ~iscell(Ac) || ~iscell(Cc)
    error('mvst:diss_dwell_cert:modes', 'AC and CC must be 1 x M cells.');
end
if ~(isscalar(Td) && isreal(Td) && Td > 0)
    error('mvst:diss_dwell_cert:dwell', 'TD (dwell-time) must be a positive scalar.');
end
M = numel(Ac); n = size(Ac{1},1); m = size(Dw,2);
if nargin < 7; Chan = {}; end
if nargin < 8; Dwc  = {}; end
Chan = normchan(Chan, M, n);
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
Rg = g^2*eye(m) - Dw.'*Dw;
[~,pch] = chol((Rg+Rg.')/2);
if pch ~= 0; ok = false; return; end
tol = 1e-7; p = n*(n+1)/2;

% per-mode stabilizing solutions (monotone limit from zero terminal data)
Salg = cell(1,M); feas0 = true;
for ii = 1:M
    [Si, oki] = bre_limit(Ac{ii}, Cc{ii}, Chan{ii});
    Salg{ii} = Si; feas0 = feas0 && oki;
end
if ~feas0; ok = false; return; end             % gamma below a per-mode gain

% The violation objective is a non-smooth max of eigenvalue terms over M*p
% variables; a single Nelder-Mead run from one start can stall in a basin
% far from a feasible certificate even when one exists (fv plateaus around a
% small positive value at every gamma, so a single-start acceptance test is
% never satisfied and any bisection on gamma runs away).  We therefore
% MULTI-START.  The dominant coupling is the switch condition
% S_j(0) <= S_i(Td), which is relaxed by inflating every terminal storage by
% a common multiple of the identity; we sweep that inflation deterministically
% (scale-identity starts) and add a few reproducible random restarts, keeping
% the global RNG untouched via a local stream.  We accept as soon as any run
% reaches feasibility.
oo = optimset('Display','off','MaxFunEvals',20000,'MaxIter',20000, ...
              'TolFun',1e-10,'TolX',1e-10);
scales = [0.05 0.2 0.5 1.0 2.0];               % common identity inflations
vI = mvst_vech(eye(n));
rs = RandStream('twister','Seed',20240808);
nrand = 16;
fv = inf; ok = false;
for s = 1:(numel(scales)+nrand)
    if s <= numel(scales)
        x0 = repmat(scales(s)*vI, M, 1);        % same inflation on every mode
    else
        x0 = 0.5*randn(rs, M*p, 1);
    end
    [~, fvs] = fminsearch(@viol, x0, oo);
    fv = min(fv, fvs);
    if fv <= tol; ok = true; break; end        % feasible certificate found
end

    function v = viol(x)
        St = cell(1,M); S0 = cell(1,M); v = -inf;
        for i = 1:M
            D = mvst_unvech(x((i-1)*p+1:i*p), n);
            St{i} = Salg{i} + (D+D.')/2;                % S_i(Td)
            v = max(v, maxeig(bri(Ac{i}, Cc{i}, St{i}, Chan{i})));   % (iii)
            [S0{i}, tr] = back(Ac{i}, Cc{i}, St{i}, Chan{i});    % (i)
            v = max(v, tr);                             % (iv) + escape
        end
        for i = 1:M
            for j = 1:M
                if j ~= i
                    v = max(v, maxeig(S0{j} - St{i}));  % (ii)
                end
            end
        end
    end

    function G = bri(A, C, S, chi)
        if ~hasNoise
            w = S*Bw + C.'*Dw;
            G = A.'*S + S*A + C.'*C + w*(Rg\w.'); G = (G+G.')/2;
        else
            mult = zeros(n); w = S*Bw + C.'*Dw; Th = Rg;
            for l = 1:qd
                mult = mult + chi{l}.'*S*chi{l};
                w  = w  + chi{l}.'*S*Dwc{l};
                Th = Th - Dwc{l}.'*S*Dwc{l};
            end
            [~,pc] = chol((Th+Th.')/2);
            if pc ~= 0                       % pivot lost: dominant positive term
                G = A.'*S + S*A + C.'*C + mult + 1e6*eye(n);
            else
                G = A.'*S + S*A + C.'*C + mult + w*(Th\w.');
            end
            G = (G+G.')/2;
        end
    end

    function [S0, tr] = back(A, C, ST, chi)
        % backward BRE on [0,Td] from S(Td)=ST; tr = max over the grid of
        % (negativity of S, escape) as a violation measure
        opts = odeset('RelTol',1e-7,'AbsTol',1e-9, ...
                      'Events',@(s,q) mvst_odecap(q,1e6));
        [ss, Y] = ode45(@(s,q) rhsb(q,A,C,chi), linspace(0,Td,25), ST(:), opts);
        if ss(end) < Td - 1e-9
            S0 = ST; tr = 1e3; return                   % escaped: infeasible
        end
        tr = -inf;
        for k = 1:size(Y,1)
            Sk = reshape(Y(k,:).',n,n); Sk = (Sk+Sk.')/2;
            tr = max(tr, maxeig(-Sk));                  % (iv)
        end
        S0 = reshape(Y(end,:).',n,n); S0 = (S0+S0.')/2;
    end

    function dq = rhsb(q, A, C, chi)
        % s = Td - tau (backward in the timer): dS/ds = +BRI RHS
        S = reshape(q,n,n); S = (S+S.')/2;
        G = bri(A, C, S, chi); dq = G(:);
    end

    function [S, okl] = bre_limit(A, C, chi)
        % stabilizing solution by growing-horizon backward integration
        S = zeros(n); okl = false;
        opts = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
                      'Events',@(s,q) mvst_odecap(q,1e6));
        [ss, Y] = ode45(@(s,q) rhsb(q,A,C,chi), [0 80], zeros(n*n,1), opts);
        if ss(end) < 80 - 1e-6; return; end
        S = reshape(Y(end,:).',n,n); S = (S+S.')/2;
        okl = (maxeig(bri(A,C,S,chi)) < 1e-5);          % settled
    end

    function e = maxeig(X)
        e = max(eig((X+X.')/2));
    end
end

% ------------------------------------------------------------------------
function C = normchan(Chan, M, ~)
%NORMCHAN  Normalise per-mode diffusion channels to a 1xM cell of 1xq cells.
if isempty(Chan); C = repmat({{}}, 1, M); return; end
if ~iscell(Chan)
    error('mvst:diss_dwell_cert:chan', 'CHAN must be a 1 x M cell (per mode).');
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
