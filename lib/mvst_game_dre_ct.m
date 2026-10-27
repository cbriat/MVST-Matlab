function [t, Pt, Kut, Kwt] = mvst_game_dre_ct(Afun, Bu, Bw, Q, Ru, gam, tspan, Pterm, reltol, Cc, Dwc)
%MVST_GAME_DRE_CT  Backward game (HJI) differential Riccati, continuous time.
%   [t,Pt,Kut,Kwt] = MVST_GAME_DRE_CT(AFUN,BU,BW,Q,RU,GAM,TSPAN) integrates
%   the state-feedback Hamilton-Jacobi-Isaacs differential Riccati equation
%   of the zero-sum game (minimizer control u, maximizer disturbance w)
%
%     -Pdot = A(t)'P + P A(t) + Q
%             - P Bu Ru^{-1} Bu' P + (1/gam^2) P Bw Bw' P,
%
%   backward on TSPAN=[t0 T] from P(T)=PTERM (default 0).  This is the I-HJI
%   GRE after eliminating the (indefinite) Isaacs block diag(Ru,-gam^2 I).
%   AFUN may be a constant matrix or a function handle @(t) (time-varying
%   drift), evaluated via MVST_EVALC.  Returns increasing times t, the dual
%   Pt(:,:,k)=P(t(k)), and the saddle-point gains
%       Ku = -Ru^{-1} Bu' P     (optimal control u = Ku x),
%       Kw =  gam^{-2} Bw' P     (worst-case disturbance w = Kw x);
%   with the diffusion channels below, Kw = Theta^{-1}(Bw' P + sum_l Dwc{l}' P Cc{l})
%   (the w-block of the saddle set eq:game:Ks of the Game chapter).
%   The infinite-horizon game value is obtained by taking T large (P settles
%   to the stabilizing solution when gam exceeds the critical value gam*).
%
%   [...] = MVST_GAME_DRE_CT(...,CC,DWC) adds the Brownian diffusion channels of
%   the multiplicative-noise game
%     dx = (A x + Bu u + Bw w)dt + sum_l (Cc{l} x + Dwc{l} w) dW_l:
%   the game Riccati gains the state-multiplicative fold sum_l Cc{l}'P Cc{l} and,
%   when disturbance-multiplicative channels DWC are present, the maximizer's
%   quadratic form uses the storage-dependent Hessian
%   Theta = gam^2 I - sum_l Dwc{l}'P Dwc{l} in place of gam^2 I: the disturbance
%   term becomes (P Bw + sum_l Cc{l}'P Dwc{l}) Theta^{-1} (.)'.  Loss of
%   concavity (Theta not positive definite) returns a non-finite derivative so
%   the boundedness bisection of MVST_GAME_GAMMA_CT marks the level infeasible.
%   With CC and DWC absent or empty the noise-free equation above is integrated
%   unchanged (bit-identical to the pre-extension call).
%
%   See also MVST_GAME_GAMMA_CT, MVST_GAME_DRE_CT_POISSON, MVST_DRE_CT, MVST_EVALC.
if size(Q,1) ~= size(Q,2)
    error('mvst:game_dre_ct:notSquare', 'Q must be square.');
end
if size(Bu,1) ~= size(Q,1) || size(Bw,1) ~= size(Q,1)
    error('mvst:game_dre_ct:dim', 'BU and BW must have %d rows to match Q.', size(Q,1));
end
if isnumeric(Afun) && ~isequal(size(Afun), [size(Q,1) size(Q,1)])
    error('mvst:game_dre_ct:dim', 'A must be %d x %d to match Q.', size(Q,1), size(Q,1));
end
if nargin < 9 || isempty(reltol); reltol = 1e-8; end
if nargin < 10; Cc  = {}; end
if nargin < 11; Dwc = {}; end
if ~iscell(Cc);  Cc  = tocell3(Cc);  end
if ~iscell(Dwc); Dwc = tocell3(Dwc); end
n = size(Q,1);
if nargin < 8 || isempty(Pterm); Pterm = zeros(n); end
% Structurally all-zero channels stay on the deterministic exact path.
hasNoise = anynz(Cc) || anynz(Dwc);
hasDwc   = anynz(Dwc);
if hasNoise
    qd = max(numel(Cc), numel(Dwc)); m = size(Bw,2);
    if numel(Cc)  < qd; Cc(end+1:qd)  = {zeros(n)};   end
    if numel(Dwc) < qd; Dwc(end+1:qd) = {zeros(n,m)}; end
end
BuRu = Bu*(Ru\Bu.');                 % Bu Ru^{-1} Bu'
BwBw = Bw*Bw.';
odef = @(tt,pv) reshape(-gres(tt,pv,Afun,BuRu,BwBw,Q,gam,n,Bw,Cc,Dwc,hasNoise,hasDwc), [], 1);
opts = odeset('RelTol',reltol,'AbsTol',reltol*1e-2);
[ts, Pv] = ode45(odef, [tspan(2) tspan(1)], Pterm(:), opts);
[t, idx] = sort(ts); Pv = Pv(idx,:);
nt = numel(t);
Pt = zeros(n,n,nt); Kut = zeros(size(Bu,2),n,nt); Kwt = zeros(size(Bw,2),n,nt);
for k = 1:nt
    Pk = reshape(Pv(k,:), n, n); Pk = (Pk+Pk.')/2; Pt(:,:,k) = Pk;
    Kut(:,:,k) = -(Ru\(Bu.'*Pk));
    if hasDwc
        % Maximizer gain from the Isaacs block: the Ito term of the w-channels
        % adds sum_l Dwc{l}'P Cc{l} to the cross block and loads the pivot
        % Theta = gam^2 I - sum_l Dwc{l}'P Dwc{l} (fix 2026-10-04: the
        % previous gam^{-2} Bw'P ignored both).
        m = size(Bw,2); vk = Bw.'*Pk; Th = gam^2*eye(m);
        for l = 1:numel(Cc)
            vk = vk + Dwc{l}.'*Pk*Cc{l};
            Th = Th - Dwc{l}.'*Pk*Dwc{l};
        end
        Kwt(:,:,k) = Th\vk;
    else
        Kwt(:,:,k) = (Bw.'*Pk)/gam^2;
    end
end
end

function dP = gres(t, pv, Afun, BuRu, BwBw, Q, gam, n, Bw, Cc, Dwc, hasNoise, hasDwc)
%GRES  +RHS (i.e. -Pdot) of the game differential Riccati at time t.
A = mvst_evalc(Afun, t);
P = reshape(pv, n, n); P = (P+P.')/2;
if ~hasNoise
    dP = A.'*P + P*A + Q - P*BuRu*P + (P*BwBw*P)/gam^2;
else
    mult = zeros(n);
    for l = 1:numel(Cc); mult = mult + Cc{l}.'*P*Cc{l}; end
    if ~hasDwc
        dP = A.'*P + P*A + Q + mult - P*BuRu*P + (P*BwBw*P)/gam^2;
    else
        m = size(Bw,2); v = P*Bw; Th = gam^2*eye(m);
        for l = 1:numel(Cc)
            v  = v  + Cc{l}.'*P*Dwc{l};
            Th = Th - Dwc{l}.'*P*Dwc{l};
        end
        [~,pc] = chol((Th+Th.')/2);
        if pc ~= 0; dP = nan(n); return; end   % maximizer wins at this level
        dP = A.'*P + P*A + Q + mult - P*BuRu*P + v*(Th\v.');
    end
end
dP = (dP + dP.')/2;
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
