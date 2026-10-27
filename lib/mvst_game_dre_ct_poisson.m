function [t, Pt, Kut, Kwt] = mvst_game_dre_ct_poisson(A, Bu, Bw, Cc, Jc, lam, Q, Ru, gam, tspan, Pterm, reltol, BuJc, BwJc, Rdc)
%MVST_GAME_DRE_CT_POISSON  Backward game (HJI) Riccati for GBM+Poisson systems.
%   Same as MVST_GAME_DRE_CT but with the multiplicative-noise term
%   sum_l Cc{l}'P Cc{l} and the Poissonian jump-rate term added to the flow.
%
%   STATE-ONLY JUMP (default, BuJc/BwJc/Rdc omitted or empty), reset x -> J_r x:
%     -Pdot = A'P + P A + sum_l Cc{l}'P Cc{l} + sum_r lam_r (Jc{r}'P Jc{r} - P)
%             + Q - P Bu Ru^{-1} Bu' P + (1/gam^2) P Bw Bw' P.
%
%   JUMP WITH BOTH PLAYERS' CHANNELS (BuJc,BwJc,Rdc supplied as cells, one
%   entry per jump map r), reset x -> J_r x + BuJ_r u_J + BwJ_r w_J: at each
%   Poisson event the players solve a STATIC ISAACS SADDLE, minimizer over the
%   jump control u_J (weight Rd_r > 0), maximizer over the jump disturbance w_J
%   (weight -gam^2 I).  With B_r = [BuJ_r, BwJ_r] and D_r = blkdiag(Rd_r,
%   -gam^2 I), the saddle value uses the reduced post-jump operator
%     Phat_r = P - P B_r (B_r'P B_r + D_r)^{-1} B_r'P ,
%   and the jump-rate term becomes sum_r lam_r (Jc{r}'Phat_r Jc{r} - P).  The
%   saddle is well posed only while the maximizer subproblem stays concave
%   (BwJ_r'P BwJ_r - gam^2 I < 0); when it fails the maximizer wins at the jump
%   and the routine returns a non-finite derivative so the boundedness
%   bisection of MVST_GAME_GAMMA_CT_POISSON marks the level infeasible.
%
%   Returns increasing times t, Pt(:,:,k)=P(t(k)), and the FLOW saddle gains
%   Ku = -Ru^{-1}Bu'P and Kw = gam^{-2}Bw'P.
%
%   See also MVST_GAME_DRE_CT, MVST_GAME_GAMMA_CT_POISSON, MVST_DRE_CT_POISSON.
if size(A,1) ~= size(A,2)
    error('mvst:game_dre_ct_poisson:notSquare', 'A must be square.');
end
if size(Q,1) ~= size(A,1)
    error('mvst:game_dre_ct_poisson:dim', 'Q must be %d x %d to match A.', size(A,1), size(A,1));
end
if size(Bu,1) ~= size(A,1) || size(Bw,1) ~= size(A,1)
    error('mvst:game_dre_ct_poisson:dim', 'BU and BW must have %d rows to match A.', size(A,1));
end
if nargin < 12 || isempty(reltol); reltol = 1e-8; end
n = size(Q,1);
p = numel(Jc);
if nargin < 11 || isempty(Pterm); Pterm = zeros(n); end
if nargin < 13; BuJc = {}; end
if nargin < 14; BwJc = {}; end
if nargin < 15; Rdc  = {}; end
hasIO = ~isempty(BuJc) || ~isempty(BwJc);
if hasIO
    if isempty(BuJc); BuJc = repmat({zeros(n,0)},1,p); end
    if isempty(BwJc); BwJc = repmat({zeros(n,0)},1,p); end
    if isempty(Rdc);  Rdc  = repmat({[]},1,p);         end
end
BuRu = Bu*(Ru\Bu.'); BwBw = Bw*Bw.';
opts = odeset('RelTol',reltol,'AbsTol',reltol*1e-2);
[ts, Pv] = ode45(@(tt,pv) -rhs(pv), [tspan(2) tspan(1)], Pterm(:), opts);
[t, idx] = sort(ts); Pv = Pv(idx,:);
nt = numel(t);
Pt = zeros(n,n,nt); Kut = zeros(size(Bu,2),n,nt); Kwt = zeros(size(Bw,2),n,nt);
for k = 1:nt
    P = reshape(Pv(k,:).',n,n); P = (P+P.')/2;
    Pt(:,:,k) = P;
    Kut(:,:,k) = -(Ru\(Bu.'*P));
    Kwt(:,:,k) = (1/gam^2)*(Bw.'*P);
end

    function g = rhs(pv)
        P = reshape(pv,n,n); P = (P+P.')/2;
        mult = zeros(n);
        for l = 1:numel(Cc); mult = mult + Cc{l}.'*P*Cc{l}; end
        pois = zeros(n);
        for r = 1:p
            Phat = P;
            if hasIO && (~isempty(BuJc{r}) || ~isempty(BwJc{r}))
                BuJ = BuJc{r}; BwJ = BwJc{r};
                mwJ = size(BwJ,2);
                Rd  = Rdc{r}; if isempty(Rd); Rd = eye(size(BuJ,2)); end
                % concavity guard: maximizer subproblem must stay bounded
                if mwJ > 0
                    Ww = BwJ.'*P*BwJ - gam^2*eye(mwJ);
                    if max(real(eig(Ww))) >= -1e-9
                        g = nan(n*n,1); return;   % maximizer wins at the jump
                    end
                end
                Br = [BuJ, BwJ];
                Dr = blkdiag(Rd, -gam^2*eye(mwJ));
                Mr = Br.'*P*Br + Dr;
                if rcond(Mr) < 1e-12; g = nan(n*n,1); return; end
                Phat = P - P*Br*(Mr\(Br.'*P)); Phat = (Phat+Phat.')/2;
            end
            pois = pois + lam(r)*(Jc{r}.'*Phat*Jc{r} - P);
        end
        G = A.'*P + P*A + mult + pois + Q - P*BuRu*P + (1/gam^2)*P*BwBw*P;
        G = (G+G.')/2; g = G(:);
    end
end
