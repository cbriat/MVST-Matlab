function [tg, P, K, ok] = mvst_dre_ct_poisson(A, B, Cc, Dc, Jc, lam, Q, R, T, PT, nT)
%MVST_DRE_CT_POISSON  Backward LQ Riccati for GBM + Poissonian-jump systems.
%   [TG,P,K,OK] = MVST_DRE_CT_POISSON(A,B,CC,DC,JC,LAM,Q,R,T,PT,NT) integrates
%   the generalized differential Riccati equation
%       -Pdot = A'P + P A + sum_l Cc{l}'P Cc{l}
%               + sum_r LAM(r) (Jc{r}'P Jc{r} - P)
%               + Q - S' Theta^{-1} S,
%       Theta = R + sum_l Dc{l}'P Dc{l},  S = B'P + sum_l Dc{l}'P Cc{l},
%   backward on [0,T] from P(T)=PT (substitution s = T-t, forward in s).
%   CC,DC,JC are cell arrays (possibly empty), LAM a vector of rates.
%   TG is a 1 x NT time grid on [0,T]; P,K are n x n x NT and m x n x NT
%   sampled at TG (in forward time t).  No toolbox is required.
%
%   OK is false when the integration escaped in finite time (norm cap hit,
%   e.g. for an indefinite terminal PT during a Newton line search); the
%   missing samples are then filled by holding the last computed value, so
%   callers should treat the output as unusable and reject the step.
%
%   See also MVST_CARE_CT_POISSON, MVST_DRE_CT, MVST_GAME_DRE_CT_POISSON.
if size(A,1) ~= size(A,2)
    error('mvst:dre_ct_poisson:notSquare', ...
        'A must be square; got %d x %d.', size(A,1), size(A,2));
end
if size(B,1) ~= size(A,1)
    error('mvst:dre_ct_poisson:dim', 'B must have %d rows to match A.', size(A,1));
end
if ~isequal(size(Q), [size(A,1) size(A,1)])
    error('mvst:dre_ct_poisson:dim', 'Q must be %d x %d to match A.', size(A,1), size(A,1));
end
n = size(A,1); m = size(B,2);
tg = linspace(0, T, nT);
cap = 1e8;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10, ...
              'Events',@(s,p) escape(p,cap));
ok = true;
try
    [ss, Y] = ode45(@(s,p) rhs(p), T - tg(end:-1:1), PT(:), opts);
catch
    ss = 0; Y = PT(:).';
end
if numel(ss) < nT || ss(end) < T - 1e-9
    ok = false;
    Y = [Y; repmat(Y(end,:), nT - size(Y,1), 1)];   % hold the last value
end
P = zeros(n,n,nT); K = zeros(m,n,nT);
for k = 1:nT
    Pk = reshape(Y(nT-k+1,:).', n, n); Pk = (Pk+Pk.')/2;
    P(:,:,k) = Pk;
    [~, Kk] = gain(Pk);
    K(:,:,k) = Kk;
end

    function dp = rhs(p)
        Pm = reshape(p,n,n); Pm = (Pm+Pm.')/2;
        [G,~] = gain(Pm);
        dp = G(:);
    end

    function [G, Kk] = gain(Pm)
        Th = R; S = B.'*Pm;
        mult = zeros(n);
        for l = 1:numel(Cc)
            Th = Th + Dc{l}.'*Pm*Dc{l};
            S  = S  + Dc{l}.'*Pm*Cc{l};
            mult = mult + Cc{l}.'*Pm*Cc{l};
        end
        pois = zeros(n);
        for r = 1:numel(Jc)
            pois = pois + lam(r)*(Jc{r}.'*Pm*Jc{r} - Pm);
        end
        Kk = -Th\S;
        G  = A.'*Pm + Pm*A + mult + pois + Q - S.'*(Th\S);
        G  = (G+G.')/2;
    end

    function [val, terminate, dir] = escape(p, capv)
        val = capv - norm(p); terminate = 1; dir = -1;
    end
end
