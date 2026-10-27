function gam = mvst_l2gain_ct_poisson(A, Bw, C, Dw, Cc, Dwc, Jc, lam, Dr, Thoriz, rtol)
%MVST_L2GAIN_CT_POISSON  Mean-square L2-gain of a GBM+Poisson system.
%   GAM = MVST_L2GAIN_CT_POISSON(A,BW,C,DW,CC,DWC,JC,LAM,DR,THORIZ) returns
%   the mean-square L2-gain by bisection for the Ito system
%       dx = (A x + Bw w)dt + sum_l (Cc{l} x + Dwc{l} w) dW_l
%            + sum_r ((Jc{r} - I) x + Dr{r} w) dEta_r,
%   with output z = C x + Dw w, Poisson rates LAM(r).  The disturbance may
%   drive the spontaneous resets through the jump-disturbance channels Dr{r}
%   (a cell matching Jc; empty or omitted for state-only resets).  Feasibility
%   of a trial gamma is decided by integrating the stochastic bounded-real
%   Riccati equation
%       -Sdot = A'S + S A + sum_l Cc{l}'S Cc{l}
%               + sum_r LAM(r)(Jc{r}'S Jc{r} - S) + C'C
%               + v(S) Theta(S)^{-1} v(S)',
%       v(S)     = S Bw + sum_l Cc{l}'S Dwc{l}
%                        + sum_r LAM(r) Jc{r}'S Dr{r} + C'Dw,
%       Theta(S) = gamma^2 I - Dw'Dw - sum_l Dwc{l}'S Dwc{l}
%                        - sum_r LAM(r) Dr{r}'S Dr{r},
%   backward from S(T)=0 over the horizon THORIZ (default 60): the gain is
%   <= gamma iff the solution stays bounded with Theta(S) positive
%   definite (no finite escape).  There is no Hamiltonian eigenvalue test
%   in the stochastic case, so this shooting test replaces it.  DWC may be
%   an empty cell for drift-only disturbance.  The finite-horizon escape test
%   returns a slight UNDER-estimate of the gain (escape times grow without
%   bound as gamma approaches the gain from below); a long THORIZ and a small
%   optional RTOL (11th argument, relative bisection tolerance, default 1e-3)
%   tighten it.
%
%   See also MVST_L2GAIN_CT, MVST_GAME_GAMMA_CT_POISSON, MVST_CARE_CT_POISSON.
if size(A,1) ~= size(A,2)
    error('mvst:l2gain_ct_poisson:notSquare', 'A must be square.');
end
if size(Bw,1) ~= size(A,1)
    error('mvst:l2gain_ct_poisson:dim', 'BW must have %d rows to match A.', size(A,1));
end
if size(C,2) ~= size(A,1)
    error('mvst:l2gain_ct_poisson:dim', 'C must have %d columns to match A.', size(A,1));
end
if nargin < 10 || isempty(Thoriz); Thoriz = 60; end
% optional relative bisection tolerance RTOL (default 1e-3, added 2026-10-04)
if nargin < 11 || isempty(rtol); rtol = 1e-3; end
n = size(A,1); m = size(Dw,2);
if nargin < 9 || isempty(Dr); Dr = repmat({zeros(n,m)}, 1, numel(Jc)); end
dnorm = norm(Dw);
lo = dnorm + 1e-6;
hi = max(2*lo, 1.0);
for kk = 1:60
    if feasible(hi); break; end
    hi = hi*2;
    if kk == 60
        warning('mvst:l2gain_ct_poisson:nohi','No feasible gamma bracketed.');
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
        opts = odeset('RelTol',1e-6,'AbsTol',1e-8, ...
                      'Events',@(s,p) escape(p,Rg,big));
        try
            [ss,Y] = ode45(@(s,p) rhs(p,Rg), [0 Thoriz], zeros(n*n,1), opts);
            % The escape() event stops before Theta degenerates; the isfinite
            % backstop guards the corner where ode45 returns ss(end)=Thoriz with
            % a NaN state (a NaN derivative does not truncate the time span).
            ok = (ss(end) >= Thoriz - 1e-6) && all(isfinite(Y(end,:)));
        catch
            ok = false;                       % integrator failure = escape
        end
    end

    function dp = rhs(p, Rg)
        S = reshape(p,n,n); S = (S+S.')/2;
        mult = zeros(n); v = S*Bw + C.'*Dw; Th = Rg;
        for l = 1:numel(Cc)
            mult = mult + Cc{l}.'*S*Cc{l};
            v  = v  + Cc{l}.'*S*Dwc{l};
            Th = Th - Dwc{l}.'*S*Dwc{l};
        end
        pois = zeros(n);
        for r = 1:numel(Jc)
            pois = pois + lam(r)*(Jc{r}.'*S*Jc{r} - S);
            v  = v  + lam(r)*(Jc{r}.'*S*Dr{r});
            Th = Th - lam(r)*(Dr{r}.'*S*Dr{r});
        end
        G = A.'*S + S*A + mult + pois + C.'*C + v*(Th\v.');
        G = (G+G.')/2; dp = G(:);
    end

    function [val, terminate, dir] = escape(p, Rg, big)
        % stop on norm blow-up or loss of pivot positive-definiteness
        S = reshape(p,n,n); S = (S+S.')/2;
        Th = Rg;
        for l = 1:numel(Cc); Th = Th - Dwc{l}.'*S*Dwc{l}; end
        for r = 1:numel(Jc); Th = Th - lam(r)*(Dr{r}.'*S*Dr{r}); end
        val = [big - norm(p); min(eig((Th+Th.')/2)) - 1e-9];
        terminate = [1; 1]; dir = [-1; -1];
    end
end
