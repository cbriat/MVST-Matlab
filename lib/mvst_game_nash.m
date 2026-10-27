function res = mvst_game_nash(sys, opts)
%MVST_GAME_NASH  Infinite-horizon N-player Nash equilibrium via the coupled
%   I-GRE system (Chapter "Game", prop:game:nash:verify, eq:game:nash:GRE).
%   Each player i solves its best-response optimal-control problem given the
%   opponents' gains K_{-i}; imposing best response for all players is the
%   coupled system.  Solved by a Gauss-Seidel iteration over the players that
%   mirrors the MJLS coupling of MVST_SOLVE_INFINITE: at each sweep, player i's
%   opponents' loops are closed at their current gains (the congruence
%   T_i(K_{-i}) of eq:game:nash:Ti) and player i's reduced CARE is solved.
%
%   RES = MVST_GAME_NASH(SYS, OPTS).  SYS fields:
%     .A        common drift (n x n); .Ac (n x n x L) state-multiplicative
%               Brownian channels (optional);
%     .Bu       1 x N cell, player i's control input Bu{i} (n x m_i);
%     .Q        1 x N cell, player i's state cost Q{i} (n x n, symmetric);
%     .R        N x N cell (or 1 x N cell of 1 x N cells): R{i}{j} is player i's
%               weight on control j; R{i}{i} is player i's own control weight;
%     .X0       initial state second moment (optional; default I) for the costs.
%   For a ZERO-SUM two-player game (H-infinity) set, with control u and
%   disturbance w at level gamma: Bu={Bu,Bw}, Q={Q,-Q},
%   R={{Ru,-gamma^2 I},{-Ru,gamma^2 I}}; the Nash equilibrium is the saddle,
%   P{1} is the H-infinity dual, P{2}=-P{1}, K{1}=K, K{2}=L.
%
%   OPTS: .maxit (200), .tol (1e-9), plus CARE options forwarded per player.
%   RES: .P (1 x N cell), .K (1 x N cell), .residual (max player I-GRE residual),
%        .iters, .cost (1 x N, J^{(i)}=<P_i,X0>), .N.
%
%   Validation anchors: N=1 reduces to MVST_SOLVE_INFINITE; a zero-sum two-player
%   instance reproduces MVST_GAME_GAMMA's saddle; the per-player I-GRE residuals
%   (MVST_RIC_CT of each reduced problem) vanish at the equilibrium.
%
%   See also MVST_SOLVE_INFINITE, MVST_GAME_GAMMA, MVST_CARE_CT, MVST_RIC_CT.

    if nargin < 2; opts = struct(); end
    A  = mvst_getf(sys,'A', mvst_getf(sys,'A0',[]));
    n  = size(A,1);
    Ac = mvst_getf(sys,'Ac', zeros(n,n,0));
    Bu = mvst_getf(sys,'Bu',{}); if ~iscell(Bu); Bu={Bu}; end
    Q  = mvst_getf(sys,'Q',{});  if ~iscell(Q);  Q={Q};   end
    R  = mvst_getf(sys,'R',{});
    X0 = mvst_getf(sys,'X0', eye(n));
    N  = numel(Bu);
    maxit = mvst_getf(opts,'maxit',200); tol = mvst_getf(opts,'tol',1e-9);

    K = cell(1,N); for i=1:N; K{i} = zeros(size(Bu{i},2), n); end
    P = repmat({zeros(n)},1,N);
    it = 0; delta = inf;
    while it < maxit && delta > tol
        it = it + 1; Pold = P;
        for i = 1:N
            md = reduced_mode(A,Ac,Bu,Q,R,K,i,n);
            [P{i}, K{i}] = mvst_care_ct(md, opts);
        end
        delta = 0; for i=1:N; delta = max(delta, norm(P{i}-Pold{i},'fro')); end
    end

    % verification residuals: each reduced I-GRE residual at the equilibrium
    resid = 0;
    for i = 1:N
        md = reduced_mode(A,Ac,Bu,Q,R,K,i,n);
        G = mvst_ric_ct(P{i}, md);
        resid = max(resid, norm(G,'fro'));
    end
    cost = zeros(1,N); for i=1:N; cost(i) = trace(P{i}*X0); end
    res = struct('P',{P},'K',{K},'residual',resid,'iters',it,'cost',cost,'N',N);
end

% ----------------------------------------------------------------------- %
function md = reduced_mode(A,Ac,Bu,Q,R,K,i,n)
% Player i's best-response mode: opponents' loops closed at K_{-i}
% (congruence T_i(K_{-i}), eq:game:nash:Ti); drift A+sum_{j~=i} Bu_j K_j,
% state cost Q_i+sum_{j~=i} K_j' R_ij K_j, own control weight R_ii.
    N = numel(Bu);
    Acl = A; Qcl = Q{i};
    for j = 1:N
        if j == i; continue; end
        Acl = Acl + Bu{j}*K{j};
        Qcl = Qcl + K{j}.'*Rij(R,i,j)*K{j};
    end
    Qcl = (Qcl + Qcl.')/2;
    md = mvst_mode(Acl, Bu{i}, Qcl, Rij(R,i,i));
    if ~isempty(Ac) && size(Ac,3) > 0
        md = mvst_mode(Acl, Bu{i}, Qcl, Rij(R,i,i), 'Ac', Ac);
    end
end

function W = Rij(R,i,j)
% Player i's weight on control j from the R cell (R{i}{j}), robust to the
% 1 x N cell-of-cells or N x N cell layout.
    Ri = R{i};
    if iscell(Ri); W = Ri{j}; else; W = Ri; end   % scalar/matrix R{i} => own only
end
