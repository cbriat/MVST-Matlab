function res = mvst_solve_infinite(sys, opts)
%MVST_SOLVE_INFINITE  Infinite-horizon linear-cost optimal control (CT).
%   RES = MVST_SOLVE_INFINITE(SYS) solves the infinite-horizon I-GRE for a
%   continuous-time system SYS with one or several modes:
%     SYS.modes : 1xM cell array of mode structs (MVST_MODE)
%     SYS.Pi    : M x M Markov generator coupling the modes (optional;
%                 default 0 => independent / single mode)
%     SYS.pi0   : 1xM initial mode distribution (optional; default e_1)
%   For M=1 this is the stochastic algebraic Riccati equation; for M>1 it is
%   the coupled Markov-jump Riccati system, solved by a Gauss-Seidel
%   iteration in which the off-diagonal coupling sum_{j~=i} pi_ij P_j is
%   folded into the cost and the rate pi_ii into the drift.
%   RES has fields .P (1xM cell), .K (1xM cell), .cost, .iters, .residual.
%
%   Example:
%       S = mvst.System.ct([-1 0;0 -2], [1;1], 'Q', eye(2), 'R', 1);
%       r = mvst_solve_infinite(S.thin());   r.cost
%
%   See also MVST_SOLVE_FINITE, MVST_SOLVE_AVERAGE, MVST_CARE_CT, MVST_RIC_CT.
if nargin < 2; opts = struct(); end
if ~isstruct(sys) || ~isfield(sys, 'modes')
    error('mvst:solve_infinite:sys', ...
        'SYS must be a thin system struct with a .modes cell (see MVST_SYSTEM).');
end
M  = numel(sys.modes);
n  = sys.modes{1}.n;
Pi = mvst_getf(sys,'Pi',zeros(M));
pi0= mvst_getf(sys,'pi0',[1 zeros(1,M-1)]);
maxit = mvst_getf(opts,'maxit',300);
tol   = mvst_getf(opts,'tol',1e-8);

P = repmat({zeros(n)}, 1, M);
it = 0; delta = inf;
while it < maxit && delta > tol
    it = it + 1; Pold = P;
    for i = 1:M
        md  = sys.modes{i};
        mdi = md;
        mdi.A = md.A + 0.5*Pi(i,i)*eye(n);     % fold pi_ii into the drift
        Qadd = zeros(n);
        for j = 1:M
            if j ~= i; Qadd = Qadd + Pi(i,j)*P{j}; end
        end
        Qadd = (Qadd + Qadd.')/2;
        mdi.Q = md.Q + Qadd;                    % fold coupling into the cost
        P{i} = mvst_care_ct(mdi, opts);
    end
    delta = 0;
    for i = 1:M; delta = max(delta, norm(P{i}-Pold{i},'fro')); end
end

K = cell(1,M); res_i = 0;
for i = 1:M
    [G, K{i}] = mvst_ric_ct(P{i}, sys.modes{i}, coupling(Pi,P,i));
    res_i = max(res_i, norm(G,'fro'));
end
cost = 0;
for i = 1:M; cost = cost + pi0(i)*trace(P{i}*sys.modes{i}.X0); end

res.P = P; res.K = K; res.cost = cost; res.iters = it; res.residual = res_i;
end

function C = coupling(Pi, P, i)
M = numel(P); C = zeros(size(P{1}));
for j = 1:M; C = C + Pi(i,j)*P{j}; end
end
