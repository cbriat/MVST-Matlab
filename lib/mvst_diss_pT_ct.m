function [Tgrid, P0] = mvst_diss_pT_ct(A, Bw, C, Dw, gam, Tmax, nT)
%MVST_DISS_PT_CT  Finite-horizon storage matrix versus horizon length.
%   [TGRID,P0] = MVST_DISS_PT_CT(A,BW,C,DW,GAM,TMAX,NT) illustrates the
%   convergence of the finite-horizon available-storage matrix to its
%   infinite-horizon limit as the horizon T grows.
%
%   For each T in TGRID = linspace(0.1,TMAX,NT) the bounded-real Riccati
%   differential equation
%       -Pdot = A'P + P A + C'C
%               + (P Bw + C' Dw) (gam^2 I - Dw'Dw)^{-1} (Bw'P + Dw'C),
%   is integrated BACKWARD on [0,T] from the terminal condition P(T)=0.
%   Using the substitution s = T - t (so Pdot_s = +RHS) the integration
%   runs FORWARD in s from 0 to T; P is symmetrised after the solve.
%
%   P0(:,:,k) is P(0) for horizon TGRID(k).  In the monograph sign
%   convention the available-storage matrix is -P0; this function RETURNS
%   P0 = +P(0) (so callers plot -P0 for the storage matrix).
%
%   TGRID is 1 x NT and P0 is n x n x NT.  No toolbox is required.
%
%   See also MVST_L2GAIN_CT, MVST_CARE_CT, MVST_SOLVE_FINITE.
if size(A,1) ~= size(A,2)
    error('mvst:diss_pT_ct:notSquare', 'A must be square.');
end
if size(Bw,1) ~= size(A,1)
    error('mvst:diss_pT_ct:dim', 'BW must have %d rows to match A.', size(A,1));
end
if size(C,2) ~= size(A,1)
    error('mvst:diss_pT_ct:dim', 'C must have %d columns to match A.', size(A,1));
end
n  = size(A,1);
nz = size(C,1);
R  = gam^2*eye(nz) - Dw.'*Dw;       % must be positive definite (gam > norm(Dw))
Ri = inv(R);                        % small (nz x nz); precompute once
CtC = C.'*C;
CtD = C.'*Dw;                       % n x nz

Tgrid = linspace(0.1, Tmax, nT);
P0 = zeros(n, n, nT);

opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

for k = 1:nT
    T = Tgrid(k);
    % Integrate forward in s = T - t from 0 to T, starting at P(s=0)=P(T)=0.
    [~, Y] = ode45(@(s,p) bre_rhs(p, A, Bw, Ri, CtC, CtD, n), [0 T], ...
                   zeros(n*n,1), opts);
    Pend = reshape(Y(end,:).', n, n);     % P at s=T  <=>  P at t=0
    Pend = (Pend + Pend.')/2;             % symmetrise
    P0(:,:,k) = Pend;
end
end

function dp = bre_rhs(p, A, Bw, Ri, CtC, CtD, n)
%BRE_RHS  +RHS of the bounded-real Riccati ODE in reverse time s = T - t.
P = reshape(p, n, n);
P = (P + P.')/2;                          % symmetrise the state each step
PBw_CtD = P*Bw + CtD;                     % n x nz
dP = A.'*P + P*A + CtC + PBw_CtD*(Ri*PBw_CtD.');
dP = (dP + dP.')/2;
dp = dP(:);
end
