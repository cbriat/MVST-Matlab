function [P, K, info] = mvst_care_ct_poisson(A, B, Cc, Dc, Jc, lam, Q, R)
%MVST_CARE_CT_POISSON  Algebraic generalized Riccati for GBM+Poisson systems.
%   [P,K,INFO] = MVST_CARE_CT_POISSON(A,B,CC,DC,JC,LAM,Q,R) computes the
%   stabilizing solution of the algebraic Riccati equation of
%   MVST_DRE_CT_POISSON by the monotone-limit construction: integrate the
%   differential equation backward from P(T)=0 over a growing horizon until
%   the residual is negligible.  INFO returns .residual and .horizon.
%
%   Inputs: A (n x n drift), B (n x m control), CC/DC (1 x q cells of the
%   diffusion state/input channels), JC (1 x p cell of jump maps), LAM
%   (1 x p rate vector), Q (n x n state weight), R (m x m input weight).
%
%   See also MVST_DRE_CT_POISSON, MVST_CARE_CT, MVST_L2GAIN_CT_POISSON.
if size(A,1) ~= size(A,2)
    error('mvst:care_ct_poisson:notSquare', ...
        'A must be square; got %d x %d.', size(A,1), size(A,2));
end
if size(B,1) ~= size(A,1)
    error('mvst:care_ct_poisson:dim', ...
        'B must have %d rows to match A.', size(A,1));
end
if ~isequal(size(Q), [size(A,1) size(A,1)])
    error('mvst:care_ct_poisson:dim', ...
        'Q must be %d x %d to match A.', size(A,1), size(A,1));
end
n = size(A,1);
% tol sits above the ode-integration accuracy floor (RelTol 1e-8 accumulated
% over a long horizon), so slow-but-converged solves do not raise warnings.
tol = 1e-6; T = 5; Tmax = 800; res = inf; P = zeros(n); K = zeros(size(B,2),n);
while res > tol && T <= Tmax
    [~, Pt, Kt, ok] = mvst_dre_ct_poisson(A,B,Cc,Dc,Jc,lam,Q,R,T,zeros(n),3);
    if ~ok
        warning('mvst:care_ct_poisson:escape', ...
            'Backward Riccati escaped in finite time: not mean-square stabilizable.');
        break
    end
    P = Pt(:,:,1); K = Kt(:,:,1);
    % residual of the algebraic equation at P
    Th = R; S = B.'*P; mult = zeros(n);
    for l = 1:numel(Cc)
        Th = Th + Dc{l}.'*P*Dc{l}; S = S + Dc{l}.'*P*Cc{l};
        mult = mult + Cc{l}.'*P*Cc{l};
    end
    pois = zeros(n);
    for r = 1:numel(Jc); pois = pois + lam(r)*(Jc{r}.'*P*Jc{r} - P); end
    G = A.'*P + P*A + mult + pois + Q - S.'*(Th\S);
    res = norm((G+G.')/2,'fro');
    T = T*1.6;
end
info = struct('residual',res,'horizon',T/1.6);
if res > tol
    warning('mvst:care_ct_poisson:converge', ...
        'Riccati residual %.2e exceeds tol %.1e; check mean-square stabilizability.', res, tol);
end
end
