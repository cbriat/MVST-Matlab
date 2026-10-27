function [G, K, Theta, S] = mvst_ric_ct(P, md, Pcoupling)
%MVST_RIC_CT  Continuous-time I-GRE residual and optimal gain for one mode.
%   [G,K,Theta,S] = MVST_RIC_CT(P,MD) returns, for a symmetric P and a mode
%   MD built by MVST_MODE, the residual of the continuous-time I-GRE
%
%       G(P) = He[A'P] + sum_l Ac_l'P Ac_l + Q - S'*Theta^{-1}*S,
%       Theta = R + sum_l Bc_l'P Bc_l,
%       S     = B'P + Scross' + sum_l Bc_l'P Ac_l,
%
%   together with the optimal gain K = -Theta^{-1} S (control law u = K x),
%   where He[M] = M + M'. The differential I-GRE is Pdot = -G(P) with
%   P(T) = ZT; the infinite-horizon dual variable solves G(P) = 0.
%
%   [..] = MVST_RIC_CT(P,MD,PCOUPLING) adds the Markov-coupling term
%   PCOUPLING (an n x n matrix, e.g. sum_j pi_ij P_j) to G; this is used by
%   the coupled Markov-jump solver and defaults to zero.
%
%   See also MVST_DRE_CT, MVST_CARE_CT, MVST_RIC_CT_TV, MVST_MODE.
if ~isstruct(md) || ~isfield(md,'A')
    error('mvst:ric_ct:input', ...
        'MD must be a mode struct carrying at least .A (see MVST_MODE).');
end
if size(P,1) ~= size(md.A,1) || size(P,2) ~= size(md.A,1)
    error('mvst:ric_ct:dim', ...
        'P must be %d x %d to match MD.A.', size(md.A,1), size(md.A,1));
end
A = md.A; B = md.B; Q = md.Q; R = md.R; Scr = md.S;
Theta = R;
S     = B.'*P + Scr.';
acc   = A.'*P + P*A + Q;
for l = 1:md.Nc
    Acl = md.Ac(:,:,l);  Bcl = md.Bc(:,:,l);
    Theta = Theta + Bcl.'*P*Bcl;
    S     = S     + Bcl.'*P*Acl;
    acc   = acc   + Acl.'*P*Acl;
end
Theta = (Theta + Theta.')/2;
K = -(Theta \ S);
G = acc - S.'*(Theta \ S);
if nargin >= 3 && ~isempty(Pcoupling)
    G = G + Pcoupling;
end
G = (G + G.')/2;
end
