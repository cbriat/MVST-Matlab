function [G, K, Theta, S] = mvst_ric_ct_tv(t, P, md, Pcoupling)
%MVST_RIC_CT_TV  Time-varying continuous-time I-GRE residual and gain.
%   [G,K,Theta,S] = MVST_RIC_CT_TV(T,P,MD) is the time-varying counterpart
%   of MVST_RIC_CT: the mode coefficients (A,B,Q,R,S,Ac,Bc) may be constant
%   matrices OR function handles of time, and are evaluated at time T via
%   MVST_EVALC.  The formulas are identical to MVST_RIC_CT:
%
%       G(t,P) = He[A(t)'P] + sum_l Ac_l(t)'P Ac_l(t) + Q(t) - S'*Theta^{-1}*S,
%       Theta  = R(t) + sum_l Bc_l(t)'P Bc_l(t),
%       S      = B(t)'P + Scross(t)' + sum_l Bc_l(t)'P Ac_l(t),
%       K      = -Theta^{-1} S.
%
%   The differential I-GRE is Pdot = -G(t,P).  MVST_RIC_CT_TV(T,P,MD,PCOUP)
%   adds the Markov-coupling term PCOUP (e.g. sum_j pi_ij P_j) to G.
%
%   See also MVST_RIC_CT, MVST_DRE_CT_TV, MVST_EVALC, MVST_MODE.
if size(P,1) ~= size(P,2)
    error('mvst:ric_ct_tv:dim', 'P must be square.');
end
A = mvst_evalc(md.A,t);  B = mvst_evalc(md.B,t);
Q = mvst_evalc(md.Q,t);  R = mvst_evalc(md.R,t);  Scr = mvst_evalc(md.S,t);
Q = (Q+Q.')/2;  R = (R+R.')/2;

Theta = R;
S     = B.'*P + Scr.';
acc   = A.'*P + P*A + Q;

Ac = mvst_evalc(md.Ac,t);  Bc = mvst_evalc(md.Bc,t);
for l = 1:size(Ac,3)
    Acl = Ac(:,:,l);  Bcl = Bc(:,:,l);
    Theta = Theta + Bcl.'*P*Bcl;
    S     = S     + Bcl.'*P*Acl;
    acc   = acc   + Acl.'*P*Acl;
end
Theta = (Theta + Theta.')/2;
K = -(Theta \ S);
G = acc - S.'*(Theta \ S);
if nargin >= 4 && ~isempty(Pcoupling); G = G + Pcoupling; end
G = (G + G.')/2;
end
