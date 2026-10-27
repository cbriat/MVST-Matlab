function [Q, S, R] = mvst_qsr_supply(supply, p, m, where)
%MVST_QSR_SUPPLY  Validate and unpack a quadratic (Q,S,R) supply contract.
%   [Q,S,R] = MVST_QSR_SUPPLY(SUPPLY, P, M) returns the weights of the supply
%   rate
%       s(z,w) = [z; w]' * [Q S; S' R] * [z; w]
%   for a performance output z in R^P and a disturbance w in R^M (Ch08 sign
%   convention: dissipativity is  d/dt <P,x> <= s ; the l2-gain supply is
%   Q = -I, S = 0, R = gamma^2 I, and the passivity supply with indices
%   (nu,rho) is Q = -rho I, S = I, R = -nu I).
%
%   CONTRACT  supply = struct('type','qsr','Q',Q,'S',S,'R',R)
%     Q  P x P symmetric,  S  P x M,  R  M x M symmetric.
%   Optional JUMP weights  .Qk .Sk .Rk  (used when WHERE = 'jump').  If ANY of
%   them is given, ALL THREE are required (error mvst:qsr:jumpweights) and
%   they are used.  When none is given the flow weights are reused provided
%   the jump dimensions (pj,md) equal (P,M); otherwise they are required.
%
%   [Q,S,R] = MVST_QSR_SUPPLY(SUPPLY, P, M, 'jump') selects the jump weights.
%   This is the single authoritative validation of the contract; the SOS
%   engines (diss, diss_modes, dissyn) and MVST.SYSTEM.DISSIPATIVITY use it.
%
%   See also MVST_DISS_SOS, MVST_DISS_SOS_MODES, MVST_DISSYN_SOS.
if nargin < 4 || isempty(where); where = 'flow'; end
sfx = ''; if strcmp(where,'jump'); sfx = 'k'; end
if ~isstruct(supply) || ~isfield(supply,'type') || ~strcmp(supply.type,'qsr')
    error('mvst:qsr:type','supply.type must be ''qsr''.');
end
Qn = ['Q' sfx]; Sn = ['S' sfx]; Rn = ['R' sfx];
if strcmp(where,'jump')
    hv = [isfield(supply,'Qk') isfield(supply,'Sk') isfield(supply,'Rk')];
    if any(hv) && ~all(hv)
        error('mvst:qsr:jumpweights', ...
            'partial jump weights: supply.Qk, Sk and Rk must all be given (or none).');
    end
end
if strcmp(where,'jump') && ~(isfield(supply,Qn) && isfield(supply,Sn) && isfield(supply,Rn))
    % reuse the flow weights iff the jump dimensions match
    if isfield(supply,'Q') && isfield(supply,'S') && isfield(supply,'R') && ...
            isequal(size(supply.Q),[p p]) && isequal(size(supply.R),[m m]) && ...
            isequal(size(supply.S),[p m])
        Qn = 'Q'; Sn = 'S'; Rn = 'R';
    else
        error('mvst:qsr:jumpweights', ['jump dimensions (pj=%d, md=%d) differ from ' ...
            'the flow ones: supply.Qk/Sk/Rk are required.'], p, m);
    end
end
for f = {Qn,Sn,Rn}
    if ~isfield(supply,f{1})
        error('mvst:qsr:missing','supply.%s is missing (type ''qsr'').', f{1});
    end
end
Q = supply.(Qn); S = supply.(Sn); R = supply.(Rn);
if p == 0 && isempty(Q); Q = zeros(0); end
if p == 0 && isempty(S); S = zeros(0,m); end
if m == 0 && isempty(R); R = zeros(0); end
if m == 0 && isempty(S); S = zeros(p,0); end
if ~isnumeric(Q) || ~isnumeric(S) || ~isnumeric(R) || ~isreal(Q) || ~isreal(S) || ~isreal(R)
    error('mvst:qsr:real','Q,S,R must be real numeric arrays.');
end
if ~isequal(size(Q),[p p]); error('mvst:qsr:sizeQ','%s must be %d x %d (output dim), got %s.',Qn,p,p,mat2str(size(Q))); end
if ~isequal(size(S),[p m]); error('mvst:qsr:sizeS','%s must be %d x %d, got %s.',Sn,p,m,mat2str(size(S))); end
if ~isequal(size(R),[m m]); error('mvst:qsr:sizeR','%s must be %d x %d (disturbance dim), got %s.',Rn,m,m,mat2str(size(R))); end
tolsym = 1e-10*max([1, norm(Q,'fro'), norm(R,'fro')]);
if norm(Q-Q','fro') > tolsym; error('mvst:qsr:symQ','%s must be symmetric.',Qn); end
if norm(R-R','fro') > tolsym; error('mvst:qsr:symR','%s must be symmetric.',Rn); end
if ~isequal(Q,Q'); Q = (Q+Q')/2; end      % exact data kept when already symmetric
if ~isequal(R,R'); R = (R+R')/2; end
end
