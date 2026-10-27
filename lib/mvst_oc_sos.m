function res = mvst_oc_sos(sys, cost, dwell, opts)
%MVST_OC_SOS  Optimal control by the dual (storage-function) SOS/LMI program
%  for a SINGLE linear Hermitian matrix-valued mode, covering the continuous-
%  time (GBM), Poissonian-jump, and discrete-time classes of the Linear-Cost
%  Optimal Control chapter.  This is the optimal-control counterpart of the
%  SOS family MVST_STAB_SOS / MVST_STABZ_SOS / MVST_DISS_SOS: the certificate
%  P is the DUAL variable of the chapter, and the routine solves the dual LMI
%
%     J^*(X0) = max_P  <P,X0>   s.t.   C_Z(P) >= 0     (flow, CT),
%
%  the infinite-horizon dual LMI of Theorem th:IMP:LQ:inf:LMI (eq
%  eq:IMP:LQ:LMI:inf), whose continuous-time specialization is eq:CT:LQ:LMI
%  (Corollary of th:IMP:LQ:LMI) and whose discrete-time specialization is
%  eq:DT:LQ:LMI (Corollary cor:DT:LQ:fh).  The flow operator C_Z(P) is the
%  augmented KYP / Schur matrix of the I-GRE (Definition def:IMP:GRE, blocks
%  eq:IMP:GRE:fPSD/fC/fR; sufficiency Schur decomposition pf:suff:schur:flow),
%
%     C_Z(P) = [ (C_Z)_1  (C_Z)_2 ;  (C_Z)_2^*  (C_Z)_3 ] ,
%     (C_Z)_1 = He[A0^*P] + Q + sum_l Ac_l^*P Ac_l
%                 + sum_r lamf_r (Jf_r^*P Jf_r - P) ,
%     (C_Z)_2 = P B0 + S + sum_l Ac_l^*P Bc_l + sum_r lamf_r Jf_r^*P BJf_r ,
%     (C_Z)_3 = R + sum_l Bc_l^*P Bc_l + sum_r lamf_r BJf_r^*P BJf_r ,
%
%  so C_Z(P) >= 0 is the pre-Schur bounded-real form of the I-GRE and the
%  maximal feasible P is the stabilizing I-GRE solution P^star.  SOSTOOLS
%  poses C_Z(P) >= 0 as the SOS condition [x;u]^*C_Z(P)[x;u] in Sigma over the
%  indeterminates z = [x;u] (pvar mode), and the optimum recovers P^star and
%  the gain K^star = -(C_Z)_3^{-1}(C_Z)_2^* to solver accuracy; complementary
%  slackness C_Z(P^star)[I;K^star] = 0 then holds.  The Poisson terms
%  lamf_r (Jf_r^*P Jf_r - P) (eq:ex:poisson:riccati) enter P LINEARLY, exactly
%  as the multiplicative-noise terms Ac_l^*P Ac_l do, so the storage LMI stays
%  affine in the decision variable.  Solver: SOSTOOLS + SeDuMi.
%
%  RES = MVST_OC_SOS(SYS) uses zero cross-cost, R = I, X0 = I, Q = 0 unless the
%  input carries its own cost data.  RES = MVST_OC_SOS(SYS,COST,DWELL,OPTS).
%
%  SYS  (rich, or a THIN MVST_MODE struct routed through MVST_TORICH):
%    .A0  n x n         flow drift            (thin .A)
%    .B0  n x m         flow control input    (thin .B)
%    .Ac  n x n x q     flow state-multiplicative (Brownian/GBM) channels
%    .Bc  n x m x q     input-multiplicative channels (paired with Ac)
%    .Jf  n x n x pf    SPONTANEOUS (Poisson) jump maps
%    .lamf 1 x pf       Poisson rates (paired with Jf)
%    .BJf n x m x pf    spontaneous-jump input selectors (paired with Jf)
%  A thin mode built by MVST_MODE carries the cost blocks (.Q/.R/.S/.X0); they
%  are used when COST is omitted or empty.
%
%  COST  (struct; each field falls back to the SYS field, then to a default):
%    .Q   n x n   state cost               (default 0)
%    .R   m x m   control cost             (default I)
%    .S   n x m   cross cost               (default 0)
%    .X0  n x n   initial second moment    (default I); objective weight <P,X0>
%
%  DWELL  (struct; .type default 'arbitrary'):
%    'arbitrary' = continuous-time infinite-horizon dual LMI, C_Z(P) >= 0 with
%                  a constant certificate P (the GBM / Poisson class; the
%                  Poisson jumps live in the flow through Jf/lamf/BJf);
%    'discrete'  = discrete-time dual LMI eq:DT:LQ:LMI, D_Z(P,P) >= 0 with the
%                  transition A0, input B0, noise Ac/Bc.
%
%  OPTS  (struct or MVST.OPTIONS, optional; normalized by MVST_OPTS):
%    .solver  SDP backend 'sedumi' (default) | 'sdpt3' | 'mosek' (errors
%             'mvst:opts:solverMissing' if absent),  .fid  verbosity (default 0)
%    .psd     also impose P >= 0 as in the book dual LMI (default false: the
%             maximal-C_Z(P) formulation used in the chapter examples; P >= 0 is
%             inactive at the stabilizing optimum, so the two agree, but leaving
%             it off keeps the SDP identical to the plain storage program).
%
%  RES:
%    .P     n x n   certificate (the stabilizing I-GRE solution P^star)
%    .K     m x n   optimal gain K^star = -(C_Z)_3^{-1}(C_Z)_2^*
%    .cost  struct  .primal  = <P,X0>                        (the optimal value)
%                   .dual    = trace([Q S;S' R] Z^star), Z^star=[I;K]Xocc[I;K]'
%                              the moment (complementary-slackness) dual cost
%                   .slackres= ||C_Z(P^star)[I;K^star]||_F   (compl. slackness)
%    .feas  clean-solve flag (pinf==0 && dinf==0 && numerr<2)
%    .M     the KYP / Schur matrix C_Z(P^star)
%    .Theta the (C_Z)_3 block ;  .hist  SeDuMi solver-info struct
%
%  See also MVST_OC_SOS_MODES (coupled Markov-jump), MVST_CARE_CT,
%  MVST_CARE_CT_POISSON, MVST_STABZ_SOS, MVST_DISS_SOS.

if nargin < 2; cost  = []; end
if nargin < 3 || isempty(dwell); dwell = struct('type','arbitrary'); end
if nargin < 4; opts  = struct(); end
if ~isfield(dwell,'type') || isempty(dwell.type); dwell.type = 'arbitrary'; end
opts = mvst_opts(opts);   % normalize + validate .solver ('sedumi' default)

raw = sys;
if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
n = size(sys.A0,1);
if ~isfield(sys,'B0') || isempty(sys.B0); sys.B0 = zeros(n,0); end
m = size(sys.B0,2);
sys = fill_defaults(sys, n, m);
cst = extract_cost(cost, raw, n, m);

% -------- indeterminates z = [x;u] and the constant certificate P ----------
V = mpvar('z', [n+m, 1]);
prog = sosprogram(V);
[prog, P] = sospolymatrixvar(prog, monomials(V,0), [n n], 'symmetric');

% -------- dual LMI: max <P,X0> s.t. C_Z(P) >= 0 (Schur / KYP form) ----------
M = oc_kyp(P, sys, cst, zeros(n), dwell.type);
prog = sosineq(prog, V.'*M*V);                 % [x;u]^* C_Z(P) [x;u] in Sigma
if getf(opts,'psd',false)
    prog = sosineq(prog, V(1:n).'*P*V(1:n));    % book dual LMI also has P >= 0
end
prog = sossetobj(prog, -trace(cst.X0*P));       % maximize <P,X0>

o = mvst_sosopts(opts);
o.params.fid = getf(opts,'fid',0);
[prog, info] = mvst_sossolve(prog, o);

P = double(sosgetsol(prog, P)); P = (P + P.')/2;
[M, K, Theta] = oc_kyp(P, sys, cst, zeros(n), dwell.type);

res.P     = P;
res.K     = K;
res.cost  = dual_cost(P, K, M, sys, cst, dwell.type);
res.feas  = (info.pinf == 0 && info.dinf == 0 && info.numerr < 2);
res.M     = M;
res.Theta = Theta;
res.hist  = info;
end

% =========================================================================
function [M, K, Theta] = oc_kyp(P, sys, cst, coup, type)
%OC_KYP  The augmented KYP / Schur matrix C_Z(P) = [acc S'; S Theta] of the
%  I-GRE (def:IMP:GRE), with the Markov-coupling term COUP added to the (1,1)
%  block (COUP = sum_j pi_ij P_j for the coupled solver; 0 in the single case).
%  TYPE 'arbitrary' -> continuous-time C_Z; 'discrete' -> discrete-time D_Z.
q = size(sys.Ac,3);  pf = size(sys.Jf,3);
if strcmp(type,'discrete')
    Theta = cst.R;                 S = sys.B0.'*P*sys.A0 + cst.S.';
    acc   = sys.A0.'*P*sys.A0 - P + cst.Q + coup;
    for l = 1:q
        Ac = sys.Ac(:,:,l); Bc = sys.Bc(:,:,l);
        Theta = Theta + Bc.'*P*Bc;  S = S + Bc.'*P*sys.A0;  acc = acc + Ac.'*P*Ac;
    end
else                                % 'arbitrary' (continuous-time)
    Theta = cst.R;                 S = sys.B0.'*P + cst.S.';
    acc   = sys.A0.'*P + P*sys.A0 + cst.Q + coup;
    for l = 1:q
        Ac = sys.Ac(:,:,l); Bc = sys.Bc(:,:,l);
        Theta = Theta + Bc.'*P*Bc;  S = S + Bc.'*P*Ac;      acc = acc + Ac.'*P*Ac;
    end
    for r = 1:pf                     % spontaneous (Poisson) jumps
        Jf = sys.Jf(:,:,r); BJf = sys.BJf(:,:,r); lr = sys.lamf(r);
        acc   = acc   + lr*(Jf.'*P*Jf - P);
        S     = S     + lr*(BJf.'*P*Jf);
        Theta = Theta + lr*(BJf.'*P*BJf);
    end
end
Theta = (Theta + Theta.')/2;
M = [acc, S.'; S, Theta];  M = (M + M.')/2;
K = [];
if nargout >= 2; K = -(Theta \ S); end   % numeric P only (recovery), not pvar
end

% =========================================================================
function cost = dual_cost(P, K, M, sys, cst, type)
%DUAL_COST  Primal value <P,X0> and the dual moment cost trace([Q S;S' R] Z^*)
%  with Z^* = [I;K] Xocc [I;K]^* supported on u = K x (complementary slackness).
%  Xocc is the mean-square occupation second moment: the closed-loop generator
%  of the (jump-)diffusion applied to the second moment, solved as a (discrete
%  or continuous) generalized Lyapunov equation.
n = size(P,1);
q = size(sys.Ac,3);  pf = size(sys.Jf,3);
Acl = sys.A0 + sys.B0*K;
if strcmp(type,'discrete')
    Lyap = kron(Acl, Acl) - eye(n*n);
    for l = 1:q
        Accl = sys.Ac(:,:,l) + sys.Bc(:,:,l)*K;  Lyap = Lyap + kron(Accl, Accl);
    end
else
    Lyap = kron(eye(n), Acl) + kron(Acl, eye(n));      % vec(Acl X + X Acl')
    for l = 1:q
        Accl = sys.Ac(:,:,l) + sys.Bc(:,:,l)*K;  Lyap = Lyap + kron(Accl, Accl);
    end
    for r = 1:pf
        Jcl = sys.Jf(:,:,r) + sys.BJf(:,:,r)*K;
        Lyap = Lyap + sys.lamf(r)*(kron(Jcl, Jcl) - eye(n*n));
    end
end
xocc = -Lyap \ cst.X0(:);  Xocc = reshape(xocc, n, n);  Xocc = (Xocc + Xocc.')/2;
QSR = [cst.Q, cst.S; cst.S.', cst.R];
Z   = [eye(n); K]*Xocc*[eye(n); K].';
cost.primal   = trace(cst.X0*P);
cost.dual     = trace(QSR*Z);
cost.slackres = norm(M*[eye(n); K], 'fro');    % C_Z(P^*)[I;K] = 0
end

% =========================================================================
function cst = extract_cost(cost, raw, n, m)
%EXTRACT_COST  Cost blocks from COST, falling back to the SYS fields, then to
%  Q = 0, R = I, S = 0, X0 = I.
cst.Q  = pick(cost,'Q',  raw,'Q',  zeros(n));
cst.R  = pick(cost,'R',  raw,'R',  eye(m));
cst.S  = pick(cost,'S',  raw,'S',  zeros(n,m));
cst.X0 = pick(cost,'X0', raw,'X0', eye(n));
cst.Q  = (cst.Q  + cst.Q.') /2;
cst.R  = (cst.R  + cst.R.') /2;
cst.X0 = (cst.X0 + cst.X0.')/2;
end

function v = pick(a, fa, b, fb, d)
if isstruct(a) && isfield(a,fa) && ~isempty(a.(fa)); v = a.(fa);
elseif isstruct(b) && isfield(b,fb) && ~isempty(b.(fb)); v = b.(fb);
else; v = d; end
end

function sys = fill_defaults(sys, n, m)
if ~isfield(sys,'Ac') ||isempty(sys.Ac);  sys.Ac  = zeros(n,n,0); end
if ~isfield(sys,'Bc') ||isempty(sys.Bc)||size(sys.Bc,3)~=size(sys.Ac,3)
    sys.Bc = zeros(n,m,size(sys.Ac,3)); end
if ~isfield(sys,'Jf') ||isempty(sys.Jf);  sys.Jf  = zeros(n,n,0); end
if ~isfield(sys,'BJf')||isempty(sys.BJf)||size(sys.BJf,3)~=size(sys.Jf,3)
    sys.BJf = zeros(n,m,size(sys.Jf,3)); end
if ~isfield(sys,'lamf')||isempty(sys.lamf); sys.lamf = ones(1,size(sys.Jf,3)); end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
