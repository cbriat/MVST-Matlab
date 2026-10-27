function res = mvst_oc_sos_modes(sys, cost, dwell, opts)
%MVST_OC_SOS_MODES  Optimal control by the dual (storage-function) SOS/LMI
%  program for a COUPLED Markov-jump system, the multi-mode counterpart of
%  MVST_OC_SOS.  One SOS program over the mode certificates {P_i} solves all
%  modes at once: the coupled dual LMI of the Markov-jump case (I-GRE example
%  eq:ex:mjls:riccati, infinite-horizon dual LMI th:IMP:LQ:inf:LMI /
%  eq:IMP:LQ:LMI:inf) couples the flow operators across modes through the
%  generator term  sum_j pi_ij P_j  added to the (1,1) block of each C_Z:
%
%     max_{P_1..P_N}  sum_i <P_i,X0_i>
%     s.t.  C_Z^i(P_i) + [sum_j pi_ij P_j , 0 ; 0 , 0] >= 0,   i = 1..N,
%
%  where C_Z^i is the per-mode augmented KYP / Schur matrix of MVST_OC_SOS
%  (Definition def:IMP:GRE, Schur decomposition pf:suff:schur:flow).  SOSTOOLS
%  poses each C_Z^i(P_i) >= 0 as [x;u]^* C_Z^i [x;u] in Sigma over z = [x;u]
%  (pvar mode); the optimum recovers the coupled stabilizing I-GRE solution and
%  the per-mode gains K_i = -(C_Z^i)_3^{-1}(C_Z^i)_2^* to solver accuracy.
%  Solver: SOSTOOLS + SeDuMi.
%
%  SYS  a multi-mode system (MVST_SYSTEM):
%    .modes  1 x N cell of thin modes (MVST_MODE): drift .A, input .B, noise
%            .Ac/.Bc, and cost blocks .Q/.R/.S/.X0
%    .Pi     N x N continuous-time Markov generator (rows sum to 0)
%    .pi0    1 x N initial mode distribution (default uniform); the cost is
%            weighted  sum_i pi0_i <P_i,X0_i>
%
%  COST  (optional) 1 x N cell of cost structs (.Q/.R/.S/.X0) overriding the
%    per-mode cost blocks; omit or [] to use the modes' own cost data.
%
%  DWELL  (struct; .type default 'arbitrary'):  'arbitrary' = continuous-time,
%    'discrete' = discrete-time Markov jump (per-mode transition A0, noise Ac).
%
%  OPTS  (struct or MVST.OPTIONS, optional; normalized by MVST_OPTS):
%    .solver   SDP backend 'sedumi' (default) | 'sdpt3' | 'mosek',  .fid
%              verbosity (default 0)
%    .riccati  also return the Gauss-Seidel coupled-I-GRE benchmark solution
%              (default false): .Pric, .Kric filled when true.
%
%  RES:
%    .P     1 x N cell of certificates {P_i}  (coupled I-GRE solution)
%    .K     1 x N cell of optimal gains {K_i}
%    .cost  struct  .primal  = sum_i pi0_i <P_i,X0_i>   (the optimal value)
%                   .dual    = sum_i trace([Q_i S_i;S_i' R_i] Z_i^*), the
%                              mode-coupled moment (complementary-slackness) cost
%                   .riccati = sum_i pi0_i <P_i,X0_i> read on the SOS solution
%    .feas  clean-solve flag ;  .hist  SeDuMi solver-info struct
%    .Pric,.Kric  1 x N cells: Gauss-Seidel coupled-I-GRE benchmark (OPTS.riccati)
%
%  See also MVST_OC_SOS, MVST_CARE_CT, MVST_STAB_SOS_MODES, MVST_DISS_SOS_MODES.

if nargin < 2; cost  = []; end
if nargin < 3 || isempty(dwell); dwell = struct('type','arbitrary'); end
if nargin < 4; opts  = struct(); end
if ~isfield(dwell,'type') || isempty(dwell.type); dwell.type = 'arbitrary'; end
opts = mvst_opts(opts);   % normalize + validate .solver ('sedumi' default)
if ~isfield(sys,'modes')
    error('mvst:oc_sos_modes:input','SYS must be a multi-mode MVST_SYSTEM struct.');
end

modes = sys.modes;  N = numel(modes);
Pi = getf(sys,'Pi', zeros(N));
pi0 = getf(sys,'pi0', ones(1,N)/N);
n = modes{1}.n;  m = modes{1}.m;

% Per-mode rich dynamics + cost blocks.
dyn = cell(1,N);  cst = cell(1,N);
for i = 1:N
    dyn{i} = to_dyn(modes{i}, n, m);
    cst{i} = mode_cost(cost, modes{i}, i, n, m);
end

% -------- one SOS program over {P_i} ---------------------------------------
V = mpvar('z', [n+m, 1]);
prog = sosprogram(V);
P = cell(1,N);
for i = 1:N
    [prog, P{i}] = sospolymatrixvar(prog, monomials(V,0), [n n], 'symmetric');
end
obj = 0;
for i = 1:N
    coup = zeros(n);  for j = 1:N; coup = coup + Pi(i,j)*P{j}; end
    M = oc_kyp(P{i}, dyn{i}, cst{i}, coup, dwell.type);
    prog = sosineq(prog, V.'*M*V);
    obj  = obj - trace(cst{i}.X0*P{i});          % maximize sum_i <P_i,X0_i>
end
prog = sossetobj(prog, obj);
o = mvst_sosopts(opts);
o.params.fid = getf(opts,'fid',0);
[prog, info] = mvst_sossolve(prog, o);

Pc = cell(1,N);  Kc = cell(1,N);
for i = 1:N; Pc{i} = double(sosgetsol(prog,P{i})); Pc{i} = (Pc{i}+Pc{i}.')/2; end
primal = 0;  ric = 0;
for i = 1:N
    coup = zeros(n);  for j = 1:N; coup = coup + Pi(i,j)*Pc{j}; end
    [~, Kc{i}] = oc_kyp(Pc{i}, dyn{i}, cst{i}, coup, dwell.type);
    primal = primal + pi0(i)*trace(cst{i}.X0*Pc{i});
    ric    = ric    + pi0(i)*trace(cst{i}.X0*Pc{i});
end
% Coupled dual: mode-coupled occupation second moments on u = K_i x.
Xocc = coupled_occupation(dyn, cst, Kc, Pi, pi0, dwell.type);
dual = 0;
for i = 1:N
    QSR  = [cst{i}.Q, cst{i}.S; cst{i}.S.', cst{i}.R];
    dual = dual + trace(QSR*[eye(n);Kc{i}]*Xocc{i}*[eye(n);Kc{i}].');
end

res.P    = Pc;
res.K    = Kc;
res.cost = struct('primal',primal,'dual',dual,'riccati',ric);
res.feas = (info.pinf == 0 && info.dinf == 0 && info.numerr < 2);
res.hist = info;

if getf(opts,'riccati',false)
    [res.Pric, res.Kric] = coupled_riccati(modes, dyn, cst, Pi, dwell.type, n);
end
end

% =========================================================================
function [M, K, Theta] = oc_kyp(P, dyn, cst, coup, type)
%OC_KYP  Per-mode augmented KYP / Schur matrix C_Z(P) with Markov coupling COUP
%  added to the (1,1) block (see MVST_OC_SOS).
q = size(dyn.Ac,3);  pf = size(dyn.Jf,3);
if strcmp(type,'discrete')
    Theta = cst.R;                 S = dyn.B0.'*P*dyn.A0 + cst.S.';
    acc   = dyn.A0.'*P*dyn.A0 - P + cst.Q + coup;
    for l = 1:q
        Ac = dyn.Ac(:,:,l); Bc = dyn.Bc(:,:,l);
        Theta = Theta + Bc.'*P*Bc;  S = S + Bc.'*P*dyn.A0;  acc = acc + Ac.'*P*Ac;
    end
else
    Theta = cst.R;                 S = dyn.B0.'*P + cst.S.';
    acc   = dyn.A0.'*P + P*dyn.A0 + cst.Q + coup;
    for l = 1:q
        Ac = dyn.Ac(:,:,l); Bc = dyn.Bc(:,:,l);
        Theta = Theta + Bc.'*P*Bc;  S = S + Bc.'*P*Ac;      acc = acc + Ac.'*P*Ac;
    end
    for r = 1:pf
        Jf = dyn.Jf(:,:,r); BJf = dyn.BJf(:,:,r); lr = dyn.lamf(r);
        acc   = acc   + lr*(Jf.'*P*Jf - P);
        S     = S     + lr*(BJf.'*P*Jf);
        Theta = Theta + lr*(BJf.'*P*BJf);
    end
end
Theta = (Theta + Theta.')/2;
M = [acc, S.'; S, Theta];  M = (M + M.')/2;
K = [];
if nargout >= 2; K = -(Theta \ S); end
end

% =========================================================================
function Xocc = coupled_occupation(dyn, cst, K, Pi, pi0, type)
%COUPLED_OCCUPATION  Mode-coupled occupation second moments {X_i} solving
%  (CT)   Acl_i X_i + X_i Acl_i' + sum_l Acl_il X_i Acl_il'
%             + sum_j Pi(j,i) X_j + pi0_i X0_i = 0,
%  (DT)   Acl_i X_i Acl_i' + sum_l Acl_il X_i Acl_il' - X_i
%             + sum_j Pi(j,i) X_j + pi0_i X0_i = 0.
N = numel(dyn); n = size(dyn{1}.A0,1); nn = n*n;
L = zeros(N*nn); rhs = zeros(N*nn,1);
for i = 1:N
    Acl = dyn{i}.A0 + dyn{i}.B0*K{i};
    if strcmp(type,'discrete')
        Di = kron(Acl, Acl) - eye(nn);
        for l = 1:size(dyn{i}.Ac,3)
            Accl = dyn{i}.Ac(:,:,l) + dyn{i}.Bc(:,:,l)*K{i}; Di = Di + kron(Accl,Accl);
        end
    else
        Di = kron(eye(n),Acl) + kron(Acl,eye(n));
        for l = 1:size(dyn{i}.Ac,3)
            Accl = dyn{i}.Ac(:,:,l) + dyn{i}.Bc(:,:,l)*K{i}; Di = Di + kron(Accl,Accl);
        end
    end
    Di = Di + Pi(i,i)*eye(nn);                          % j = i coupling
    L((i-1)*nn+(1:nn),(i-1)*nn+(1:nn)) = Di;
    for j = 1:N
        if j ~= i
            L((i-1)*nn+(1:nn),(j-1)*nn+(1:nn)) = Pi(j,i)*eye(nn);
        end
    end
    rhs((i-1)*nn+(1:nn)) = -pi0(i)*cst{i}.X0(:);
end
x = L\rhs;  Xocc = cell(1,N);
for i = 1:N; Xi = reshape(x((i-1)*nn+(1:nn)),n,n); Xocc{i} = (Xi+Xi.')/2; end
end

% =========================================================================
function [Pc, Kc] = coupled_riccati(modes, dyn, cst, Pi, type, n)
%COUPLED_RICCATI  Gauss-Seidel benchmark for the coupled Markov-jump I-GRE.
%  Fold pi_ii P_i into A_i -> A_i + (pi_ii/2)I and the off-diagonal coupling
%  sum_{j~=i} pi_ij P_j into Q_i; each mode is then a standard single-mode ARE
%  (MVST_CARE_CT).  Continuous-time only (the benchmark for OPTS.riccati).
if strcmp(type,'discrete')
    error('mvst:oc_sos_modes:riccati','coupled Riccati benchmark is continuous-time.');
end
N = numel(modes);  Pc = cell(1,N);  Kc = cell(1,N);
for i = 1:N; Pc{i} = zeros(n); end
for it = 1:400
    Pold = Pc;
    for i = 1:N
        coupoff = zeros(n);
        for j = 1:N; if j ~= i; coupoff = coupoff + Pi(i,j)*Pc{j}; end; end
        mdi = modes{i};
        mdi.A = mdi.A + 0.5*Pi(i,i)*eye(n);
        mdi.Q = (mdi.Q + coupoff + (mdi.Q+coupoff).')/2;
        [Pc{i}, Kc{i}] = mvst_care_ct(mdi);
    end
    d = 0;  for i = 1:N; d = max(d, norm(Pc{i}-Pold{i},'fro')); end
    if d < 1e-11; break; end
end
% Report the gains with the true (undamped) mode data and full coupling.
for i = 1:N
    coup = zeros(n);  for j = 1:N; coup = coup + Pi(i,j)*Pc{j}; end
    [~, Kc{i}] = oc_kyp(Pc{i}, dyn{i}, cst{i}, coup, type);
end
end

% =========================================================================
function dyn = to_dyn(md, n, m)
%TO_DYN  Thin mode -> rich per-mode dynamics struct for OC_KYP.
dyn.A0 = md.A;  dyn.B0 = md.B;
if isfield(md,'Ac') && ~isempty(md.Ac); dyn.Ac = md.Ac; else; dyn.Ac = zeros(n,n,0); end
if isfield(md,'Bc') && ~isempty(md.Bc) && size(md.Bc,3)==size(dyn.Ac,3)
    dyn.Bc = md.Bc; else; dyn.Bc = zeros(n,m,size(dyn.Ac,3)); end
dyn.Jf = zeros(n,n,0);  dyn.BJf = zeros(n,m,0);  dyn.lamf = zeros(1,0);
end

function cst = mode_cost(cost, md, i, n, m)
%MODE_COST  Per-mode cost blocks from the COST cell override, else the mode.
if iscell(cost) && numel(cost) >= i && ~isempty(cost{i}); c = cost{i}; else; c = []; end
cst.Q  = pick(c,'Q',  md,'Q',  zeros(n));
cst.R  = pick(c,'R',  md,'R',  eye(m));
cst.S  = pick(c,'S',  md,'S',  zeros(n,m));
cst.X0 = pick(c,'X0', md,'X0', eye(n));
cst.Q  = (cst.Q  + cst.Q.') /2;
cst.R  = (cst.R  + cst.R.') /2;
cst.X0 = (cst.X0 + cst.X0.')/2;
end

function v = pick(a, fa, b, fb, d)
if isstruct(a) && isfield(a,fa) && ~isempty(a.(fa)); v = a.(fa);
elseif isstruct(b) && isfield(b,fb) && ~isempty(b.(fb)); v = b.(fb);
else; v = d; end
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
