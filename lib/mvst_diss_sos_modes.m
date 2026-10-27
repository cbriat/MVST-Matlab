function res = mvst_diss_sos_modes(sys, supply, dwell, mode, opts)
%MVST_DISS_SOS_MODES  Generic SOSTOOLS dissipativity-ANALYSIS program on a
%  PRODUCT CONE (H^n)^M for mode-indexed linear Hermitian matrix-valued
%  systems: the coupled storage certificates {P_i} of the Dissipativity
%  chapter's Markov-jump and switched dwell-time results.  This is the
%  product-cone counterpart of MVST_DISS_SOS (single certificate) and the
%  dissipativity counterpart of MVST_STAB_SOS_MODES: here every mode carries
%  its own storage certificate P_i, coupled through the Markov generator
%  (CT-MJLS) or the all-pairs switch conditions (minimum dwell-time switched).
%  As in MVST_DISS_SOS the L2-gain level is found by bisection and dissipativity
%  is decided by common-slack maximization.
%
%  SYS fields:
%    .A0  n x n x M          per-mode flow drift
%    .Ac  n x n x q x M      per-mode flow state-multiplicative channels
%                            (optional; q may be 0)
%    .Dwc n x m x q          SHARED flow disturbance-multiplicative channels
%                            D_{w,l} (optional; paired with .Ac, padded with
%                            zeros to match; loads the flow pivot Theta as in
%                            the exact backends MVST_L2GAIN_CT_MJLS/_CT_DWELL)
%    .Bw  n x m              common disturbance input
%    .C   p x n x M          per-mode flow output (or p x n common)
%    .Dw  p x m              common output feedthrough
%    .Pig M x M              Markov GENERATOR (CT-MJLS), rows sum to 0
%    .J   n x n              switch reset (switched; default identity)
%
%  SUPPLY:  .type 'l2gain'  (.gamma), .nu/.rho  (see MVST_DISS_SOS).
%  DWELL:
%    .type 'ctmjls' | 'switched'
%       ('ctmjls'   = CT Markov jump diffusion, constant {P_i}, generator
%                     coupling sum_j Pig(i,j) P_j added to the flow C_M block;
%                     feasible iff the coupled bounded-real Riccati has a
%                     bounded solution (cross-check mvst_l2gain_ct_mjls).
%        'switched' = minimum-dwell-time SWITCHED L2-gain: mode-indexed timer
%                     certificates P_i(theta), theta = tau/Tbar in [0,1], with
%                     per-mode flow C_M >= 0 on [0,1], the frozen-tail flow
%                     C_M^0(P_i(1)) >= 0, and the all-pairs switch conditions
%                     J^*P_j(0)J - P_i(1) >= 0 (Theorem th:diss:MDT); the SOS
%                     value is a certified UPPER bound on the switched gain
%                     (cross-check mvst_l2gain_ct_dwell).)
%    .Tbar scalar (switched)          .deg  timer degree (switched; default 4).
%
%  MODE: 'feas' | 'gamma'.
%  OPTS (struct or MVST.OPTIONS, optional): .solver selects the SDP backend
%  ('sedumi' (default) | 'sdpt3' | 'mosek' | 'nosdp'), normalized by MVST_OPTS.
%  RES: .feas, .slack, .gamma, .Pi (cell of certificates), .info.

if nargin < 4 || isempty(mode); mode = 'feas'; end
if nargin < 5; opts = struct(); end
oo = mvst_opts(opts); solver = oo;   % normalized options travel to the solve (reltol/verbosity)
% SYS may be a THIN multi-mode system (MVST_SYSTEM): if it lacks .A0 it is
% routed through MVST_TORICH, which stacks the modes and exposes the
% generator under .Pig.  (A thin system carries no per-mode output data, so
% .C/.Dw default to zero unless attached to the struct before the call.)
if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
n = size(sys.A0,1);  M = size(sys.A0,3);
m = size(sys.Bw,2);
sys = fill_defaults(sys, n, m, M);

switch mode
    case 'gamma'
        if ~strcmp(supply.type,'l2gain')
            error('mvst:diss_sos_modes:gammamode','gamma mode needs l2gain.');
        end
        lo = norm(sys.Dw) + 1e-4; hi = max(2*lo,1); Rbest = []; grow = 0;
        while true
            r = solve_once(sys, supply, dwell, n, m, M, hi, true, solver);
            if r.feas; Rbest = r; break; end
            hi = 2*hi; grow = grow+1; if grow>30; break; end
        end
        for it = 1:24
            mid = 0.5*(lo+hi);
            r = solve_once(sys, supply, dwell, n, m, M, mid, true, solver);
            if r.feas; hi = mid; Rbest = r; else; lo = mid; end
        end
        res = Rbest; if isempty(res); res = solve_once(sys,supply,dwell,n,m,M,hi,true,solver); end
        res.gamma = hi;
    case 'feas'
        g = getf(supply,'gamma',1);
        res = solve_once(sys, supply, dwell, n, m, M, g, true, solver);
    otherwise
        error('mvst:diss_sos_modes:mode','unknown mode "%s"', mode);
end
end

% =========================================================================
function res = solve_once(sys, supply, dwell, n, m, M, gam, useSlackObj, solver)
tol = 1e-7; Rcap = 1e4;
switch dwell.type
    case 'ctmjls';   res = solve_ctmjls(sys, supply, n, m, M, gam, useSlackObj, tol, Rcap, solver);
    case 'switched'; res = solve_switched(sys, supply, dwell, n, m, M, gam, useSlackObj, tol, Rcap, solver);
    otherwise; error('mvst:diss_sos_modes:type','unknown dwell "%s"', dwell.type);
end
end

% =========================================================================
function res = solve_ctmjls(sys, supply, n, m, M, gam, useSlackObj, tol, Rcap, solver)
% Constant coupled storages {P_i}: per-mode flow C_M(P_i) + generator
% coupling sum_j Pig(i,j) P_j >= slack I, and 0 <= -P_i <= Rcap I.
pvar th; sl = dpvar('sl');
prog = sosprogram(th); prog = sosdecvar(prog, sl);
P = cell(1,M);
for i=1:M; [prog,P{i}] = sospolymatrixvar(prog, monomials(th,0), [n n], 'symmetric'); end
if useSlackObj; slk = sl; else; slk = 1e-7; end

for i = 1:M
    Mb = supply_blocks(supply, sys, gam, m, i);
    C = flow_cM_mode(sys, P, i, M, Mb, m);            % augmented (n+m)
    prog = sosmatrixineq(prog, C - slk*eye(n+m));
    prog = sosmatrixineq(prog, -P{i});
    prog = sosmatrixineq(prog, Rcap*eye(n) + P{i});
end
if useSlackObj; prog = sossetobj(prog, -sl); end
o = mvst_sosopts(solver); prog = mvst_sossolve(prog, o);
sval = double(sosgetsol(prog, sl));
Pv = cell(1,M); for i=1:M; Pv{i} = double(sosgetsol(prog,P{i})); end
res.type='ctmjls'; res.slack=sval; res.feas=sval>tol; res.gamma=gam;
res.Pi=Pv; res.info=struct('P',{Pv});
end

% -------------------------------------------------------------------------
function C = flow_cM_mode(sys, P, i, M, Mb, m)
n = size(sys.A0,1); q = size(sys.Ac,3);
A = sys.A0(:,:,i);
G11 = A'*P{i} + P{i}*A + Mb.M11;
for j = 1:M; G11 = G11 + sys.Pig(i,j)*P{j}; end       % generator coupling
G12 = P{i}*sys.Bw + Mb.M12;
G22 = Mb.M22;
for l = 1:q
    Cl = sys.Ac(:,:,l,i);  Dwcl = sys.Dwc(:,:,l);
    G11 = G11 + Cl'*P{i}*Cl;
    G12 = G12 + Cl'*P{i}*Dwcl;             % state x disturbance channel cross-term
    G22 = G22 + Dwcl'*P{i}*Dwcl;           % disturbance channel pivot loading (Theta)
end
C = [ G11, G12 ; G12', G22 ]; C = (C+C')/2;
end

% =========================================================================
function res = solve_switched(sys, supply, dwell, n, m, M, gam, useSlackObj, tol, Rcap, solver)
% Minimum-dwell-time switched L2-gain: mode-indexed timer certificates.
T = dwell.Tbar; deg = getf(dwell,'deg',4);
pvar th; sl = dpvar('sl'); gt = th*(1-th);
prog = sosprogram(th); prog = sosdecvar(prog, sl);
P = cell(1,M);
for i=1:M; [prog,P{i}] = sospolymatrixvar(prog, monomials(th,0:deg), [n n], 'symmetric'); end
if useSlackObj; slk = sl; else; slk = 1e-7; end

for i = 1:M
    Mb = supply_blocks(supply, sys, gam, m, i);
    Pdot = diff(P{i},th)/T;
    Cf = flow_cM_single(sys, P{i}, Pdot, i, Mb, m);      % flow on [0,1]
    [prog,S1] = sospolymatrixvar(prog, monomials(th,0:2*deg), [n+m n+m], 'symmetric');
    prog = sosmatrixineq(prog, Cf - slk*eye(n+m) - gt*S1); prog = sosmatrixineq(prog, S1);
    Cf0 = flow_cM_single(sys, subs(P{i},th,1), zeros(n), i, Mb, m);  % frozen tail
    prog = sosmatrixineq(prog, Cf0 - slk*eye(n+m));
    % storage nonnegativity 0 <= -P_i <= Rcap I on [0,1]
    [prog,Xlo]=sospolymatrixvar(prog,monomials(th,0:2*deg),[n n],'symmetric');
    prog = sosmatrixineq(prog, -P{i} - gt*Xlo); prog = sosmatrixineq(prog, Xlo);
    [prog,Xhi]=sospolymatrixvar(prog,monomials(th,0:2*deg),[n n],'symmetric');
    prog = sosmatrixineq(prog, Rcap*eye(n)+P{i} - gt*Xhi); prog = sosmatrixineq(prog, Xhi);
end
% switch conditions, all ordered pairs, reset J, no jump supply:
% J^*P_j(0)J - P_i(1) >= slack I
for i = 1:M
    Pi1 = subs(P{i},th,1);
    for j = 1:M
        if j ~= i
            Pj0 = subs(P{j},th,0);
            prog = sosmatrixineq(prog, sys.J'*Pj0*sys.J - Pi1 - slk*eye(n));
        end
    end
end
if useSlackObj; prog = sossetobj(prog, -sl); end
o = mvst_sosopts(solver); prog = mvst_sossolve(prog, o);
sval = double(sosgetsol(prog, sl));
Pv = cell(1,M); for i=1:M; Pv{i}=sosgetsol(prog,P{i}); end
res.type='switched'; res.slack=sval; res.feas=sval>tol; res.gamma=gam;
res.Pi=Pv; res.info=struct('T',T,'deg',deg);
end

% -------------------------------------------------------------------------
function C = flow_cM_single(sys, Pi, Pdot, i, Mb, m)
% flow C_M for a single mode i (no generator coupling): augmented (n+m).
n = size(sys.A0,1); q = size(sys.Ac,3); A = sys.A0(:,:,i);
G11 = Pdot + A'*Pi + Pi*A + Mb.M11;
G12 = Pi*sys.Bw + Mb.M12;
G22 = Mb.M22;
for l=1:q
    Cl=sys.Ac(:,:,l,i); Dwcl=sys.Dwc(:,:,l);
    G11 = G11 + Cl'*Pi*Cl;
    G12 = G12 + Cl'*Pi*Dwcl;               % state x disturbance channel cross-term
    G22 = G22 + Dwcl'*Pi*Dwcl;             % disturbance channel pivot loading (Theta)
end
C = [ G11, G12 ; G12', G22 ]; C = (C+C')/2;
end

% =========================================================================
function Mb = supply_blocks(supply, sys, gam, m, i)
n = size(sys.A0,1);
C = sys.C(:,:,i); Dw = sys.Dw; p = size(C,1);
nu = getf(supply,'nu',0); rho = getf(supply,'rho',0);
switch supply.type
    case 'l2gain'
        Mb.M11 = -(C'*C); Mb.M12 = -(C'*Dw); Mb.M22 = gam^2*eye(m) - Dw'*Dw;
    case 'passivity'
        Mb.M11 = -rho*(C'*C); Mb.M12 = C' - rho*(C'*Dw);
        Mb.M22 = (Dw+Dw') - nu*eye(m) - rho*(Dw'*Dw);
    case 'qsr'
        [Qs,Ss,Rs] = mvst_qsr_supply(supply, p, m, 'flow');
        Mb.M11 = C'*Qs*C; Mb.M12 = C'*(Qs*Dw + Ss);
        Mb.M22 = Dw'*Qs*Dw + Dw'*Ss + Ss'*Dw + Rs;
    otherwise; error('mvst:diss_sos_modes:supply','unknown supply.');
end
end

function sys = fill_defaults(sys, n, m, M)
if ~isfield(sys,'Ac')||isempty(sys.Ac); sys.Ac = zeros(n,n,0,M); end
if ndims(sys.Ac) < 4 && size(sys.Ac,3) >= 1; sys.Ac = reshape(sys.Ac,n,n,size(sys.Ac,3),M); end
if ~isfield(sys,'Bw')||isempty(sys.Bw); sys.Bw = zeros(n,m); end
if ~isfield(sys,'Dw')||isempty(sys.Dw); sys.Dw = zeros(size(sys.C,1),m); end
if size(sys.C,3)==1 && M>1; sys.C = repmat(sys.C,1,1,M); end
if ~isfield(sys,'Pig')||isempty(sys.Pig); sys.Pig = zeros(M); end
if ~isfield(sys,'J')||isempty(sys.J); sys.J = eye(n); end
% Disturbance-multiplicative (Brownian) channels D_{w,l}: SHARED across modes,
% n x m x q, paired with the per-mode state channels sys.Ac (n x n x q x M).
% They load the flow pivot exactly as in the exact backends' Theta.
if ~isfield(sys,'Dwc')||isempty(sys.Dwc); sys.Dwc = zeros(n,m,0); end
if iscell(sys.Dwc); sys.Dwc = cat(3, sys.Dwc{:}); end
qAc = size(sys.Ac,3); qD = size(sys.Dwc,3); qmax = max(qAc,qD);
if qAc < qmax; sys.Ac(:,:,qAc+1:qmax,:) = 0; end          % pad state channels
if qD  < qmax; sys.Dwc(:,:,qD+1:qmax)   = 0; end          % pad disturbance channels
end

function v = getf(s, f, d); if isfield(s,f)&&~isempty(s.(f)); v=s.(f); else; v=d; end; end
