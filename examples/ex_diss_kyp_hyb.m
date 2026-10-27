function ex_diss_kyp_hyb()
%EX_DISS_KYP_HYB  Non-exhaustion of the hybrid-harmonic sweep (Chapter
%   "Dissipativity", Example ex:diss:kyp:hyb:cex, Proposition
%   prop:diss:kyp:hyb:nonexh).
%
%   Scalar impulsive system (eq:diss:kyp:hyb:det), n=m=1, no jump disturbance:
%       dx/dt = a x + b w on [0,Tp),   x+ = J x at the period boundary,
%   a=-0.3, b=1, J=0.5, Tp=1, with the Hermitian supply matrices M (flow) and
%   Mk (jump) on (x,w) printed in the example.
%
%   Computes, with no SDP solver (closed-form profiles + adaptive quadrature):
%     (1) the monodromy Jt = J e^{a Tp} (Schur check);
%     (2) the rank-one hybrid-harmonic sweep Sigma(Y_{omega,theta}) of
%         prop:diss:kyp:hyb:realize with |w|=1: minimum over a coarse and a
%         refined (omega,theta) grid, then local refinement (fminsearch);
%     (3) a sweep-positive, section-negative witness built as the comb
%         superposition of lem:diss:kyp:hyb:comb: the 2K+1 single-line points
%         on the theta-comb omega_r=(theta+2 pi r)/Tp, r=-K..K, all with the
%         same discrete phase theta (so their sum stays in K_hyb).  For
%         amplitudes c_r, Sigma(Y) = c' (Gf + Gj) c with Gram matrices
%             Gf(r,r') = int_0^Tp zeta_r' M zeta_r' dtau,
%             Gj(r,r') = zeta_r(Tp)' Mk zeta_r'(Tp),
%         c the unit eigenvector of the least eigenvalue, theta refined by
%         fminbnd.  Reported: flow and jump supplies, the sum of single-line
%         values, the jump cross-term defect Delta_jp (off-diagonal part of
%         Gj) and the flow cross-term (off-diagonal part of Gf), and the
%         minimum of the disturbance profile |w(tau)|^2 on [0,Tp];
%     (4) fewer-line superpositions: for every subset of 2 or 3 lines from
%         r=-4..4 and of 4 lines from r=-2..2, the minimum over theta (1-degree
%         grid + fminbnd) of the least eigenvalue of the sub-Gram.  Result
%         (2026-10-04): 2 lines 0.177 > 0; 3 lines -0.029 < 0 (r=-4,-2,-1);
%         4 lines -0.199 < 0.
%     (5) numerical check of the restated Lemma lem:diss:kyp:hyb:comb
%         (2026-10-04): membership of the superposition in K_hyb, the
%         decomposition Sigma = sum single + Delta_jp + Delta_fl, the closed
%         form of the transient amplitudes delta_r, Delta_fl from the
%         transients, and the pure-comb identity (c) on an m=2 variant with
%         a jump disturbance channel F, the jump disturbance being tied to
%         the arc (vhat = e^{j om Tp} what) and the second disturbance
%         component tuned so that delta_r=0.  Prints residuals; asserts they
%         are below 1e-8 (ODE residual below 1e-7).
%
%   Macros: figures/data/diss_kyp_hyb.tex.  No figure.
%
%   See also EX_KYP, MVST_KYP.

here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)

% ---- data (Example ex:diss:kyp:hyb:cex) ------------------------------------
a = -0.3; b = 1; J = 0.5; Tp = 1;
M  = [0.448, 0.841-0.374i; 0.841+0.374i, 0.821];
Mk = [-0.451, 0.494-0.501i; 0.494+0.501i, 0.025];
assert(norm(M-M')==0 && norm(Mk-Mk')==0, 'supply matrices must be Hermitian');

% ---- (1) monodromy ----------------------------------------------------------
Jt = J*exp(a*Tp);
fprintf('Hybrid KYP cex: monodromy Jt = J e^{a Tp} = %.6f (Schur: %d)\n', Jt, abs(Jt)<1);
assert(abs(Jt) < 1, 'monodromy must be Schur');

% ---- (2) rank-one sweep ------------------------------------------------------
% Single-line profile with |w|=1: zeta(tau) = u(tau) xi0 + v(tau),
%   u = [e^{a tau}; 0],  v = [f(tau); e^{j om tau}],
%   f(tau) = b (e^{j om tau} - e^{a tau})/(j om - a),
%   xi0 = (e^{j theta} - Jt)^{-1} J f(Tp).
% Sigma is quadratic in xi0: Sigma = q0 + 2 Re(conj(xi0) q1) + |xi0|^2 q2.
sweep = @(om,th) sweepval(om,th,a,b,J,Tp,M,Mk);
thg = linspace(0, 2*pi, 361); thg(end) = [];
grids = {linspace(-30,30,601), linspace(-30,30,2401)};
gmin = zeros(1,2); gargs = zeros(2,2);
for gi = 1:2
    best = inf;
    for om = grids{gi}
        [q0,q1,q2] = sweepcoef(om,a,b,J,Tp,M,Mk);
        xi0 = J*ffun(Tp,om,a,b) ./ (exp(1i*thg) - Jt);
        s = q0 + 2*real(conj(xi0)*q1) + abs(xi0).^2*q2;
        [v,k] = min(s);
        if v < best; best = v; gargs(gi,:) = [om, thg(k)]; end
    end
    gmin(gi) = best;
end
opt = optimset('TolX',1e-10,'TolFun',1e-12,'MaxFunEvals',4000,'MaxIter',4000);
[p, smin] = fminsearch(@(p) sweep(p(1),p(2)), gargs(2,:), opt);
omMin = p(1); thMin = mod(p(2), 2*pi);
fprintf('  sweep min: coarse grid %.6f, fine grid %.6f, refined %.6f at (omega,theta)=(%.4f,%.4f)\n', ...
    gmin(1), gmin(2), smin, omMin, thMin);
% high-frequency limit: xi -> 0, Sigma -> Tp M(2,2) + Mk(2,2)
fprintf('  sweep high-frequency limit Tp*M22+Mk22 = %.6f\n', Tp*real(M(2,2))+real(Mk(2,2)));
assert(smin > 0, 'sweep must be strictly positive');

% ---- (3) comb-superposition witness -----------------------------------------
K = 2; rr = -K:K;
lam = @(th) combmin(th,rr,a,b,J,Tp,M,Mk);
thc = linspace(0, 2*pi, 181); thc(end) = [];
lv = arrayfun(lam, thc); [~,k] = min(lv);
dth = thc(2)-thc(1);
thW = fminbnd(lam, thc(k)-dth, thc(k)+dth, optimset('TolX',1e-10));
[lmin, c, Gf, Gj] = combmin(thW,rr,a,b,J,Tp,M,Mk);
flowW = real(c'*Gf*c);  jumpW = real(c'*Gj*c);
single = real(sum(abs(c).^2 .* (diag(Gf)+diag(Gj))));
Djp   = real(c'*(Gj - diag(diag(Gj)))*c);       % jump cross-term defect
Dfl   = real(c'*(Gf - diag(diag(Gf)))*c);       % flow cross-term
% disturbance profile w(tau) = sum_r c_r e^{j omega_r tau}
tt = linspace(0,Tp,20001); omr = (thW + 2*pi*rr)/Tp;
wt = sum(c(:).*exp(1i*omr(:)*tt), 1);
wmin = min(abs(wt).^2);
fprintf(['  comb witness (r=%d..%d, theta=%.4f): Sigma=%.6f = flow %.6f + jump %.6f;\n', ...
    '    sum of single-line values %.6f, jump cross-term defect %.6f, flow cross-term %.6f;\n', ...
    '    min |w(tau)|^2 = %.4e\n'], rr(1), rr(end), thW, lmin, flowW, jumpW, single, Djp, Dfl, wmin);
assert(abs(flowW + jumpW - lmin) < 1e-9, 'witness split inconsistent');
assert(abs(single + Djp + Dfl - lmin) < 1e-9, 'defect decomposition inconsistent');
assert(lmin < 0, 'comb witness must be negative (section-negative)');

% normalize the global phase so that the r=0 amplitude is real positive
k0 = find(rr==0); c = c*exp(-1i*angle(c(k0)));
fprintf('  comb frequencies omega_r = (theta+2 pi r)/Tp: %s\n', mat2str(omr,4));
fprintf('  amplitudes c_r (r=%d..%d, unit norm, c_0>0):\n', rr(1), rr(end));
for i = 1:numel(rr); fprintf('    r=%+d: %+.4f %+.4fj\n', rr(i), real(c(i)), imag(c(i))); end

% ---- (4) fewer-line superpositions on the same comb ----------------------------
% For every subset S of comb lines (sizes 2,3 from the pool r=-4..4, size 4
% from the witness lines r=-2..2), minimize over theta in [0,2 pi) the least
% eigenvalue of the sub-Gram (Gf+Gj)(S,S): this is the minimum of Sigma over
% all unit-amplitude superpositions of exactly the lines S on one theta-comb.
pool = -4:4; nth = 360; thg2 = (0:nth-1)*2*pi/nth;
G = cell(1,nth);
for i = 1:nth; [~,~,Gf_,Gj_] = combmin(thg2(i),pool,a,b,J,Tp,M,Mk); G{i} = Gf_+Gj_; end
sizes = [2 3 4]; subMin = zeros(1,3); subArg = cell(1,3); subTh = zeros(1,3);
for si = 1:3
    if sizes(si) < 4; cand = pool; else; cand = rr; end
    S = nchoosek(cand, sizes(si)); best = inf;
    for q = 1:size(S,1)
        idx = arrayfun(@(r) find(pool==r), S(q,:));
        lv = cellfun(@(Gm) min(real(eig(Gm(idx,idx)))), G);
        [v,k] = min(lv);
        f = @(t) subeig(t,S(q,:),a,b,J,Tp,M,Mk);
        dt = thg2(2)-thg2(1);
        [tt_,v2] = fminbnd(f, thg2(k)-dt, thg2(k)+dt, optimset('TolX',1e-10));
        v = min(v,v2);
        if v < best; best = v; subArg{si} = S(q,:); subTh(si) = mod(tt_,2*pi); end
    end
    subMin(si) = best;
    fprintf('  min over all %d-line superpositions (lines from r=%d..%d, theta in [0,2pi)): %.6f, lines %s, theta=%.4f\n', ...
        sizes(si), cand(1), cand(end), best, mat2str(subArg{si}), subTh(si));
end

% ---- (5) verification of Lemma lem:diss:kyp:hyb:comb (restated 2026-10-04) ----
% (a) the comb superposition of hybrid single-line profiles lies in K_hyb and
%     Sigma(Y) = sum_r Sigma(Y_r) + Delta_jp + Delta_fl;
% (b) xi_r(tau) = e^{j om_r tau} xihat_r + e^{a tau} delta_r with
%     delta_r = e^{j th} (e^{j th} - Jt)^{-1} ((J-1) xihat_r + F what_r), and
%     Delta_fl is carried by the transients delta_r;
% (c) pure comb (delta_r = 0, needs a jump disturbance F what_r here since
%     J ~= 1): Delta_fl = 0 and Sigma(Y) = Tp sum_r <M,zh_r zh_r^*> + <Mk,zh_S zh_S^*>.
% The jump disturbance is TIED (2026-10-04 v-hat tie): the jump acts on the
% arc's own disturbance block at Tp^-, i.e. vhat = e^{j om Tp} what = e^{j th} what
% on the comb.  Since (1-J) xihat_r = F what_r with a common F cannot hold on
% several lines when m=1, check (c) uses an m=2 variant B=[b 0], F=[0 1].
verifyComb(thW, rr, c, a, b, J, Tp, M, Mk, single, Djp, Dfl, lmin);

% ---- macros -------------------------------------------------------------------
w = mvst_texdata(fullfile(mvst_figdir(),'data','diss_kyp_hyb.tex'));
w.num('DissKypHybJt', Jt, 3);
w.num('DissKypHybSweepMin', smin, 3);
w.num('DissKypHybSweepOm', omMin, 1);
w.num('DissKypHybSweepTh', thMin, 2);
w.num('DissKypHybCombTh', thW, 2);
w.num('DissKypHybWitTot', lmin, 3);
w.num('DissKypHybWitFlow', flowW, 3);
w.num('DissKypHybWitJump', jumpW, 3);
w.num('DissKypHybWitSingle', single, 3);
w.num('DissKypHybWitDjp', Djp, 3);
w.num('DissKypHybWitDfl', Dfl, 3);
w.num('DissKypHybWitWmin', wmin, 3);
w.num('DissKypHybCombThFull', thW, 4);
w.raw('DissKypHybWitAmp', cvec(c,3));
w.num('DissKypHybTwoLineMin', subMin(1), 3);
w.num('DissKypHybThreeLineMin', subMin(2), 3);
w.num('DissKypHybFourLineMin', subMin(3), 3);
w.close();
end

% =============================================================================
function s = cvec(c,p)
% column vector of complex entries a+bj as a LaTeX bmatrix
f = ['%.' num2str(p) 'f'];
e = cell(1,numel(c));
for i = 1:numel(c)
    re = real(c(i)); im = imag(c(i));
    if abs(im) < 0.5*10^(-p); e{i} = sprintf(f,re);
    elseif im >= 0; e{i} = [sprintf(f,re) '+' sprintf(f,im) 'j'];
    else; e{i} = [sprintf(f,re) '-' sprintf(f,-im) 'j']; end
end
s = ['\begin{bmatrix}' strjoin(e, ' \\\\ ') '\end{bmatrix}'];
end

function v = subeig(th,S,a,b,J,Tp,M,Mk)
[~,~,Gf,Gj] = combmin(th,S,a,b,J,Tp,M,Mk);
v = min(real(eig(Gf+Gj)));
end

function f = ffun(t,om,a,b)
f = b*(exp(1i*om*t) - exp(a*t))/(1i*om - a);
end

function [q0,q1,q2] = sweepcoef(om,a,b,J,Tp,M,Mk) %#ok<INUSD>
u = @(t) [exp(a*t); 0*t];
v = @(t) [ffun(t,om,a,b); exp(1i*om*t)];
iq = @(g) integral(g, 0, Tp, 'AbsTol',1e-13,'RelTol',1e-11);
q2 = iq(@(t) real(sum(conj(u(t)).*(M*u(t)),1))) + real(u(Tp)'*Mk*u(Tp));
q1 = iq(@(t) sum(conj(u(t)).*(M*v(t)),1)) + u(Tp)'*Mk*v(Tp);
q0 = iq(@(t) real(sum(conj(v(t)).*(M*v(t)),1))) + real(v(Tp)'*Mk*v(Tp));
end

function s = sweepval(om,th,a,b,J,Tp,M,Mk)
[q0,q1,q2] = sweepcoef(om,a,b,J,Tp,M,Mk);
xi0 = J*ffun(Tp,om,a,b)/(exp(1i*th) - J*exp(a*Tp));
s = q0 + 2*real(conj(xi0)*q1) + abs(xi0)^2*q2;
end

function z = lineprof(t,om,th,a,b,J,Tp)
% single-line hybrid-harmonic profile zeta(tau) with unit amplitude
xi0 = J*ffun(Tp,om,a,b)/(exp(1i*th) - J*exp(a*Tp));
z = [exp(a*t)*xi0 + ffun(t,om,a,b); exp(1i*om*t)];
end

function z = lineprofv(t,om,th,a,b,J,Tp,F,wh,vh)
% single-line hybrid-harmonic profile (eq:diss:kyp:hyb:profile) with
% disturbance amplitude wh and jump-disturbance amplitude vh (channel F)
xi0 = (J*ffun(Tp,om,a,b)*wh + F*vh)/(exp(1i*th) - J*exp(a*Tp));
z = [exp(a*t)*xi0 + ffun(t,om,a,b)*wh; exp(1i*om*t)*wh];
end

function z = lineprofT(t,om,th,a,B,J,Tp,F,wh)
% single-line hybrid-harmonic profile (eq:diss:kyp:hyb:profile), n=1, m=numel(wh),
% with the TIED jump disturbance e^{j om Tp} wh (taken from the arc at Tp^-):
%   xi(0) = (e^{j th} - J e^{a Tp})^{-1} (J g_om + e^{j om Tp} F wh)
bw  = B*wh;
xi0 = (J*ffun(Tp,om,a,bw) + exp(1i*om*Tp)*(F*wh))/(exp(1i*th) - J*exp(a*Tp));
z = [exp(a*t)*xi0 + ffun(t,om,a,bw); exp(1i*om*t)*wh];
end

function verifyComb(th,rr,c,a,b,J,Tp,M,Mk,single,Djp,Dfl,lmin)
% Numerical check of the restated Lemma lem:diss:kyp:hyb:comb.
iq  = @(g) integral(g, 0, Tp, 'AbsTol',1e-13,'RelTol',1e-11,'ArrayValued',true);
om  = (th + 2*pi*rr)/Tp;  R = numel(rr);  e = exp(1i*th);  Jt = J*exp(a*Tp);
c   = c(:);
% ---- (a) membership of the sum in K_hyb (no jump disturbance, F absent) ----
zs  = @(t) sum(cell2mat(arrayfun(@(i) c(i)*lineprofv(t,om(i),th,a,b,J,Tp,0,1,0), ...
           (1:R)', 'UniformOutput', false)'), 2);
z0 = zs(0); zT = zs(Tp);
resBnd = abs(J*zT(1) - e*z0(1));                       % J xi(Tp) = e^{j th} xi(0)
resBndCov = abs(abs(J*zT(1))^2 - abs(z0(1))^2);        % E Y(0) E^* = E J(Y(Tp^-)) E^*
wfun = @(t) sum(c.*exp(1i*om(:)*t), 1);
opts = odeset('RelTol',1e-11,'AbsTol',1e-13);
[tt,xx] = ode45(@(t,x) a*x + b*wfun(t), linspace(0,Tp,201), z0(1), opts);
zz = zeros(numel(tt),1); for k = 1:numel(tt); q = zs(tt(k)); zz(k) = q(1); end
resFlow = max(abs(xx - zz));                           % xi' = a xi + b w
% direct value of Sigma(Y) and decomposition (a)
SigDirect = real(iq(@(t) zs(t)'*M*zs(t))) + real(zT'*Mk*zT);
resDecomp = abs(SigDirect - (single + Djp + Dfl));
resDirect = abs(SigDirect - lmin);
% ---- (b) transient formula and Delta_fl from the transients ----
resDelta = 0; resProf = 0; tg = linspace(0,Tp,401);
zh = zeros(2,R); dl = zeros(1,R);
for i = 1:R
    xh = b*c(i)/(1i*om(i) - a);  zh(:,i) = [xh; c(i)];
    p0 = c(i)*lineprofv(0,om(i),th,a,b,J,Tp,0,1,0);
    dl(i) = p0(1) - xh;
    dform = (e*(J-1)*xh)/(e - Jt);
    resDelta = max(resDelta, abs(dl(i) - dform));
    for t = tg
        p = c(i)*lineprofv(t,om(i),th,a,b,J,Tp,0,1,0);
        resProf = max(resProf, abs(p(1) - (exp(1i*om(i)*t)*xh + exp(a*t)*dl(i))));
    end
end
DflTrans = 0;
for i = 1:R
    for k = 1:R
        if i == k; continue; end
        hi  = @(t) exp(1i*om(i)*t)*zh(:,i);  hk = @(t) exp(1i*om(k)*t)*zh(:,k);
        ei  = @(t) [exp(a*t)*dl(i); 0];      ek = @(t) [exp(a*t)*dl(k); 0];
        % <M, h_i eta_k^* + eta_i h_k^* + eta_i eta_k^*> = eta_k^* M h_i + h_k^* M eta_i + eta_k^* M eta_i
        DflTrans = DflTrans + iq(@(t) ek(t)'*M*hi(t) + hk(t)'*M*ei(t) + ek(t)'*M*ei(t));
    end
end
resDflTrans = abs(real(DflTrans) - Dfl) + abs(imag(DflTrans));
harmCross = 0;                                         % harmonic-harmonic cross terms
for i = 1:R; for k = 1:R
    if i ~= k; harmCross = max(harmCross, abs(iq(@(t) exp(1i*(om(i)-om(k))*t)))); end
end; end
% ---- (b') transient formula with a jump disturbance F=1, tied vhat ----
% tied jump disturbance: vhat_r = e^{j om_r Tp} what_r = e^{j th} c_r
F = 1; resDeltaV = 0;
for i = 1:R
    xh = b*c(i)/(1i*om(i) - a);
    p0 = lineprofv(0,om(i),th,a,b,J,Tp,F,c(i),exp(1i*om(i)*Tp)*c(i));
    resDeltaV = max(resDeltaV, abs((p0(1) - xh) - e*((J-1)*xh + F*c(i))/(e - Jt)));
end
% ---- (c) pure comb, m=2 variant with tied jump disturbance ----
% B2 = [b 0], F2 = [0 1]; line amplitude what_r = [c_r; (1-J) xihat_r] gives
% (1-J) xihat_r = F2 what_r, i.e. delta_r = 0.  Random Hermitian 3x3 supplies.
rng(1);
B2 = [b 0]; F2 = [0 1];
G = randn(3) + 1i*randn(3); M3  = (G + G')/2;
G = randn(3) + 1i*randn(3); Mk3 = (G + G')/2;
wp = zeros(2,R); zh3 = zeros(3,R);
for i = 1:R
    xh = b*c(i)/(1i*om(i) - a);
    wp(:,i) = [c(i); (1-J)*xh];
    zh3(:,i) = [xh; wp(:,i)];
end
prof = @(t,i) lineprofT(t,om(i),th,a,B2,J,Tp,F2,wp(:,i));
zp  = @(t) sum(cell2mat(arrayfun(@(i) prof(t,i), (1:R)', 'UniformOutput', false)'), 2);
resPureDelta = 0;
for i = 1:R; p0 = prof(0,i); resPureDelta = max(resPureDelta, abs(p0(1) - zh3(1,i))); end
DflPure = 0; DjpPure = 0; singlePure = 0;
for i = 1:R
    for k = 1:R
        val = iq(@(t) prof(t,k)'*M3*prof(t,i));
        jv  = prof(Tp,k)'*Mk3*prof(Tp,i);
        if i == k; singlePure = singlePure + real(val + jv);
        else; DflPure = DflPure + val; DjpPure = DjpPure + jv; end
    end
end
% boundary: the jump acts on the arc's own pre-jump disturbance block zpT(2:3)
zpT = zp(Tp); zp0 = zp(0);
resPureBnd = abs(J*zpT(1) + F2*zpT(2:3) - e*zp0(1));
SigPure = real(iq(@(t) zp(t)'*M3*zp(t))) + real(zpT'*Mk3*zpT);
zhS = sum(zh3,2);
SigComb = Tp*real(sum(sum(conj(zh3).*(M3*zh3),1))) + real(zhS'*Mk3*zhS);
resComb = abs(SigPure - SigComb);
DjpComb = real(zhS'*Mk3*zhS) - real(sum(sum(conj(zh3).*(Mk3*zh3),1)));
resDjp = abs(real(DjpPure) - DjpComb);
fprintf(['  Lemma lem:diss:kyp:hyb:comb check (witness, theta=%.4f):\n', ...
    '    (a) boundary |J xi(Tp)-e^{j th} xi(0)| = %.2e, state-block boundary %.2e, flow ODE residual %.2e\n', ...
    '        Sigma direct %.6f; |direct - (single+Djp+Dfl)| = %.2e; |direct - witness value| = %.2e\n', ...
    '    (b) |delta_r - formula| = %.2e (F absent), %.2e (F=1, tied vhat=e^{j om Tp}what); profile split residual %.2e\n', ...
    '        Dfl from transients %.6f vs Gram %.6f, residual %.2e; max |int e^{j(om_r-om_r'')t}| = %.2e\n', ...
    '        transients |delta_r| = %s\n', ...
    '    (c) pure comb (m=2, B=[b 0], F=[0 1], tied vhat, what_r=[c_r;(1-J)xihat_r]): max |xi_r(0)-xihat_r| = %.2e, boundary %.2e,\n', ...
    '        Dfl = %.2e, |Sigma - (Tp sum <M,zh zh^*> + <Mk,zhS zhS^*>)| = %.2e, |Djp - comb Djp| = %.2e,\n', ...
    '        Sigma = %.6f = single %.6f + Djp %.6f\n'], ...
    th, resBnd, resBndCov, resFlow, SigDirect, resDecomp, resDirect, resDelta, resDeltaV, resProf, ...
    real(DflTrans), Dfl, resDflTrans, harmCross, mat2str(abs(dl),3), resPureDelta, resPureBnd, ...
    abs(DflPure), resComb, resDjp, SigPure, singlePure, real(DjpPure));
assert(max([resBnd resBndCov resDecomp resDirect resDelta resDeltaV resProf resDflTrans ...
    resPureDelta resPureBnd abs(DflPure) resComb resDjp]) < 1e-8 && resFlow < 1e-7, ...
    'Lemma lem:diss:kyp:hyb:comb verification failed');
end

function [lmin, c, Gf, Gj] = combmin(th,rr,a,b,J,Tp,M,Mk)
R = numel(rr); om = (th + 2*pi*rr)/Tp;
Gf = zeros(R); Gj = zeros(R);
for i = 1:R
    for k = 1:R
        gi = @(t) sum(conj(lineprof(t,om(i),th,a,b,J,Tp)).*(M*lineprof(t,om(k),th,a,b,J,Tp)),1);
        Gf(i,k) = integral(gi, 0, Tp, 'AbsTol',1e-13,'RelTol',1e-11);
        Gj(i,k) = lineprof(Tp,om(i),th,a,b,J,Tp)'*Mk*lineprof(Tp,om(k),th,a,b,J,Tp);
    end
end
Gf = (Gf+Gf')/2; Gj = (Gj+Gj')/2;
[V,D] = eig(Gf+Gj); [lmin,k] = min(real(diag(D))); c = V(:,k);
end
