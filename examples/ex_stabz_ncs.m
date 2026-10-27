function ex_stabz_ncs()
%EX_STABZ_NCS  Networked control with packet dropouts: mean-square
%  stabilization by design and certify.  Synthesis mirror of ex:stab:ncs.
%  Continuous plant, discrete controller over a lossy network: at each sample
%  t_k = k Tbar the input is refreshed to K x(t_k) with probability p (packet
%  delivered) and held with probability 1-p (dropout).  On the augmented state
%  z = [x;u] the flow is dz = Atil z dt + Actil z dW1 and the mean jump is
%     J(Z) = p Jtil Z Jtil' + (1-p) Z,   Jtil = [I 0; K 0].
%  The gain K sits in the refresh map (non-convex), so it is DESIGNED from the
%  continuous flow synthesis condition on the plant (solved as a matrix
%  sum-of-squares program with arbitrary dwell-time through mvst_stabz_sos) and
%  CERTIFIED by the periodic design condition (lifted
%  one-period monodromy radius < 1).  Open/closed measure: that radius.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

n = 2; m = 1;
A = [0 1; 0 0]; B = [0; 1]; Ac = zeros(n); Tbar = 0.2; p = 0.8;   % deterministic master double integrator
% Design the transmitted gain K on the SAMPLED (discretized) delivered map,
% then CERTIFY it on the periodic two-outcome dropout monodromy.  (A continuous
% flow gain is too aggressive for the ZOH period; the sampled design is the
% correct one for the held-input dropout loop.)
Ad = expm(A*Tbar);
Bd = intexpm(A, Tbar) * B;                       % zero-order-hold input matrix
rc = mvst_stabz_sos(struct('A0',Ad,'B0',Bd), struct('type','discrete','rho',0.7));
Kc = rc.Kfun(0);
[radOpen,  Mfo, MJo] = ncsmono(A, B, Ac, zeros(m,n), Tbar, p, n, m);
[radClosed,Mfc, MJc] = ncsmono(A, B, Ac, Kc,          Tbar, p, n, m);

T = 8*Tbar; dt = 1e-3; nreal = 16; x0 = [1; -1];
[t, Xo, ~ ] = sim_ncs(A, B, Ac, zeros(m,n), Tbar, p, x0, T, dt, nreal, 91);
[~, Xc, Uc] = sim_ncs(A, B, Ac, Kc,          Tbar, p, x0, T, dt, nreal, 91);

n2 = n + m;
[tc, trO] = ncstr(Mfo, MJo, Tbar, round(T/Tbar), n2, n);
[~,  trC] = ncstr(Mfc, MJc, Tbar, round(T/Tbar), n2, n);

mvst_stabfig(figdir, 'ncs', t, Xo, Xc, Uc, tc, trO, trC, false, gc);
mvst_stabtable(figdir, 'ncs');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_ncs.tex'));
w.num('StabzNcsOpen', radOpen); w.num('StabzNcsClosed', radClosed);
w.mat('StabzNcsK', Kc); w.num('StabzNcsTbar', Tbar); w.num('StabzNcsP', p); w.close();
fprintf('NCS stabilization: open monodromy %.4f, closed %.4f\n', radOpen, radClosed);
end

% ===========================================================================
function M = intexpm(A, T)
% integral_0^T expm(A s) ds via the augmented-matrix trick.
n = size(A,1); E = expm([A eye(n); zeros(n) zeros(n)]*T); M = E(1:n, n+1:end);
end

function [rad, Mf, MJ] = ncsmono(A, B, Ac, K, Tbar, p, n, m)
Atil  = [A B; zeros(m,n) zeros(m)];
Actil = [Ac zeros(n,m); zeros(m,n) zeros(m)];
Jtil  = [eye(n) zeros(n,m); K zeros(m)];
Mf = mvst_liftsym(Atil, {Actil});
Ud = mvst_svecbasis(n+m); d = size(Ud,2);
MJ = p*(Ud'*kron(Jtil,Jtil)*Ud) + (1-p)*eye(d);      % mean dropout jump
rad = max(abs(eig(MJ*expm(Tbar*Mf))));
end

function [t, X, U] = sim_ncs(A, B, Ac, K, Tbar, p, x0, T, dt, nreal, seed)
rng(seed); n = numel(x0); m = size(K,1); nt = round(T/dt)+1; t = (0:nt-1)*dt;
X = zeros(n,nt,nreal); U = zeros(m,nt,nreal); jt = Tbar*(1:ceil(T/Tbar));
for r = 1:nreal
    x = x0; u = zeros(m,1); nextj = 1;
    for k = 1:nt
        X(:,k,r) = x; U(:,k,r) = u;
        if k < nt
            x = x + (A*x + B*u)*dt + (Ac*x)*(sqrt(dt)*randn);
            if nextj <= numel(jt) && t(k+1) >= jt(nextj) - dt/2
                if rand < p; u = K*x; end          % delivered, else held
                nextj = nextj + 1;
            end
        end
    end
end
end

function [t, tr] = ncstr(Mf, MJ, Tbar, Ncyc, n2, n)
Ud = mvst_svecbasis(n2); Z0 = zeros(n2); Z0(1:n,1:n) = eye(n);
s = Ud'*reshape(Z0,[],1); t = []; tr = []; t0 = 0; npt = 30;
for c = 1:Ncyc
    tt = linspace(0,Tbar,npt);
    for j = 1:npt
        Z = reshape(Ud*(expm(tt(j)*Mf)*s), n2, n2);
        t(end+1) = t0+tt(j); tr(end+1) = trace(Z(1:n,1:n)); %#ok<AGROW>
    end
    s = MJ*(expm(Tbar*Mf)*s); t0 = t0 + Tbar;
end
end
