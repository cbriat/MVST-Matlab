function ex_stabz_sd()
%EX_STABZ_SD  Sampled-data system: mean-square stabilization by design and
%  certify.  Synthesis mirror of ex:stab:sd.  The plant with multiplicative
%  noise
%     dx = (A x + B u) dt + Ac x dW1,   u = K x(t_k) held on (t_k,t_{k+1}],
%  is sampled at the fixed period Tbar.  The sampled gain sits inside the
%  non-invertible refresh map Jtil = [I 0; K 0], so the synthesis is NOT
%  convex; the gain is DESIGNED from the continuous flow synthesis condition on
%  the plant (solved as a matrix sum-of-squares program with arbitrary
%  dwell-time through mvst_stabz_sos) and the closed
%  loop is then CERTIFIED by the periodic-sampling design condition (lifted
%  augmented one-period monodromy radius < 1).  Open/closed measure: that
%  monodromy radius.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

n = 2; m = 1;
A = [0 1; 2 0.1]; B = [0; 1]; Ac = 0.2*eye(n); Tbar = 0.25;
B1 = 0.1*[0; 1];        % input selector actuating the diffusion channel B_{u,1} (eq SD data)

rc = mvst_stabz_sos(struct('A0',A,'B0',B,'Ac',Ac,'Bc',B1), struct('type','arbitrary'));
Kc = rc.Kfun(0);                               % continuous design (flow SOS)
[radOpen,  Mfo, MJo] = sdmono(A, B, Ac, B1, zeros(m,n), Tbar, n, m);
[radClosed,Mfc, MJc] = sdmono(A, B, Ac, B1, Kc,          Tbar, n, m);

T = 8*Tbar; dt = 1e-3; nreal = 16; x0 = [1; -1];
[t, Xo, ~ ] = sim_sd(A, B, Ac, B1, zeros(m,n), Tbar, x0, T, dt, nreal, 71);
[~, Xc, Uc] = sim_sd(A, B, Ac, B1, Kc,          Tbar, x0, T, dt, nreal, 71);

n2 = n + m;
[tc, trO] = sdtr(Mfo, MJo, Tbar, round(T/Tbar), n2, n);
[~,  trC] = sdtr(Mfc, MJc, Tbar, round(T/Tbar), n2, n);

mvst_stabfig(figdir, 'sd', t, Xo, Xc, Uc, tc, trO, trC, false, gc);
mvst_stabtable(figdir, 'sd');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_sd.tex'));
w.num('StabzSdOpen', radOpen); w.num('StabzSdClosed', radClosed);
w.mat('StabzSdK', Kc); w.num('StabzSdTbar', Tbar); w.close();
fprintf('Sampled-data stabilization: open monodromy %.4f, closed %.4f\n', radOpen, radClosed);
end

% ===========================================================================
function [rad, Mf, MJ] = sdmono(A, B, Ac, B1, K, Tbar, n, m)
Atil  = [A B; zeros(m,n) zeros(m)];
Actil = [Ac B1; zeros(m,n) zeros(m)];   % held input actuates the Brownian channel
Jtil  = [eye(n) zeros(n,m); K zeros(m)];
Mf = mvst_liftsym(Atil, {Actil});
Ud = mvst_svecbasis(n+m); MJ = Ud' * kron(Jtil,Jtil) * Ud;
rad = max(abs(eig(MJ*expm(Tbar*Mf))));
end

function [t, X, U] = sim_sd(A, B, Ac, B1, K, Tbar, x0, T, dt, nreal, seed)
rng(seed); n = numel(x0); m = size(K,1); nt = round(T/dt)+1; t = (0:nt-1)*dt;
X = zeros(n,nt,nreal); U = zeros(m,nt,nreal); jt = Tbar*(1:ceil(T/Tbar));
for r = 1:nreal
    x = x0; u = zeros(m,1); nextj = 1;
    for k = 1:nt
        X(:,k,r) = x; U(:,k,r) = u;
        if k < nt
            x = x + (A*x + B*u)*dt + (Ac*x + B1*u)*(sqrt(dt)*randn);
            if nextj <= numel(jt) && t(k+1) >= jt(nextj) - dt/2
                u = K*x; nextj = nextj + 1;
            end
        end
    end
end
end

function [t, tr] = sdtr(Mf, MJ, Tbar, Ncyc, n2, n)
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
