function ex_stabz_switched()
%EX_STABZ_SWITCHED  Switched system: mean-square stabilization by mode-dependent
%  state feedback under a dwell-time.  Synthesis mirror of ex:stab:switched.
%  Two unstable controlled modes (A_i,B), state continuous at the switches,
%  each carrying its own gain K_i obtained by solving the continuous-time flow
%  synthesis condition (flow channel Ac) as a matrix sum-of-squares program
%  (arbitrary dwell-time) through mvst_stabz_sos, verified through the
%  one-cycle switched lifted monodromy.
%  Periodic switching of dwell-time Tsw.  Open loop (K=0): unstable modes.

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

n = 2; m = 1;
A1 = [0 1; -1 0.1]; A2 = [0 1; -3 0.2]; B = [0;1]; Ac = 0.2*eye(n);
B1 = 0.1*[0; 1];        % input selector actuating the diffusion channel (eq:stabz:sw:data)
Tsw = 0.5;

Mo = swmono({A1,A2}, {zeros(m,n),zeros(m,n)}, B, Ac, B1, Tsw, n);
r1 = mvst_stabz_sos(struct('A0',A1,'B0',B,'Ac',Ac,'Bc',B1), struct('type','arbitrary'));
r2 = mvst_stabz_sos(struct('A0',A2,'B0',B,'Ac',Ac,'Bc',B1), struct('type','arbitrary'));
Mc_K1 = r1.Kfun(0);  Mc_K2 = r2.Kfun(0);
Mc = swmono({A1,A2}, {Mc_K1,Mc_K2}, B, Ac, B1, Tsw, n);
radOpen = Mo; radClosed = Mc;   % one-cycle monodromy radii

T = 8*Tsw; dt = 2e-3; nreal = 16; x0 = [1;-1];
[t, Xo, ~ ] = sim_sw({A1,A2}, {zeros(m,n),zeros(m,n)}, B, Ac, B1, Tsw, x0, T, dt, nreal, 61);
[~, Xc, Uc] = sim_sw({A1,A2}, {Mc_K1,Mc_K2},          B, Ac, B1, Tsw, x0, T, dt, nreal, 61);

[tc, trO] = swtr({A1,A2}, {zeros(m,n),zeros(m,n)}, B, Ac, B1, Tsw, round(T/Tsw), n);
[~,  trC] = swtr({A1,A2}, {Mc_K1,Mc_K2},          B, Ac, B1, Tsw, round(T/Tsw), n);

mvst_stabfig(figdir, 'switched', t, Xo, Xc, Uc, tc, trO, trC, false, gc);
mvst_stabtable(figdir, 'switched');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_switched.tex'));
w.num('StabzSwOpen', radOpen); w.num('StabzSwClosed', radClosed);
w.mat('StabzSwKone', Mc_K1); w.mat('StabzSwKtwo', Mc_K2); w.num('StabzSwTsw', Tsw); w.close();
fprintf('Switched stabilization: open monodromy %.4f, closed %.4f\n', radOpen, radClosed);
end

% ===========================================================================
function rad = swmono(As, Ks, B, Ac, B1, Tsw, n)
d = size(mvst_svecbasis(n),2); Phi = eye(d);
for i = 1:numel(As)
    Mi = mvst_liftsym(As{i}+B*Ks{i}, {Ac+B1*Ks{i}});   % channel actuated by B1
    Phi = expm(Tsw*Mi) * Phi;
end
rad = max(abs(eig(Phi)));
end

function [t, X, U] = sim_sw(As, Ks, B, Ac, B1, Tsw, x0, T, dt, nreal, seed)
rng(seed); n = numel(x0); m = size(Ks{1},1); nt = round(T/dt)+1; t = (0:nt-1)*dt;
X = zeros(n,nt,nreal); U = zeros(m,nt,nreal);
for r = 1:nreal
    x = x0;
    for k = 1:nt
        mode = mod(floor(t(k)/Tsw), numel(As)) + 1;
        X(:,k,r) = x; U(:,k,r) = Ks{mode}*x;
        if k < nt
            Acl = As{mode} + B*Ks{mode};
            x = x + Acl*x*dt + ((Ac+B1*Ks{mode})*x)*(sqrt(dt)*randn);
        end
    end
end
end

function [t, tr] = swtr(As, Ks, B, Ac, B1, Tsw, Ncyc, n)
Us = mvst_svecbasis(n); s = Us'*reshape(eye(n),[],1);
t = []; tr = []; t0 = 0; npt = 30; nm = numel(As);
for c = 1:Ncyc
    mode = mod(c-1, nm) + 1; Mi = mvst_liftsym(As{mode}+B*Ks{mode}, {Ac+B1*Ks{mode}});
    tt = linspace(0, Tsw, npt);
    for j = 1:npt
        t(end+1) = t0+tt(j);                                      %#ok<AGROW>
        tr(end+1) = trace(reshape(Us*(expm(tt(j)*Mi)*s), n, n));  %#ok<AGROW>
    end
    s = expm(Tsw*Mi)*s; t0 = t0 + Tsw;
end
end
