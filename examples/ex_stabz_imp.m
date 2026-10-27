function ex_stabz_imp()
%EX_STABZ_IMP  Impulsive system with deterministic (periodic) impulses:
%  mean-square state-feedback stabilization.  Synthesis mirror of the analysis
%  example ex:stab:imp (ex:stab:impstabz), on the MASTER impulsive datum:
%     flow:  dx = (A0 x + B0 u) dt + (A1 x + B1 u) dW1,   A0 saddle,  u = K(tau) x,
%     jump:  x(t_k^+) = J0 x + B_{J,0} u_J + J1 x xi,  t_k = k Tbar,  u_J = K_J x,
%  with A0 = [0 1; 2 0.1] the master saddle, A1 = 0.3 I2 the (input-actuated,
%  selector B1) Brownian channel, and the reset data
%  J0 = [0.7 0.3; 0 0.7] (non-normal shear, spectrum {0.7,0.7}), J1 = [0.2 0.1;
%  -0.1 0.2] (rotation-scaled stochastic reset channel, spectrum 0.2 +- 0.1i,
%  radius 0.2236), B_{J,0} = [0;1].  The
%  flow gain K(tau) and the jump gain K_J are obtained together from the
%  periodic synthesis condition (matrix SOS in the timer, degree 2) through
%  mvst_stabz_sos with the .Bc / .BJ selectors and the two-channel reset .J.
%  Open/closed measure: one-period lifted monodromy radius (closed-loop flow
%  integrated with the timer gain over a fine expm-step product, closed reset
%  J0 + B_{J,0} K_J plus the stochastic channel J1).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

% ------------------- system (master impulsive datum) -------------------
n = 2; m = 1;
A0 = [0 1; 2 0.1];         % master saddle (open-loop unstable flow)
B0 = [0; 1];
A1 = 0.3*eye(n);           % Brownian channel
B1 = 0.1*[0; 1];           % input selector actuating the Brownian channel
J0 = [0.7 0.3; 0 0.7];     % deterministic reset, non-normal shear
J1 = [0.2 0.1; -0.1 0.2];  % stochastic reset channel, rotation-scaled
BJ0 = [0; 1];              % jump input selector B_{J,0}
Tbar = 0.6;

nstep = 40;                                     % fine expm steps over [0,Tbar]
Mfo = mvst_liftsym(A0, {A1});                   % open-loop lifted flow operator
MJo = mvst_liftsym(zeros(n), {J0, J1});         % open-loop lifted jump (2 channels)
monoOpen = max(abs(eig(MJo*expm(Tbar*Mfo))));   % > 1 : unstable

% ------------------- periodic timer synthesis (flow gain K(tau) + jump K_J) --
sys = struct('A0', A0, 'B0', B0, 'Ac', A1, 'Bc', B1, ...
             'J', cat(3, J0, J1), 'BJ', cat(3, BJ0, zeros(n,m)));
res = mvst_stabz_sos(sys, struct('type','periodic','Tbar',Tbar,'deg',2));
Kfun = @(tau) res.Kfun(tau);                    % timer-dependent flow gain
KJ   = res.KJ;                                   % scheduled-jump gain
J0cl = J0 + BJ0*KJ;                              % closed-loop deterministic reset
MJc  = mvst_liftsym(zeros(n), {J0cl, J1});       % closed-loop lifted jump
Phic = flowtrans(A0, B0, A1, B1, Kfun, Tbar, n, nstep);
monoClosed = max(abs(eig(MJc*Phic)));           % < 1 : stable

% ------------------- sample paths (flow Euler + periodic stochastic jumps) ---
T = 6*Tbar; dt = 2e-3; nreal = 16; x0 = [1; -1];
[t, Xo, ~ ] = sim_imp(A0, B0, A1, B1, J0, J1, zeros(m,n), @(tau) zeros(m,n), x0, T, dt, Tbar, nreal, 7);
[~, Xc, Uc] = sim_imp(A0, B0, A1, B1, J0, J1, KJ,          Kfun,              x0, T, dt, Tbar, nreal, 7);

% ------------------- covariance trace (exact saw) --------------------------
[tc, trO] = saw_tv(A0, B0, A1, B1, @(tau) zeros(m,n), MJo, Tbar, round(T/Tbar), n, nstep);
[~,  trC] = saw_tv(A0, B0, A1, B1, Kfun,              MJc, Tbar, round(T/Tbar), n, nstep);

mvst_stabfig(figdir, 'imp', t, Xo, Xc, Uc, tc, trO, trC, false, gc);
mvst_stabtable(figdir, 'imp');   %% generated summary-table row

w = mvst_texdata(fullfile(figdir,'data','ex_stabz_imp.tex'));
w.num('StabzImpOpen',   monoOpen);
w.num('StabzImpClosed', monoClosed);
w.num('StabzImpTbar',   Tbar);
w.mat('StabzImpKzero', res.Kfun(0));
w.mat('StabzImpKTbar', res.Kfun(Tbar));
w.mat('StabzImpKJ',    KJ);
w.close();
fprintf('Impulsive stabilization: open monodromy %.4f, closed %.4f, KJ=[%.3f %.3f]\n', ...
        monoOpen, monoClosed, KJ(1), KJ(2));
end

% ===========================================================================
function Phi = flowtrans(A0, B0, A1, B1, Kfun, Tbar, n, nstep)
% Lifted closed-loop flow transition over one period [0,Tbar] with the
% timer-dependent gain, as a fine expm-step product; the control is routed
% into the drift (B0) and into the Brownian channel (B1).
d = size(mvst_svecbasis(n),2); Phi = eye(d);
tt = linspace(0, Tbar, nstep+1); dtau = Tbar/nstep;
for j = 1:nstep
    Kt = Kfun(tt(j));
    Mf = mvst_liftsym(A0 + B0*Kt, {A1 + B1*Kt});
    Phi = expm(dtau*Mf)*Phi;
end
end

function [t, X, U] = sim_imp(A0, B0, A1, B1, J0, J1, KJ, Kfun, x0, T, dt, Tbar, nreal, seed)
rng(seed);
n = numel(x0); m = size(Kfun(0),1); BJ0 = [0; 1];
nt = round(T/dt) + 1; t = (0:nt-1)*dt;
X = zeros(n,nt,nreal); U = zeros(m,nt,nreal);
jumpt = Tbar*(1:ceil(T/Tbar));            % scheduled reset times
for r = 1:nreal
    x = x0; nextj = 1; tlast = 0;         % tlast: time of the last jump
    for k = 1:nt
        Kk = Kfun(t(k) - tlast);          % gain at the current timer value
        X(:,k,r) = x; U(:,k,r) = Kk*x;
        if k < nt
            x = x + (A0*x + B0*(Kk*x))*dt + ((A1 + B1*Kk)*x)*(sqrt(dt)*randn);
            if nextj <= numel(jumpt) && t(k+1) >= jumpt(nextj) - dt/2
                x = J0*x + BJ0*(KJ*x) + (J1*x)*randn;   % reset + jump control + stochastic channel
                nextj = nextj + 1; tlast = t(k+1);
            end
        end
    end
end
end

function [t, tr] = saw_tv(A0, B0, A1, B1, Kfun, MJ, Tbar, Ncyc, n, nstep)
% Covariance saw with a timer-dependent flow gain and the closed jump operator
% MJ (already carrying J0 + B_{J,0} K_J and the stochastic channel J1).
Us = mvst_svecbasis(n); s = Us'*reshape(eye(n),[],1);
t = []; tr = []; t0 = 0; tt = linspace(0, Tbar, nstep+1); dtau = Tbar/nstep;
for c = 1:Ncyc
    sk = s;
    for j = 1:nstep+1
        t(end+1)  = t0 + tt(j);                                %#ok<AGROW>
        tr(end+1) = trace(reshape(Us*sk, n, n));               %#ok<AGROW>
        if j <= nstep
            Kt = Kfun(tt(j));
            Mf = mvst_liftsym(A0 + B0*Kt, {A1 + B1*Kt});
            sk = expm(dtau*Mf)*sk;
        end
    end
    s = MJ*sk; t0 = t0 + Tbar;             % jump at period end (closed operator)
end
end
