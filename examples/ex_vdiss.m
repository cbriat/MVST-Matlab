%% EX_VDISS  Figures and checks for the vector dissipativity appendix.
%   Part 1: analytic validations of the worked examples (GBM mean-square
%           gain 2/|2a+sigma^2|, DT delay-with-noise gain (1-s^2)^(-1/2),
%           DT one-pole filter gain 1/(1-a) via LMI feasibility window).
%   Part 2: mean-square L2-gain of the stochastic oscillator (appendix
%           example) versus noise intensity sigma, by bisection over the
%           monotone fixed-point iteration on the generalized Riccati
%           equation, each step a generalized Lyapunov solve in the
%           reduced (svec) coordinates; sigma = 0 cross-checked against
%           Hamiltonian bisection and a frequency sweep.  Writes
%           figures/vdiss_gain.pdf.
here = fileparts(mfilename('fullpath'));
figdir = mvst_figdir();

%% ------------------------------------------------------------------
%% Part 1: analytic validations
% (i) GBM with input: dx = (a x + w)dt + s x dW, z = x.
a = -1; s = 0.6; mu = -(2*a + s^2);
gstar = 2/mu;                                    % analytic gain
gnum  = local_msgain(a, s, 1, 1, 1e-10);         % numeric (scalar lift)
fprintf('vdiss GBM: analytic gain %.6f, numeric %.6f (err %.1e)\n', ...
    gstar, gnum, abs(gstar-gnum));

% (ii) DT delay with multiplicative noise: x+ = w + s x nu, z = x.
sdt = 0.5; gdt = 1/sqrt(1 - sdt^2);
% smallest gamma with feasible p: p in [1/(1-s^2), gamma^2]
fprintf('vdiss DT noise: analytic gain %.6f, window nonempty at gamma: %d\n', ...
    gdt, 1/(1-sdt^2) <= gdt^2 + 1e-12);

% (iii) DT one-pole filter: x+ = a x + w, z = x; gain 1/(1-a).
adt = 0.7; gpole = 1/(1-adt);
feas = @(g) local_dtfeas(adt, g);                % LMI feasibility scan in p
fprintf('vdiss DT pole: gain %.4f; feasible at 1.01*gain: %d, at 0.99*gain: %d\n', ...
    gpole, feas(1.01*gpole), feas(0.99*gpole));

%% ------------------------------------------------------------------
%% Part 2: stochastic oscillator, gain versus noise intensity
A0 = [0 1; -2 -1]; E0 = [0; 1]; Cz = [1 0]; A1s = [0 0; 1 0];
% mean-square stability margin: largest sigma with lifted flow Hurwitz
f = @(sg) max(real(eig(mvst_liftsym(A0, {sg*A1s}))));
slo = 0; shi = 2; for it = 1:60, sm = (slo+shi)/2; if f(sm) < 0, slo = sm; else, shi = sm; end, end
sstar = (slo+shi)/2;
fprintf('vdiss oscillator: mean-square stability margin sigma* = %.4f\n', sstar);

sigmas = linspace(0, 0.965*sstar, 25);
gains = zeros(size(sigmas));
for i = 1:numel(sigmas)
    gains(i) = local_msgain(A0, sigmas(i)*A1s, E0, Cz, 1e-9);
end

% sigma = 0 cross-checks: Hamiltonian bisection and frequency sweep
gham = local_hinf_ham(A0, E0, Cz);
om = logspace(-2, 2, 4000);
gfreq = max(arrayfun(@(w) abs(Cz*((1i*w*eye(2)-A0)\E0)), om));
fprintf('vdiss oscillator: gain(0) lifted %.8f, Hamiltonian %.8f, freq sweep %.8f\n', ...
    gains(1), gham, gfreq);

fig = figure('Units','centimeters','Position',[2 2 12 8]); hold on
plot(sigmas, gains, 'LineWidth', 1.4);
plot(0, gham, 's', 'MarkerSize', 7, 'LineWidth', 1.2);
xline(sstar, 'k--');
xlabel('noise intensity \sigma');
ylabel('mean-square L_2-gain \gamma^*(\sigma)');
legend({'bisection over the lifted iteration', ...
        'H_\infty norm at \sigma = 0 (Hamiltonian)', ...
        'mean-square stability margin \sigma^*'}, 'Location','northwest');
title('Dissipativity degradation under multiplicative noise');
mvst_savepdf(fig, fullfile(figdir, 'vdiss_gain.pdf'));

%% ------------------------------------------------------------------
function g = local_msgain(A0, A1, E0, Cz, tol)
% Mean-square L2-gain by bisection; feasibility at gamma tested by the
% monotone fixed-point iteration on the generalized Riccati equation,
% each step a generalized Lyapunov solve in svec coordinates (the
% adjoint operator is the transpose of the lifted flow matrix).
n = size(A0,1);
M = mvst_liftsym(A0, {A1});                      % lift of F (forward)
U = mvst_svecbasis(n);
svec = @(X) U'*X(:);
unsvec = @(v) reshape(U*v, n, n);
CC = Cz'*Cz; EE = E0*E0';
glo = 0; ghi = 1;
while ~local_iter(ghi), ghi = 2*ghi; if ghi > 1e12, g = inf; return; end, end
while ghi - glo > tol*max(1,ghi)
    gm = (glo+ghi)/2;
    if local_iter(gm), ghi = gm; else, glo = gm; end
end
g = ghi;
    function ok = local_iter(gam)
        P = zeros(n); ok = false;
        for k = 1:50000
            RHS = -CC - (P*EE*P)/gam^2;
            Pn = unsvec(M.'\svec(RHS));          % adjoint gen. Lyapunov
            Pn = (Pn+Pn')/2;
            if any(~isfinite(Pn(:))) || norm(Pn,'fro') > 1e10, return; end
            if norm(Pn-P,'fro') <= 1e-12*max(1,norm(Pn,'fro')), ok = true; return; end
            P = Pn;
        end
    end
end

function g = local_hinf_ham(A, E, C)
% Deterministic Hinf norm (F = 0) by Hamiltonian bisection.
n = size(A,1);
glo = 0; ghi = 1;
noimag = @(gam) all(abs(real(eig([A, (E*E')/gam^2; -C'*C, -A']))) > 1e-9);
while ~noimag(ghi), ghi = 2*ghi; end
while ghi - glo > 1e-10*max(1,ghi)
    gm = (glo+ghi)/2;
    if noimag(gm), ghi = gm; else, glo = gm; end
end
g = ghi;
end

function ok = local_dtfeas(a, gam)
% Feasibility of the DT bounded-real LMI for x+ = a x + w, z = x, by a
% scan over the storage parameter p in the window (0, gamma^2).
ps = linspace(1e-6, gam^2*(1-1e-9), 20000);
lhs = (a^2-1)*ps + 1 + (a^2*ps.^2)./(gam^2 - ps);
ok = any(lhs < 0);
end
