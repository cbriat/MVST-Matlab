function ex_cones_quat()
%EX_CONES_QUAT  Stability of a quaternionic Hermitian matrix-valued system
%  through the complex adjoint, for the chapter "Generalization to Proper
%  Cones" (section "Quaternionic Hermitian matrices").
%
%  A quaternion matrix M = A + B j (A,B complex) is represented by its
%  complex adjoint
%        chi(M) = [ A  B ; -conj(B)  conj(A) ]  in  C^{2n x 2n},
%  an injective *-homomorphism with  M >= 0  iff  chi(M) >= 0. A quaternionic
%  system on H^n(Q) therefore reduces to a complex system on H^{2n}, and a
%  stability certificate is computed there as an ordinary SDP/Lyapunov
%  problem. The certificate is returned to the quaternions by AVERAGING over
%  the Kramers involution  sigma(Y) = Omega conj(Y) Omega',  Omega=[0 -I;I 0]:
%  the averaged certificate P = (Ptilde + sigma(Ptilde))/2 is sigma-fixed,
%  hence chi of a quaternionic Hermitian matrix.
%
%  System (n = 2, one diffusion channel), completely positive flow
%        F(X) = A X + X A* + A1 X A1*,   X in H^2(Q),
%  with quaternionic A, A1. The dual Lyapunov inequality is
%        F*(P) = A* P + P A + A1* P A1  <  0,   P in H^2(Q), P > 0.
%
%  Verifications performed:
%   (a) chi is a *-homomorphism on the chosen data (chi(A A1)=chi(A)chi(A1)).
%   (b) the doubled complex flow F_chi is exponentially stable
%       (spectral abscissa of the lifted generator < 0).
%   (c) a complex certificate Ptilde > 0 with F_chi*(Ptilde) < 0
%       (closed-form Lyapunov solve).
%   (d) the averaged certificate P is quaternionic: the Kramers residual
%       ||sigma(P) - P|| is at machine precision, and P has the chi block
%       structure  P = [P11 P12; -conj(P12) conj(P11)].
%   (e) P certifies the quaternionic system: max eig(F_chi*(P)) < 0.
%   (f) the covariance value t -> <P, X(t)> decays monotonically.

here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

n = 2;
I2 = eye(2); Om = [zeros(n) -eye(n); eye(n) zeros(n)];   % 4x4 symplectic
chi   = @(A,B) [A B; -conj(B) conj(A)];                  % complex adjoint
sigma = @(Y) Om*conj(Y)*Om';                             % Kramers involution

% ---- quaternionic coefficients A = Aa + Ab j, A1 = A1a + A1b j ----
Aa = [-1.6+0.3i, 0.4-0.2i; -0.3+0.1i, -1.9-0.4i];
Ab = [ 0.2+0.1i, 0.15;      0.1-0.05i, 0.25+0.2i];
A1a= [ 0.25,      0.1+0.1i;  0.0,       0.2-0.1i];
A1b= [ 0.1i,      0.05;      0.05-0.05i, 0.1];
Ac  = chi(Aa,Ab);      % chi(A)   in C^{4x4}
A1c = chi(A1a,A1b);    % chi(A1)  in C^{4x4}

% ---- (a) chi is a *-homomorphism (product law check) ----
% chi(A A1): compute quaternionic product then adjoint, compare with Ac*A1c
Pa = Aa*A1a - Ab*conj(A1b);
Pb = Aa*A1b + Ab*conj(A1a);
homres = norm(chi(Pa,Pb) - Ac*A1c,'fro');
fprintf('(a) chi homomorphism residual ||chi(A A1)-chi(A)chi(A1)|| = %.2e\n', homres);

% ---- (b) lifted generator of the complex flow and its spectral abscissa ----
% vec(Xdot) = L vec(X),  L = I(x)Ac + conj(Ac)(x)I + conj(A1c)(x)A1c
d = 2*n;
L = kron(eye(d),Ac) + kron(conj(Ac),eye(d)) + kron(conj(A1c),A1c);
absc = max(real(eig(L)));
fprintf('(b) spectral abscissa of lifted flow = %.4f (stable: %d)\n', absc, absc<0);

% ---- (c) complex certificate: F_chi*(P) < 0, P > 0 ----
Ptilde = solve_lyap(Ac,A1c,d);

% ---- (d) average over the Kramers involution -> quaternionic certificate ----
P = 0.5*(Ptilde + sigma(Ptilde));
P = 0.5*(P + P');                       % enforce Hermitian symmetry numerically
kres = norm(sigma(P) - P,'fro');
P11 = P(1:n,1:n); P12 = P(1:n,n+1:end); P21 = P(n+1:end,1:n); P22 = P(n+1:end,n+1:end);
struct_res = norm(P21 + conj(P12),'fro') + norm(P22 - conj(P11),'fro');
fprintf('(d) Kramers residual ||sigma(P)-P||    = %.2e\n', kres);
fprintf('    chi-block residual                 = %.2e\n', struct_res);

% ---- (e) P certifies the quaternionic system ----
Fstar = Ac'*P + P*Ac + A1c'*P*A1c;      % F_chi*(P)
mev = max(real(eig(0.5*(Fstar+Fstar'))));
pev = min(real(eig(0.5*(P+P'))));
fprintf('(e) max eig F_chi*(P) = %.4f (<0: %d) ; min eig P = %.4f (>0: %d)\n', ...
        mev, mev<0, pev, pev>0);

% ---- (f) covariance value decay <P, X(t)> ----
tt = linspace(0,6,300);
X0 = chi([1 0;0 0.5],[0 0.2;-0.2 0]); X0 = 0.5*(X0+X0');   % X0 in H^2(Q), >=0
X0 = X0 + 1e-6*eye(d);
vX0 = X0(:);
val = zeros(size(tt));
for kk = 1:numel(tt)
    Xt = reshape(expm(L*tt(kk))*vX0, d, d);
    Xt = 0.5*(Xt+Xt');
    val(kk) = real(trace(P*Xt));        % <P,X(t)> = (1/2) Re tr(chi(P)chi(X))
end
fprintf('(f) value decay: <P,X0>=%.4f  ->  <P,X(T)>=%.4f (monotone: %d)\n', ...
        val(1), val(end), all(diff(val)<=1e-9));

% ============================ figure ============================
figure('Color','w','Units','centimeters','Position',[2 2 12 7]);
plot(tt, val, '-', 'Color', gc.sapphire, 'LineWidth', 1.8); hold on;
xlabel('$t$','Interpreter','latex');
ylabel('$\langle P, X(t)\rangle$','Interpreter','latex');
title('Quaternionic covariance decay under the averaged certificate', ...
      'Interpreter','latex','FontWeight','normal');
grid on; box on;
% inset annotation: the Kramers residual certifies quaternionic structure
text(0.55, 0.85, sprintf('$\\|\\sigma(P)-P\\| = %.1e$', kres), ...
     'Units','normalized','Interpreter','latex','Color',gc.ruby);
try
    mvst_savepdf(gcf, fullfile(figdir,'cones_quat.pdf'));
catch, end

% ---------------------- LaTeX data macros ----------------------
try
    wt = mvst_texdata(fullfile(figdir,'data','cones_quat.tex'));
    wt.num('ConeQuatAbsc', absc);
    wt.num('ConeQuatHomRes', homres);
    wt.num('ConeQuatKramers', kres);
    wt.num('ConeQuatMaxEigFstar', mev);
    wt.num('ConeQuatMinEigP', pev);
    wt.mat('ConeQuatP', P, 3);
    wt.close();
catch, end
end

% ======================================================================
function P = solve_lyap(Ac,A1c,d)
%SOLVE_LYAP  P > 0 with  Ac'P + P Ac + A1c' P A1c < 0.
%  Solve the linear equation
%  Ac'P + P Ac + A1c' P A1c = -I  (a stable Lyapunov/Stein-type system,
%  solved by vectorization) and verify positivity.
% vectorized solve of the adjoint Lyapunov equation with rhs = -I,
% using vec(AXB) = kron(B.',A) vec(X):
%   F*(P) = Ac' P + P Ac + A1c' P A1c
Ladj = kron(eye(d),Ac') + kron(Ac.',eye(d)) + kron(A1c.',A1c');
p = -Ladj \ reshape(eye(d),[],1);
P = reshape(p, d, d);
P = 0.5*(P + P');
end
