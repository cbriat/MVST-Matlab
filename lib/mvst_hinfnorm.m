function res = mvst_hinfnorm(sys, opts)
%MVST_HINFNORM  H-infinity norm / (mean-square) L2-gain with a uniform result.
%   RES = MVST_HINFNORM(SYS) returns the H-infinity norm of the system SYS (a
%   struct in the MVST_L2GAIN convention) through the class-routed exact
%   backends, as a UNIFORM struct:
%     RES.gamma        the gain;
%     RES.certificate  storage certificate where available, [] otherwise:
%                        struct with .S (storage matrix, V(x)=x'S x, S >= 0, a
%                        cell of mode storages for dtmjls), .gamma (the level at
%                        which S certifies the bounded-real inequality) and
%                        .residual (Riccati residual, ct only);
%     RES.class        the resolved class tag;
%     RES.method       the backend that produced the number (see table);
%     RES.diagnostics  struct: .tol (bisection tolerance), backend info
%                      (sd: .Nper, .increment, .converged), .exact (true when
%                      the value is exact up to .tol, false for shooting
%                      lower bounds).
%
%   ROUTING TABLE (class -> backend, nature of the number)
%     ct (no noise)  MVST_L2GAIN_CT   Hamiltonian bisection, exact; certificate
%     ct (noise)     MVST_L2GAIN_CT_POISSON mean-square shooting (lower bound)
%     dt             MVST_L2GAIN_DT   bilinear map + Hamiltonian, exact; cert.
%     dtmjls         MVST_L2GAIN_DT_MJLS coupled Riccati iteration; cert. S_i
%     mjls           MVST_L2GAIN_CT_MJLS shooting (lower bound)
%     dwell          MVST_L2GAIN_CT_DWELL certified dwell-time bound
%     poisson        MVST_L2GAIN_CT_POISSON shooting (lower bound)
%     sd             MVST_L2GAIN_SD_HYBRID periodic Riccati, adaptive Nper
%                    with convergence check (see that help)
%   Class inference and the SYS fields are those of MVST_L2GAIN (to which OPTS
%   is passed; OPTS may be a struct or an MVST.OPTIONS).
%
%   Example:
%       r = mvst_hinfnorm(struct('A0',[-1 1;0 -2],'Bw',[1;0],'C',[1 0]));
%       r.gamma, r.class, r.certificate.S
%
%   See also MVST_L2GAIN, MVST_PASSIVITY, MVST.SYSTEM/HINFNORM.
if nargin < 2 || isempty(opts); opts = struct(); end
opts = mvst_opts(opts, false);
r = mvst_l2gain(sys, opts);
cls = r.method;
gam = r.gamma;
diag = struct('tol', NaN, 'exact', false);
cert = [];
switch cls
    case 'ct'
        noisy = anyfield(sys, 'Ac') || anyfield(sys, 'Dwc');
        if ~noisy
            diag.tol = 1e-4; diag.exact = true; method = 'mvst_l2gain_ct';
            cert = ct_certificate(sys, gam);
        else
            method = 'mvst_l2gain_ct_poisson'; diag.tol = 1e-3;
        end
    case 'dt'
        diag.tol = 1e-4; diag.exact = true; method = 'mvst_l2gain_dt';
        cert = dt_certificate(sys, gam);
    case 'dtmjls'
        diag.tol = 1e-4; method = 'mvst_l2gain_dt_mjls';
        if isfield(r.info, 'S') && ~isempty(r.info.S)
            cert = struct('S', {r.info.S}, 'gamma', 1.001*gam, 'residual', NaN);
        end
    case 'mjls';    method = 'mvst_l2gain_ct_mjls';    diag.tol = 1e-3;
    case 'dwell';   method = 'mvst_l2gain_ct_dwell';   diag.tol = 1e-3;
    case 'poisson'; method = 'mvst_l2gain_ct_poisson'; diag.tol = 1e-3;
    case 'sd'
        method = 'mvst_l2gain_sd_hybrid'; diag.tol = 5e-4;
        diag.exact = isfield(r.info, 'converged') && isequal(r.info.converged, true);
        if isfield(r.info, 'Nper');      diag.Nper = r.info.Nper;           end
        if isfield(r.info, 'increment'); diag.increment = r.info.increment; end
        if isfield(r.info, 'converged'); diag.converged = r.info.converged; end
    otherwise
        method = ['mvst_l2gain/' cls];
end
res = struct('gamma', gam, 'certificate', cert, 'class', cls, ...
             'method', method, 'diagnostics', diag);
end

% ------------------------------------------------------------------------
function tf = anyfield(s, f)
tf = isfield(s, f) && ~isempty(s.(f)) && any(cellfun(@(x) nnz(x) > 0, num2cell_(s.(f))));
end

function c = num2cell_(x)
if iscell(x); c = x; else; c = num2cell(x); end
end

% ------------------------------------------------------------------------
function cert = ct_certificate(sys, gam)
%CT_CERTIFICATE  Stabilizing solution of the bounded-real Riccati equation.
%   At gamma_c = 1.001*gam the stable invariant subspace of the Hamiltonian
%   gives S with  S A + A'S + (S B + C'D) R^{-1} (B'S + D'C) + C'C = 0,
%   R = gamma_c^2 I - D'D, so that V(x) = x'S x is a storage for the supply
%   gamma_c^2 |w|^2 - |z|^2.  Returns [] if the subspace is not well posed.
cert = [];
try
    A = pick(sys, {'A0','A'}); A = A(:,:,1);
    B = sys.Bw; C = pick(sys, {'C'}); C = C(:,:,1);
    n = size(A,1); nw = size(B,2);
    if isfield(sys,'Dw') && ~isempty(sys.Dw); D = sys.Dw; else; D = zeros(size(C,1), nw); end
    g = 1.001*gam;
    R = g^2*eye(nw) - D.'*D;
    RiDtC = R \ (D.'*C);
    F = A + B*RiDtC;
    Q = C.'*C + C.'*D*RiDtC;
    H = [F, B*(R\B.'); -Q, -F.'];
    [U, T] = schur(H, 'real'); [U, T] = ordschur(U, T, real(ordeig(T)) < 0); %#ok<ASGLU>
    X1 = U(1:n,1:n); X2 = U(n+1:end,1:n);
    S = X2 / X1; S = (S + S.')/2;
    res = norm(S*A + A.'*S + (S*B + C.'*D)*(R\(B.'*S + D.'*C)) + C.'*C, 'fro');
    if min(eig(S)) > -1e-8*(1 + norm(S))
        cert = struct('S', S, 'gamma', g, 'residual', res);
    end
catch
    cert = [];
end
end

% ------------------------------------------------------------------------
function cert = dt_certificate(sys, ~)
%DT_CERTIFICATE  Fixed-point storage of the DT bounded-real recursion at 1.001*gam.
cert = [];
try
    A = pick(sys, {'A0','A'}); A = A(:,:,1);
    B = sys.Bw; C = pick(sys, {'C'}); C = C(:,:,1);
    if isfield(sys,'Dw') && ~isempty(sys.Dw); D = sys.Dw; else; D = zeros(size(C,1), size(B,2)); end
    [g2, S] = mvst_l2gain_dt_mjls({A}, B, {C}, D, 1);
    if ~isempty(S)
        cert = struct('S', S{1}, 'gamma', 1.001*g2, 'residual', NaN);
    end
catch
    cert = [];
end
end

function v = pick(s, names)
v = [];
for k = 1:numel(names)
    if isfield(s, names{k}) && ~isempty(s.(names{k})); v = s.(names{k}); return; end
end
end
