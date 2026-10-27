function ex_kyp()
%EX_KYP  Operator-domain KYP illustration (Chapter "Dissipativity", Section
%   sec:diss:kyp): the general (arbitrary-sign) supply on the master CT data,
%   with BOTH the LMI verdict and the frequency-domain verdict reported, and
%   the hybrid COLLAPSE case on the sampled-data master (monodromy DT-KYP).
%
%   Master CT data:  dx=(A x+Bw w)dt, z=C x+Dw w,  A=[0 1;-1 -1] (Hurwitz),
%   Bw=[1;0], C=[1 0], Dw=0.5.
%     * bounded-real supply: the H-infinity level gamma* (least gamma with the
%       KYP LMI feasible), LMI and FDI verdicts agreeing;
%     * a general INDEFINITE supply (passivity, M=[0 C';C Dw+Dw']): the LMI
%       feasibility and the frequency inequality are reported and must agree
%       (Theorem th:diss:kypCT).
%   Hybrid collapse (prop:diss:kyp:hyb:collapse): the periodic impulsive system
%   with reset J=0.9 I and period Tp=Tbar=0.3 collapses to the monodromy DT-KYP
%   (Theorem th:diss:kypDT) with monodromy Jt=J e^{A Tp} and one-period
%   disturbance operator BwT=int_0^Tp e^{A t}Bw dt; its bounded-real level is the
%   hybrid gamma reported here.
%
%   No figure: the numbers are a short display / table (macros in
%   figures/data/kyp.tex).
%
%   See also MVST_KYP, EX_DATA.

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
figdir = mvst_figdir();

A=[0 1;-1 -1]; Bw=[1;0]; C=[1 0]; Dw=0.5;
sys = struct('A',A,'Bw',Bw,'C',C,'Dw',Dw);

% ---- bounded-real level (CT) ----------------------------------------------
br = mvst_kyp(sys, struct('domain','ct','supply','boundedreal'));

% ---- general indefinite supply: passivity (LMI + frequency verdicts) -------
pv = mvst_kyp(sys, struct('domain','ct','supply','passivity'));

% ---- hybrid collapse: monodromy DT-KYP level ------------------------------
Tp = 0.3; J = 0.9*eye(2);
Jt  = J*expm(A*Tp);
BwT = integral(@(t) expm(A*t)*Bw, 0, Tp, 'ArrayValued', true);
mono = mvst_kyp(struct('A',Jt,'Bw',BwT,'C',C,'Dw',Dw), ...
                struct('domain','dt','supply','boundedreal'));
% cross-check: the collapsed hybrid mode reproduces the DT verdict at that level
hyb = mvst_kyp(sys, struct('domain','hybrid','supply','boundedreal', ...
               'gamma',mono.gamma*1.02,'Tp',Tp,'J',J,'collapsed',true));

fprintf('[ex_kyp] CT bounded-real gamma*=%.5f (%s)\n', br.gamma, br.certifies);
fprintf('         passivity: LMI feas=%d, FDI holds=%d (%s)\n', pv.feas, pv.fdi, pv.certifies);
fprintf('         hybrid monodromy DT-KYP gamma=%.5f (%s); collapse verdict=%s\n', ...
        mono.gamma, mono.certifies, hyb.certifies);

% ---- data macros ----------------------------------------------------------
try
    w = mvst_texdata(fullfile(figdir,'data','kyp.tex'));
    w.num('KypGamCt',br.gamma);
    w.raw('KypPassLmi', tf(pv.feas)); w.raw('KypPassFdi', tf(pv.fdi));
    w.num('KypGamMono',mono.gamma);
    w.num('KypTp',Tp);
    w.close();
catch, end
end

function s = tf(b); if b; s='yes'; else; s='no'; end; end
