classdef tkyp < matlab.unittest.TestCase
%TKYP  Unit tests for the operator-domain KYP engine MVST_KYP (Chapter
%   "Dissipativity", Section sec:diss:kyp): the general matrix-valued KYP lemma
%   for an arbitrary-sign supply, in continuous time (th:diss:kypCT +
%   lem:diss:kyp:dc), discrete time (th:diss:kypDT + lem:diss:kypDT:dc), and the
%   hybrid periodic setting (th:diss:kyp:hyb:nec/loss, prop:diss:kyp:hyb:collapse,
%   def:diss:kyp:hyb:reach).  The genuinely non-collapsed case is solved directly
%   as the periodic differential-difference LMI (eq:diss:kyp:hyb:lmi) with a
%   timer-varying storage P(tau) via the SOS-timer route; the
%   rank-one hybrid-harmonic sweep is NOT used for certification (it is
%   non-exhaustive, prop:diss:kyp:hyb:nonexh) but is used to CHECK necessity.
%%   Default engine: SOSTOOLS + SeDuMi (skipped when missing; see MVST_SETUP).  The
%   whole suite can be rerun on the other configurations with the environment
%   variable MVST_TKYP_CFG = 'yalmip' (legacy engine, needs YALMIP) or 'nosdp'
%   (SOSNOSDP backend); enginesAgree and nosdpAgrees compare them directly.
%
%   Validation anchors:
%     * bounded-real supply reproduces MVST_L2GAIN's H-infinity norm (CT and DT);
%     * passivity supply reproduces MVST_PASSIVITY_CT's verdict;
%     * indefinite supply: LMI feasibility <=> zero/unit-frequency inequality,
%       with a frequency-domain WITNESS when infeasible;
%     * regularity (R)/(Rd): a reachability-violating instance is REFUSED;
%     * hybrid collapse == direct monodromy DT-KYP; an unstable monodromy is
%       refused;
%     * NON-COLLAPSED periodic LMI: the timer-varying certificate SOLVES and
%       dominates the constant-P slice; the constant-P certificate implies
%       monodromy feasibility (prop:collapse direction); a feasible periodic
%       LMI satisfies the necessity inequality Sigma(Y)>=0 on hybrid-harmonic
%       arcs (th:diss:kyp:hyb:nec).

    properties (Constant)
        tol = 1e-3;
    end

    methods (TestClassSetup)
        function addpaths(tc)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
            cfg = getenv('MVST_TKYP_CFG');       % '' | 'yalmip' | 'nosdp' (alternate full-suite configurations)
            if strcmp(cfg,'yalmip')
                tc.assumeTrue(exist('sdpvar','file')==2, 'YALMIP (sdpvar) not on the path');
            end
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS (sossolve) not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods (Test)

        function ctBoundedRealMatchesL2gain(tc)
            A=[-1 0.5;0 -2]; Bw=[1;0.5]; C=[1 1]; Dw=0.2;
            g = mvst_l2gain_ct(A,Bw,C,Dw,{},{});
            r = kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                    struct('domain','ct','supply','boundedreal'));
            tc.verifyTrue(r.regular);
            tc.verifyEqual(r.gamma, g, 'AbsTol', tc.tol);
            tc.verifyEqual(r.certifies, 'both');
        end

        function ctFdiBothWaysWithWitness(tc)
            A=[-1 0.5;0 -2]; Bw=[1;0.5]; C=[1 1]; Dw=0.2;
            g = mvst_l2gain_ct(A,Bw,C,Dw,{},{});
            ra = kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                    struct('domain','ct','supply','boundedreal','gamma',g+0.5));
            rb = kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                    struct('domain','ct','supply','boundedreal','gamma',g-0.3));
            tc.verifyTrue(ra.feas && ra.fdi);            % above norm: both hold
            tc.verifyFalse(rb.feas); tc.verifyFalse(rb.fdi);  % below: both fail
            tc.verifyNotEmpty(rb.witness);               % frequency-domain witness
        end

        function enginesAgree(tc)
            % SOSTOOLS (default) vs legacy YALMIP: same verdicts, gamma, and witness objective.
            tc.assumeTrue(exist('sdpvar','file')==2, 'YALMIP (sdpvar) not on the path');
            S = struct('A',[-1 0.5;0 -2],'Bw',[1;0.5],'C',[1 1],'Dw',0.2);
            D = struct('A',[0.4 0.2;0 0.5],'Bw',[1;0.3],'C',[1 1],'Dw',0.1);
            for dm = {'ct','dt'}
                sy = S; if strcmp(dm{1},'dt'); sy = D; end
                a = mvst_kyp(sy, struct('domain',dm{1},'supply','boundedreal'));
                b = mvst_kyp(sy, struct('domain',dm{1},'supply','boundedreal','engine','yalmip'));
                tc.verifyEqual(a.gamma, b.gamma, 'AbsTol', 1e-4);
                lo = a.gamma - 0.3;
                a = mvst_kyp(sy, struct('domain',dm{1},'supply','boundedreal','gamma',lo));
                b = mvst_kyp(sy, struct('domain',dm{1},'supply','boundedreal','gamma',lo,'engine','yalmip'));
                tc.verifyEqual([a.feas a.fdi], [b.feas b.fdi]);
                tc.verifyEqual(a.info.val, b.info.val, 'AbsTol', 1e-6);
            end
        end

        function nosdpAgrees(tc)
            % SOSNOSDP backend vs SeDuMi: same gamma and verdicts.
            tc.assumeTrue(exist('sossolve_nosdp','file')==2, 'sosnosdp not on the path');
            S = struct('A',[-1 0.5;0 -2],'Bw',[1;0.5],'C',[1 1],'Dw',0.2);
            a = mvst_kyp(S, struct('domain','ct','supply','boundedreal'));
            b = mvst_kyp(S, struct('domain','ct','supply','boundedreal','solver','nosdp'));
            tc.verifyEqual(b.gamma, a.gamma, 'AbsTol', 1e-3);
            lo = a.gamma - 0.3;
            b = mvst_kyp(S, struct('domain','ct','supply','boundedreal','gamma',lo,'solver','nosdp'));
            tc.verifyFalse(b.feas); tc.verifyFalse(b.fdi); tc.verifyNotEmpty(b.witness);
        end

        function ctPassivityMatchesEngine(tc)
            A=[-1 1;-1 -1]; B=eye(2); C=eye(2); D=0.5*eye(2);
            isp = mvst_passivity_ct(A,B,C,D);
            r = kyp(struct('A',A,'Bw',B,'C',C,'Dw',D), ...
                    struct('domain','ct','supply','passivity'));
            tc.verifyEqual(logical(r.feas), logical(isp));
        end

        function dtBoundedRealMatchesFreqSweep(tc)
            A=[0.4 0.2;0 0.5]; Bw=[1;0.3]; C=[1 1]; Dw=0.1;
            r = kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                    struct('domain','dt','supply','boundedreal'));
            th = linspace(0,pi,4000); gmax = 0;
            for w = th; G = C*((exp(1i*w)*eye(2)-A)\Bw)+Dw; gmax = max(gmax,norm(G)); end
            tc.verifyTrue(r.regular);
            tc.verifyEqual(r.gamma, gmax, 'AbsTol', 5e-3);
        end

        function regularityViolationRefused(tc)
            % Bw = 0: the noise channel is unreachable, (R) fails, so the KYP
            % theorem does not apply and the engine refuses to conclude.
            A=[-1 0.5;0 -2];
            r = kyp(struct('A',A,'Bw',zeros(2,1),'C',[1 1],'Dw',0), ...
                    struct('domain','ct','supply','boundedreal','gamma',1));
            tc.verifyFalse(r.regular);
            tc.verifyEqual(r.certifies, 'refused');
            tc.verifyEmpty(r.feas);
        end

        function hybridCollapseMatchesMonodromyDTKYP(tc)
            A=[-1 0.5;0 -2]; Bw=[1;0.5]; C=[1 1]; Dw=0.2; J=0.9*eye(2); Tp=1; gam=3.0;
            Jt = J*expm(A*Tp);
            BwT = integral(@(t) expm(A*t)*Bw, 0, Tp, 'ArrayValued', true);
            rh = kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                    struct('domain','hybrid','Tp',Tp,'J',J,'supply','boundedreal','gamma',gam,'collapsed',true));
            rd = kyp(struct('A',Jt,'Bw',BwT,'C',C,'Dw',Dw), ...
                    struct('domain','dt','supply','boundedreal','gamma',gam));
            tc.verifyTrue(rh.monoSchur);
            tc.verifyEqual(logical(rh.feas), logical(rd.feas));    % collapse == monodromy DT-KYP
        end

        function hybridNonCollapsedSolvesAndDominates(tc)
            % The non-collapsed periodic differential-difference
            % LMI is SOLVED (no longer refused); the timer-varying deg-4
            % certificate dominates the constant-P deg-0 slice (smaller BRL gain).
            A=[0 1;-2 -0.5]; Bw=[0;1]; C=[1 0]; Dw=0; J=0.5*eye(2);
            base=struct('domain','hybrid','Tp',1,'J',J,'supply','boundedreal','collapsed',false);
            r0=kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), setfield(base,'deg',0));
            r4=kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), setfield(base,'deg',4));
            tc.verifyTrue(r4.feas);
            tc.verifyEqual(r4.certifies,'lmi');
            tc.verifyLessThanOrEqual(r4.gamma, r0.gamma + 1e-3);   % dominance
            tc.verifyLessThan(r4.gamma, r0.gamma);                 % strictly non-collapsed
        end

        function hybridConstantPImpliesMonodromy(tc)
            % prop:diss:kyp:hyb:collapse (always-true direction): whenever the
            % constant-P periodic LMI (deg 0) is feasible, the monodromy DT-KYP
            % is feasible at the same gamma.
            A=[0 1;-2 -0.5]; Bw=[0;1]; C=[1 0]; Dw=0; J=0.5*eye(2); Tp=1;
            sys=struct('A',A,'Bw',Bw,'C',C,'Dw',Dw);
            for g=[1.5 2.0 3.0]
                r0=kyp(sys, struct('domain','hybrid','Tp',Tp,'J',J, ...
                        'supply','boundedreal','gamma',g,'collapsed',false,'deg',0));
                rm=kyp(sys, struct('domain','hybrid','Tp',Tp,'J',J, ...
                        'supply','boundedreal','gamma',g,'collapsed',true));
                if ~isempty(r0.feas) && r0.feas
                    tc.verifyTrue(logical(rm.feas));   % constant-P => monodromy
                end
            end
        end

        function hybridNonCollapsedNecessity(tc)
            % th:diss:kyp:hyb:nec: a feasible periodic LMI implies the hybrid-
            % frequency inequality Sigma(Y_{w,theta})>=0 on hybrid-harmonic arcs.
            A=[0 1;-2 -0.5]; Bw=[0;1]; C=[1 0]; Dw=0; J=0.5*eye(2); Tp=1;
            Jt=J*expm(A*Tp);
            r4=kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                    struct('domain','hybrid','Tp',Tp,'J',J,'supply','boundedreal', ...
                           'collapsed',false,'deg',4));
            g=r4.gamma+0.15; M=[-C.'*C -C.'*Dw; -Dw.'*C g^2-Dw.'*Dw];
            worst=inf; ns=300;
            for w=[0 0.5 1 2 5 -1.5]
              for th=linspace(0,2*pi,7)
                gw=intBexpL(A,Bw,Tp,w); xi0=(exp(1i*th)*eye(2)-Jt)\(J*gw);
                S=0; h=Tp/ns;
                for k=0:ns
                  t=k*h; xi=expm(A*t)*xi0+intBexpL(A,Bw,t,w);
                  v=[xi; exp(1i*w*t)]; wt=1-0.5*(k==0||k==ns);
                  S=S+wt*real(v'*M*v)*h;
                end
                worst=min(worst,S);
              end
            end
            tc.verifyGreaterThanOrEqual(worst, -1e-3);
        end

        function hybridUnstableMonodromyRefused(tc)
            % J = 3 I: monodromy Jt = 3 e^{A} has an eigenvalue > 1 (not Schur),
            % so the hybrid KYP hypothesis fails and the engine refuses.
            A=[-1 0.5;0 -2]; Bw=[1;0.5]; C=[1 1]; Dw=0.2;
            r = kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                    struct('domain','hybrid','Tp',1,'J',3*eye(2),'supply','boundedreal','gamma',3,'collapsed',true));
            tc.verifyFalse(r.monoSchur);
            tc.verifyEqual(r.certifies, 'refused');
        end

    end
end

% ======================================================================= %
function r = kyp(sys,o)
% mvst_kyp with the alternate full-suite configuration of MVST_TKYP_CFG merged in.
    cfg = getenv('MVST_TKYP_CFG');
    if strcmp(cfg,'yalmip'); o.engine = 'yalmip'; elseif strcmp(cfg,'nosdp'); o.solver = 'nosdp'; end
    r = mvst_kyp(sys,o);
end

function G = intBexpL(A,B,T,w)
% int_0^T e^{A(T-r)} B e^{j w r} dr  (hybrid-harmonic forcing integral).
    n=size(A,1); ns=200; h=T/max(ns,1); G=zeros(n,1); if T<=0; return; end
    for k=0:ns; r=k*h; wt=1-0.5*(k==0||k==ns); G=G+wt*expm(A*(T-r))*B*exp(1i*w*r)*h; end
end
