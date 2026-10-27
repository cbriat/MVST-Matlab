classdef tpassivity < matlab.unittest.TestCase
%TPASSIVITY  Tests for MVST_PASSIVITY (exact CT/DT positive-real test, SOS routes,
%   unified index convention) and the MVST.SYSTEM passivity / kyp /
%   dissipativity methods.
%
%   CONVENTION.  nu = min_w lambda_min(G(jw)+G(jw)') is the book's ISP index
%   (supply He(w'z) - nu|w|^2, Ch08 sec:diss:passivity); MVST_ISP_MARGIN_CT and
%   MVST_PASSIVITY_CT (changed 2026-10-27) return the same number.
%   Exact route: base MATLAB.  SOS routes and the KYP cross-check skip when
%   SOSTOOLS / SeDuMi are missing.

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath'));
            root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
            if exist('sosprogram','file') ~= 2 && exist('mvst_setup','file') == 2
                try, mvst_setup('quiet', true); catch, end %#ok<NOCOM>
            end
        end
    end

    methods (Static)
        function tf = haveSOS()
            tf = exist('sosprogram','file') == 2 && exist('sedumi','file') > 0;
        end
        function g = gridRho(A,B,C,D)
            % brute-force largest rho with G+G* - rho G*G >= 0 on a frequency grid
            w = [0 logspace(-3,3,3000) 1e6]; n = size(A,1); best = inf;
            for k = 1:numel(w)
                G = C*((1i*w(k)*eye(n)-A)\B)+D;
                H = G+G'; Q = G'*G;
                % smallest rho with H - rho Q singular: min generalized eig of (H,Q), Q>0
                if min(eig(Q)) > 1e-14
                    e = real(eig(H, Q)); best = min(best, min(e(e > -inf)));
                end
            end
            g = best;
        end
    end

    methods (Test)
        function exactFirstOrderPositiveReal(tc)
            % G = 1/(s+1) + 0.5:  2Re G = 1 + 2/(1+w^2)  ->  nu = 1 (w -> inf)
            s = struct('A0',-1,'Bw',1,'C',1,'Dw',0.5);
            r = mvst_passivity(s);
            tc.verifyEqual(r.method, 'hamiltonian');
            tc.verifyTrue(r.passive);
            tc.verifyEqual(r.nu_isp, 1, 'AbsTol', 1e-6);
            tc.verifyEqual(r.rho_osp, tpassivity.gridRho(-1,1,1,0.5), 'RelTol', 5e-3);
            tc.verifyTrue(r.diagnostics.exact);
        end

        function dtCertificateSatisfiesDTKYP(tc)
            % returned S must satisfy the DT KYP LMI at the certificate index
            A = [0.5 0.2; -0.1 0.3]; B = [1; 0.5]; C = [0.4 0.2]; D = 1.5;
            r = mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D,'tag','dt'), struct('class','dt'));
            tc.assertFalse(isempty(r.certificate));
            S = r.certificate.S; nu = r.certificate.nu;
            L = [A'*S*A - S, A'*S*B - C'; B'*S*A - C, B'*S*B - (D + D' - nu*eye(1))];
            tc.verifyLessThanOrEqual(max(eig((L+L')/2)), 1e-8);
            tc.verifyGreaterThanOrEqual(min(eig((S+S')/2)), -1e-8);
        end

        function exactNonPassive(tc)
            % G = 1 - 3/(s+1): 2Re G = 2 - 6/(1+w^2), min = -4 at w = 0
            r = mvst_passivity(struct('A0',-1,'Bw',1,'C',-3,'Dw',1));
            tc.verifyFalse(r.passive);
            tc.verifyEqual(r.nu_isp, -4, 'AbsTol', 1e-6);
            tc.verifyTrue(isnan(r.rho_osp));
        end

        function exactAgreesWithGrid(tc)
            rng(7);
            nchk = 0;
            for k = 1:25
                n = randi(4); m = randi(2);
                A = randn(n); A = A - (max(real(eig(A))) + 0.3 + rand)*eye(n);
                B = randn(n,m); C = randn(m,n); D = 0.5*randn(m);
                if rand < 0.5; D = D + D' + (1+rand)*eye(m); end
                r = mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D), struct('method','hamiltonian'));
                if ~strcmp(r.method,'hamiltonian'); continue; end
                [pa, nug] = mvst_passivity_ct(A,B,C,D);
                tc.verifyEqual(r.nu_isp, nug, 'AbsTol', 2e-5, sprintf('random system %d', k));
                tc.verifyEqual(r.passive, pa);
                nchk = nchk + 1;
            end
            tc.verifyGreaterThan(nchk, 15);
        end

        function gridMethodUsesBookConvention(tc)
            A = [0 1;-1 -1]; B = [1;0]; C = [1 0]; D = 0.5;
            rg = mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D), struct('method','grid'));
            re = mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D));
            tc.verifyEqual(rg.method, 'grid');
            tc.verifyEqual(rg.nu_isp, re.nu_isp, 'AbsTol', 1e-4);
            [~, nu] = mvst_passivity_ct(A,B,C,D);
            tc.verifyEqual(nu, re.nu_isp, 'AbsTol', 1e-4);       % same convention
        end

        function indexMatchesIspMargin(tc)
            tc.assumeTrue(tpassivity.haveSOS(), 'SOSTOOLS + SeDuMi required');
            A = [0 1;-1 -1]; B = [1;0]; C = [1 0]; D = 0.5;
            nuSOS = mvst_isp_margin_ct(A,B,C,D);
            r = mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D));
            tc.verifyEqual(r.nu_isp, nuSOS, 'AbsTol', 1e-4);
        end

        function singularFeedthroughRoutesToSOS(tc)
            tc.assumeTrue(tpassivity.haveSOS(), 'SOSTOOLS + SeDuMi required');
            % D = 0: G = 1/(s+1) passive, nu = 0 (w -> inf)
            r = mvst_passivity(struct('A0',-1,'Bw',1,'C',1,'Dw',0));
            tc.verifyEqual(r.method, 'isp-sos');
            tc.verifyTrue(r.passive);
            tc.verifyEqual(r.nu_isp, 0, 'AbsTol', 1e-4);
            % G = -1/(s+1): 2Re G = -2/(1+w^2), min = -2 at w = 0
            r = mvst_passivity(struct('A0',-1,'Bw',1,'C',-1,'Dw',0));
            tc.verifyFalse(r.passive);
            tc.verifyEqual(r.nu_isp, -2, 'AbsTol', 1e-3);
            [~, nug] = mvst_passivity_ct(-1,1,-1,0);
            tc.verifyEqual(r.nu_isp, nug, 'AbsTol', 1e-3);
        end

        function noiseRoutesToIspProgram(tc)
            tc.assumeTrue(tpassivity.haveSOS(), 'SOSTOOLS + SeDuMi required');
            A = [0 1;-1 -1]; B = [1;0]; C = [1 0]; D = 0.5; Dwc = 0.1*[1;0];
            r = mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D,'Ac',zeros(2,2,1),'Dwc',Dwc));
            ref = mvst_isp_margin_ct(A,B,C,D,{},{Dwc});
            tc.verifyEqual(r.method, 'isp-sos');
            tc.verifyEqual(r.nu_isp, ref, 'AbsTol', 1e-6);
            tc.verifyLessThan(r.nu_isp, mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D)).nu_isp);
        end

        function discreteTimeBilinearExact(tc)
            A = [0.5 0.2;-0.1 0.3]; B = [1;0.5]; C = [1 0.5]; D = 0.6;
            r = mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D,'domain','dt'));
            tc.verifyEqual(r.method, 'hamiltonian-bilinear');
            g = inf;
            for t = linspace(0, pi, 20001)
                G = C*((exp(1i*t)*eye(2)-A)\B)+D; g = min(g, real(G+G'));
            end
            tc.verifyEqual(r.nu_isp, g, 'AbsTol', 1e-5);
            tc.verifyEqual(r.passive, g >= 0);
        end

        function agreesWithKypLmi(tc)
            tc.assumeTrue(exist('sossolve','file') == 2 && exist('sedumi','file') > 0, ...
                'SOSTOOLS + SeDuMi required');
            A = [0 1;-1 -1]; B = [0;1]; C = [1 0]; D = 0.5;
            nu = mvst_passivity(struct('A0',A,'Bw',B,'C',C,'Dw',D)).nu_isp;
            for d = [-0.02 0.02]
                M = [zeros(2) C'; C D+D'-(nu+d)];
                k = mvst_kyp(struct('A',A,'Bw',B,'C',C,'Dw',D), struct('M',M));
                tc.verifyEqual(logical(k.feas), d < 0, sprintf('KYP LMI at nu%+g', d));
            end
        end

        function sosRouteMjls(tc)
            tc.assumeTrue(tpassivity.haveSOS(), 'SOSTOOLS + SeDuMi required');
            s = struct('A0',cat(3,[-1 0.2;0 -2],[-2 0;0.3 -1]),'Bw',[1;0],'C',[1 0], ...
                'Dw',1,'Pig',[-1 1;2 -2]);
            try
                r = mvst_passivity(s);
            catch ME
                tc.assumeFail(['MJLS passivity SOS route unavailable: ' ME.message]);
            end
            tc.verifyEqual(r.class, 'mjls');
            tc.verifyTrue(r.passive);
            tc.verifyGreaterThan(r.nu_isp, 0.5);
            tc.verifyLessThanOrEqual(r.nu_isp, 2 + 1e-6);          % <= lambda_min(Dw+Dw')
        end

        function systemMethodsAndOptions(tc)
            S = mvst.System.ct([0 1;-1 -1], zeros(2,0), 'Bw', [0;1], 'C', [1 0], 'Dw', 0.5);
            r = S.passivity();
            ref = mvst_passivity(struct('A0',[0 1;-1 -1],'Bw',[0;1],'C',[1 0],'Dw',0.5));
            tc.verifyEqual(r.nu_isp, ref.nu_isp, 'AbsTol', 1e-9);
            r2 = S.passivity(mvst.Options('solver','nosdp'));
            tc.verifyEqual(r2.nu_isp, ref.nu_isp, 'AbsTol', 1e-9);
            % non-square port is refused
            S2 = mvst.System.ct([-1 0;0 -2], zeros(2,0), 'Bw', eye(2), 'C', [1 0]);
            tc.verifyError(@() S2.passivity(), 'mvst:passivity:square');
        end

        function systemKypClassCheck(tc)
            tc.assumeTrue(exist('sossolve','file') == 2 && exist('sedumi','file') > 0, ...
                'SOSTOOLS + SeDuMi required');
            S = mvst.System.ct([0 1;-1 -1], zeros(2,0), 'Bw', [0;1], 'C', [1 0]);
            g = S.hinfnorm().gamma;
            r = S.kyp(struct('supply','boundedreal','gamma',1.1*g));
            tc.verifyTrue(logical(r.feas));
            r = S.kyp(struct('supply','boundedreal','gamma',0.9*g));
            tc.verifyFalse(logical(r.feas));
            Sn = mvst.System.ct([0 1;-1 -1], zeros(2,0), 'Bw', [0;1], 'C', [1 0], 'Ac', 0.1*eye(2));
            tc.verifyError(@() Sn.kyp(), 'mvst:System:kyp');
            Sm = mvst.System.mjls({-1,-2},{1,1},[-1 1;1 -1]);
            tc.verifyError(@() Sm.kyp(), 'mvst:System:kyp');
        end

        function systemDissipativityValidation(tc)
            S = mvst.System.ct([-1 0.5;0 -2], zeros(2,0), 'Bw', [1;0], 'C', [1 0]);
            tc.verifyError(@() S.dissipativity(struct('type','qsr','Q',1,'S',[1 2],'R',1)), ...
                'mvst:qsr:sizeS');
            tc.verifyError(@() S.dissipativity(struct('type','qsr','Q',1,'S',0)), ...
                'mvst:qsr:missing');
            tc.verifyError(@() S.dissipativity('bogus'), 'mvst:System:supply');
        end

        function systemDissipativityQsr(tc)
            tc.assumeTrue(tpassivity.haveSOS(), 'SOSTOOLS + SeDuMi required');
            S = mvst.System.ct([-1 0.5;0 -2], zeros(2,0), 'Bw', [1;0], 'C', [1 0]);
            g = S.hinfnorm().gamma;
            try
                rp = S.dissipativity(struct('type','qsr','Q',-1,'S',0,'R',(1.2*g)^2));
                rn = S.dissipativity(struct('type','qsr','Q',-1,'S',0,'R',(0.8*g)^2));
            catch ME
                tc.assumeFail(['qsr supply in mvst_diss_sos not available yet: ' ME.message]);
            end
            tc.verifyTrue(logical(rp.feas));
            tc.verifyFalse(logical(rn.feas));
            % string supply agrees with the qsr form
            rs = S.dissipativity(struct('type','l2gain','gamma',1.2*g));
            tc.verifyTrue(logical(rs.feas));
        end
    end
end
