classdef tstabdelay < matlab.unittest.TestCase
%TSTABDELAY  Unit tests for MVST_STAB_DELAY: delay-independent
%   stability of matrix-valued and positive time-delay systems (Chapter
%   "Stability" th:stab:delay / th:stab:delayDT, and Chapter "Linear Positive
%   Systems" prop:pos:classncs).  Default engine: SOSTOOLS + SeDuMi (skipped when missing; see MVST_SETUP);
%   the legacy YALMIP engine is compared in the engine tests (skipped without YALMIP).
%
%   Validation anchors:
%     * CT: the delay-independent LMI verdict (eq:stab:delayLMI) equals the
%       zero-delay operator spectral abscissa test (F0+sum Hi on S^n), both ways;
%     * DT: the LMI verdict (eq:stab:delayLMIDT) equals the zero-delay log
%       spectral radius test (J0+sum Ti), both ways;
%     * positive slice: the verdict equals MVST_POS_STABILITY on the zero-delay
%       Metzler matrix A0+A1 (prop:pos:classncs), and the copositive
%       Lyapunov-Krasovskii vector satisfies (A0+A1)'p<0, p>0.

    methods (TestClassSetup)
        function addpaths(tc)
            here = fileparts(mfilename('fullpath')); root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
            tc.assumeTrue(exist('sossolve','file')==2, 'SOSTOOLS (sossolve) not on the path: run mvst_setup');
            tc.assumeTrue(exist('sedumi','file')>0, 'SeDuMi not on the path: run mvst_setup');
        end
    end

    methods (Test)

        function ctLMIequalsSpectral(tc)
            Cs = struct('A0',diag([-2 -3]),'Ad',{{diag([0.5 0.5])}});
            Cu = struct('A0',diag([-1 -1]),'Ad',{{diag([1.5 1.5])}});
            rs = mvst_stab_delay(Cs, struct('domain','ct'));
            ru = mvst_stab_delay(Cu, struct('domain','ct'));
            tc.verifyTrue(rs.stable);  tc.verifyLessThan(rs.margin,0);
            tc.verifyFalse(ru.stable); tc.verifyGreaterThan(ru.margin,0);
            tc.verifyEqual(rs.stable, rs.margin<0);   % LMI <=> spectral
            tc.verifyEqual(ru.stable, ru.margin<0);
        end

        function dtLMIequalsSpectral(tc)
            Ds = struct('Jl',{{0.4*eye(2)}},'Hd',{{0.3*eye(2)}});
            Du = struct('Jl',{{0.7*eye(2)}},'Hd',{{0.9*eye(2)}});
            ds = mvst_stab_delay(Ds, struct('domain','dt'));
            du = mvst_stab_delay(Du, struct('domain','dt'));
            tc.verifyTrue(ds.stable);  tc.verifyLessThan(ds.margin,0);
            tc.verifyFalse(du.stable); tc.verifyGreaterThan(du.margin,0);
            tc.verifyEqual(ds.stable, ds.margin<0);
            tc.verifyEqual(du.stable, du.margin<0);
        end

        function positiveSliceMatchesPosStability(tc)
            A0=[-2 0.5;0.3 -3];
            Ps=struct('A0',A0,'A1',[0.2 0.1;0.1 0.2]);
            Pu=struct('A0',A0,'A1',[2 1;1 2]);
            ps=mvst_stab_delay(Ps,struct('domain','positive'));
            pu=mvst_stab_delay(Pu,struct('domain','positive'));
            lps=mvst_pos_stability(struct('A',A0+Ps.A1));
            lpu=mvst_pos_stability(struct('A',A0+Pu.A1));
            tc.verifyEqual(ps.stable, lps.stable);
            tc.verifyEqual(pu.stable, lpu.stable);
            tc.verifyEqual(ps.stable, all(real(eig(A0+Ps.A1))<0));
            tc.verifyEqual(pu.stable, all(real(eig(A0+Pu.A1))<0));
        end

        function enginesAgree(tc)
            % New default engine ('sostools') vs legacy 'yalmip' on CT and DT cases,
            % and the SOS certificate satisfies the LMIs it claims.
            tc.assumeTrue(exist('sdpvar','file')==2, 'YALMIP (sdpvar) not on the path: run mvst_setup');
            ct = {struct('A0',-1,'Ad',{{1.35}}), struct('A0',-1,'Ad',{{1.45}}), ...
                  struct('A0',diag([-2 -3]),'Ad',{{diag([0.5 0.5])}}), ...
                  struct('A0',[-2 1;0 -3],'Al',{{0.3*eye(2)}},'Ad',{{0.4*eye(2),[0 0.3;0.2 0]}})};
            for k = 1:numel(ct)
                a = mvst_stab_delay(ct{k}, struct('domain','ct','engine','sostools'));
                b = mvst_stab_delay(ct{k}, struct('domain','ct','engine','yalmip'));
                tc.verifyEqual(a.stable, b.stable, sprintf('CT case %d',k));
                if a.stable
                    P = a.P; Q = a.Q; n = size(P,1);
                    F = ct{k}.A0.'*P + P*ct{k}.A0;
                    if isfield(ct{k},'Al'); for l=1:numel(ct{k}.Al); F = F + ct{k}.Al{l}.'*P*ct{k}.Al{l}; end; end
                    for i=1:numel(Q); F = F + Q{i}; end
                    tc.verifyLessThan(max(eig(F)), 1e-7);            % F0*(P)+sum Qi < 0
                    tc.verifyGreaterThan(min(eig(P)), 0);
                end
            end
            dt = {struct('Jl',{{0.6}},'Hd',{{0.79}}), struct('Jl',{{0.6}},'Hd',{{0.81}}), ...
                  struct('Jl',{{0.4*eye(2)}},'Hd',{{0.3*eye(2)}}), struct('Jl',{{0.7*eye(2)}},'Hd',{{0.9*eye(2)}})};
            for k = 1:numel(dt)
                a = mvst_stab_delay(dt{k}, struct('domain','dt','engine','sostools'));
                b = mvst_stab_delay(dt{k}, struct('domain','dt','engine','yalmip'));
                tc.verifyEqual(a.stable, b.stable, sprintf('DT case %d',k));
            end
        end

        function nosdpBackend(tc)
            % The SOS engine also runs on the SDP-free backend and agrees on clear cases.
            tc.assumeTrue(exist('sossolve_nosdp','file')>0, 'SOSNOSDP not on the path: run mvst_setup');
            rs = mvst_stab_delay(struct('A0',diag([-2 -3]),'Ad',{{diag([0.5 0.5])}}), struct('domain','ct','solver','nosdp'));
            ru = mvst_stab_delay(struct('A0',diag([-1 -1]),'Ad',{{diag([1.5 1.5])}}), struct('domain','ct','solver','nosdp'));
            tc.verifyTrue(rs.stable); tc.verifyFalse(ru.stable);
            ds = mvst_stab_delay(struct('Jl',{{0.4*eye(2)}},'Hd',{{0.3*eye(2)}}), struct('domain','dt','solver','nosdp'));
            du = mvst_stab_delay(struct('Jl',{{0.7*eye(2)}},'Hd',{{0.9*eye(2)}}), struct('domain','dt','solver','nosdp'));
            tc.verifyTrue(ds.stable); tc.verifyFalse(du.stable);
        end

        function scalarBlocksBothBackends(tc)
            % n = 1: all Gram blocks are 1x1 (SeDuMi pretransfo needs them merged);
            % SeDuMi and SOSNOSDP agree with the zero-delay spectral test.
            tc.assumeTrue(exist('sossolve_nosdp','file')>0, 'SOSNOSDP not on the path: run mvst_setup');
            for s = {'sedumi','nosdp'}
                o = struct('solver',s{1});
                a = mvst_stab_delay(struct('A0',-1,'Ad',{{0.5}}), setfield(o,'domain','ct'));
                b = mvst_stab_delay(struct('A0',-1,'Ad',{{1.45}}), setfield(o,'domain','ct'));
                c = mvst_stab_delay(struct('Jl',{{0.4}},'Hd',{{0.3}}), setfield(o,'domain','dt'));
                tc.verifyTrue(a.stable, s{1}); tc.verifyFalse(b.stable, s{1}); tc.verifyTrue(c.stable, s{1});
            end
        end

        function positiveLKFvector(tc)
            % Copositive Lyapunov-Krasovskii vector: (A0+A1)'p<0, p>0.
            A0=[-2 0.5;0.3 -3]; A1=[0.2 0.1;0.1 0.2];
            ps=mvst_stab_delay(struct('A0',A0,'A1',A1),struct('domain','positive'));
            tc.verifyTrue(ps.stable);
            if ~isempty(ps.p)
                tc.verifyTrue(all(ps.p>0));
                tc.verifyLessThan(max((A0+A1).'*ps.p), 1e-7);
            end
        end

    end
end
