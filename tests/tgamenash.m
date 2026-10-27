classdef tgamenash < matlab.unittest.TestCase
%TGAMENASH  Unit tests for MVST_GAME_NASH: infinite-horizon
%   N-player Nash equilibrium via the coupled I-GRE system (Chapter "Game",
%   prop:game:nash:verify, eq:game:nash:GRE), solved by a Gauss-Seidel sweep
%   over per-player best-response CAREs.  Requires only ODE integration.
%
%   Validation anchors:
%     * N=1 reduces to the single-player infinite-horizon solver
%       MVST_SOLVE_INFINITE (identical dual and gain);
%     * a two-player ZERO-SUM instance reproduces MVST_GAME_GAMMA's H-infinity
%       saddle (dual P, control gain K, disturbance gain L), with P2=-P1;
%     * the per-player I-GRE verification residuals vanish at the equilibrium
%       (prop:game:nash:verify).

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath')); root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
        end
    end

    methods (Test)

        function nEquals1MatchesSolveInfinite(tc)
            A=[-1 0.5;0 -2]; Bu=[1;0.5]; Q=eye(2); Ru=1;
            si = mvst_solve_infinite(mvst.System.ct(A,Bu,'Q',Q,'R',Ru).thin());
            r  = mvst_game_nash(struct('A',A,'Bu',{{Bu}},'Q',{{Q}},'R',{{ {Ru} }}));
            tc.verifyLessThan(norm(r.P{1}-si.P{1}), 1e-6);
            tc.verifyLessThan(norm(r.K{1}-si.K{1}), 1e-6);
            tc.verifyLessThan(r.residual, 1e-6);
        end

        function zeroSumReproducesGameGamma(tc)
            A=[-1 0.5;0 -2]; Bu=[1;0.5]; Bw=[0.3;0.6]; Q=eye(2); Ru=1;
            gg = mvst_game_gamma(struct('A',A,'Bu',Bu,'Bw',Bw,'Q',Q,'Ru',Ru));
            g  = 1.3*gg.gammastar;
            [~,Pt]=mvst_game_dre_ct(A,Bu,Bw,Q,Ru,g,[0 40],[],1e-8,{},{});
            Pref=Pt(:,:,1); Kref=-Ru\(Bu.'*Pref); Lref=(1/g^2)*(Bw.'*Pref);
            zs = mvst_game_nash(struct('A',A,'Bu',{{Bu,Bw}},'Q',{{Q,-Q}}, ...
                                       'R',{{ {Ru,-g^2}, {-Ru,g^2} }}));
            tc.verifyLessThan(norm(zs.P{1}-Pref), 1e-6);
            tc.verifyLessThan(norm(zs.K{1}-Kref), 1e-6);   % control saddle gain
            tc.verifyLessThan(norm(zs.K{2}-Lref), 1e-6);   % disturbance saddle gain
            tc.verifyLessThan(norm(zs.P{2}+zs.P{1}), 1e-6);% zero-sum: P2=-P1
            tc.verifyLessThan(zs.residual, 1e-6);
        end

        function nonzeroSumVerificationResiduals(tc)
            A=[-1 0.5;0 -2]; Bu=[1;0.5]; Bw=[0.3;0.6];
            nz = mvst_game_nash(struct('A',A,'Bu',{{Bu,Bw}}, ...
                     'Q',{{eye(2),2*eye(2)}}, 'R',{{ {1,0}, {0,2} }}));
            tc.verifyLessThan(nz.residual, 1e-6);
            tc.verifyEqual(nz.N, 2);
        end

    end
end
