classdef tgamemixed < matlab.unittest.TestCase
%TGAMEMIXED  Unit tests for MVST_GAME_MIXED: infinite-horizon
%   mixed H2/H-infinity design as a two-player nonzero-sum Nash game (Chapter
%   "Game", prop:game:mixed:verify, eq:game:mixed:GRE).  The H-infinity leg
%   (control closed) and H2 leg (disturbance closed) are solved by a
%   Gauss-Seidel iteration.  Requires only ODE integration.
%
%   Validation anchors:
%     * gamma -> infinity recovers the pure H2 solution (MVST_SOLVE_INFINITE),
%       the worst-case disturbance gain vanishing;
%     * the H-infinity-critical level from MVST_GAME_GAMMA is the feasibility
%       boundary: the mixed game solves for gamma above it and blows up below it;
%     * at the equilibrium the closed-loop w->z H-infinity norm is at most gamma.

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath')); root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
            warning('off','MATLAB:ode45:IntegrationTolNotMet');
        end
    end

    methods (Test)

        function gammaInfRecoversH2(tc)
            A=[-1 0.5;0 -2]; Bu=[1;0.5]; Bw=[0.3;0.6]; Q=eye(2); Ru=1;
            base=struct('A',A,'Bu',Bu,'Bw',Bw,'Q',Q,'Ru',Ru,'Su',zeros(2,1));
            H2 = mvst_solve_infinite(mvst.System.ct(A,Bu,'Q',Q,'R',Ru).thin());
            big = mvst_game_mixed(setfield(base,'gamma',1e4)); %#ok<*SFLD>
            tc.verifyLessThan(norm(big.P2-H2.P{1}), 1e-4);
            tc.verifyLessThan(norm(big.Ku-H2.K{1}), 1e-4);
            tc.verifyLessThan(norm(big.Kw), 1e-4);
            tc.verifyTrue(big.feasible);
        end

        function criticalLevelIsFeasibilityBoundary(tc)
            A=[-1 0.5;0 -2]; Bu=[1;0.5]; Bw=[0.3;0.6]; Q=eye(2); Ru=1;
            base=struct('A',A,'Bu',Bu,'Bw',Bw,'Q',Q,'Ru',Ru,'Su',zeros(2,1));
            gstar = mvst_game_gamma(struct('A',A,'Bu',Bu,'Bw',Bw,'Q',Q,'Ru',Ru)).gammastar;
            above = mvst_game_mixed(setfield(base,'gamma',4*gstar));
            below = mvst_game_mixed(setfield(base,'gamma',0.6*gstar));
            tc.verifyTrue(above.feasible);
            tc.verifyLessThan(above.residual, 1e-5);
            tc.verifyFalse(below.feasible);        % below the game value: blow-up
        end

        function hinfConstraintHolds(tc)
            A=[-1 0.5;0 -2]; Bu=[1;0.5]; Bw=[0.3;0.6]; C=eye(2); Du=[0;0]; Ru=1;
            base=struct('A',A,'Bu',Bu,'Bw',Bw,'Q',C.'*C,'Ru',Ru,'Su',zeros(2,1));
            gstar = mvst_game_gamma(struct('A',A,'Bu',Bu,'Bw',Bw,'Q',C.'*C,'Ru',Ru)).gammastar;
            r = mvst_game_mixed(setfield(base,'gamma',5*gstar));
            Acl = A + Bu*r.Ku; Cz = C + Du*r.Ku; g=0;
            for w=logspace(-2,3,3000); g=max(g,norm(Cz*((1i*w*eye(2)-Acl)\Bw))); end
            tc.verifyLessThanOrEqual(g, r.gamma + 1e-2);
            tc.verifyTrue(r.feasible);
        end

    end
end
