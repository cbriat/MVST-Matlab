classdef tgameaverage < matlab.unittest.TestCase
%TGAMEAVERAGE  Unit tests for MVST_GAME_AVERAGE: infinite-horizon
%   average (ergodic) zero-sum game value (Chapter "Game", th:game:AC:finite /
%   th:game:AC:suff, cor:game:AC:equiv).  The average value is the persistent
%   forcing paired against the HOMOGENEOUS game dual P_inf (eq:game:AC:rhostar).
%   Requires only ODE integration.
%
%   Validation anchors:
%     * homogeneous consistency: P_inf solves the game GARE that MVST_GAME_GAMMA
%       settles (residual ~ 0), the settle has converged, and the critical level
%       gamma* agrees with MVST_GAME_GAMMA;
%     * periodic pairing exactness: for periodic forcing, the per-period pairing
%       (1/Tp)int<P_inf,W>dt equals the long-horizon time-average (the limsup of
%       eq:game:AC:rhostar);
%     * rho is linear (hence monotone) in the forcing intensity, and monotone
%       decreasing in gamma (weaker disturbance => smaller value).

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath')); root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
            warning('off','MATLAB:ode45:IntegrationTolNotMet');
        end
    end

    methods (Test)

        function homogeneousConsistency(tc)
            base=struct('A',[-1 0.5;0 -2],'Bu',[1;0.5],'Bw',[0.3;0.6],'Q',eye(2),'Ru',1);
            gg=mvst_game_gamma(base); g=1.3*gg.gammastar;
            r  = mvst_game_average(setfield(base,'gamma',g), struct()); %#ok<*SFLD>
            r80= mvst_game_average(setfield(base,'gamma',g), struct('Tlong',80));
            tc.verifyLessThan(r.residual, 1e-6);            % solves the game GARE
            tc.verifyLessThan(norm(r.Pinf-r80.Pinf), 1e-6); % settle converged
            tc.verifyEqual(r.gammastar, gg.gammastar, 'AbsTol', 1e-9);
        end

        function constantForcingLinear(tc)
            base=struct('A',[-1 0.5;0 -2],'Bu',[1;0.5],'Bw',[0.3;0.6],'Q',eye(2),'Ru',1);
            g=1.3*mvst_game_gamma(base).gammastar; W=[1 0.2;0.2 0.5];
            r1=mvst_game_average(setmany(base,'gamma',g,'W',W));
            r2=mvst_game_average(setmany(base,'gamma',g,'W',3*W));
            tc.verifyEqual(r1.rho, trace(r1.Pinf*W), 'AbsTol', 1e-9);   % rho=<Pinf,W>
            tc.verifyEqual(r2.rho/r1.rho, 3, 'AbsTol', 1e-6);           % linear
        end

        function periodicPairingExactness(tc)
            base=struct('A',[-1 0.5;0 -2],'Bu',[1;0.5],'Bw',[0.3;0.6],'Q',eye(2),'Ru',1);
            g=1.3*mvst_game_gamma(base).gammastar; W=[1 0.2;0.2 0.5]; Tp=1.7;
            Wf=@(t) W*(1+0.6*sin(2*pi*t/Tp));
            rp=mvst_game_average(setmany(base,'gamma',g,'W',Wf,'Tp',Tp));
            Pinf=rp.Pinf; Tbig=Tp*400; s=linspace(0,Tbig,200000);
            avg=trapz(s, arrayfun(@(t) trace(Pinf*Wf(t)), s))/Tbig;
            tc.verifyEqual(rp.rho, avg, 'AbsTol', 1e-4);                % per-period == long-run
            tc.verifyEqual(rp.rho, trace(Pinf*W), 'AbsTol', 1e-4);
        end

        function monotoneInGamma(tc)
            base=struct('A',[-1 0.5;0 -2],'Bu',[1;0.5],'Bw',[0.3;0.6],'Q',eye(2),'Ru',1);
            gstar=mvst_game_gamma(base).gammastar; W=[1 0.2;0.2 0.5];
            ra=mvst_game_average(setmany(base,'gamma',1.3*gstar,'W',W));
            rb=mvst_game_average(setmany(base,'gamma',3.0*gstar,'W',W));
            tc.verifyGreaterThanOrEqual(ra.rho, rb.rho - 1e-9);
        end

    end
end

function s=setmany(s,varargin); for i=1:2:numel(varargin); s.(varargin{i})=varargin{i+1}; end; end
