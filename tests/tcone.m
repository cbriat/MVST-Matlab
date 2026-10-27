classdef tcone < matlab.unittest.TestCase
%TCONE  Unit tests for the 'cone' class of MVST_POS_STABILITY (Chapter
%   "Generalization to Proper Cones"): stability of dot x = A x on a polyhedral
%   proper cone Omega = cone(V) via the Metzler generator factorization
%   A V = V H (Vandergraft / Berman-Plemmons) and the copositive-Lyapunov LP
%   exists p>0 with H'p<0.  Requires only LINPROG (Optimization Toolbox).
%
%   Validation anchors:
%     * ORTHANT COLLAPSE: V=I_n reproduces the 'ct' class exactly (same verdict,
%       same certificate p, H=A) -- the lem:pos:metzstab(d) instance;
%     * SIMPLICIAL non-orthant case: A=[-1 -1;0 -2] (non-Metzler) on the shear
%       cone V=[1 1;0 1] gives H=V^{-1}A V=diag(-1,-2), stable;
%     * NON-SIMPLICIAL cone (m>n): a non-Metzler rotation-contraction on the
%       square-pyramid cone in R^3 (4 generators) is certified, A V = V H holds,
%       H is Metzler, and e^{At} keeps Omega invariant (direct sampling);
%     * a cone-NON-invariant instance is correctly refused (coneInvariant=false).

    methods (TestClassSetup)
        function addpaths(~)
            here = fileparts(mfilename('fullpath')); root = fileparts(here);
            addpath(root); addpath(fullfile(root,'lib'));
        end
    end

    methods (Test)

        function orthantCollapseMatchesCt(tc)
            A = [-2 1;0.5 -3];                       % Metzler Hurwitz
            rc = mvst_pos_stability(struct('A',A), struct('class','ct'));
            rk = mvst_pos_stability(struct('A',A,'V',eye(2)), struct('class','cone'));
            tc.verifyEqual(rk.stable, rc.stable);
            tc.verifyLessThan(norm(rk.p-rc.p), 1e-9);      % same certificate
            tc.verifyLessThan(norm(rk.H-A,'fro'), 1e-9);   % H = A
        end

        function simplicialNonOrthant(tc)
            A = [-1 -1;0 -2]; V = [1 1;0 1];          % non-Metzler; shear cone
            r = mvst_pos_stability(struct('A',A,'V',V), struct('class','cone'));
            tc.verifyTrue(r.coneInvariant);
            tc.verifyTrue(r.stable);
            tc.verifyLessThan(norm(r.H-diag([-1 -2]),'fro'), 1e-9);
        end

        function nonSimplicialSquareCone(tc)
            v1=[1;1;1]; v2=[1;-1;1]; v3=[-1;-1;1]; v4=[-1;1;1];
            V = [v1 v2 v3 v4];                        % m=4 > n=3
            P4 = [0 0 0 1;1 0 0 0;0 1 0 0;0 0 1 0];
            H0 = -eye(4) + 0.3*P4;                    % Metzler, cyclic
            A  = V*H0*pinv(V);                        % non-Metzler (rotation block)
            tc.verifyLessThan(min(A(~eye(3))), -1e-3);% A is genuinely non-Metzler
            r = mvst_pos_stability(struct('A',A,'V',V), struct('class','cone'));
            tc.verifyTrue(r.coneInvariant);
            tc.verifyTrue(r.stable);
            tc.verifyLessThan(norm(A*V-V*r.H,'fro'), 1e-8);        % A V = V H
            H = r.H; tc.verifyGreaterThanOrEqual(min(H(~eye(4))), -1e-8); % Metzler
            tc.verifyTrue(all(r.p>0));
            tc.verifyLessThan(max(r.H.'*r.p), 1e-7);              % H'p < 0
            % e^{At} keeps Omega invariant (direct sampling cross-check)
            rng(0); ok = true;
            for trn=1:60
                x0 = V*rand(4,1);
                for t=[0.2 1 3]
                    if ~inCone(V, expm(A*t)*x0); ok=false; end
                end
            end
            tc.verifyTrue(ok);
        end

        function nonInvariantRefused(tc)
            v1=[1;1;1]; v2=[1;-1;1]; v3=[-1;-1;1]; v4=[-1;1;1]; V=[v1 v2 v3 v4];
            Abad = [0 1 0;-1 0 0;0 0 -1];             % rotation: no square-cone factorization
            r = mvst_pos_stability(struct('A',Abad,'V',V), struct('class','cone'));
            tc.verifyFalse(r.coneInvariant);
            tc.verifyFalse(r.stable);
        end

    end
end

% ======================================================================= %
function tf = inCone(V, x)
% x in cone(V) iff exists xi>=0 with V xi = x (LP feasibility).
    m = size(V,2);
    opt = optimoptions('linprog','Display','off');
    [~,~,ef] = linprog(zeros(m,1), [], [], V, x, zeros(m,1), [], opt);
    tf = (ef == 1);
end
