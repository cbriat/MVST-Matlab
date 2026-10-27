function res = mvst_game_average(sys, opts)
%MVST_GAME_AVERAGE  Infinite-horizon average (ergodic) zero-sum game value
%   (Chapter "Game", th:game:AC:finite / th:game:AC:suff, cor:game:AC:equiv).
%   A bounded solution P_inf of the HOMOGENEOUS infinite-horizon I-HJI-GRE,
%   together with Isaacs' condition, determines the average game value as the
%   forcing data paired against that solution (eq:game:AC:rhostar):
%       rho = limsup_{T->inf} (1/T)[ int_0^T <P_inf,W(t)> dt + sum_k <P_inf,W_k> ].
%   The homogeneous part is exactly the H-infinity game dual (MVST_GAME_GAMMA);
%   the average value is linear in the persistent forcing W.
%
%   RES = MVST_GAME_AVERAGE(SYS, OPTS).  SYS fields:
%     .A,.Bu,.Bw,.Q,.Ru   the zero-sum game data (as in MVST_GAME_GAMMA);
%     .gamma   disturbance level (must exceed gamma*; default 1.3*gamma*);
%     .W       persistent forcing: a constant covariance matrix, OR a function
%              handle @(t) with .Tp the forcing period (periodic forcing);
%     .Wk,.Tp  optional per-period jump forcing W_k with dwell period Tp
%              (adds (1/Tp)<P_inf,Wk>, th:game:AC:dwell:periodic).
%
%   OPTS: .Tlong (homogeneous DRE settle horizon, default 40), .nq (periodic
%   quadrature nodes, default 400).
%   RES: .rho (average game value), .Pinf (homogeneous game dual), .Ku,.Kw
%   (saddle gains), .gamma, .gammastar, .residual (game GARE residual).
%
%   Validation anchors: the homogeneous part P_inf matches MVST_GAME_GAMMA at
%   gamma*; for periodic forcing the finite-horizon average converges to the
%   exact per-period pairing (1/Tp)int_0^Tp<P_inf,W>dt; rho is linear (hence
%   monotone) in the forcing intensity.
%
%   See also MVST_GAME_GAMMA, MVST_GAME_DRE_CT, MVST_SOLVE_AVERAGE.

    if nargin < 2; opts = struct(); end
    A=mvst_getf(sys,'A',mvst_getf(sys,'A0',[])); Bu=mvst_getf(sys,'Bu',[]);
    Bw=mvst_getf(sys,'Bw',[]); Q=mvst_getf(sys,'Q',[]); Ru=mvst_getf(sys,'Ru',[]);
    n=size(A,1);
    Tlong=mvst_getf(opts,'Tlong',40); nq=mvst_getf(opts,'nq',400);

    gg = mvst_game_gamma(struct('A',A,'Bu',Bu,'Bw',Bw,'Q',Q,'Ru',Ru));
    gstar = gg.gammastar;
    g = mvst_getf(sys,'gamma', mvst_getf(opts,'gamma', 1.3*gstar));

    % homogeneous infinite-horizon game dual: settle the game DRE backward
    [~,Pt] = mvst_game_dre_ct(A,Bu,Bw,Q,Ru,g,[0 Tlong],[],1e-8,{},{});
    Pinf = Pt(:,:,1);
    Ku = -Ru\(Bu.'*Pinf); Kw = (1/g^2)*(Bw.'*Pinf);
    resid = norm(A.'*Pinf+Pinf*A+Q - Pinf*Bu*(Ru\(Bu.'*Pinf)) ...
                 + (1/g^2)*Pinf*Bw*Bw.'*Pinf, 'fro');

    % forcing pairing rho = time-average of <P_inf, W(t)> (+ jump forcing)
    W = mvst_getf(sys,'W', zeros(n));
    Tp = mvst_getf(sys,'Tp', []);
    if isa(W,'function_handle')
        if isempty(Tp); error('mvst:gameavg:Tp','periodic W needs SYS.Tp.'); end
        s = linspace(0,Tp,nq+1); vals = zeros(1,nq+1);
        for k=1:nq+1; vals(k) = trace(Pinf*W(s(k))); end
        rho = trapz(s, vals)/Tp;                       % (1/Tp) int_0^Tp <Pinf,W> dt
    else
        rho = trace(Pinf*W);                           % constant forcing
    end
    Wk = mvst_getf(sys,'Wk', []);
    if ~isempty(Wk) && ~isempty(Tp); rho = rho + trace(Pinf*Wk)/Tp; end

    res = struct('rho',rho,'Pinf',Pinf,'Ku',Ku,'Kw',Kw,'gamma',g, ...
                 'gammastar',gstar,'residual',resid);
end
