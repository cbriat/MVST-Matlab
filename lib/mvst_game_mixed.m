function res = mvst_game_mixed(sys, opts)
%MVST_GAME_MIXED  Infinite-horizon mixed H2/H-infinity design as a two-player
%   nonzero-sum Nash game (Chapter "Game", prop:game:mixed:verify,
%   eq:game:mixed:GRE).  Player 1 (designer) carries the H2 control cost;
%   player 2 (fictitious worst-case disturbance) carries the bounded-real cost
%   at level gamma.  The two legs (eq:game:mixed:GRE:hinf / :h2) are solved by a
%   Gauss-Seidel iteration: the H-infinity leg (control closed at Ku) yields the
%   worst-case disturbance gain Kw and the bounded-real dual P_inf; the H2 leg
%   (disturbance closed at Kw) yields the optimal control Ku and the H2 dual P2.
%
%   RES = MVST_GAME_MIXED(SYS, OPTS).  SYS fields:
%     .A     drift (n x n); .Ac (n x n x L) state Brownian channels (optional);
%     .Bu    control input; .Bw disturbance input;
%     .C     performance output (z = C x + Du u); OR give .Q=C'C directly;
%     .Du    control feedthrough (Ru=Du'Du, Su=C'Du); OR give .Ru (,.Su);
%     .gamma H-infinity level (required); .X0 initial 2nd moment (default I).
%
%   Two legs (min-form, eq:game:mixed:GRE):
%     H-infinity (control at Ku):  bounded-real dual for A+Bu Ku, disturbance
%       Bw, output energy C'C, maximizing weight -gamma^2 I;  Kw=(1/gamma^2)Bw'P_inf.
%     H2 (disturbance at Kw):  optimal-control I-GRE for A+Bw Kw, control Bu,
%       cost (Q=C'C, Ru, Su);  Ku=-Ru^{-1}(Bu'P2+Su').
%
%   OPTS: .maxit (200), .tol (1e-9), CARE options forwarded.
%   RES: .Pinf, .P2, .Ku, .Kw, .gamma, .residual (max leg I-GRE residual),
%        .iters, .h2cost (<P2,X0>), .N=2.
%
%   Validation anchors: gamma->inf recovers the pure H2 solution
%   (MVST_SOLVE_INFINITE, Kw->0); the H-infinity-critical level from
%   MVST_GAME_GAMMA is the feasibility boundary of the mixed game.
%
%   See also MVST_GAME_NASH, MVST_GAME_GAMMA, MVST_SOLVE_INFINITE, MVST_CARE_CT.

    if nargin < 2; opts = struct(); end
    A  = mvst_getf(sys,'A', mvst_getf(sys,'A0',[])); n = size(A,1);
    Ac = mvst_getf(sys,'Ac', zeros(n,n,0));
    Bu = mvst_getf(sys,'Bu',[]); Bw = mvst_getf(sys,'Bw',[]);
    C  = mvst_getf(sys,'C',[]);
    Q  = mvst_getf(sys,'Q', ternary(~isempty(C), @() C.'*C, @() []));
    Du = mvst_getf(sys,'Du',[]);
    Ru = mvst_getf(sys,'Ru', ternary(~isempty(Du), @() Du.'*Du, @() []));
    Su = mvst_getf(sys,'Su', ternary(~isempty(C)&&~isempty(Du), @() C.'*Du, @() zeros(n,size(Bu,2))));
    g  = mvst_getf(sys,'gamma', mvst_getf(opts,'gamma',[]));
    X0 = mvst_getf(sys,'X0', eye(n));
    if isempty(g); error('mvst:mixed:gamma','need SYS.gamma.'); end
    if isempty(Q)||isempty(Ru); error('mvst:mixed:cost','need C/Q and Du/Ru.'); end
    maxit = mvst_getf(opts,'maxit',200); tol = mvst_getf(opts,'tol',1e-9);
    hasAc = ~isempty(Ac) && size(Ac,3)>0;

    careopts = opts; if ~isfield(careopts,'Tmax'); careopts.Tmax = 80; end
    Ku = zeros(size(Bu,2), n); Kw = zeros(size(Bw,2), n);
    Pinf = zeros(n); P2 = zeros(n); it = 0; delta = inf; feasible = true;
    while it < maxit && delta > tol
        it = it + 1; P2old = P2;
        % --- H2 leg (disturbance closed at Kw): optimal control with Q,Ru,Su
        %     (done first so the control is warm-started at the H2 design before
        %     the bounded-real leg evaluates the closed-loop disturbance gain)
        mdH2 = mkmode(A + Bw*Kw, Bu, Q, Ru, Su, Ac, hasAc);
        [P2, Ku] = mvst_care_ct(mdH2, careopts);
        if any(~isfinite(P2(:))) || norm(P2,'fro') > 1e8
            feasible = false; break;
        end
        % --- H-infinity leg (control closed at Ku): min -||z||^2+gamma^2||w||^2
        mdH = mkmode(A + Bu*Ku, Bw, -Q, (g^2)*eye(size(Bw,2)), [], Ac, hasAc);
        [PhinfMin, KwMin] = mvst_care_ct(mdH, careopts);
        Pinf = -PhinfMin;  Kw = KwMin;      % bounded-real dual, worst-case gain
        if any(~isfinite(Pinf(:))) || norm(Pinf,'fro') > 1e8
            feasible = false; break;        % below the Hinf-critical level: blow-up
        end
        delta = norm(P2 - P2old,'fro');
    end
    if ~feasible
        res = struct('Pinf',Pinf,'P2',P2,'Ku',Ku,'Kw',Kw,'gamma',g, ...
                     'residual',inf,'iters',it,'h2cost',inf,'N',2,'feasible',false);
        return;
    end

    % verification residuals (each leg's I-GRE at the equilibrium)
    mdH  = mkmode(A + Bu*Ku, Bw, -Q, (g^2)*eye(size(Bw,2)), [], Ac, hasAc);
    mdH2 = mkmode(A + Bw*Kw, Bu, Q, Ru, Su, Ac, hasAc);
    rH  = norm(mvst_ric_ct(-Pinf, mdH),'fro');
    rH2 = norm(mvst_ric_ct(P2, mdH2),'fro');
    res = struct('Pinf',Pinf,'P2',P2,'Ku',Ku,'Kw',Kw,'gamma',g, ...
                 'residual',max(rH,rH2),'iters',it,'h2cost',trace(P2*X0),'N',2, ...
                 'feasible',true);
end

% ----------------------------------------------------------------------- %
function md = mkmode(A,B,Q,R,S,Ac,hasAc)
    Q = (Q+Q.')/2;
    args = {}; if ~isempty(S); args=[args,{'S',S}]; end
    if hasAc; args=[args,{'Ac',Ac}]; end
    md = mvst_mode(A, B, Q, R, args{:});
end

function y = ternary(c,f,g)
    if c; y=f(); else; y=g(); end
end
