%% EX_GAME_GAMMA_ITER  Worked gamma-iteration (bisection) for the H-infinity game.
%  The critical disturbance level gamma* is the smallest gamma for which the
%  game (HJI) Riccati admits a bounded stabilizing solution.  Feasibility at a
%  trial gamma is decided by integrating the game differential Riccati backward
%  over a long horizon and testing boundedness (the finite-escape / imaginary-
%  axis test).  This script records the bisection bracket [lo,hi] at each step
%  and plots its geometric convergence to gamma*.
% Each boundedness test uses the LONG-HORIZON criterion of the game library
% (horizon 1000, finite solution with ||P|| < 1e6, and the uniform Isaacs
% margin lambda_min(Theta(P(t))) >= 1e-6 gamma^2 on the whole horizon,
% Theta = gamma^2 I - sum_l Dwc{l}'P Dwc{l}).
% The bracket tolerance of the recorded iteration is 1e-4.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();

% ----- game instance: mode 1 of the Markov-jump example (same w-actuated
%       diffusion channel: state-multiplicative C_1 = 0.2 I_2 and shared
%       disturbance-multiplicative column D_{w,1} = 0.1[1;0]) -----
A = [0 1;-2 -1];  Bu = [0;1];  Bw = [1;0];  Q = eye(2);  Ru = 1;
Cc = {0.2*eye(2)};  Dwc = {0.1*[1;0]};
Tlong = 1000;  cap = 1e6;  threl = 1e-6;
feas = @(g) gamefeas(A,Bu,Bw,Q,Ru,g,Tlong,cap,Cc,Dwc,threl);

% ----- bracket: grow hi until feasible, keep lo infeasible -----
lo = 0.05;  hi = 1.0;
while ~feas(hi);  hi = hi*1.6;  end

% ----- bisection, recording the bracket history -----
loh = lo;  hih = hi;                          % history (step 0)
for it = 1:40
    mid = 0.5*(lo+hi);
    if feas(mid);  hi = mid;  else;  lo = mid;  end
    loh(end+1) = lo;  hih(end+1) = hi;        %#ok<SAGROW>
    if (hi-lo) < 1e-4;  break;  end
end
gstar = hi;  niter = numel(hih)-1;
fprintf('gamma-iteration: gstar = %.5f in %d bisection steps\n', gstar, niter);

% ----- export -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','game_gamma.tex'));
% GGamIterN is a step count, written with 0 decimals.
w.num('GGamIterStar', gstar);  w.num('GGamIterN', niter, 0);  w.close();

% ----- figure: bracket [lo,hi] converging to gamma* -----
fig = figure('Color','w','Units','centimeters','Position',[2 2 13 8]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
it = 0:niter;
h1 = plot(ax, it, hih, '-o','LineWidth',1.8,'Color',gc.order(1,:),'MarkerSize',4);
h2 = plot(ax, it, loh, '-s','LineWidth',1.8,'Color',gc.order(2,:),'MarkerSize',4);
yline(ax, gstar, '--','Color',gc.onyx,'LineWidth',1.0,'HandleVisibility','off');
xlabel(ax,'bisection step $k$','Interpreter','latex');
ylabel(ax,'disturbance level $\gamma$','Interpreter','latex');
lg = legend(ax,[h1 h2], ...
      {'feasible bound $\gamma_{\mathrm{hi}}$','infeasible bound $\gamma_{\mathrm{lo}}$'}, ...
      'Interpreter','latex','Orientation','horizontal','Location','southoutside');
set(lg,'Box','off');
text(ax, niter*0.5, gstar, sprintf('$\\gamma^\\star\\approx%.3f$',gstar), ...
     'Interpreter','latex','VerticalAlignment','bottom','FontSize',11);
set(ax,'FontSize',11);
mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_game_gamma_iter.pdf'));

% ========================================================================
function ok = gamefeas(A,Bu,Bw,Q,Ru,g,Tlong,cap,Cc,Dwc,threl)
%GAMEFEAS  Feasibility of gamma: boundedness of the game DRE over [0,TLONG]
%   AND the uniform Isaacs margin lambda_min(Theta(P(t))) >= THREL*g^2 at every
%   output time, folding the state-multiplicative CC and
%   disturbance-multiplicative DWC Brownian channels into the mean-square game
%   Riccati (same criterion as MVST_GAME_GAMMA_CT with RTOL supplied).
ok = false;
try
    [t,Pt] = mvst_game_dre_ct(A,Bu,Bw,Q,Ru,g,[0 Tlong],[],1e-7,Cc,Dwc);
catch
    return
end
if t(1) > 1e-6; return; end                 % integration did not reach t = 0
m = size(Bw,2);
for k = 1:size(Pt,3)
    P = Pt(:,:,k);
    if ~all(isfinite(P(:))) || norm(P) >= cap; return; end
    Th = g^2*eye(m);
    for l = 1:numel(Dwc); Th = Th - Dwc{l}.'*P*Dwc{l}; end
    if min(eig((Th+Th.')/2)) < threl*g^2; return; end
end
ok = true;
end
