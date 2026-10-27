function [gstar, lo] = mvst_game_gamma_sd_held(A, Bu, Bw, Q, Ru, Tp, Dwc, relgap, Nmax)
%MVST_GAME_GAMMA_SD_HELD  Critical level of the held-input sampled-data game.
%   [GSTAR,LO] = MVST_GAME_GAMMA_SD_HELD(A,BU,BW,Q,RU,TP,DWC) returns the
%   smallest disturbance level GAMMA for which the backward iteration of the
%   periodic held-input HJI-GRE (MVST_GAME_SD_HELD) from a zero terminal value
%   converges with Isaacs' condition in force (Theta > 0 on the flow, Puu > 0
%   at the samples).  Bisection: GSTAR is the feasible upper end and LO the
%   infeasible lower end of the final bracket, (GSTAR-LO) <= RELGAP*GSTAR
%   (default 1e-6).  NMAX is passed to MVST_GAME_SD_HELD (default 5000).
%
%   See also MVST_GAME_SD_HELD, MVST_GAME_GAMMA_CT.
if nargin < 7; Dwc = []; end
if nargin < 8 || isempty(relgap); relgap = 1e-6; end
if nargin < 9 || isempty(Nmax); Nmax = 5000; end
feas = @(g) mvst_game_sd_held(A, Bu, Bw, Q, Ru, g, Tp, Dwc, Nmax);
lo = 0.05; hi = 1.0;
for it = 1:60
    o = feas(hi); if o.ok; break; end
    lo = hi; hi = 1.6*hi;
    if it == 60; error('mvst:game_gamma_sd_held:nohi', 'no feasible level found'); end
end
o = feas(lo);
if o.ok; error('mvst:game_gamma_sd_held:lo', 'lower bracket end is feasible'); end
while hi - lo > relgap*hi
    mid = 0.5*(lo + hi);
    o = feas(mid);
    if o.ok; hi = mid; else; lo = mid; end
end
gstar = hi;
end
