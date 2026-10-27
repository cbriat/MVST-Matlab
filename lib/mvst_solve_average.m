function res = mvst_solve_average(sys, opts)
%MVST_SOLVE_AVERAGE  Long-term average-cost optimal control (CT).
%   RES = MVST_SOLVE_AVERAGE(SYS) solves the infinite-horizon I-GRE and
%   then returns the optimal long-term average cost rate
%       rho* = sum_i pi_stat(i)*<P_i,W_i> = sum_i pi_stat(i)*trace(P_i W_i),
%   where W_i = SYS.modes{i}.W is the additive process-noise covariance of
%   mode i and pi_stat is the STATIONARY distribution of the Markov chain
%   (pi_stat*Pi = 0, sum = 1).  The long-run average cost is independent of
%   the initial distribution pi0, so the stationary distribution -- not pi0 --
%   is the correct weighting; this is the value to which the closed-loop
%   running-average cost of MVST_SIMULATE_MJLS converges.  RES extends the
%   output of MVST_SOLVE_INFINITE with fields .rho and .piss.
%
%   See also MVST_SOLVE_INFINITE, MVST_SIMULATE_MJLS, MVST_CARE_CT.
if nargin < 2; opts = struct(); end
if ~isstruct(sys) || ~isfield(sys, 'modes')
    error('mvst:solve_average:sys', ...
        'SYS must be a thin system struct with a .modes cell (see MVST_SYSTEM).');
end
res = mvst_solve_infinite(sys, opts);
M   = numel(sys.modes);
if M == 1
    piss = 1;                                   % single mode: trivial stationary law
else
    Pi  = mvst_getf(sys,'Pi',zeros(M));
    Aeq = [Pi.'; ones(1,M)];  beq = [zeros(M,1); 1];   % pi*Pi = 0, sum(pi) = 1
    piss = (Aeq\beq).';
end
rho = 0;
for i = 1:M
    rho = rho + piss(i)*trace(res.P{i}*sys.modes{i}.W);
end
res.rho = rho;  res.piss = piss;
end
