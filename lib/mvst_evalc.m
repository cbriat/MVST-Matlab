function M = mvst_evalc(c, t)
%MVST_EVALC  Evaluate a (possibly time-varying) coefficient at time t.
%   M = MVST_EVALC(C,T) returns C(T) if C is a function handle, and C
%   otherwise.  This lets the time-varying solvers treat constant and
%   time-varying mode data uniformly: a mode field may be a constant matrix
%   or a handle @(t) returning a matrix of the appropriate size.
%
%   INPUTS
%     C   constant matrix, or a function handle @(t) returning a matrix;
%     T   time at which to evaluate C when it is a handle.
%   OUTPUT
%     M   C(T) if C is a function handle, otherwise C unchanged.
%
%   Example:
%       A = mvst_evalc(@(t) [0 1; -1 -t], 2);   % [0 1; -1 -2]
%
%   See also MVST_MODE_TV, MVST_RIC_CT_TV.
if isa(c,'function_handle'); M = c(t); else; M = c; end
end
