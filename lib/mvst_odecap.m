function [val, terminate, direction] = mvst_odecap(y, cap)
%MVST_ODECAP  Terminal ODE event: stop when the state norm reaches CAP.
%   Shared event function for the backward Riccati integrators: a finite
%   escape (infeasible gamma, indefinite terminal condition from a Newton
%   trial step) then ends the integration cleanly instead of letting the
%   solver grind its step size down and emit tolerance warnings.  Callers
%   detect the truncated time span and treat the run as infeasible.
%
%   INPUTS
%     Y     current ODE state vector;
%     CAP   norm threshold at which to terminate the integration.
%   OUTPUTS (ODEEVENT signature)
%     VAL         event value CAP - norm(Y) (zero crossing triggers the stop);
%     TERMINATE   1 (stop integration on the event);
%     DIRECTION   -1 (detect the decreasing crossing).
%
%   Usage:
%       opts = odeset('Events', @(t,y) mvst_odecap(y, 1e6));
%
%   See also ODESET, ODE45.
val = cap - norm(y);
terminate = 1;
direction = -1;
end
