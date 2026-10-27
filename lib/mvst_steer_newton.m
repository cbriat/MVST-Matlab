function [M, hist, iters] = mvst_steer_newton(gapfn, M0, n, opts)
%MVST_STEER_NEWTON  Damped Newton for the steering multiplier.
%   [M,HIST,ITERS] = MVST_STEER_NEWTON(GAPFN,M0,N,OPTS) solves GAPFN(vech(M))=0
%   for the symmetric N-by-N multiplier M, where GAPFN(mv) returns the
%   vectorized terminal gap vech(Sigma(T;M)-X1f).  Uses a numerical Jacobian
%   of M -> Sigma(T;M) and a backtracking line search on the gap norm, giving
%   the quadratically convergent forward-backward variant of the chapter.
%   HIST is the gap-norm history.  This is the shared core used by the
%   continuous-time, Markov-jump, sampled-data, and switched steering solvers.
%
%   See also MVST_STEER, MVST_VECH, MVST_UNVECH.
if ~isa(gapfn, 'function_handle')
    error('mvst:steer_newton:gapfn', 'GAPFN must be a function handle.');
end
if size(M0,1) ~= n || size(M0,2) ~= n
    error('mvst:steer_newton:dim', 'M0 must be %d x %d.', n, n);
end
if nargin < 4; opts = struct(); end
tol   = mvst_getf(opts,'tol',1e-7);
maxit = mvst_getf(opts,'maxit',60);
p = n*(n+1)/2;
mv = mvst_vech(M0); hist = zeros(1,0); iters = 0;
while iters < maxit
    iters = iters + 1;
    r = gapfn(mv); nr = norm(r); hist(end+1) = nr;        %#ok<AGROW>
    if nr < tol; break; end
    h = 1e-5*max(1, norm(mv)); J = zeros(p);
    for j = 1:p
        e = zeros(p,1); e(j) = h; J(:,j) = (gapfn(mv+e) - r)/h;
    end
    dm = -(J \ r);
    if any(~isfinite(dm)); dm = -r; end                  % safeguard
    step = 1; accepted = false;
    for ls = 1:30
        rn = gapfn(mv + step*dm);
        if all(isfinite(rn)) && norm(rn) < nr; accepted = true; break; end
        step = step/2;
    end
    if ~accepted; break; end
    mv = mv + step*dm;
end
M = mvst_unvech(mv, n);
end
