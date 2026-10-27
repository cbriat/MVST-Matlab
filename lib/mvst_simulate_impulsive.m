function res = mvst_simulate_impulsive(intervals, tk, jumpfun, X0, tf, opts)
%MVST_SIMULATE_IMPULSIVE  Closed-loop covariance of an impulsive system.
%   RES = MVST_SIMULATE_IMPULSIVE(INTERVALS,TK,JUMPFUN,X0,TF) simulates the
%   matrix-valued covariance X(t) of an impulsive/switched closed loop on
%   [0,TF], with jump (event) times TK (strictly increasing in (0,TF)).
%   On interval i (between consecutive events, with interval 1 starting at
%   t=0) the flow is
%       dX/dt = Acl X + X Acl' + sum_l Lc(:,:,l) X Lc(:,:,l)' + W,
%   with data INTERVALS{i} = struct('Acl',..,'Lc',..,'W',..,'L',..), where
%   L is the running-cost weight (cost rate tr(L X)).  At each event time
%   TK(k) the state is reset by X <- JUMPFUN(X,k).  INTERVALS must have
%   numel(TK)+1 entries.
%
%   RES has fields .t, .X (d x d x numel(t)), .inst (cost rate), .cost
%   (accumulated), .runavg (running average rate).
%
%   See also MVST_SIMULATE_CT, MVST_SIMULATE_MJLS, MVST_GETF.
if nargin < 6; opts = struct(); end
rt   = mvst_getf(opts,'reltol',1e-8);
npts = mvst_getf(opts,'npts',60);
d    = size(X0,1);
edges = [0, tk(:).', tf];
assert(numel(intervals) == numel(edges)-1, ...
    'Need numel(tk)+1 interval data structs.');

t = []; Xall = zeros(d,d,0); inst = []; gacc = 0; g = [];
X = (X0+X0.')/2;
for i = 1:numel(intervals)
    seg = intervals{i};
    Lc = mvst_getf(seg,'Lc',zeros(d,d,0));
    W  = mvst_getf(seg,'W', zeros(d));
    L  = mvst_getf(seg,'L', zeros(d));
    ti = linspace(edges(i), edges(i+1), npts);
    odef = @(tt,y) segderiv(y, seg.Acl, Lc, W, L, d, size(Lc,3));
    oo = odeset('RelTol',rt,'AbsTol',rt*1e-2);
    [ts, Y] = ode45(odef, ti, [X(:); 0], oo);
    for k = 1:numel(ts)
        Xk = reshape(Y(k,1:d*d), d, d); Xk = (Xk+Xk.')/2;
        Xall(:,:,end+1) = Xk;             %#ok<AGROW>
        inst(end+1,1) = trace(L*Xk);      %#ok<AGROW>
        g(end+1,1) = gacc + Y(k,end);     %#ok<AGROW>
    end
    t = [t; ts];                          %#ok<AGROW>
    gacc = g(end);
    Xend = reshape(Y(end,1:d*d), d, d); X = (Xend+Xend.')/2;
    if i <= numel(tk)
        X = jumpfun(X, i); X = (X+X.')/2; % apply reset at event i
    end
end
runavg = g ./ max(t, eps); if ~isempty(inst); runavg(1) = inst(1); end
% Blank the burn-in: dividing by a near-zero elapsed time makes the first
% samples spike; NaNs are skipped by plot so the curve starts cleanly.
runavg(t - t(1) < 0.02*(t(end) - t(1))) = NaN;
res.t = t; res.X = Xall; res.inst = inst; res.cost = g; res.runavg = runavg;
end

function dy = segderiv(y, Acl, Lc, W, L, d, Nc)
X = reshape(y(1:d*d), d, d); X = (X+X.')/2;
dX = Acl*X + X*Acl.' + W;
for l = 1:Nc; dX = dX + Lc(:,:,l)*X*Lc(:,:,l).'; end
dy = [dX(:); trace(L*X)];
end
