function mvst_steer_ellipse(Sig, X0, Xf, ttl, fname) %#ok<INUSD>
%MVST_STEER_ELLIPSE  Single-panel covariance-steering figure.
%   Plots the 1-sigma ellipses of the covariance trajectory Sig(:,:,k)
%   (faint, shaded from the initial colour to the target colour) morphing from
%   the initial covariance X0 (bold, gem blue = gc.order(1,:)) to the target
%   Xf (bold, gem red = gc.order(2,:)).  No legend is drawn: the color key is
%   carried by line samples in the book caption.  Saves a white-background
%   vector PDF to FNAME.
%
%   See also MVST_STEER, MVST_STEER_COMPONENTS, MVST_SAVEPDF.
gc = mvst_gemcolors();
fig = figure('Color','w','Units','centimeters','Position',[2 2 11 10]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); axis(ax,'equal');
th = linspace(0,2*pi,180);
idx = unique(round(linspace(1,size(Sig,3),14)));      % more intermediate ellipses
EkAll = zeros(2,numel(th),numel(idx));
for q = 1:numel(idx)
    a = (q-1)/max(numel(idx)-1,1);                     % 0 (initial) -> 1 (target)
    col = (1-a)*gc.order(1,:) + a*gc.order(2,:);
    Ek = ellpts(Sig(:,:,idx(q)), th); EkAll(:,:,q) = Ek;
    plot(ax, Ek(1,:), Ek(2,:), '-','LineWidth',0.7,'Color',[col 0.35]);
end
E0 = ellpts(X0, th); Ef = ellpts(Xf, th);
plot(ax, E0(1,:), E0(2,:), '-','LineWidth',2.6,'Color',gc.order(1,:));
plot(ax, Ef(1,:), Ef(2,:), '-','LineWidth',2.6,'Color',gc.order(2,:));
% ~5 percent padding (of the larger data span, on every side, so the gap is
% the same physical length on both axes): no ellipse stroke touches the box.
P = [E0, Ef, reshape(EkAll,2,[])];
pad = 0.05*max(max(P(1,:))-min(P(1,:)), max(P(2,:))-min(P(2,:)));
xlim(ax, [min(P(1,:))-pad, max(P(1,:))+pad]);
ylim(ax, [min(P(2,:))-pad, max(P(2,:))+pad]);
xlabel(ax,'$x_1$','Interpreter','latex'); ylabel(ax,'$x_2$','Interpreter','latex');
set(ax,'FontSize',11);
mvst_savepdf(fig, fname);
end

function P = ellpts(S, th)
%ELLPTS  1-sigma ellipse points of a (symmetrised, PSD-projected) covariance.
%   Projecting onto the PSD cone avoids degenerate complex square roots that
%   would otherwise collapse a near-singular covariance to a straight line.
S = (S+S.')/2; [V,D] = eig(S); d = max(real(diag(D)),0);
P = real(V*diag(sqrt(d))*V.') * [cos(th); sin(th)];
end
