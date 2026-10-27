%% EX_TV_STEER_COMBINED  Consolidated time-varying covariance-steering.
%  Native replacement for the former three-across minipage tiling
%  (fig:steer:tv:mjls): the 1-sigma covariance-steering ellipses under the
%  time-varying drift A(t) = A0 + dA*exp(-t), drawn in ONE figure with
%  readable fonts: (a) and (b) side by side on top, (c) centered below.
%    (a) Markov-jump, (b) sampled-data, (c) switched.
%  Each class is recomputed by the generic steering wrapper MVST_STEER exactly
%  as in the per-class EX_TV_STEER_* scripts.  Figure only: no data macros.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors(); D = ex_data();

dA = [0 0; 1.5 0]; Xf = [1.0 0;0 1.6];

% ---- (a) Markov-jump (time-varying modes) ------------------------------
base = D.mjls; n = base.n; m = base.m; Mn = base.M;
Xa0 = base.modes{1}.X0;
md = cell(1,Mn);
for i = 1:Mn
    b = base.modes{i};
    md{i} = mvst_mode_tv(n, m, 'A', @(t) b.A + dA*exp(-t), 'B', b.B, ...
        'Q', b.Q, 'R', b.R, 'Ac', b.Ac, 'Bc', b.Bc, 'W', b.W, 'X0', b.X0);
end
sysA = mvst_system(md, 'Pi', base.Pi, 'pi0', base.pi0);
ra = mvst_steer(sysA, Xf, 2.0, struct('class','mjls','X0',Xa0));

% ---- (b) sampled-data (time-varying flow) ------------------------------
s = D.sd; Xb0 = [3 0;0 0.5]; s.A = @(t) s.A + dA*exp(-t);
rb = mvst_steer(s, Xf, 3.0, struct('class','sd','X0',Xb0));

% ---- (c) switched (time-varying modes) ---------------------------------
modes = D.sw.modes; Td = D.sw.Tdwell; nc = modes{1}.n; mc = modes{1}.m; Mc = numel(modes);
Xc0 = [3 0;0 0.5];
mdc = cell(1,Mc);
for i = 1:Mc
    b = modes{i};
    mdc{i} = mvst_mode_tv(nc, mc, 'A', @(t) b.A + dA*exp(-t), 'B', b.B, ...
        'Q', b.Q, 'R', b.R, 'Ac', b.Ac, 'Bc', b.Bc, 'W', b.W, 'X0', b.X0);
end
sysC = struct('Tdwell', Td); sysC.modes = mdc;
rc = mvst_steer(sysC, Xf, 4.0, struct('class','switched','X0',Xc0));

% ---- combined figure: 2+1 layout, no legend -----------------------------
% tiledlayout(2,4): (a) spans tiles 1-2, (b) tiles 3-4, (c) tiles 6-7
% (centered in the bottom row).  A two-column span of the 2x4 grid has the
% same width as one quadrant of the former 17 x 16 cm 2x2 grid, so every
% panel keeps the old quadrant geometry.
% The color key (X_1^0 thick gem blue = gc.order(1,:), X_1^f thick gem red =
% gc.order(2,:), thin intermediate ellipses shaded blue -> red) is carried by
% line samples in the book caption, not by an in-figure legend.
% Each panel keeps equal axis scaling (1-sigma ellipses undistorted) and its
% own padded data window (set by ELL_AX).  The three tiles differ slightly
% in size (tick-label decorations), so one common plot-box size is taken as
% the smallest tile width and height; each panel's shorter data direction is
% then widened so that, at equal scaling, its plot box has that common aspect
% (FILL_TILE below): the three panels come out the same size.
fig = figure('Color','w','Units','centimeters','Position',[2 2 17 16]);
tl  = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
ax = gobjects(1,3);
ax(1) = nexttile(tl,1,[1 2]); ell_ax(ax(1), ra.Sigma, Xa0, Xf, gc); title(ax(1),'(a) Markov-jump','Interpreter','latex','FontWeight','normal');
ax(2) = nexttile(tl,3,[1 2]); ell_ax(ax(2), rb.Sigma, Xb0, Xf, gc); title(ax(2),'(b) sampled-data','Interpreter','latex','FontWeight','normal');
ax(3) = nexttile(tl,6,[1 2]); ell_ax(ax(3), rc.Sigma, Xc0, Xf, gc); title(ax(3),'(c) switched','Interpreter','latex','FontWeight','normal');
% Release the equal-aspect constraint to measure the full tile rectangles
% (with DataAspectRatio manual, the reported position is the shrunk box).
set(ax,'DataAspectRatioMode','auto','PlotBoxAspectRatioMode','auto'); drawnow;
pp = cell2mat(arrayfun(@getpixelposition, ax(:), 'UniformOutput', false));
r  = min(pp(:,3))/min(pp(:,4));                      % common aspect (width/height)
for k = 1:3, fill_tile(ax(k), r); end
mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_tv_steer_combined.pdf'));

% ========================================================================
function [h0,hf] = ell_ax(ax, Sig, X0, Xf, gc)
%ELL_AX  Draw the 1-sigma steering ellipses into AX (verbatim from
%  MVST_STEER_ELLIPSE, minus the figure/legend creation).
hold(ax,'on'); box(ax,'on'); grid(ax,'on'); axis(ax,'equal');
th = linspace(0,2*pi,180);
idx = unique(round(linspace(1,size(Sig,3),14)));
EkAll = zeros(2,numel(th),numel(idx));
for q = 1:numel(idx)
    a = (q-1)/max(numel(idx)-1,1);
    col = (1-a)*gc.order(1,:) + a*gc.order(2,:);
    Ek = ellpts(Sig(:,:,idx(q)), th); EkAll(:,:,q) = Ek;
    plot(ax, Ek(1,:), Ek(2,:), '-','LineWidth',0.7,'Color',[col 0.35]);
end
E0 = ellpts(X0, th); Ef = ellpts(Xf, th);
h0 = plot(ax, E0(1,:), E0(2,:), '-','LineWidth',2.6,'Color',gc.order(1,:));
hf = plot(ax, Ef(1,:), Ef(2,:), '-','LineWidth',2.6,'Color',gc.order(2,:));
% ~5 percent padding (of the larger data span, on every side, so the gap is
% the same physical length on both axes): no ellipse stroke touches the box.
P = [E0, Ef, reshape(EkAll,2,[])];
pad = 0.05*max(max(P(1,:))-min(P(1,:)), max(P(2,:))-min(P(2,:)));
xlim(ax, [min(P(1,:))-pad, max(P(1,:))+pad]);
ylim(ax, [min(P(2,:))-pad, max(P(2,:))+pad]);
xlabel(ax,'$x_1$','Interpreter','latex'); ylabel(ax,'$x_2$','Interpreter','latex');
set(ax,'FontSize',12);
end

function fill_tile(ax, r)
%FILL_TILE  Widen the limits of the shorter direction of AX (about their
%  center) so that, at equal axis scaling, the plot box has aspect R
%  (width/height).  The padded data window set by ELL_AX is kept; only empty
%  margin is added.
xl = xlim(ax); yl = ylim(ax); w = diff(xl); h = diff(yl);
if w/h < r
    xl = mean(xl) + [-1 1]*(r*h)/2;
else
    yl = mean(yl) + [-1 1]*(w/r)/2;
end
set(ax,'XLim',xl,'YLim',yl,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[r 1 1]);
end

function P = ellpts(S, th)
S = (S+S.')/2; [V,Dg] = eig(S); d = max(real(diag(Dg)),0);
P = real(V*diag(sqrt(d))*V.') * [cos(th); sin(th)];
end
