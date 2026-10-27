%% EX_STEER_COMBINED  Consolidated constant-coefficient covariance-steering.
%  Native replacement for the former four-across 2x2 minipage tiling
%  (fig:steer:mjls): the 1-sigma covariance-steering ellipses for the four
%  example classes drawn in ONE figure with readable fonts.
%    (a) Markov-jump, (b) sampled-data, (c) switched, (d) Poissonian-jump.
%  Each class is recomputed by the generic steering wrapper MVST_STEER exactly
%  as in the per-class EX_STEER_* scripts.  Figure only: writes no data macros.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors(); D = ex_data();

% ---- (a) Markov-jump --------------------------------------------------
Xf = [1.0 0;0 1.6];
Xa0 = D.mjls.modes{1}.X0;
ra = mvst_steer(D.mjls, Xf, 2.0, struct('class','mjls','X0',Xa0));

% ---- (b) sampled-data -------------------------------------------------
Xb0 = [3 0;0 0.5];
rb = mvst_steer(D.sd, Xf, 3.0, struct('class','sd','X0',Xb0,'M0',-eye(2)));

% ---- (c) switched -----------------------------------------------------
Xc0 = [3 0;0 0.5];
rc = mvst_steer(D.sw, Xf, 4.0, struct('class','switched','X0',Xc0));

% ---- (d) Poissonian-jump ----------------------------------------------
A = [0 1; 2 -0.5]; B = [0;1]; C1 = 0.2*eye(2); D1 = [0;0];
J1 = 1.25*eye(2);  lam = 0.5; Q = eye(2); R = 1;
Xd0 = [2 0.3;0.3 1]; Xdf = [0.8 0;0 0.8];
sysd = struct('A',A,'B',B,'Cc',{{C1}},'Dc',{{D1}},'Jc',{{J1}},'lam',lam,'Q',Q,'R',R);
rd = mvst_steer(sysd, Xdf, 2.0, struct('class','poisson','X0',Xd0,'M0',eye(2)));

% ---- combined 2x2 figure, no legend -----------------------------------
% The color key (X_1^0 thick gem blue = gc.order(1,:), X_1^f thick gem red =
% gc.order(2,:), thin intermediate ellipses shaded blue -> red) is carried by
% line samples in the book caption, not by an in-figure legend.
fig = figure('Color','w','Units','centimeters','Position',[2 2 17 16]);
tl  = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
ell_ax(nexttile(tl), ra.Sigma, Xa0, Xf,  gc); title('(a) Markov-jump','Interpreter','latex','FontWeight','normal');
ell_ax(nexttile(tl), rb.Sigma, Xb0, Xf,  gc); title('(b) sampled-data','Interpreter','latex','FontWeight','normal');
ell_ax(nexttile(tl), rc.Sigma, Xc0, Xf,  gc); title('(c) switched','Interpreter','latex','FontWeight','normal');
ell_ax(nexttile(tl), rd.Sigma, Xd0, Xdf, gc); title('(d) Poissonian-jump','Interpreter','latex','FontWeight','normal');
mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_steer_combined.pdf'));

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

function P = ellpts(S, th)
S = (S+S.')/2; [V,Dg] = eig(S); d = max(real(diag(Dg)),0);
P = real(V*diag(sqrt(d))*V.') * [cos(th); sin(th)];
end
