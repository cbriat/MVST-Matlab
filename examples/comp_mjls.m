%% COMP_MJLS  Relaxation (Gauss-Seidel) solution of the coupled Markov-jump ARE.
%  Illustrates Section "Coupled and switched algebraic equations": the coupled
%  infinite-horizon I-GRE of the two-mode Markov-jump example is solved by a
%  Gauss-Seidel sweep over the modes, each inner step solving one single-mode
%  ARE with the off-diagonal coupling frozen into the cost and the diagonal
%  rate folded into the drift.  The per-mode residual of the FULL coupled
%  equation is recorded at each sweep and shown to decrease monotonically.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
D = ex_data(); sys = D.mjls; M = numel(sys.modes); n = sys.modes{1}.n;
Pi = sys.Pi;

nsweep = 8;
P = repmat({zeros(n)}, 1, M);
res = zeros(nsweep, M);
for s = 1:nsweep
    for i = 1:M                                  % Gauss-Seidel sweep
        md = sys.modes{i}; mdi = md;
        mdi.A = md.A + 0.5*Pi(i,i)*eye(n);       % diagonal rate -> drift
        Qadd = zeros(n);
        for j = 1:M; if j ~= i; Qadd = Qadd + Pi(i,j)*P{j}; end; end
        mdi.Q = md.Q + (Qadd+Qadd.')/2;          % off-diagonal coupling -> cost
        P{i} = mvst_care_ct(mdi, struct());
    end
    for i = 1:M                                  % full coupled residual
        Cpl = zeros(n);
        for j = 1:M; Cpl = Cpl + Pi(i,j)*P{j}; end
        G = mvst_ric_ct(P{i}, sys.modes{i}, Cpl);
        res(s,i) = norm(G, 'fro') + eps;         % +eps so log scale is finite
    end
end
fprintf('coupled MJLS ARE: final residuals = [%.2e %.2e]\n', res(end,1), res(end,2));

fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
ax = axes(fig); hold(ax,'on'); box(ax,'on'); grid(ax,'on'); set(ax,'YScale','log');
h(1) = plot(ax, 1:nsweep, res(:,1), '-o','LineWidth',1.8,'MarkerSize',4,'Color',gc.order(1,:));
h(2) = plot(ax, 1:nsweep, res(:,2), '-s','LineWidth',1.8,'MarkerSize',4,'Color',gc.order(2,:));
xlabel(ax,'sweep','Interpreter','latex');
ylabel(ax,'coupled ARE residual $\|G_i(P_i)\|_F$','Interpreter','latex');
lg = legend(ax, h, {'mode 1','mode 2'}, 'Interpreter','latex', ...
            'Location','northeast'); set(lg,'Box','off'); set(ax,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'comp_mjls.pdf'));
