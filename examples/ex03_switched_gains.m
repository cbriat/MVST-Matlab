%% EX03_SWITCHED_GAINS  Switched system: P(t), gains, and the three regimes.
%  Infinite-horizon per-mode gains are constant (exported).  The
%  finite-horizon dual variable P(t), integrated backward through a
%  switching schedule, is continuous at switches (identity jump) but its
%  derivative and the active gain K(t) change at each switch; both are
%  plotted.  The long-term average rate rho* comes from a forced simulation.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
modes = D.sw.modes; Td = D.sw.Tdwell; M = numel(modes); n = modes{1}.n;

% per-mode infinite-horizon gains (constant)
Kinf = cell(1,M);
for i = 1:M; [~, Kinf{i}] = mvst_care_ct(modes{i}); end
fprintf('K1_inf = [% .4f % .4f],  K2_inf = [% .4f % .4f]\n', Kinf{1}, Kinf{2});

% ---- finite-horizon backward Riccati through a switching schedule ----
T = 4; sw = Td:Td:(T-1e-9); edges = [0 sw T];
seq = 1 + mod(numel(edges)-2:-1:0, M);      % active mode per interval (time order)
rt = 1e-8; oo = odeset('RelTol',rt,'AbsTol',rt*1e-2);
t = []; Pt = zeros(n,n,0); Kt = zeros(1,n,0);
Pcur = zeros(n);                            % terminal P(T)=0
for i = numel(edges)-1:-1:1                 % integrate intervals backward
    md = modes{seq(i)};
    f = @(tt,pv) reshape(-mvst_ric_ct(reshape(pv,n,n), md), [], 1);
    [ts, Pv] = ode45(f, [edges(i+1) edges(i)], Pcur(:), oo);
    [ts, ix] = sort(ts); Pv = Pv(ix,:);
    for k = 1:numel(ts)
        Pk = reshape(Pv(k,:), n, n); Pk = (Pk+Pk.')/2;
        Pt(:,:,end+1) = Pk;                 %#ok<AGROW>
        [~, Kk] = mvst_ric_ct(Pk, md); Kt(:,:,end+1) = Kk; %#ok<AGROW>
    end
    t = [t; ts];                            %#ok<AGROW>
    Pcur = reshape(Pv(1,:), n, n); Pcur = (Pcur+Pcur.')/2;  % continuous across switch
end
[t, ix] = sort(t); Pt = Pt(:,:,ix); Kt = Kt(:,:,ix);

% ---- long-term average rate via forced switched simulation ----
tf = 30; tks = Td:Td:(tf-1e-9); sq = 1 + mod(0:numel(tks), M);
Acl = cell(1,M); Lc = cell(1,M); L = cell(1,M);
for i = 1:M
    md = modes{i}; Acl{i} = md.A + md.B*Kinf{i};
    Lc{i} = zeros(n,n,md.Nc);
    for l = 1:md.Nc; Lc{i}(:,:,l) = md.Ac(:,:,l) + md.Bc(:,:,l)*Kinf{i}; end
    L{i} = md.Q + Kinf{i}.'*md.R*Kinf{i};
end
iv = cell(1,numel(tks)+1);
for i = 1:numel(iv); mi = sq(i); iv{i} = struct('Acl',Acl{mi},'Lc',Lc{mi},'W',modes{mi}.W,'L',L{mi}); end
av = mvst_simulate_impulsive(iv, tks, @(X,k) X, modes{1}.X0, tf, struct('npts',16));
rho = av.runavg(end);
fprintf('rho* (forced steady state) = %.4f\n', rho);

% ---- export ----
w = mvst_texdata(fullfile(mvst_figdir(),'data','ex_switched.tex'));
w.mat('AswOne', modes{1}.A); w.mat('AswTwo', modes{2}.A);
w.mat('BswOne', modes{1}.B); w.mat('BswTwo', modes{2}.B);
w.mat('QswOne', modes{1}.Q); w.mat('QswTwo', modes{2}.Q);
w.num('RswOne', modes{1}.R); w.num('RswTwo', modes{2}.R);
w.mat('KswInfOne', Kinf{1}); w.mat('KswInfTwo', Kinf{2});
w.num('Tdwellsw', Td); w.num('RhoSw', rho);
w.close();

% ---- figure: P(t) (continuous, kinks) and K(t) (jumps at switches) ----
fig = figure('Color','w','Units','centimeters','Position',[2 2 16 11]);
tl = tiledlayout(fig,2,1,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
plot(ax1, t, squeeze(Pt(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(1,:));
plot(ax1, t, squeeze(Pt(2,2,:)), '-','LineWidth',2.0,'Color',gc.order(2,:));
plot(ax1, t, squeeze(Pt(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(3,:));
for x = sw; xline(ax1, x, ':','Color',[0.6 0.6 0.6]); end
ylabel(ax1,'dual $P(t)$','Interpreter','latex');
lg=legend(ax1,{'$P_{11}$','$P_{22}$','$P_{12}$'},'Interpreter','latex','Location','northeastoutside');
set(lg,'Box','off'); set(ax1,'FontSize',11);

ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, t, squeeze(Kt(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(1,:));
plot(ax2, t, squeeze(Kt(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(2,:));
for x = sw; xline(ax2, x, ':','Color',[0.6 0.6 0.6]); end
xlabel(ax2,'time $t$','Interpreter','latex'); ylabel(ax2,'gain $K(t)$','Interpreter','latex');
lg2=legend(ax2,{'$K_1(t)$','$K_2(t)$'},'Interpreter','latex','Location','northeastoutside');
set(lg2,'Box','off'); set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_switched_gains.pdf'));
