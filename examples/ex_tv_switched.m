%% EX_TV_SWITCHED  Time-varying switched optimal control: P(t) and gains K_i(t).
%  Time-varying counterpart of EX03_SWITCHED_GAINS: each mode's flow drift
%  decays in time, A_i(t) = A_i + dA*exp(-t), so the system is genuinely
%  time-varying but settles to the LTI switched system of EX03 as t grows.
%  The finite-horizon dual variable P(t) is integrated backward through the
%  switching schedule (continuous across switches, identity jump) using the
%  time-varying differential I-GRE; the active gain K(t) changes at each
%  switch.  Both are plotted.  This illustrates the main contribution of the
%  paper for the switched case over a time-varying horizon.
here = fileparts(mfilename('fullpath'));
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults();
D = ex_data(); gc = mvst_gemcolors();
modes = D.sw.modes; Td = D.sw.Tdwell; M = numel(modes); n = modes{1}.n;
m = modes{1}.m;
dA = [0 0; 1.5 0];                          % decaying drift perturbation

% ----- build the time-varying modes A_i(t) = A_i + dA*exp(-t) -----
md = cell(1,M);
for i = 1:M
    b = modes{i};
    md{i} = mvst_mode_tv(n, m, ...
        'A', @(t) b.A + dA*exp(-t), 'B', b.B, 'Q', b.Q, 'R', b.R, ...
        'Ac', b.Ac, 'Bc', b.Bc, 'W', b.W, 'X0', b.X0);
end

% ----- finite-horizon backward time-varying I-GRE through the schedule ----
T = 4; sw = Td:Td:(T-1e-9); edges = [0 sw T];
seq = 1 + mod(numel(edges)-2:-1:0, M);      % active mode per interval (time order)
rt = 1e-8; oo = odeset('RelTol',rt,'AbsTol',rt*1e-2);
t = []; Pt = zeros(n,n,0);
Kall = cell(1,M); for i = 1:M; Kall{i} = zeros(m,n,0); end
Pcur = zeros(n);                            % terminal P(T)=0
for i = numel(edges)-1:-1:1                 % integrate intervals backward
    am = md{seq(i)};
    f = @(tt,pv) reshape(-mvst_ric_ct_tv(tt, sym(pv,n), am), [], 1);
    [ts, Pv] = ode45(f, [edges(i+1) edges(i)], Pcur(:), oo);
    [ts, ix] = sort(ts); Pv = Pv(ix,:);
    for k = 1:numel(ts)
        Pk = sym(Pv(k,:).', n);
        Pt(:,:,end+1) = Pk;                 %#ok<AGROW>
        for q = 1:M
            [~, Kq] = mvst_ric_ct_tv(ts(k), Pk, md{q});
            Kall{q}(:,:,end+1) = Kq;        %#ok<AGROW>
        end
    end
    t = [t; ts];                            %#ok<AGROW>
    Pcur = sym(Pv(1,:).', n);               % continuous across switch
end
[t, ix] = sort(t); Pt = Pt(:,:,ix);
for q = 1:M; Kall{q} = Kall{q}(:,:,ix); end

% ----- export data macros -----
w = mvst_texdata(fullfile(mvst_figdir(),'data','tv_sw.tex'));
w.mat('TVdAsw', dA); w.num('TVTdwellsw', Td); w.close();

% ----- figure: P(t) entries (top) and gains K1(t),K2(t) (bottom) -----
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

% The two modes share the control-relevant data (B, Bc, Ac, R), so the
% optimal gain K=-Theta^{-1}S is mode-independent: only the two ENTRIES of
% the active gain vector K(t)=[K_1(t) K_2(t)] are distinct (matching the LTI
% counterpart EX03_SWITCHED_GAINS), not four per-mode curves.
ax2 = nexttile(tl); hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on');
plot(ax2, t, squeeze(Kall{1}(1,1,:)), '-','LineWidth',2.0,'Color',gc.order(1,:));
plot(ax2, t, squeeze(Kall{1}(1,2,:)), '-','LineWidth',2.0,'Color',gc.order(2,:));
for x = sw; xline(ax2, x, ':','Color',[0.6 0.6 0.6]); end
xlabel(ax2,'time $t$','Interpreter','latex'); ylabel(ax2,'gain $K(t)$','Interpreter','latex');
lg2=legend(ax2,{'$K_1(t)$','$K_2(t)$'}, ...
      'Interpreter','latex','Location','northeastoutside');
set(lg2,'Box','off'); set(ax2,'FontSize',11);

mvst_savepdf(fig, fullfile(mvst_figdir(),'ex_tv_switched.pdf'));

% ========================================================================
function S = sym(v, n)
S = reshape(v, n, n);  S = (S + S.')/2;
end
