function ex_pos_gain_dwell()
%EX_POS_GAIN_DWELL  L-infinity (peak) gain versus the minimum dwell-time for a
%   POSITIVE impulsive system (Chapter "Linear Positive Systems"), the
%   copositive dwell-time sweep of the companion NAHS study.  The certified peak
%   gain gamma(T_min) decreases from the arbitrary-dwell level (T_min -> 0) to
%   the jump-free continuous-time level (T_min -> infinity) as the guaranteed
%   dwell grows, both marked.
%
%   Positive impulsive instance (Metzler flow, nonnegative channels/reset):
%     flow  dx=A x+E_c w, z=C_c x+F_c w,  A=[-1 0;1 -2] (Metzler, Hurwitz),
%           E_c=[0.1;1.1], C_c=[0 1], F_c=0.3;
%     jump  x+=J x+E_d w, z_d=C_d x+F_d w, J=[0.1 1;0.1 0.1]>=0,
%           E_d=[0.3;0.3], C_d=[0 1], F_d=0.03.
%   gamma(T_min) is the minimum-dwell L-infinity gain (forward timer certificate,
%   frozen tail eq:pos:pktail); the arbitrary-dwell level is the common-certificate
%   copositive LP; the asymptotic level is the jump-free continuous L-infinity gain
%   ||C_c A^{-1} E_c + F_c||.
%
%   Outputs: figures/pos_gain_dwell.pdf (gamma vs T_min, two reference levels),
%   macros in figures/data/pos_gain_dwell.tex.
%
%   See also MVST_POS_GAIN.

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

s = struct('A',[-1 0;1 -2],'Ec',[0.1;1.1],'Cc',[0 1],'Fc',0.3, ...
           'J',[0.1 1;0.1 0.1],'Ed',[0.3;0.3],'Cd',[0 1],'Fd',0.03);

% ---- minimum-dwell L-infinity gain sweep ----------------------------------
Tg = 0.4:0.2:2.4; gI = nan(size(Tg));
for i = 1:numel(Tg)
    r = mvst_pos_gain(s, struct('norm','Linf','class','impulsive', ...
                                'dwell','min','Tmin',Tg(i),'deg',4));
    if r.feas; gI(i) = r.gamma; end
end

% ---- reference levels ------------------------------------------------------
arb = mvst_pos_gain(s, struct('norm','Linf','class','impulsive','dwell','arbitrary'));
ct  = mvst_pos_gain(struct('A',s.A,'Ec',s.Ec,'Cc',s.Cc,'Fc',s.Fc), ...
                    struct('norm','Linf','class','ct'));
gArb = arb.gamma; gAsym = ct.gamma;
% representative L1 (trace) gains for the text
l1a = mvst_pos_gain(s, struct('norm','L1','class','impulsive','dwell','min','Tmin',1,'deg',4));
l1b = mvst_pos_gain(s, struct('norm','L1','class','impulsive','dwell','min','Tmin',2,'deg',4));

fprintf('[ex_pos_gain_dwell] Linf: arb=%.4f, gamma(1)=%.4f, gamma(2)=%.4f, asym=%.4f\n', ...
        gArb, interp1(Tg,gI,1), interp1(Tg,gI,2), gAsym);

% ---- figure ---------------------------------------------------------------
figure('Color','w','Units','centimeters','Position',[2 2 13 7.5]);
plot(Tg, gI, '-o', 'Color', gc.order(1,:), 'LineWidth', 1.8, 'MarkerSize',4, ...
     'MarkerFaceColor', gc.order(1,:)); hold on
yline(gArb, '--', 'Color', gc.order(2,:), 'LineWidth', 1.3);
yline(gAsym, ':', 'Color', gc.order(3,:), 'LineWidth', 1.5);
grid on; box on
xlabel('minimum dwell-time $T_{\min}$','Interpreter','latex');
ylabel('certified $L_\infty$ (peak) gain $\gamma$','Interpreter','latex');
% Legend entries are padded on the right (\hspace) so the frame clears the longest LaTeX entry.
legend({'minimum-dwell-time $\gamma(T_{\min})$\hspace{0.8em}', ...
        'arbitrary-dwell-time level\hspace{0.8em}', ...
        'jump-free level\hspace{0.8em}'},'Interpreter','latex','Location','east');
mvst_savepdf(gcf, fullfile(figdir,'pos_gain_dwell.pdf'));

% ---- data macros ----------------------------------------------------------
try
    w = mvst_texdata(fullfile(figdir,'data','pos_gain_dwell.tex'));
    w.num('PosGdArb',gArb); w.num('PosGdAsym',gAsym);
    w.num('PosGdOne',interp1(Tg,gI,1)); w.num('PosGdTwo',interp1(Tg,gI,2));
    w.num('PosGdLoneOne',l1a.gamma); w.num('PosGdLoneTwo',l1b.gamma);
    w.close();
catch, end
end
