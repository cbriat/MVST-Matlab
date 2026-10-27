function ex_oc_dwell()
%EX_OC_DWELL  Minimum-dwell-time optimal control on an impulsive instance
%   (Chapter "Linear Cost", Section sec:OC:dwell; Theorem th:OC:MDT).  A stable
%   continuous flow is destabilized by EXPANSIVE resets; a minimum dwell-time
%   makes the guaranteed cost finite, and longer guaranteed dwell is strictly
%   valuable (the certified cost decreases toward the jump-free CARE value).
%
%   Instance (stated fully in the book text):
%     flow    dx = (A x + B u) dt + sigma x dW,  A=[-1 1;-1 -1] (Hurwitz),
%             B=[0;1],  sigma=0.2 (single multiplicative channel Ac=sigma I),
%             running cost x'Q x + u'R u,  Q=I, R=1;
%     jump    x+ = J x + Bd ud,  J=1.2 I (EXPANSIVE),  jump cost x'Qd x+ud'Rd ud,
%             Qd=0.05 I, Rd=1.
%   The certified cost bound is J_T(X*,X10) <= <P_s(0),X10> (Theorem th:OC:MDT),
%   an UPPER bound realized by the recovered dwell policy; X10=I.
%
%   Outputs: the cost bound versus the minimum dwell-time T_min (single-panel
%   figure figures/oc_dwell.pdf, the critical dwell T* marked), and the numbers
%   as macros in figures/data/oc_dwell.tex.
%
%   See also MVST_OC_DWELL, MVST_CARE_CT, EX_DATA.

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

% ---- instance --------------------------------------------------------------
A = [-1 1;-1 -1]; B = [0;1]; sigma = 0.2; Q = eye(2); R = 1;
J = 1.2*eye(2); Bd = [0;1]; Qd = 0.05*eye(2); Rd = 1; X10 = eye(2);
md   = mvst_mode(A, B, Q, R, 'Ac', sigma*eye(2));
jump = struct('J',J,'Bd',Bd,'Qd',Qd,'Rd',Rd);

% ---- critical dwell T* (feasibility onset, Theorem th:OC:MDT LMI (iii)) -----
lo = 0.1; hi = 0.5;
for it = 1:30
    m = 0.5*(lo+hi);
    r = mvst_oc_dwell(md, jump, struct('dwell','min','Tmin',m,'X10',X10));
    if r.feas; hi = m; else; lo = m; end
end
Tstar = hi;

% ---- cost bound versus minimum dwell-time ---------------------------------
Tg = linspace(Tstar+0.01, 3, 40); Jb = nan(size(Tg));
for i = 1:numel(Tg)
    r = mvst_oc_dwell(md, jump, struct('dwell','min','Tmin',Tg(i),'X10',X10));
    if r.feas; Jb(i) = r.Jbound; end
end
Pc  = mvst_care_ct(md);  Jcare = trace(Pc*X10);   % jump-free CARE limit
% representative points
r1 = mvst_oc_dwell(md, jump, struct('dwell','min','Tmin',0.5,'X10',X10));
r2 = mvst_oc_dwell(md, jump, struct('dwell','min','Tmin',1.0,'X10',X10));

fprintf('[ex_oc_dwell] T*=%.4f  J(0.5)=%.4f  J(1)=%.4f  CARE=%.4f\n', ...
        Tstar, r1.Jbound, r2.Jbound, Jcare);

% ---- figure ---------------------------------------------------------------
figure('Color','w','Units','centimeters','Position',[2 2 13 7.5]);
plot(Tg, Jb, '-', 'Color', gc.order(1,:), 'LineWidth', 1.8); hold on
yline(Jcare, '--', 'Color', gc.order(2,:), 'LineWidth', 1.3);
xline(Tstar, ':', 'Color', gc.order(3,:), 'LineWidth', 1.5);
grid on; box on
xlabel('minimum dwell-time $T_{\min}$','Interpreter','latex');
ylabel('certified cost bound $\langle P_s(0),X_1^0\rangle$','Interpreter','latex');
% Legend entries are padded on the right (\hspace) so the frame clears the longest LaTeX entry.
legend({'certified cost bound\hspace{0.8em}','jump-free CARE limit\hspace{0.8em}', ...
        'critical dwell-time $T^\star$\hspace{0.8em}'}, ...
       'Interpreter','latex','Location','northeast');
mvst_savepdf(gcf, fullfile(figdir,'oc_dwell.pdf'));

% ---- data macros ----------------------------------------------------------
try
    w = mvst_texdata(fullfile(figdir,'data','oc_dwell.tex'));
    w.num('OcDwTstar',Tstar); w.num('OcDwJhalf',r1.Jbound);
    w.num('OcDwJone',r2.Jbound); w.num('OcDwCare',Jcare);
    w.num('OcDwSigma',sigma);
    w.close();
catch, end
end
