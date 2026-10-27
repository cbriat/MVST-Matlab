function ex_pos_smallnoise()
%EX_POS_SMALLNOISE  Small-noise mean-square stability threshold for a stochastic
%   positive system (Chapter "Linear Positive Systems", Proposition
%   prop:pos:dnndiag(b)): the canonical diagonal certificate D=diag(e_i/d_i)
%   certifies mean-square stability iff the (isotropic) noise level s satisfies
%       s < s* = min_i (d_i + e_i)/(d_i e_i),
%   where d=-A0^{-1} 1 and e=-A0^{-T} 1 (the flow and adjoint occupation vectors).
%   The certified verdict flips exactly at s* as s is swept.
%
%   Instance:  Metzler flow A0=[-1 0;1 -2], isotropic Brownian channel
%   Ac=sqrt(s) I (so the covariance drift gains the diagonal term
%   Sigma_i=sum_l (Ac_l)_ii^2 = s); the certificate is the diagonal quadratic of
%   prop:pos:dnndiag, canonical mode.
%
%   Outputs: figures/pos_smallnoise.pdf (certified verdict versus s, s* marked),
%   macros in figures/data/pos_smallnoise.tex.
%
%   See also MVST_POS_MS.

here = fileparts(mfilename('fullpath')); root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
gc = mvst_gemcolors(); mvst_latexdefaults();
figdir = mvst_figdir();

A0 = [-1 0;1 -2];
d = -A0\ones(2,1);  e = -A0.'\ones(2,1);
sstar = min((d+e)./(d.*e));

% ---- sweep the noise level, record the certified verdict -------------------
sg = linspace(0.5, 2.5, 41); feas = zeros(size(sg));
for i = 1:numel(sg)
    r = mvst_pos_ms(struct('A0',A0,'Ac',sqrt(sg(i))*eye(2)), ...
                    struct('what','stability','cert','diag','canonical',true));
    feas(i) = double(r.feas);
end
% empirical flip location (last feasible sweep point)
kflip = find(feas>0,1,'last'); sflip = sg(kflip);

fprintf('[ex_pos_smallnoise] s*=%.5f (closed form)  empirical flip in (%.3f,%.3f]\n', ...
        sstar, sflip, sg(min(kflip+1,end)));

% ---- figure ---------------------------------------------------------------
figure('Color','w','Units','centimeters','Position',[2 2 13 7.5]);
stairs(sg, feas, '-', 'Color', gc.order(1,:), 'LineWidth', 1.8); hold on
xline(sstar, '--', 'Color', gc.order(2,:), 'LineWidth', 1.5);
ylim([-0.1 1.1]); yticks([0 1]); yticklabels({'infeasible','certified'});
grid on; box on
xlabel('noise level $s$ (isotropic channel $A_c=\sqrt{s}\,I$)','Interpreter','latex');
ylabel('canonical certificate verdict','Interpreter','latex');
legend({'certified verdict','$s^\star=\min_i (d_i+e_i)/(d_i e_i)$'}, ...
       'Interpreter','latex','Location','east');
mvst_savepdf(gcf, fullfile(figdir,'pos_smallnoise.pdf'));

% ---- data macros ----------------------------------------------------------
try
    w = mvst_texdata(fullfile(figdir,'data','pos_smallnoise.tex'));
    w.num('PosSnStar',sstar); w.mat('PosSnD',d.'); w.mat('PosSnE',e.');
    w.close();
catch, end
end
