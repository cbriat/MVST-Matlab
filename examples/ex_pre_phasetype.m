function ex_pre_phasetype()
%EX_PRE_PHASETYPE  Phase-type density shapes for Chapter 2 (preliminaries),
%  figure fig:pre:cov:phasetype:densities.  Plots three densities on a common
%  s-grid: an exponential, an Erlang, and a BIMODAL MIXTURE OF TWO ERLANG
%  densities whose two modes sit at distant values.  The densities are the
%  closed forms of the phase-type specializations of eq:pre:cov:sdph:density:
%     exponential:  f(s) = lam e^{-lam s},
%     Erlang(nu,lam):  f(s) = lam^nu s^{nu-1} e^{-lam s} / (nu-1)!.
%  Pure figure: no data macros are written (as for the ex_vstab_* scripts).

here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..');
% requires mvst_setup (library and solvers on the path)
mvst_latexdefaults(); gc = mvst_gemcolors();
figdir = mvst_figdir();

% ------------------------------------------------------------- parameters
% Exponential(lam): mean 1/lam.  Take lam = 1  ->  mean 1.
lamE = 1;
% Erlang(nu,lam): mean nu/lam.  Take nu = 5, lam = 5  ->  mean 1 (matched
% to the exponential for comparability).
nuA = 5;  lamA = 5;
% Bimodal mixture  w*Erlang(nu1,lam1) + (1-w)*Erlang(nu2,lam2).  The mode of
% Erlang(nu,lam) is (nu-1)/lam:  Erlang(5,4) has mode 1, Erlang(12,2) has
% mode 5.5, so the two modes are well separated.
w   = 0.40;
nu1 = 5;   lam1 = 4;
nu2 = 12;  lam2 = 2;

erlpdf = @(s,nu,lam) lam.^nu .* s.^(nu-1) .* exp(-lam.*s) ./ factorial(nu-1);
exppdf = @(s,lam)    lam .* exp(-lam.*s);

s  = linspace(0, 10, 2001);
fE = exppdf(s, lamE);
fA = erlpdf(s, nuA, lamA);
fM = w*erlpdf(s, nu1, lam1) + (1-w)*erlpdf(s, nu2, lam2);

% ---------------------------- verify bimodality: two interior local maxima
dfM   = diff(fM);
sg    = sign(dfM);
ismax = [false, (sg(1:end-1) > 0) & (sg(2:end) < 0), false];   % rise then fall
midx  = find(ismax);
modes = s(midx);
assert(numel(midx) == 2, ...
    'mixture is not bimodal: found %d interior local maxima', numel(midx));
assert(abs(modes(1) - 1) < 0.3 && abs(modes(2) - 5.5) < 0.5, ...
    'mixture modes not at the expected locations: %s', mat2str(modes,3));

fprintf('exponential : lam = %.2f  (mean %.2f)\n', lamE, 1/lamE);
fprintf('Erlang      : nu = %d, lam = %.2f  (mean %.2f)\n', nuA, lamA, nuA/lamA);
fprintf('mixture     : w = %.2f  Erlang(%d,%.2f) mode %.2f  +  (1-w) Erlang(%d,%.2f) mode %.2f\n', ...
        w, nu1, lam1, (nu1-1)/lam1, nu2, lam2, (nu2-1)/lam2);
fprintf('mixture interior local maxima at s = %s  (peaks %.3f, %.3f)\n', ...
        mat2str(modes,3), fM(midx(1)), fM(midx(2)));

% ----------------------------------------- figure: single full-width axes
fig = figure('Color','w','Units','centimeters','Position',[2 2 14 8.5]);
ax  = axes(fig); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
plot(ax, s, fE, '-', 'Color', gc.sapphire, 'LineWidth', 1.4);
plot(ax, s, fA, '-', 'Color', gc.ruby,     'LineWidth', 1.4);
plot(ax, s, fM, '-', 'Color', gc.topaz,    'LineWidth', 1.4);
xlabel(ax, '$s$', 'Interpreter','latex');
ylabel(ax, 'density $f(s)$', 'Interpreter','latex');
legend(ax, {'exponential', ...
            sprintf('Erlang($\\nu=%d$)', nuA), ...
            'mixture of two Erlang'}, ...
       'Location','northeast','Interpreter','latex');
set(ax,'FontSize',11);
mvst_savepdf(fig, fullfile(figdir,'pre_phasetype_densities.pdf'));
fprintf('  wrote pre_phasetype_densities.pdf\n');
end
