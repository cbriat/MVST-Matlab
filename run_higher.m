%% RUN_HIGHER  Figures and data for the Higher-order Systems chapter.
%  Runs the tensor-power-lift illustrations:
%    ex_high_lift        lifted certificates for an impulsive system
%                        (Theorem th:high:lifted), data only, no figure;
%    ex_high_hierarchy   minimum dwell-time hierarchy T_m decreasing toward
%                        the exact value (Proposition prop:high:monotone,
%                        Theorem th:high:approach), needs the MathWorks LMI
%                        Lab (feasp); writes high_hierarchy.pdf;
%    ex_high_transient   transient-growth instance of sec:high:emergence:
%                        no Loewner greatest element of the lifted Gramian
%                        set on a finite grid of S_min[Tbar], m = 1..6
%                        (numerical evidence), needs LMI Lab (feasp);
%                        no figure, no data macro;
% ex_high_compound (exterior/compound material) is not part of the batch; the script stays in examples/.
%    ex_high_compound    second additive compound and 2-positivity
%                        (Theorem th:high:kpos); writes high_compound.pdf.
%  Outputs: vector PDFs to figures/, data macros to figures/data/.
close all
% requires mvst_setup (library and solvers on the path)
% ex_high_compound (exterior/compound material) is not part of the batch; the script stays in examples/.
mvst_runscripts({ ...
    'ex_high_lift', ...
    'ex_high_hierarchy', ...
    'ex_high_transient'});
%    'ex_high_compound'  not part of the batch.
