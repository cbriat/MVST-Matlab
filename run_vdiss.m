%% RUN_VDISS  Figure and checks for the vector dissipativity appendix.
%  Runs the vector-dissipativity example (ex_vdiss).  Outputs: vector PDF to
%  figures/, data macros to figures/data/.  Runtime: order of seconds
%  (no SDP solver required).
close all
% requires mvst_setup (library and solvers on the path)
mvst_runscripts({'ex_vdiss'});
