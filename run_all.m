%% RUN_ALL  Reproduce every figure and data file of the monograph.
%  Thin wrapper over the per-chapter runners; to regenerate a single
%  chapter, run the corresponding script instead:
%    run_stability       Stability Analysis
%    run_oc              Linear-Cost Optimal Control
%    run_steering        Finite-Time Optimal Steering
%    run_dissipativity   Dissipativity
%    run_game            Minimax Optimal Control (games)
%    run_computational   Computational Considerations
%    run_higher          Higher-order Systems
%  Vector PDFs (white background, black text) are written to figures/
%  and data macros to figures/data/ under the CURRENT folder, or under the
%  folder named by the environment variable MVST_FIGDIR when it is set.
%  Run MVST_SETUP once first.
%
%  Runtime: several minutes end to end (the SOS/SDP examples dominate);
%  requires an SDP solver (SeDuMi/SDPT3/MOSEK) on the path for those.
run_stability
run_oc
run_steering
run_dissipativity
run_game
run_computational
run_higher
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_cones.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_cones_soc.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_cones_invariance.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_cones_quat.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_pos.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_perf.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_perf_classes.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_pos_switched.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_pos_linf.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_pos_gain_dwell.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_pos_smallnoise.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_perf_dwell.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_oc_sos_poisson.m'))
run(fullfile(fileparts(mfilename('fullpath')),'examples','ex_oc_sos_sd.m'))
