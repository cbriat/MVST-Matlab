% Matrix-Valued Systems Toolbox (MVST)
% Version 1.0.0-rc1
%
% Setup and diagnostics
%   mvst_setup        - Put the library on the path and run the solver doctor.
%   mvst.version      - Toolbox version string.
%   mvst.solvercheck  - Report which SDP / SOS backends resolve and solve.
%   mvst.Options      - Solver and algorithm options.
%   mvst.System       - Unified system object (factories ct, mjls, sampled, poisson, switched, impulsive).
%
% Sum-of-squares engines (SOSTOOLS)
%   mvst_dof_sos                - Generic SOSTOOLS timer-dependent output-feedback bracket SYNTHESIS for the scalar dwell-time output-feedback problem.
%   mvst_oc_sos                 - Optimal control by the dual (storage-function) SOS/LMI program for a SINGLE linear Hermitian matrix-valued mode.
%   mvst_oc_sos_modes           - Optimal control by the dual (storage-function) SOS/LMI program for a COUPLED Markov-jump system.
%   mvst_sosmatrvar             - Free (sign-indefinite) polynomial matrix variable.
%   mvst_sosposmatrvar          - SOS (positive-semidefinite) polynomial matrix variable.
%   mvst_sossolve               - Solve a SOSTOOLS program with the backend selected by OPTS.SOLVER.
%
% Stability and stabilization
%   mvst_stab_critical          - Critical parameter value on the mean-square stability margin.
%   mvst_stab_delay             - Delay-independent stability of matrix-valued (and positive) time-delay systems.
%   mvst_stab_sos               - Generic SOSTOOLS stability-ANALYSIS test for linear impulsive Hermitian matrix-valued systems.
%   mvst_stab_sos_modes         - Generic SOSTOOLS stability-ANALYSIS test on a PRODUCT CONE (S^n)^M for mode-indexed linear Hermitian matrix-valued systems.
%   mvst_stability              - Exact mean-square stability analysis for the book's system classes.
%   mvst_stabilize              - Exact state-feedback stabilization for the book's system classes.
%   mvst_stabz_sd_range         - Aperiodic (range) sampled-data refresh-gain synthesis for a held-input matrix-valued plant.
%   mvst_stabz_sos              - Generic SOSTOOLS state-feedback stabilization synthesis for linear impulsive Hermitian matrix-valued systems.
%
% Riccati equations and linear-cost optimal control
%   mvst_care_ct                - Infinite-horizon CT I-GRE (stochastic ARE) for one mode.
%   mvst_care_ct_poisson        - Algebraic generalized Riccati for GBM+Poisson systems.
%   mvst_dre_ct                 - Backward solution of the continuous-time differential I-GRE.
%   mvst_dre_ct_poisson         - Backward LQ Riccati for GBM + Poissonian-jump systems.
%   mvst_dre_ct_tv              - Backward time-varying differential I-GRE.
%   mvst_igre_periodic_sd       - Periodic-sampling sampled-data optimal control.
%   mvst_oc_dwell               - Causal minimum / range dwell-time optimal control of a single-cone linear impulsive system.
%   mvst_ric_ct                 - Continuous-time I-GRE residual and optimal gain for one mode.
%   mvst_ric_ct_tv              - Time-varying continuous-time I-GRE residual and gain.
%   mvst_solve_average          - Long-term average-cost optimal control (CT).
%   mvst_solve_finite           - Finite-horizon linear-cost optimal control (CT).
%   mvst_solve_finite_tv        - Finite-horizon linear-cost optimal control, TIME-VARYING.
%   mvst_solve_infinite         - Infinite-horizon linear-cost optimal control (CT).
%   mvst_solve_infinite_tv      - Infinite-horizon dual variable, TIME-VARYING.
%
% Covariance steering
%   mvst_steer                  - Class-dispatched finite-time optimal covariance steering.
%   mvst_steer_newton           - Damped Newton for the steering multiplier.
%
% Dissipativity, KYP, passivity, and H-infinity norm
%   mvst_diss_dwell_cert        - MDT dissipativity certificate for a switched system.
%   mvst_diss_pT_ct             - Finite-horizon storage matrix versus horizon length.
%   mvst_diss_sos               - Generic SOSTOOLS dissipativity-ANALYSIS program for linear impulsive Hermitian matrix-valued systems.
%   mvst_diss_sos_modes         - Generic SOSTOOLS dissipativity-ANALYSIS program on a PRODUCT CONE (H^n)^M.
%   mvst_dissyn_sos             - Generic SOSTOOLS dissipative state-feedback SYNTHESIS program for linear impulsive Hermitian matrix-valued systems.
%   mvst_hinfnorm               - H-infinity norm / (mean-square) L2-gain with a uniform result.
%   mvst_isp_margin_ct          - Input-strict-passivity index of a (stochastic) LTI flow.
%   mvst_kyp                    - Operator-domain Kalman-Yakubovich-Popov engine for matrix-valued dissipativity: the general KYP lemma.
%   mvst_l2gain                 - Class-dispatched L2-gain (thin router over the exact backends).
%   mvst_l2gain_ct              - Continuous-time L2-gain (H-infinity norm) by bisection.
%   mvst_l2gain_ct_dwell        - Dwell-time-certified L2-gain bound, switched system.
%   mvst_l2gain_ct_mjls         - Exact coupled mean-square L2-gain of an MJLS.
%   mvst_l2gain_ct_poisson      - Mean-square L2-gain of a GBM+Poisson system.
%   mvst_l2gain_dt              - Discrete-time L2-gain (H-infinity norm) by exact bisection.
%   mvst_l2gain_dt_mjls         - Coupled mean-square l2-gain of a discrete-time MJLS.
%   mvst_l2gain_sd_hybrid       - Exact hybrid L2-gain, flow disturbance / sampled output.
%   mvst_passivity              - Routed passivity analysis with indices (exact CT test by default).
%   mvst_passivity_ct           - Passivity test of a square LTI system by a frequency sweep.
%   mvst_qsr_supply             - Validate and unpack a quadratic (Q,S,R) supply contract.
%
% Zero-sum and Nash games
%   mvst_game_average           - Infinite-horizon average (ergodic) zero-sum game value.
%   mvst_game_dre_ct            - Backward game (HJI) differential Riccati, continuous time.
%   mvst_game_dre_ct_poisson    - Backward game (HJI) Riccati for GBM+Poisson systems.
%   mvst_game_gamma             - Class-dispatched H-infinity game level (thin router).
%   mvst_game_gamma_ct          - Critical disturbance level gamma* of the H-infinity game.
%   mvst_game_gamma_ct_dwell    - Certified achievable H-infinity level, switched.
%   mvst_game_gamma_ct_mjls     - Exact coupled critical level of the MJLS game.
%   mvst_game_gamma_ct_poisson  - Critical level gamma* for GBM+Poisson games.
%   mvst_game_gamma_sd_held     - Critical level of the held-input sampled-data game.
%   mvst_game_mixed             - Infinite-horizon mixed H2/H-infinity design as a two-player nonzero-sum Nash game.
%   mvst_game_nash              - Infinite-horizon N-player Nash equilibrium via the coupled I-GRE system.
%   mvst_game_sd_held           - Periodic HJI-GRE of the held-input sampled-data game.
%   mvst_game_sos               - Generic SOSTOOLS state-feedback zero-sum GAME / minimax program for linear impulsive Hermitian matrix-valued systems.
%
% Further performance measures
%   mvst_perf                   - Generic performance-analysis dispatcher.
%   mvst_perf_dwell             - Dwell-time and deterministic-peak backends for MVST_PERF.
%
% Linear positive systems
%   mvst_pos_gain               - Class-dispatched induced-gain copositive programs for linear positive systems.
%   mvst_pos_gain_sos           - Timer-dependent copositive-SOS induced-gain programs.
%   mvst_pos_lp                 - Internal LINPROG wrapper for the positive-systems LP family.
%   mvst_pos_ms                 - Second-moment (mean-square) certificates for STOCHASTIC positive systems.
%   mvst_pos_ms_dwell           - Timer-dependent (minimum / range dwell-time) second-moment certificates for stochastic positive systems.
%   mvst_pos_oc                 - Linear-cost optimal control of linear positive systems.
%   mvst_pos_stability          - Copositive linear-programming stability analysis for linear positive systems.
%   mvst_pos_stabilize          - Positivity-preserving state-feedback stabilization by copositive linear programming.
%   mvst_pos_stabilize_sos      - Timer-dependent copositive-SOS state-feedback synthesis.
%   mvst_pos_supply             - General-linear-supply dissipativity by copositive linear programming for LTI positive systems.
%
% System construction, conversion, and options
%   mvst_evalc                  - Evaluate a (possibly time-varying) coefficient at time t.
%   mvst_getf                   - Field of struct S named F, or default D if absent/empty.
%   mvst_mode                   - Build a (normalized) controlled-mode data structure.
%   mvst_mode_tv                - Build a TIME-VARYING controlled-mode data structure.
%   mvst_odecap                 - Terminal ODE event: stop when the state norm reaches CAP.
%   mvst_opts                   - Normalize solver options and verify the SDP/SOS solver is present.
%   mvst_system                 - Assemble a (possibly multi-mode) system for the solvers.
%   mvst_torich                 - Convert a thin mode/system struct to the rich SOS-layer struct.
%   mvst_tothin                 - Convert a rich SOS-layer struct back to a thin mode struct.
%
% Simulation
%   mvst_sample_ct              - Euler--Maruyama sample paths of a controlled diffusion.
%   mvst_sample_mjls            - Sample paths of a Markov-jump controlled diffusion.
%   mvst_simulate_ct            - Closed-loop second-moment (covariance) simulation (CT).
%   mvst_simulate_impulsive     - Closed-loop covariance of an impulsive system.
%   mvst_simulate_mjls          - Closed-loop second-moment simulation for an MJLS.
%
% Linear-algebra helpers
%   mvst_liftsym                - Reduced matrix representation of X -> He(A0 X) + sum C_l X C_l'.
%   mvst_svecbasis              - Orthonormal basis matrix of the symmetric subspace.
%   mvst_unvech                 - Rebuild a symmetric n-by-n matrix from MVST_VECH output.
%   mvst_vech                   - Stack the independent (upper-triangular) entries of symmetric M.
%
% Output folder of the examples
%   mvst_figdir                 - Output folder of the example scripts (figures and data macros).
%
