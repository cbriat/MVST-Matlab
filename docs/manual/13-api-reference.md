# 13. API reference

[Manual index](index.md) | Prev: [Further performance measures](11-performance.md)

Every public function of the toolbox, grouped by area, with the first line of
its help (`help <name>` for the full text). The table below is generated
directly from the source files, so it stays accurate. A few summaries are the
first physical help line and read as sentence fragments; the continuation is in
the function's own help.

The recommended entry point is the `+mvst` package (`mvst.System`,
`mvst.Options`, `mvst.solvercheck`); every `mvst_*` engine remains callable
directly.

## How to read this reference

The tables group the toolbox by area, in the order the [manual map](index.md)
walks the material. Each row names a public function and the first line of its
help. To use a function you need three things, all reachable from here:

- **What it takes.** Call `help <name>` for the full input contract. Systems are
  built through the `mvst.System` factories (see [Systems](02-systems.md)); the
  low-level engines also accept the thin/rich structs described in
  `lib/MODEL_CONVENTION.md`.
- **What it returns.** Every analysis or synthesis routine returns a single
  result struct. The recurring field names follow the conventions in
  [Result-struct conventions](#result-struct-conventions) below, so once you
  know the pattern you can read any engine's output.
- **How it can fail.** Errors carry a structured identifier in the `mvst:...`
  namespace ([Error-id namespace](#error-id-namespace)); catch them by id rather
  than by message text.

The two-layer surface is intentional. The object layer (`S.stability()`,
`S.l2gain()`, `S.perf(...)`) is the release-facing API and is what the manual
uses; the `mvst_*` engines are the computational core the object methods
dispatch to, and stay callable directly for scripting and for the book's example
files. Where a name below reads as a sentence fragment, it is the first physical
help line; the continuation is in the function's own help.

### Object API (`+mvst` package)
| Function | Summary |
|---|---|
| `mvst.System` | Unified value-type wrapper for a matrix-valued dynamical system. |
| `mvst.Options` | Solver options for the MVST package. |
| `mvst.solvercheck` | Detect and probe the SDP / SOS toolchain on the MATLAB path. |

### Model builders and bridges
| Function | Summary |
|---|---|
| `mvst_mode` | Build a (normalized) controlled-mode data structure. |
| `mvst_mode_tv` | Build a TIME-VARYING controlled-mode data structure. |
| `mvst_system` | Assemble a (possibly multi-mode) system for the solvers. |
| `mvst_torich` | Convert a thin mode/system struct to the rich SOS-layer struct. |
| `mvst_tothin` | Convert a rich SOS-layer struct back to a thin mode struct. |
| `mvst_getf` | Field of struct S named F, or default D if absent/empty. |
| `mvst_evalc` | Evaluate a (possibly time-varying) coefficient at time t. |

### Symmetric-subspace algebra
| Function | Summary |
|---|---|
| `mvst_svecbasis` | Orthonormal basis matrix of the symmetric subspace. |
| `mvst_liftsym` | Reduced matrix representation of X -> He(A0 X) + sum C_l X C_l'. |
| `mvst_vech` | Stack the independent (upper-triangular) entries of symmetric M. |
| `mvst_unvech` | Rebuild a symmetric n-by-n matrix from MVST_VECH output. |

### Stability and stabilization
| Function | Summary |
|---|---|
| `mvst_stability` | Exact mean-square stability analysis for the book's system classes. |
| `mvst_stabilize` | Exact state-feedback stabilization for the book's system classes. |
| `mvst_stab_sos` | Generic SOSTOOLS stability-ANALYSIS test for linear impulsive |
| `mvst_stab_sos_modes` | Generic SOSTOOLS stability-ANALYSIS test on a PRODUCT CONE |
| `mvst_stabz_sos` | Generic SOSTOOLS state-feedback stabilization synthesis for |
| `mvst_stabz_sd_range` | Aperiodic (range) sampled-data refresh-gain synthesis for a held-input matrix-valued plant, as a matrix sum-of-squares program in the timer. |
| `mvst_stab_critical` | Critical parameter value on the mean-square stability margin. |
| `mvst_stab_delay` | Delay-independent stability of matrix-valued (and positive) time-delay systems (Chapter "Stability", Section sec:stab:delay, and Chapter "Linear Positive Systems", prop:pos:classncs). |

### Optimal control (Riccati and SOS)
| Function | Summary |
|---|---|
| `mvst_solve_finite` | Finite-horizon linear-cost optimal control (CT). |
| `mvst_solve_infinite` | Infinite-horizon linear-cost optimal control (CT). |
| `mvst_solve_average` | Long-term average-cost optimal control (CT). |
| `mvst_solve_finite_tv` | Finite-horizon linear-cost optimal control, TIME-VARYING. |
| `mvst_solve_infinite_tv` | Infinite-horizon dual variable, TIME-VARYING. |
| `mvst_ric_ct` | Continuous-time I-GRE residual and optimal gain for one mode. |
| `mvst_ric_ct_tv` | Time-varying continuous-time I-GRE residual and gain. |
| `mvst_dre_ct` | Backward solution of the continuous-time differential I-GRE. |
| `mvst_dre_ct_tv` | Backward time-varying differential I-GRE. |
| `mvst_dre_ct_poisson` | Backward LQ Riccati for GBM + Poissonian-jump systems. |
| `mvst_care_ct` | Infinite-horizon CT I-GRE (stochastic ARE) for one mode. |
| `mvst_care_ct_poisson` | Algebraic generalized Riccati for GBM+Poisson systems. |
| `mvst_igre_periodic_sd` | Periodic-sampling sampled-data optimal control. |
| `mvst_oc_sos` | Optimal control by the dual (storage-function) SOS/LMI program |
| `mvst_oc_sos_modes` | Optimal control by the dual (storage-function) SOS/LMI |
| `mvst_oc_dwell` | Causal minimum / range dwell-time optimal control of a single-cone linear impulsive system (Chapter "Linear Cost", Section sec:OC:dwell): the timer-dependent forward-Riccati / jump-GRE programs Theorem th:OC:MDT (minimum dwell-time) and Theorem th:OC:RDT (range dwell-time), with the average-cost bounds th:IMP:AC:dwell:MDT / th:IMP:AC:dwell:RDT as a secondary mode. |
| `mvst_odecap` | Terminal ODE event: stop when the state norm reaches CAP. |

### Steering
| Function | Summary |
|---|---|
| `mvst_steer` | Class-dispatched finite-time optimal covariance steering. |
| `mvst_steer_newton` | Damped Newton for the steering multiplier. |
| `mvst_steer_components` | Covariance entries versus time. |
| `mvst_steer_ellipse` | Single-panel covariance-steering figure. |

### Dissipativity, L2-gain, and games
| Function | Summary |
|---|---|
| `mvst_l2gain` | Class-dispatched L2-gain (thin router over the exact backends). |
| `mvst_l2gain_ct` | Continuous-time L2-gain (H-infinity norm) by bisection. |
| `mvst_l2gain_ct_mjls` | Exact coupled mean-square L2-gain of an MJLS. |
| `mvst_l2gain_dt` | Discrete-time L2-gain (H-infinity norm) by exact bisection. |
| `mvst_l2gain_dt_mjls` | Coupled mean-square l2-gain of a discrete-time MJLS. |
| `mvst_l2gain_ct_dwell` | Dwell-time-certified L2-gain bound, switched system. |
| `mvst_l2gain_ct_poisson` | Mean-square L2-gain of a GBM+Poisson system. |
| `mvst_l2gain_sd_hybrid` | Exact hybrid L2-gain, flow disturbance / sampled output. |
| `mvst_diss_sos` | Generic SOSTOOLS dissipativity-ANALYSIS program for linear |
| `mvst_diss_sos_modes` | Generic SOSTOOLS dissipativity-ANALYSIS program on a |
| `mvst_diss_sos_selftest` | Validation of MVST_DISS_SOS against exact ground truth. |
| `mvst_diss_dwell_cert` | MDT dissipativity certificate for a switched system. |
| `mvst_diss_pT_ct` | Finite-horizon storage matrix versus horizon length. |
| `mvst_passivity_ct` | Passivity test of a square LTI system by a frequency sweep. |
| `mvst_passivity` | Routed passivity analysis with indices (exact CT test by default). |
| `mvst_hinfnorm` | H-infinity norm / (mean-square) L2-gain with a uniform result. |
| `mvst_isp_margin_ct` | Input-strict-passivity index of a (stochastic) LTI flow. |
| `mvst_kyp` | Operator-domain Kalman-Yakubovich-Popov engine for matrix-valued dissipativity (Chapter "Dissipativity", Section sec:diss:kyp): the general KYP lemma for an arbitrary-sign supply matrix M, in continuous time, discrete time, and the hybrid (periodic impulsive) setting. |
| `mvst_qsr_supply` | Validate and unpack a quadratic (Q,S,R) supply contract. |
| `mvst_dissyn_sos` | Generic SOSTOOLS dissipative state-feedback SYNTHESIS program |
| `mvst_dissyn_sos_selftest` | Validation of MVST_DISSYN_SOS by closed-loop |
| `mvst_dof_sos` | Generic SOSTOOLS timer-dependent output-feedback bracket |
| `mvst_game_gamma` | Class-dispatched H-infinity game level (thin router). |
| `mvst_game_gamma_ct` | Critical disturbance level gamma* of the H-infinity game. |
| `mvst_game_gamma_ct_mjls` | Exact coupled critical level of the MJLS game. |
| `mvst_game_gamma_ct_dwell` | Certified achievable H-infinity level, switched. |
| `mvst_game_gamma_ct_poisson` | Critical level gamma* for GBM+Poisson games. |
| `mvst_game_dre_ct` | Backward game (HJI) differential Riccati, continuous time. |
| `mvst_game_dre_ct_poisson` | Backward game (HJI) Riccati for GBM+Poisson systems. |
| `mvst_game_average` | Infinite-horizon average (ergodic) zero-sum game value (Chapter "Game", th:game:AC:finite / th:game:AC:suff, cor:game:AC:equiv). |
| `mvst_game_mixed` | Infinite-horizon mixed H2/H-infinity design as a two-player nonzero-sum Nash game (Chapter "Game", prop:game:mixed:verify, eq:game:mixed:GRE). |
| `mvst_game_nash` | Infinite-horizon N-player Nash equilibrium via the coupled I-GRE system (Chapter "Game", prop:game:nash:verify, eq:game:nash:GRE). |
| `mvst_game_sd_held` | Periodic HJI-GRE of the held-input sampled-data game. |
| `mvst_game_gamma_sd_held` | Critical level of the held-input sampled-data game. |
| `mvst_game_sos` | Generic SOSTOOLS state-feedback zero-sum GAME / minimax program |
| `mvst_game_sos_selftest` | Validation of MVST_GAME_SOS against exact ground truth. |

### Further performance measures
| Function | Summary |
|---|---|
| `mvst_perf` | Generic performance-analysis dispatcher (Chapter "Further |
| `mvst_perf_dwell` | Dwell-time and deterministic-peak backends for MVST_PERF (Chapter "Further Performance Measures"). |

### Positive systems (copositive LP)
| Function | Summary |
|---|---|
| `mvst_pos_stability` | Copositive linear-programming stability analysis for |
| `mvst_pos_gain` | Induced-gain copositive linear programs for LTI positive |
| `mvst_pos_stabilize` | Positivity-preserving state-feedback stabilization by |
| `mvst_pos_supply` | General-linear-supply dissipativity by copositive linear |
| `mvst_pos_lp` | Internal LINPROG wrapper for the positive-systems LP family. |
| `mvst_pos_gain_sos` | Timer-dependent copositive-SOS induced-gain programs for positive impulsive / switched systems (minimum and range dwell-time), the timer-dependent branch of MVST_POS_GAIN. |
| `mvst_pos_ms` | Second-moment (mean-square) certificates for STOCHASTIC positive systems (Chapter "Linear Positive Systems", stochastic layer): the covariance-level SDP certificates that see the Brownian and Poisson noise, which the mean-level copositive-LP layer (MVST_POS_STABILITY / _GAIN) is blind to. |
| `mvst_pos_ms_dwell` | Timer-dependent (minimum / range dwell-time) second-moment certificates for stochastic positive systems, the SOS-timer backend of MVST_POS_MS. |
| `mvst_pos_oc` | Linear-cost optimal control of linear positive systems by the copositive costate programs of Chapter "Linear Positive Systems" (Section sec:pos:oc). |
| `mvst_pos_stabilize_sos` | Timer-dependent copositive-SOS state-feedback synthesis for positive impulsive systems under constant / minimum / range dwell-time, the timer-dependent branch of MVST_POS_STABILIZE. |

### Simulation
| Function | Summary |
|---|---|
| `mvst_sample_ct` | Euler--Maruyama sample paths of a controlled diffusion. |
| `mvst_sample_mjls` | Sample paths of a Markov-jump controlled diffusion. |
| `mvst_simulate_ct` | Closed-loop second-moment (covariance) simulation (CT). |
| `mvst_simulate_mjls` | Closed-loop second-moment simulation for an MJLS. |
| `mvst_simulate_impulsive` | Closed-loop covariance of an impulsive system. |

### Solver dispatch
| Function | Summary |
|---|---|
| `mvst_sossolve` | Solve a SOSTOOLS program with the backend selected by OPTS.SOLVER. |
| `mvst_sosmatrvar` | Free (sign-indefinite) polynomial matrix variable. |
| `mvst_sosposmatrvar` | SOS (positive-semidefinite) polynomial matrix variable. |
| `mvst_opts` | Normalize solver options and verify the SDP/SOS solver is present. |

### Figures and export
`mvst_figdir` (in `lib/`) returns the output folder of the examples. The other
helpers live in `bookutil/` and are used by the example scripts only.

| Function | Summary |
|---|---|
| `mvst_figdir` | Output folder of the example scripts (figures and data macros). |
| `mvst_gemcolors` | MATLAB's built-in "gem" color order (with a fallback). |
| `mvst_latexdefaults` | Use the LaTeX interpreter for all new text objects. |
| `mvst_figstyle` | Force a white background and black text on a figure. |
| `mvst_savepdf` | Apply the white style and export a figure to a vector PDF. |
| `mvst_texdata` | Writer for LaTeX data macros (gains, costs, systems). |
| `mvst_fmtmat` | Format a numeric matrix as a LaTeX bmatrix string. |
| `mvst_stabfig` | Standard open-vs-closed stabilization figures. |
| `mvst_stabtable` | Emit one row of the Stabilization-chapter summary table. |
| `mvst_runscripts` | Run a list of example scripts with a fresh workspace each. |

## Result-struct conventions

Every engine returns one struct. The field names are consistent across the
toolbox, so the same conventions apply everywhere. The recurring fields are:

- **`res.feas`** (logical): the program was feasible / a certificate was found.
  Used by the synthesis and copositive-LP routines (`mvst_pos_*`,
  `mvst_stabilize`, `mvst_dissyn_sos`, ...). When `res.feas` is false the
  solution fields (`res.K`, `res.p`, ...) are empty rather than garbage.
- **`res.stable`** (logical): the analysis certificate holds. Returned by the
  stability analyses (`mvst_stability`, `mvst_pos_stability`).
- **`res.gamma`** / **`res.gammastar`**: the computed level. `res.gamma` is a
  gain or cost (L2-gain, induced positive-system gain, performance measure);
  `res.gammastar` is the critical game level (`mvst_game_gamma`).
- **`res.K`** (and `res.Kd` for a jump gain): the synthesized state-feedback
  gain. Synthesis routines also return the certificate in `res.P` (a matrix
  storage / Lyapunov certificate) or `res.lambda` / `res.p` (a copositive
  weight).
- **`res.margin`**: the achieved strict-inequality slack of a certificate (the
  smallest row slack; positive certifies the strict inequalities).
- **`res.cost`**: the optimal value of an optimal-control problem
  (`mvst_solve_*`), with the associated `res.P` / `res.P0` matrices.
- **`res.class`** / **`res.method`** / **`res.measure`**: the resolved dispatch
  tag (which system class or which measure was actually computed), useful when
  the class was inferred rather than passed explicitly.
- **`res.conv`** / **`res.residual`** / **`res.iters`**: convergence
  diagnostics for the iterative and flow-based engines (`mvst_solve_infinite`,
  `mvst_steer`).

Trajectory-valued engines add time-series fields (`res.t`, `res.Sigma`,
`res.Pi`, `res.trSigma`, `res.curve`, `res.spectrum`) alongside the scalar
summary. Some synthesis routines also return a `res.cert` / `res.closedloop` /
`res.check` sub-struct that re-certifies the designed closed loop.

## Error-id namespace

All errors raised by the toolbox carry a structured identifier of the form
`mvst:<area>:<reason>`; catch them by identifier, not by message text. The
`<area>` segment mirrors the function family:

| Area prefix | Raised by | Typical reason |
|---|---|---|
| `mvst:System:` | the `mvst.System` factories | dimension mismatch, bad generator |
| `mvst:opts:` / `mvst:mode:` | option and model builders | missing / malformed field |
| `mvst:pos:` | the positive-systems LP family | missing data, bad class/norm, no Optimization Toolbox |
| `mvst:perf:` | `mvst_perf` | unsupported class-measure pair, missing I/O |
| `mvst:stability:` / `mvst:stabilize:` | stability / stabilization | bad class, coupling, dwell-time data |
| `mvst:l2gain:` / `mvst:game_gamma:` | gain / game routers | missing per-class data |
| `mvst:steer:` | steering | bad target, horizon, or class |
| `mvst:tests:` | `run_all_tests` | a test truly failed (skips are not failures) |

Two identifiers are worth knowing specifically:
`mvst:pos:noOptimToolbox` (the positive-systems LPs need `linprog`) and
`mvst:perf:unsupported` (a performance measure the book does not characterize for
that class; the message names the nearest supported pair).

## Troubleshooting

**A solver is missing.** The SOS layer (`mvst_*_sos`) needs SOSTOOLS and an SDP
backend (SeDuMi / SDPT3 / MOSEK); the positive-systems LPs need the Optimization
Toolbox. Run `mvst.solvercheck` to see what resolves and actually solves on your
path. `mvst_setup` adds the bundled solvers of `thirdparty/`; SeDuMi and SDPT3 may
need `install_sedumi` / `install_sdpt3` once, and `mvst.Options('solver','nosdp')`
runs the SOSTOOLS-based engines without any SDP solver (see
[Getting started](01-getting-started.md)). A missing LP
solver surfaces as `mvst:pos:noOptimToolbox`; a missing SDP backend surfaces from
the SOS routine's solver call.

**An analysis returns `stable=0` or a synthesis returns `feas=0`.** This is a
result, not an error: the certificate does not exist for the given data. Check
that the model matches the class (a Metzler drift for the positive family, a
Hurwitz drift for the continuous-time Gramian measures), and that any dwell-time
is on the feasible side of its threshold. The dwell-time routines expose
`opts.threshold` to bisect the feasibility boundary and return it in
`res.Tstar`, which tells you how far the data is from feasibility.

**Tolerance knobs.** The copositive LPs take `opts.margin` (the strict-inequality
slack) and `opts.pmin` (the lower bound normalizing the positive certificate);
`mvst_pos_stabilize` takes `opts.eps` (the copositive decrease margin, raise it
for a faster closed-loop decay). The Gramian engines take ODE
tolerances (`opts.reltol`, `opts.abstol`) and quadrature counts (`opts.npts`,
`opts.profpts`).

**An `mvst.System` call errors on construction.** The factories validate
dimensions eagerly and raise `mvst:System:dim` (or `:generator`) with the
offending block named; fix the block sizes rather than the call site. `help
mvst.System.<factory>` lists the expected shapes.

