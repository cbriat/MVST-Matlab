# MVST companion code: conventions and release surface

This document is the release-facing reference for the companion MATLAB code of
the monograph.  It subsumes the older `MODEL_CONVENTION.md` (kept for
historical continuity) and additionally describes the recommended object entry
point, the problem-to-solver dispatcher map, the error-identifier namespace,
the solver-options contract, and the figure/export layer.

The code is layered:

* `+mvst/` : the release-facing object API (the `mvst.System` class, the
  `mvst.Options` class, `mvst.solvercheck`).  This is the recommended way in.
* `lib/mvst_*.m` : the functional core (model builders, numerical engines,
  SOS/SDP programs, dispatchers) plus a self-contained figure/export layer.
* `examples/`, `run_*.m` : book-reproduction scripts and per-chapter runners.
* `tests/` : the golden regression harness (`run_regress`) and the unit-test
  suites (`tsystem`, `tstability`, `tperf`, `tpos`), all `classdef`
  `matlab.unittest.TestCase` classes auto-discovered by `run_all_tests` (each
  reports per-method pass/fail counts).

---

## 1. Model conventions: thin, rich, and the recommended entry

Two struct conventions coexist and share ONE field naming, differing only in
the drift/input pair.

* **thin** (`mvst_mode` / `mvst_system`): the OC / Riccati / game / steering
  layer.  Drift `A`, input `B`, quadratic-cost blocks `Q,R,S,W,X0,ZT`,
  Brownian noise `Ac/Bc`, Markov coupling `Pi/pi0`.
* **rich** (the SOS layer): drift `A0`, input `B0`, and the disturbance /
  performance-output / jump channels the supply-rate machinery needs.

The genuine difference is only `A`/`B` (thin) versus `A0`/`B0` (rich).  All
other shared fields carry the same name in both.

### Recommended entry: the `mvst.System` class

`mvst.System` is the recommended way to build a model.  Static factories cover
every class and validate their inputs (`mvst:System:*`):

    S = mvst.System.ct(A, B, 'Q',Q, 'R',R);        % continuous-time
    S = mvst.System.mjls(Alist, Blist, Pi);         % Markov-jump
    S = mvst.System.sampled(A, B, J, T);            % sampled-data
    S = mvst.System.poisson(A, B, J, lambda);       % spontaneous jumps
    S = mvst.System.switched(Alist, Blist);         % switched / dwell-time
    S = mvst.System.impulsive(A, B, J);             % impulsive

Bridges connect the object to the functional core:

* `mvst.System.fromstruct(sys)` : adopt an existing thin OR rich struct.
* `S.thin()` : emit the thin struct the Riccati/OC engines consume.
* `S.rich()` : emit the rich struct the SOS engines consume (via
  `mvst_torich`).

The class also carries method sugar that forwards to the dispatchers of
Section 3 (`S.l2gain`, `S.steer`, `S.gamegamma`, `S.stability`,
`S.stabilize`, `S.dissyn`, `S.pos`, `S.perf`, `S.ocfinite/ocinfinite/
ocaverage`).

### Field mapping (thin <-> rich)

| Meaning                                   | Thin (`mvst_mode`/`mvst_system`) | Rich (SOS layer)          |
|-------------------------------------------|----------------------------------|---------------------------|
| flow drift                                | `A`                              | `A0`                      |
| flow control input                        | `B`                              | `B0`                      |
| flow state-multiplicative noise channels  | `Ac`   (n x n x Nc)              | `Ac`   (n x n x q)        |
| flow input-multiplicative noise channels  | `Bc`   (n x m x Nc)              | `Bc`   (n x m x q)        |
| scheduled (dwell-time) jump / reset maps  | `J`  *(ext)* (n x n [x p])       | `J`    (n x n [x p])      |
| scheduled-jump input selector             | `BJ` *(ext)* (n x m [x p])       | `BJ`   (n x m [x p])      |
| spontaneous (Poisson) jump maps           | `Jf` *(ext)* (n x n x pf)        | `Jf`   (n x n x pf)       |
| Poisson rates                             | `lamf` *(ext)* (1 x pf)          | `lamf` (1 x pf)           |
| spontaneous-jump input selector           | `BJf` *(ext)* (n x m x pf)       | `BJf`  (n x m x pf)       |
| CT Markov generator (multi-mode)          | `Pi` (`mvst_system`)             | `Pi` (stab), `Pig` (diss) |
| DT Markov transition (multi-mode)         | `P`  *(ext, `mvst_system`)*      | `P`                       |
| switch / impulse reset (multi-mode)       | `J`  *(ext, `mvst_system`)*      | `J`                       |
| quadratic cost blocks                     | `Q,R,S,W,X0,ZT`                  | (none)                    |
| disturbance input                         | (none)                           | `Bw`                      |
| disturbance-multiplicative noise          | (none)                           | `Dwc`                     |
| performance output / feedthrough (flow)   | (none)                           | `C`, `Dw`                 |
| Poisson jump-disturbance channels         | (none)                           | `Bdf`                     |
| jump disturbance / output / feedthrough   | (none)                           | `Bd`, `CJ`, `DJ`          |

*(ext)* marks the optional extension fields of `mvst_mode` / `mvst_system`.
Every extension defaults to **absent**: a struct built without those arguments
is byte-identical to the pre-extension output, so the Riccati solvers (which
read fields by name) are unaffected.  Multi-mode rich structs stack per-mode
data along the last dimension (`A0` n x n x M, `Ac` n x n x q x M, etc.).

### Bridge helpers

* `mvst_torich(sys)` : thin (+extensions) -> rich.  Idempotent (a struct with
  `.A0` is returned unchanged) and permissive, so SOS entry points call it
  unconditionally; cost data is dropped (no SOS counterpart).
* `mvst_tothin(rich, 'Q',Q, 'R',R, ...)` : rich -> thin, the inverse on the
  shared fields, with cost data supplied as name/value pairs.  Fields with no
  thin home (mode-indexed data, dissipativity I/O, `Pig`/`P`) raise
  `mvst:tothin:unrepresentable`.

Each SOS entry point begins with `if ~isfield(sys,'A0'); sys =
mvst_torich(sys); end`, so a thin struct is accepted transparently while the
rich fast path is untouched.

---

## 2. The mode/plumbing builders (`lib`)

* `mvst_mode(A,B,Q,R,...)` : one controlled diffusion mode (thin).  Validates
  A square, B row count, Q/R square (`mvst:mode:*`).
* `mvst_mode_tv(n,m,...)` : a time-varying mode; fields may be `@(t)` handles,
  evaluated through `mvst_evalc`.
* `mvst_system(modes,...)` : assemble a (multi-mode) system; validates the mode
  cell, generator size, and initial distribution (`mvst:system:*`).
* `mvst_getf`, `mvst_evalc` : permissive struct-field access and time-varying
  coefficient evaluation.
* `mvst_svecbasis`, `mvst_liftsym`, `mvst_vech`, `mvst_unvech` : the symmetric
  subspace algebra (reduced lifts used by the exact stability engines).  These
  validate dimensions (`mvst:svecbasis:*`, `mvst:liftsym:*`, `mvst:vech:*`,
  `mvst:unvech:*`).
* `mvst_odecap` : shared terminal ODE event for the backward Riccati
  integrators.
* `mvst_sample_*`, `mvst_simulate_*` : sample-path and
  second-moment (covariance) simulators used to draw the book figures.

---

## 3. Dispatcher map: which function solves which problem

Each generic function normalizes the system and routes to the matching
per-class backend, returning a uniform result struct.  "exact" = closed-form /
Riccati / spectral engine; "SOS" = SOSTOOLS/LMI relaxation (needs an SDP
solver, see Section 5).

| Problem                                 | Generic function      | Classes handled                                             | Kind  |
|-----------------------------------------|-----------------------|-------------------------------------------------------------|-------|
| Mean-square stability analysis          | `mvst_stability`      | ct, dt, mjls, dwell/switched, sampled-data, impulsive, poisson | exact |
| State-feedback stabilization            | `mvst_stabilize`      | same classes as `mvst_stability`                            | exact |
| Stabilization synthesis (relaxed)       | `mvst_stabz_sos`      | rich struct, all classes in one program                     | SOS   |
| L2-gain                                 | `mvst_l2gain`         | ct, mjls, dwell/switched, poisson, sampled-data/hybrid      | exact |
| Dissipativity analysis                  | `mvst_diss_sos` (`_modes`) | rich struct, supply rate + dwell/mode options          | SOS   |
| Dissipative state-feedback synthesis    | `mvst_dissyn_sos`     | rich struct, supply rate + dwell/mode options               | SOS   |
| H-infinity / minimax game level         | `mvst_game_gamma`     | ct, mjls, dwell, poisson (gamma-bisection Riccati)          | exact |
| Zero-sum game synthesis (relaxed)       | `mvst_game_sos`       | rich struct, all classes in one program                     | SOS   |
| Finite-time optimal covariance steering | `mvst_steer`          | ct, mjls, sampled-data, switched, poisson (+ time-varying)  | exact |
| Further performance measures            | `mvst_perf`           | H2 / peak / energy-to-peak per class (see its help)         | exact |
| Optimal control (dual storage SOS)      | `mvst_oc_sos` (`_modes`) | rich struct, quadratic cost + dwell/mode options         | SOS   |
| Positive-systems certificates           | `mvst_pos_*`          | copositive LP: stability, stabilize, gain, supply, lp       | LP    |

Optimal-control **Riccati solvers** (thin structs), driven by the mode/system
builders: `mvst_solve_finite`, `mvst_solve_infinite`, `mvst_solve_average`
(and the `_tv` time-varying variants), themselves built on the per-class
Riccati/differential-Riccati engines (`mvst_ric_ct`, `mvst_dre_ct`,
`mvst_care_ct`, `mvst_igre_periodic_sd`, and their poisson/tv variants).

The `mvst_pos_*` family: `mvst_pos_stability` (copositive LP stability),
`mvst_pos_stabilize` (positivity-preserving synthesis), `mvst_pos_gain`
(induced-gain LPs), `mvst_pos_supply` (general-linear-supply dissipativity),
with `mvst_pos_lp` the internal LINPROG wrapper. The object front end is the
`mvst.System.pos(what, opts)` sugar method, which routes to this family.

**SOS entry-point arity.**  The SOS routines are NOT uniform in signature.
Most take a per-program argument BEFORE the dwell/mode pair:
`mvst_diss_sos(sys, supply, dwell, mode, opts)`,
`mvst_game_sos(sys, game, dwell, mode, opts)`,
`mvst_oc_sos(sys, cost, dwell, mode, opts)`.  The dynamic-output-feedback
synthesis `mvst_dof_sos(sys, dwell, opts)` deliberately has a SHORTER arity: it
carries no separate supply/cost/game argument and no `mode` selector (the
program is fixed).  This is an intentional inconsistency; do not "normalize" its
signature to match the others.

---

## 4. Error-identifier namespace

All errors raised by the package use the identifier form

    mvst:<name>:<reason>

where `<name>` is the (short) function or class name and `<reason>` is a
lowerCamel tag.  The colon namespace is uniform: there is NO `mvst_<name>:...`
underscore-prefixed variant anywhere in the tree.  Examples:
`mvst:mode:notSquare`, `mvst:system:generator`, `mvst:liftsym:dim`,
`mvst:opts:solverMissing`, `mvst:perf:unsupported`, `mvst:pos:solver`,
`mvst:System:dim` (the class uses its CamelCase name).

Input-validation errors are additive: they only fire on inputs that would
otherwise crash cryptically (a non-struct, an empty `struct()`, a missing
required field), and never on a currently-working call.  Raw-input guards use
the reason tag `input`, e.g. `mvst:steer:input`, `mvst:ric_ct:input`,
`mvst:game_sos:input`, `mvst:fmtmat:input`, `mvst:opts:input` (a char/string
`OPTS`, including `''`, is rejected: only `[]` or an omitted argument requests
the defaults).  Tests assert specific identifiers (see `tsystem`, `tperf`,
`tpos`), so treat the identifier strings as a stable contract.

---

## 5. Solver-options contract

The SDP/SOS engines take a trailing options argument, normalized by
`mvst_opts`:

* accepts `[]` (defaults), a plain struct, or an `mvst.Options` object;
* returns a plain struct with at least `.solver` (default `'sedumi'`);
* validates the requested backend against the MATLAB path and errors
  `mvst:opts:solverMissing` with NO silent fallback;
* `mvst_opts(opts, false)` skips the path probe (for routers whose backends
  are purely numerical), still validating the solver name.

Recognized solvers and their probe entry points:

    'sedumi' -> sedumi ,   'sdpt3' -> sdpt3 ,   'mosek' -> mosekopt

`mvst.Options` fixes the release-facing surface (`solver`, `reltol`,
`verbosity`) with validation, and threads the `.solver` choice through every
SDP/SOS engine.  SOSTOOLS receives the choice via `sossolve`'s `.solver`
field; the legacy `engine = 'yalmip'` paths (YALMIP is an optional external
dependency, not bundled) via `sdpsettings('solver', ...)`.  `mvst.solvercheck`
reports which backends are installed.

---

## 6. Figure and export layer

This layer reproduces the book's figures and data macros and is deliberately
separated from the numerical core: nothing in Sections 1-5 depends on it.

* `mvst_gemcolors`, `mvst_latexdefaults`, `mvst_figstyle` : the MATLAB "gem"
  palette, the LaTeX-interpreter defaults, and the white-background / black-text
  styling applied before export.
* `mvst_savepdf` : style + export a figure to a vector PDF (validates the
  target path, `mvst:savepdf:*`).
* `mvst_texdata`, `mvst_fmtmat` : write LaTeX `\def` data macros (gains, costs,
  matrices) for `\input` into the manuscript (`mvst:texdata:*`).
* `mvst_stabfig`, `mvst_stabtable` : the standard
  Stabilization-chapter figures, the generated summary-table rows
  (`mvst:stabtable:*`).
* `mvst_runscripts` and the `run_*.m` runners : run the example scripts with a
  clean workspace each; vector PDFs go to `figures/` under the current folder (or to the folder
  named by the environment variable `MVST_FIGDIR`, see `mvst_figdir`), data
  macros to its `data/` subfolder.

Reproduction: build the book first, then `run_all` (or a per-chapter runner).
The golden regression harness `tests/run_regress` compares the emitted data
macros against `tests/golden/data`.
