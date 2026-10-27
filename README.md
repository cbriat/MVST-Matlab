# Matrix-Valued Systems Toolbox (MVST): companion solver library

The Matrix-Valued Systems Toolbox (MVST) is the companion MATLAB toolbox for the monograph **_Linear Hermitian
Matrix-Valued Dynamical Systems_** by Corentin Briat. It is a reusable solver
library for the analysis and synthesis of linear matrix-valued (impulsive and
stochastic) dynamical systems, covering the book's chapters end to end:

- **stability and stabilization** (exact mean-square analysis, state-feedback
  synthesis, dwell-time-constrained design),
- **linear-cost optimal control** over the three horizons (finite, infinite,
  long-term average),
- **finite-time optimal steering** of the state covariance,
- **dissipativity** (supply rates, KYP, L2 gain) and **dissipativity-based
  synthesis**,
- **zero-sum differential games** (H-infinity level),
- **further performance measures** (H2, energy-to-peak, peak-to-peak),
- **linear positive systems** (copositive-LP analysis and synthesis),
- **higher-order systems**.

- **routed H-infinity norm and passivity** (`S.hinfnorm()`, `S.passivity()`):
  uniform result structs with certificates and indices; exact imaginary-axis
  Hamiltonian tests for CT and DT LTI systems, the DT and DT-MJLS classes, a
  converged (adaptive) sampled-data gain, and the book's single input-strict
  passivity index convention (supply `He(w'z) - nu|w|^2`),
- **no SDP solver required for the SOS engines**: `mvst.Options('solver','nosdp')`
  solves the SOS programs through SOSNOSDP (see Requirements; only the legacy
  `opts.engine = 'yalmip'` paths still need an SDP solver).

The system classes handled are the continuous-time diffusion, the Markov-jump
linear system (MJLS), the sampled-data / impulsive system, the Poissonian-jump
diffusion, and the dwell-time switched system.

The release-facing surface is the `+mvst` package (`mvst.System`,
`mvst.Options`, `mvst.solvercheck`); every legacy engine remains available as a
`mvst_*` function in `lib/`.

## Documentation

The full user manual (author-ordered, one chapter per problem, with runnable
examples) is in [`docs/manual/index.md`](docs/manual/index.md). A typeset PDF of
the same manual is at [`docs/latex/mvst_manual.pdf`](docs/latex/mvst_manual.pdf).

## Quick start

From the `matlab/` directory, put the library on the path with the setup
helper (it adds `lib/`, `bookutil/`, `tests/`, the `+mvst` package parent, and the bundled
solvers under `thirdparty/`, and prints the solver doctor report):

```matlab
run('mvst_setup.m')          % or simply:  mvst_setup
```

`mvst_setup` accepts `'examples', true` to also add `examples/`, and
`'quiet', true` to skip the doctor report. The toolbox version is available as
`mvst.version`.

Then build a continuous-time system and compute its L2 gain and stability:

```matlab
A = [-1 0.5; 0 -2];
S = mvst.System.ct(A, zeros(2,0), 'Bw', [1;0], 'C', [1 0]);  % dx = A x + Bw w, z = C x
g  = S.l2gain().gamma;      % exact (mean-square) L2 gain
ok = S.stability().stable;  % exact mean-square stability
fprintf('L2 gain = %.4f, stable = %d\n', g, ok);   % -> L2 gain = 1.0000, stable = 1

h = S.hinfnorm();           % uniform struct: gamma, certificate, class, method, diagnostics
P = mvst.System.ct([0 1;-1 -1], zeros(2,0), 'Bw', [0;1], 'C', [1 0], 'Dw', 0.5);
p = P.passivity();          % exact test: p.passive, p.nu_isp, p.rho_osp, p.method
```

The main analysis and synthesis routines are exposed as one-line methods of
`mvst.System` (`l2gain`, `steer`, `gamegamma`, `ocfinite`, `ocinfinite`,
`ocaverage`, `stability`, `stabilize`, `dissyn`, `pos`, `perf`, `hinfnorm`,
`passivity`, `kyp`, `dissipativity`), each delegating to the matching backend
of the resolved system class; every `opts` argument takes a plain struct or an
`mvst.Options`. `S.dissipativity(supply)` takes `'l2gain'`, `'passivity'` or the
quadratic supply `struct('type','qsr','Q',Q,'S',S,'R',R)`, meaning
`[z;w]'*[Q S;S' R]*[z;w]`. The remaining engines (the SOS programs `mvst_*_sos`, the game variants,
dwell-time and positive-systems backends, ...) are called directly as `mvst_*`
functions; `Contents.m` (`help` on the toolbox folder) lists all of them.

## Requirements

- **MATLAB R2020a or later.** The code uses `arguments` validation blocks
  (R2019b) and `exportgraphics` (R2020a, figure export in `mvst_savepdf`). It is
  developed and tested on **R2025a**; older releases are not tested. GNU Octave
  is **not** supported (the `arguments` blocks and the `+mvst` classes are
  MATLAB-only).
- **Base MATLAB** suffices for the Riccati-based core: stability, optimal
  control, covariance steering, games, and the L2 gain.
- **Optimization Toolbox**: the positive-systems copositive linear programs
  (`mvst_pos_*`, via `linprog`).
- **SOSTOOLS and an SDP solver**: the sum-of-squares engines (`mvst_*_sos`,
  dissipativity, games, stabilization, and optimal control by SOS). SeDuMi is
  the default backend; SDPT3 and MOSEK are alternatives, selected through
  `mvst.Options('solver', ...)`.
- **YALMIP** (optional, not bundled): needed only by the legacy `opts.engine = 'yalmip'`
  paths. The LMIs of `mvst_kyp`, `mvst_pos_ms`, `mvst_stab_delay` and the decay-rate
  branch of `mvst_perf_dwell` are SOSTOOLS programs by default (`opts.engine = 'sostools'`,
  every backend including `'nosdp'`); their former YALMIP formulation remains available
  for one release as `opts.engine = 'yalmip'` and fails with `mvst:opts:yalmipMissing` when
  YALMIP is not installed. Install it yourself (https://yalmip.github.io/); `mvst_setup`
  and `mvst.solvercheck` detect an installation already on the path and report it as optional.
- **Control System Toolbox**: not required (the Lyapunov solves in `mvst_perf` and
  `mvst_kyp` use base MATLAB: `sylvester` and a Kronecker solve).
- **sosnosdp** (optional): required for `mvst.Options('solver','nosdp')`, which
  solves the SOS programs (the `mvst_*_sos` engines, `mvst_stabz_sd_range`, `mvst_kyp`,
  `mvst_isp_margin_ct`, `mvst_stab_delay`, `mvst_pos_ms`, and the timer-varying and decay-rate
  branches of `mvst_perf_dwell`) without any SDP solver, see https://github.com/cbriat/sosnosdp. SOSTOOLS is
  still needed; the legacy YALMIP engines reject `'nosdp'`
  (`mvst:opts:nosdpYalmip`).

  Limits of `'nosdp'` (a different algorithm, not an interior-point SDP solver):
  accuracy is lower (bisected gains agree with SeDuMi to about 1e-3 to 1e-2 relative),
  quantities recovered from a non-unique optimum (for example the synthesis gains of
  `mvst_stabz_sd_range`) can differ from the SeDuMi ones by a few percent, programs with
  large Gram blocks are much slower (the SDP-range synthesis suite takes minutes instead
  of seconds), and SOSNOSDP's certification step may print benign
  `Matrix is close to singular or badly scaled` warnings near the feasibility boundary
  (the certificate is verified by explicit Gram-matrix residual and eigenvalue tests, so
  the warning does not affect soundness). A SOSNOSDP optimum that stays uncertified (status
  `uncertified_*`) is returned by the synthesis engines with `info.certified = false` and a
  warning (`mvst:uncertified`); the classical backends treat such a solve as failed.

### Bundled third-party solvers

SOSTOOLS, SeDuMi, SDPT3, and sosnosdp are shipped in `thirdparty/`,
each as an unmodified snapshot of the upstream release in its own subfolder
`thirdparty/<solver>/`, together with its own license file. `mvst_setup` adds
whatever is present there, so installation is one step. To use a newer version
or your own installation, replace or delete the subfolder: `mvst_setup` and
`mvst.solvercheck` discover what is available. The component list, versions,
and licenses are in [`THIRD-PARTY.md`](THIRD-PARTY.md) and
[`thirdparty/VERSIONS.md`](thirdparty/VERSIONS.md). SDPT3 ships without compiled
MEX files: run `install_sdpt3` once (from `thirdparty/sdpt3/`) before selecting it.
SeDuMi runs from its MATLAB source, but fast solving needs its platform-specific
MEX files (not tracked in git): run `install_sedumi` once (from `thirdparty/sedumi/`).
`mvst.solvercheck` prints this hint when a solver is found but does not solve.

YALMIP is not bundled. If you use it (legacy engine paths), please cite it (J. Löfberg, "YALMIP: A toolbox for modeling
and optimization in MATLAB", CACSD 2004), as its license requests.

- **Solver doctor.** Run `mvst.solvercheck` at any time to see which backends
  resolve and actually solve a trivial SDP:

  ```matlab
  mvst.solvercheck
  ```

## Layout

| Path              | Contents                                                        |
|-------------------|-----------------------------------------------------------------|
| `+mvst/`          | Release package: `System`, `Options`, `solvercheck`.            |
| `lib/`            | Reusable engine layer, the `mvst_*` functions (Riccati, SOS, positive-systems, simulation). |
| `bookutil/`       | Figure style, PDF export, and table helpers used by the examples. |
| `thirdparty/`     | Bundled upstream solver snapshots (see `THIRD-PARTY.md`).       |
| `docs/`           | User manual (Markdown and typeset PDF).                         |
| `examples/`       | Example systems and figure scripts (one per book example).      |
| `tests/`          | Unit-test suites (the release battery). |
| `mvst_setup.m`    | Path setup + solver doctor (top level).                         |

## Figure output

When an example script is run, vector PDFs (white background, black text, so exports are correct even in dark
mode) are written to `figures/` under the current folder (or to the folder named
by the environment variable `MVST_FIGDIR`, see `mvst_figdir`), and the LaTeX data
macros the book `\input`s are written to its `data/` subfolder. The SOS/SDP examples require an SDP solver on
the path; the rest need only base MATLAB.

## Testing

One consolidated entry point runs every unit-test suite, prints a
PASS/FAIL/SKIPPED summary, and exits nonzero only if a test fails (tests that
need a solver not on the path are reported as SKIPPED):

```matlab
cd tests; run_all_tests
```

or from the shell:

```
matlab -batch "cd tests; run_all_tests"
```

The battery has 26 classdef `matlab.unittest` suites (about 250 tests; a full run
takes a few minutes), all auto-discovered so a new `t*.m` `TestCase` is
picked up without editing the harness. The development tree additionally
carries a golden-data regression harness tied to the book's figure data; it is
not part of the release battery. Individual suites can be run with `runtests('tpos')` etc.

## License

MIT (see [`LICENSE`](LICENSE)). The bundled solvers in `thirdparty/` carry
their own licenses (SOSTOOLS GPLv3, SeDuMi GPLv2, SDPT3
GPLv2, sosnosdp MIT); the distribution is an aggregate of separately licensed
components, listed in [`THIRD-PARTY.md`](THIRD-PARTY.md).

## Citation

See [`CITATION.cff`](CITATION.cff). Release history: [`CHANGELOG.md`](CHANGELOG.md).
