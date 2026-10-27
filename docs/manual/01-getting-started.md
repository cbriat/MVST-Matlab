# 1. Getting started

[Manual index](index.md) | Next: [Systems](02-systems.md)

## Installation

The toolbox is plain MATLAB source: there is nothing to compile. The solvers it
needs (SOSTOOLS, SeDuMi, SDPT3, SOSNOSDP) ship with it, each as an
unmodified upstream snapshot in its own subfolder of `thirdparty/`. Copy or
clone the repository, start MATLAB in the `matlab/` directory, and run the
setup helper once per session:

```matlab
mvst_setup
```

`mvst_setup` puts the following on the path:

- the `+mvst` package parent (the `matlab/` folder itself);
- the engine layer `lib/`;
- the plotting and export helpers `bookutil/` (figure style, PDF export);
- the test battery `tests/`;
- every bundled solver folder found under `thirdparty/`
  (`sostools`, `sedumi`, `sdpt3`, `sosnosdp`).

It then runs the solver doctor (next section). The setup prints one line per
solver, for example:

```
MVST companion library added to the path.
  root : .../matlab
  lib  : .../matlab/lib
  sostools : added from .../matlab/thirdparty/sostools
  sedumi   : added from .../matlab/thirdparty/sedumi
  sdpt3    : added from .../matlab/thirdparty/sdpt3
  sosnosdp : added from .../matlab/thirdparty/sosnosdp
```

Options:

- `mvst_setup('examples', true)` also adds `examples/` to the path;
- `mvst_setup('quiet', true)` skips the doctor report (path only);
- `mvst_setup('thirdparty', false)` does not add the bundled solvers, so that
  your own installations are used.

### Using your own solver installation

The bundled folders are independent and replaceable:

- A solver that is already on the MATLAB path (found through its entry point,
  for instance `sosprogram` or `sedumi`) is respected: `mvst_setup` prints
  `already on the path, bundled copy not added` and leaves the bundled copy out.
- A subfolder of `thirdparty/` can be deleted or replaced by a newer upstream
  snapshot; a missing subfolder is skipped silently.
- The components, versions, and licenses are listed in `THIRD-PARTY.md` and
  `thirdparty/VERSIONS.md`.

### Compiling the SDP solvers once

SeDuMi runs from its MATLAB source, but fast solving needs its platform MEX
files, and SDPT3 ships without any compiled MEX file. Run these once, from the
solver's own folder:

```matlab
cd thirdparty/sedumi; install_sedumi
cd ../sdpt3;          install_sdpt3
```

The doctor prints a hint (`SDPT3 found but does not solve: run install_sdpt3
once`) whenever a solver is found but fails its probe. The exact engines never
need either solver.

### MATLAB version

The floor is **R2020a**: the code uses `arguments` blocks (R2019b) and
`exportgraphics` (R2020a, used for figure export). It is developed and tested on
**R2025a**; older releases are not tested, and GNU Octave is not supported. The
exact engines (Riccati, stability, steering, games, L2 gain) need only base
MATLAB: no Control System Toolbox is required (the Lyapunov solves use
`sylvester` and a Kronecker solve).

## The solver doctor

`mvst.solvercheck` detects and probes the SDP / SOS toolchain. It has six rows:
SeDuMi, SDPT3, MOSEK, SOSTOOLS, SOSNOSDP, and YALMIP (optional, not bundled). After `mvst_setup` on a
machine with the bundled solvers (SDPT3 not yet compiled) it reports:

```
MVST solver / SOS toolchain check
------------------------------------------------------------------------
  tool       found   solves   version      location
------------------------------------------------------------------------
  SeDuMi     yes     yes      1.3.7        .../thirdparty/sedumi/sedumi.m
  SDPT3      yes     no       4.0          .../thirdparty/sdpt3/sdpt3.m
  MOSEK      no      ?        -            (not on path)
  SOSTOOLS   yes     -        4.00         .../thirdparty/sostools/sosprogram.m
  SOSNOSDP   yes     -        0.1.0        .../thirdparty/sosnosdp/bridge/sossolve_nosdp.m
  YALMIP     no      -        -            (not on path; optional, legacy engines only)
------------------------------------------------------------------------
  4 of 6 tools detected; 1 SDP solver(s) verified to solve.
  hint: SDPT3 found but does not solve: run install_sdpt3 once (compiles its MEX files).
  solves column:  yes = trivial 2x2 SDP solved,  no = failed,
                  ? = not probed (no SOSTOOLS),  - = not an SDP solver.
```

The columns:

- **found**: whether the tool resolves on the MATLAB path.
- **solves**: whether a trivial 2x2 SDP was actually solved by that backend.
  The probe is a trivial SOSTOOLS program solved by that backend (YALMIP is not used). It is `?` when SOSTOOLS is absent, `yes`/`no` once probed,
  and `-` for the tools that are not themselves SDP solvers (SOSTOOLS, YALMIP,
  SOSNOSDP).
- **version** / **location**: what was found and where.

With `mvst_setup('thirdparty', false)` and nothing else on the path, all six
rows read `no`; the exact engines still work, only the SOS layer of the
Stability, Stabilization, Optimal Control, Dissipativity, and Game chapters
needs a backend.

## External solvers

The SOS layer requires **SOSTOOLS** (the SOS front end) and one backend,
selected through `mvst.Options('solver', ...)` or the `solver` field of a plain
options struct:

| `solver` | Backend | Needs |
|---|---|---|
| `'sedumi'` (default) | SeDuMi | bundled; run `install_sedumi` once for speed |
| `'sdpt3'` | SDPT3 | bundled; run `install_sdpt3` once |
| `'mosek'` | MOSEK | your own MOSEK installation |
| `'nosdp'` | SOSNOSDP | bundled; no SDP solver at all |

The `'nosdp'` option solves the SOS programs with SOSNOSDP (no semidefinite
programming). It applies to the SOSTOOLS-based engines (`mvst_*_sos` and the
engines built on them), which reach the backend through the single dispatcher
`mvst_sossolve`. The constant-certificate LMIs of `mvst_pos_ms`,
`mvst_stab_delay`, and the decay-rate branch of `mvst_perf_dwell` are also
SOSTOOLS programs by default (`opts.engine = 'sostools'`), so they accept
`'nosdp'` too; their former **YALMIP** formulation is kept for one release as
`opts.engine = 'yalmip'` (cross-checks only; YALMIP is an optional external
dependency, not bundled, and the call fails with `mvst:opts:yalmipMissing` when it is
absent). `mvst_kyp` is a SOSTOOLS program by default as well. The
**Optimization Toolbox** (`linprog`) is required for the positive-systems
copositive LPs.

```matlab
sys   = struct('A0',[0 1;-1 0], 'J', 0.5*eye(2));     % example of the Stability chapter
dwell = struct('type','periodic','Tbar',0.5,'deg',4);
o = mvst.Options('solver', 'nosdp');       % no SDP solver involved
r = mvst_stab_sos(sys, dwell, o);          % same call, other backend
fprintf('feas=%d  slack=%.3e\n', r.feas, r.slack);
```
```
feas=1  slack=5.000e-01
```

If you request a backend that is not installed, the options normalizer refuses
with no silent fallback:

```matlab
mvst_opts(struct('solver','sdpt3'));       % with SDPT3 absent
```
```
mvst:opts:solverMissing : Requested solver 'sdpt3' is not on the MATLAB path
(probe 'sdpt3' failed).  Install it or choose another solver; no silent
fallback is performed.
```

YALMIP is not bundled. If you use it in published work (legacy engine paths), please cite it (J. Lofberg, "YALMIP: A
toolbox for modeling and optimization in MATLAB", CACSD 2004), as its license
requests.

## Running the test battery

One entry point runs every unit-test suite and prints a PASS / FAIL / SKIPPED
summary:

```matlab
cd tests; run_all_tests
```

or from the shell:

```
matlab -batch "cd tests; run_all_tests"
```

Tests that need a solver which is not on the path are not failures: they are
reported as `SKIPPED` lines and counted separately. The exit code is nonzero
only if a test truly fails. Without `mvst_setup` (library only) the battery
therefore reports many skipped tests and no failure; after `mvst_setup` the
solver-dependent tests run. A full run takes about 10 minutes; the tail of the
report looks like:

```
---------------------------- summary ----------------------------
  suite                           pass   fail   status
  -------------------------------------------------------------
  tperf                              8      0   PASS
  tpos                              14      0   PASS
  tstability                        13      0   PASS
  tsystem                           31      0   PASS
  ...
  -------------------------------------------------------------
  TOTAL: <passed> passed, 0 failed, <skipped> skipped
-----------------------------------------------------------------
  ALL GREEN.
```

The suites are auto-discovered `matlab.unittest` classes, so a new `t*.m`
`TestCase` is picked up without editing the harness. Individual suites run with
`runtests('tpos')` and so on.
