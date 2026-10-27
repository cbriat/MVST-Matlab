# 12. Book reproduction

[Manual index](index.md) | Prev: [Further performance measures](11-performance.md) | Next: [API reference](13-api-reference.md)

This chapter explains how the toolbox regenerates every figure and every
printed number of the companion monograph, where the outputs are written, and
(for the development tree) how a golden regression harness freezes those
numbers.

## The figure and data pipeline

The reproduction layer rests on one idea: the numbers the book prints and the
numbers the code computes are the **same** numbers, carried across by generated
files. There are three stages.

1. An **example script** (`examples/ex_*.m`, `examples/fig_*.m`, or
   `examples/comp_*.m`) sets up a system, calls the analysis or synthesis
   engines, and produces two kinds of artifact.
2. Every plotted figure is exported as a vector PDF to the output folder.
3. Every scalar, matrix, or string the book quotes in its running text is
   written as a LaTeX macro to the `data/` subfolder of the output folder.

The output folder is returned by `mvst_figdir`, which is the same for every
example. In order of precedence it is the environment variable `MVST_FIGDIR`
when set and nonempty, and otherwise `./figures` under the current folder. The
folder and its `data/` subfolder are created when missing, so the examples run
standalone from any clone:

```matlab
setenv('MVST_FIGDIR', fullfile(tempdir, 'mvst_out'));
mvst_figdir('data')        % -> .../mvst_out/data
```

The book manuscript then `\input`s the `data/*.tex` files and refers to
each quantity by macro name. A sentence such as "the critical level is
$\gamma^\star = \GameGammaSd$" resolves, at LaTeX compile time, to the value the
MATLAB run wrote. No number is transcribed by hand, so the text cannot drift
away from the code.

Schematically, for one example:

```
ex_game_sd.m  --->  <figdir>/ex_game_sd.pdf          (a figure in the book)
              --->  <figdir>/data/game_sd.tex        (\def\GameGammaSd{...})
                         |
                         v
              Book_Springer.tex \input{data/game_sd.tex}
                         |
                         v
              "...gamma* = \GameGammaSd ..."   (typeset value)
```

The data macros are written through `mvst_texdata`, which opens a `.tex` file
and exposes `.num(name,value)` for a formatted scalar, `.mat(name,M)` for a
LaTeX `bmatrix`, and `.str(name,text)` for a string, then `.close()`. Matrices
are formatted by `mvst_fmtmat`. This writer is the single choke point for
everything the book quotes.

## Figures and the export layer

The figure and export layer lives in its own folder, `bookutil/`, separate from
the engine layer `lib/`: nothing in the analysis engines depends on it, so the
solvers can be used headless without ever loading a plotting routine.
`mvst_setup` adds `bookutil/` to the path; only `mvst_figdir` (the output
folder, see above) stays in `lib/`.

- `mvst_gemcolors`, `mvst_latexdefaults`, and `mvst_figstyle` set the MATLAB
  "gem" color order, route all text through the LaTeX interpreter, and force a
  white background with black text before export;
- `mvst_savepdf` applies the white style and writes a cropped vector PDF, so the
  exports are correct even when MATLAB itself is in dark mode;
- `mvst_texdata` and `mvst_fmtmat` write the data macros;
- `mvst_stabfig`, `mvst_stabtable`, `mvst_steer_components`, and
  `mvst_steer_ellipse` draw the recurring chapter figures and generated table
  rows.

## The runners

Each book chapter has a per-chapter runner at the top level of `matlab/`. A
runner executes its chapter's example scripts with a fresh workspace (through
`mvst_runscripts`), writing the PDFs and data macros described above.

| Runner | Regenerates (examples) | Needs SOS solver | Rough runtime |
|---|---|---|---|
| `run_stability` | dwell-scalar, Levy, rate-SOS, sampled-data-SOS, stabz, switched-SOS | yes (the `*_sos`) | minutes |
| `run_stabilization` | the `ex_stabz_*` family (CT, DT, delay, MJLS, sampled-data, NCS, ...) | one SOS script | tens of seconds |
| `run_oc` | time-varying MJLS / sampled / switched costs and trajectories, the four horizon figures | no | tens of seconds |
| `run_steering` | covariance steering for MJLS, Poisson, sampled-data, switched | no | seconds |
| `run_dissipativity` | dissipativity for MJLS, Poisson, sampled-data, switched (plus time-varying) | no | tens of seconds |
| `run_game` | H-infinity game level for the five classes plus the gamma iteration | no | tens of seconds |
| `run_computational` | ARE, MJLS, and the SOS optimal-control / dissipativity / game examples | yes | minutes |
| `run_vdiss` | vector-dissipativity example | no | seconds |

The runtimes are indicative: the base-MATLAB runners finish in seconds to a
minute, while the SOS runs (`run_computational` and the `*_sos` scripts) are
dominated by the external SDP solver and take several minutes.

`run_all` is the thin wrapper that calls every runner plus the remaining
standalone example scripts (the cones, positive-systems, and performance
examples):

```matlab
run_all
```

The SOS / SDP examples require an SDP solver on the path (see
[Getting started](01-getting-started.md)); the rest need only base MATLAB. To
regenerate just one chapter, call its runner, for example `run_dissipativity`.
Build the book first if you will typeset it: the data macros are `\input` by the
manuscript, and some example scripts cross-reference the built document.

## The golden regression harness (development tree only)

The development tree carries a golden-data regression harness tied to the
book's figure data. It is not part of the release battery (the release battery
is the 26 unit-test suites run by `run_all_tests`, see
[Getting started](01-getting-started.md)), and it is left out of the release
tree.

The idea is mechanical. The book quotes hundreds of computed numbers, and a
refactor of the numerical core must not change any of them. When the code was
known good, the whole data-macro corpus was captured into `tests/golden/data/`.
`tests/mvst_regress` re-parses every golden file and the live `data/*.tex`
files and compares them macro by macro with numeric tolerance (default relative
`1e-6`, absolute `1e-12`); on any drift it raises `mvst:regress:failed`.
`tests/run_regress` is the strict `-batch` wrapper. Drift can be inspected
without failing through report mode:

```matlab
rep = mvst_regress('Report', true);
rep.files(find(~[rep.files.pass], 1)).diffs   % the first failing file's diffs
```

A change that is intended to alter a printed number requires re-baselining on
purpose: regenerate the live data, review the report so that only the expected
macros moved, re-capture the snapshot with a refreshed `tests/golden/MANIFEST.md`,
and commit the new golden files together with the code change that caused them.
