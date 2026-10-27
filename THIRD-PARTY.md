# Third-party notice

MVST (the code in `lib/`, `+mvst/`, `bookutil/`, `examples/`, `tests/` and the runners) is released under the MIT
license (`LICENSE`).

The folder `thirdparty/` is an aggregate of independent components that are NOT part of MVST and
remain under their own licenses, unmodified and with their license files intact. Details
(versions, upstream URLs, snapshot dates) are in `thirdparty/VERSIONS.md`.

| Component | License | License file |
|---|---|---|
| SOSTOOLS | GNU GPL v3 | `thirdparty/sostools/docs/GNU Public License.txt` |
| SeDuMi | GNU GPL v2 | `thirdparty/sedumi/COPYING` |
| SDPT3 | GNU GPL v2 | `thirdparty/sdpt3/GNU_General_Public_License_v2` (notice: `thirdparty/sdpt3/Copyright`) |
| SOSNOSDP | MIT | `thirdparty/sosnosdp/LICENSE` |

Attribution (MVST code, not bundled): `lib/mvst_sosmatrvar.m` and `lib/mvst_sosposmatrvar.m`
reimplement, with the current SOSTOOLS API (`sospolymatrixvar`, `sossosvar`, `sosmatrixineq`), the
helpers `sosmatrvar.m` and `sosposmatrvar.m` / `sosposmatr.m` of M. Peet's MatrixTools v.03
(matthew.peet@inria.fr), as distributed with the scripts of the NAHS "Linf Positive Systems"
package. The MatrixTools source files carry no license statement.

YALMIP: NOT bundled (removed from the distribution). It is an OPTIONAL external dependency, needed
only by the legacy `opts.engine = 'yalmip'` paths (`mvst_kyp`, `mvst_pos_ms`, `mvst_stab_delay`,
`mvst_perf_dwell`), which fail with `mvst:opts:yalmipMissing` when it is absent. Install it yourself
from https://yalmip.github.io/; `mvst_setup` and `mvst.solvercheck` detect an installation on the
path. Its license (Johan Lofberg) requires that YALMIP be referenced when it is used: please cite
YALMIP in any work that uses the MVST routines relying on it.

SDPT3: the snapshot ships without compiled MEX files; run `install_sdpt3` once (from
`thirdparty/sdpt3/`) before using it as a backend.

SeDuMi: the compiled MEX files are platform-specific and are not tracked in git
(`thirdparty/sedumi/.gitignore` ignores `*.mex*`). SeDuMi runs from MATLAB source
without them but slowly; run `install_sedumi` once (from `thirdparty/sedumi/`) to
locate or rebuild them for your platform.

Replacement: any subfolder of `thirdparty/` can be deleted or replaced by a newer upstream version
(or by your own installation already on the MATLAB path). `mvst_setup` adds only the subfolders that
exist, and `mvst.solvercheck` reports which solvers were detected. Use `mvst_setup('thirdparty', false)`
to skip the bundled copies entirely.
