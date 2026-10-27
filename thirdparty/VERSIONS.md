# Third-party components bundled with MVST

Each subfolder is an independent, unmodified snapshot of its upstream project, kept with its own
license file. Any subfolder can be deleted or replaced by a newer upstream version; `mvst_setup`
and `mvst.solvercheck` use whatever is present (a copy already on the MATLAB path takes precedence).

| Folder | Name | Version | Upstream | Snapshot date | License (file) | Status |
|---|---|---|---|---|---|---|
| `sostools/` | SOSTOOLS | 4.00 (14 Sep 2021, from `Contents.m`) | https://github.com/oxfordcontrol/SOSTOOLS | 2026-10-27 (folder relocated; contents as vendored earlier) | GNU GPL v3 (`sostools/docs/GNU Public License.txt`) | unmodified |
| `sedumi/` | SeDuMi | 1.3.7 (2023-03-18, from `Changelog.md`) | https://github.com/sqlp/sedumi | 2026-10-27 (folder relocated; contents as vendored earlier) | GNU GPL v2 (`sedumi/COPYING`) | unmodified |
| `sdpt3/` | SDPT3 | 4.0 (upstream commit 8456f3a44ce5b0180cb8f59717d367dda1a1a3fb, 2024-04-10) | https://github.com/sqlp/sdpt3 | 2026-10-27 | GNU GPL v2 (`sdpt3/GNU_General_Public_License_v2`, notice in `sdpt3/Copyright`) | unmodified (no precompiled MEX files: run `install_sdpt3` once) |
| `sosnosdp/` | SOSNOSDP | 0.1.0 (from `CITATION.cff`; includes the nsd_gram mixed-block fix, commit f6139d7) | https://github.com/cbriat/sosnosdp | 2026-10-27 | MIT (`sosnosdp/LICENSE`) | unmodified (`.git` excluded) |

Removed: `yalmip/` (YALMIP 26-June-2025) is no longer bundled (2026-10-27). It is an optional external dependency of the legacy `opts.engine = 'yalmip'` paths only; see `THIRD-PARTY.md`.
