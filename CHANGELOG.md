# Changelog

All notable changes to the Matrix-Valued Systems Toolbox (MVST) are documented
here. The project follows [Semantic Versioning](https://semver.org/).

## 1.0.0-rc1 (2026-10-27)

First release candidate. Companion toolbox of the monograph *Linear Hermitian
Matrix-Valued Dynamical Systems* by Corentin Briat.

- Library of `mvst_*` engines in `lib/` for the analysis and synthesis of linear
  matrix-valued (impulsive and stochastic) dynamical systems: exact
  mean-square stability and stabilization, linear-cost optimal control (finite,
  infinite, and long-term average horizons), finite-time covariance steering,
  dissipativity (supply rates, KYP, L2 gain) and dissipativity-based
  synthesis, zero-sum and Nash games, further performance measures (H2,
  energy-to-peak, peak-to-peak), and linear positive systems.
- Object interface in the `+mvst` package: `mvst.System` (factories `ct`,
  `mjls`, `sampled`, `poisson`, `switched`, `impulsive`), `mvst.Options`,
  `mvst.solvercheck`, and `mvst.version`.
- Sum-of-squares engines built on SOSTOOLS, solved through the single dispatcher
  `mvst_sossolve` by an SDP solver (SeDuMi by default, SDPT3, MOSEK) or without
  one by SOSNOSDP (`solver = 'nosdp'`). The constant-certificate LMIs of
  `mvst_pos_ms`, `mvst_stab_delay`, and `mvst_perf_dwell` are SOSTOOLS programs
  by default (`opts.engine = 'sostools'`; legacy `'yalmip'` kept for one release).
  `mvst_kyp` is ported the same way (primal SDPs, witness as the primal variable). YALMIP is no longer
  bundled: it is an optional external dependency of the legacy `engine = 'yalmip'` paths only.
  The SOS slack programs of `mvst_pos_ms` (dwell-time diagonal), `mvst_perf_dwell` (peak-to-peak timer
  feasibility) and `mvst_pos_stabilize` (pure stabilization) cap the maximized slack at 1 (by a monotonicity
  argument: scaling a feasible point by c>1 stays feasible, so the slack is unbounded without the cap), and `mvst_pos_stabilize` no longer reports a program the
  solver flagged infeasible as feasible.
- Routed H-infinity norm and passivity (`mvst_hinfnorm`, `mvst_passivity`, and
  the `mvst.System` methods `hinfnorm`, `passivity`, `kyp`, `dissipativity`),
  exact discrete-time gains (`mvst_l2gain_dt`, `mvst_l2gain_dt_mjls`), and the
  general quadratic supply `struct('type','qsr','Q',Q,'S',S,'R',R)`
  (`mvst_qsr_supply`). The passivity index follows the book's convention
  (supply `He(w'z) - nu|w|^2`) in every routine; `mvst_passivity_ct` formerly
  returned `nu/2`.
- Test battery in `tests/` (26 suites). The development tree also carries a
  golden-data regression harness tied to the book's figure data; it is not part
  of the release battery.
- Chaptered user manual in `docs/`.
- Estimation and data-driven routines are not part of this release.
