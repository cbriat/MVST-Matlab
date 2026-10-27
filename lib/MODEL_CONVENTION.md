# Unified model convention

The companion code historically carried two incompatible system-struct
conventions:

* the **thin** convention (`mvst_mode` / `mvst_system`) consumed by the
  Riccati / optimal-control / game / steering layer, with quadratic-cost data
  and Brownian + Markov-jump dynamics only;
* the **rich** convention consumed by the SOS layer (`mvst_stab_sos`,
  `mvst_stab_sos_modes`, `mvst_stabz_sos`, `mvst_diss_sos`,
  `mvst_diss_sos_modes`), covering CT / DT / Brownian / Levy / Poisson /
  switched / sampled-data / impulsive / MJLS in one call.

The two now share ONE convention. The thin structs are extended with the
optional jump/reset fields of the rich naming, and two helpers translate
between them. The only genuine naming difference is the drift/input pair:
the thin structs call them `A`, `B`; the rich structs call them `A0`, `B0`.

## Field mapping table

| Meaning                                   | Thin (`mvst_mode`/`mvst_system`) | Rich (SOS layer)        |
|-------------------------------------------|----------------------------------|-------------------------|
| flow drift                                | `A`                              | `A0`                    |
| flow control input                        | `B`                              | `B0`                    |
| flow state-multiplicative noise channels  | `Ac`   (n x n x Nc)              | `Ac`   (n x n x q)      |
| flow input-multiplicative noise channels  | `Bc`   (n x m x Nc)              | `Bc`   (n x m x q)      |
| scheduled (dwell-time) jump / reset maps  | `J`  *(ext)* (n x n [x p])       | `J`    (n x n [x p])    |
| scheduled-jump input selector             | `BJ` *(ext)* (n x m [x p])       | `BJ`   (n x m [x p])    |
| spontaneous (Poisson) jump maps           | `Jf` *(ext)* (n x n x pf)        | `Jf`   (n x n x pf)     |
| Poisson rates                             | `lamf` *(ext)* (1 x pf)          | `lamf` (1 x pf)         |
| spontaneous-jump input selector           | `BJf` *(ext)* (n x m x pf)       | `BJf`  (n x m x pf)     |
| CT Markov generator (multi-mode)          | `Pi` (`mvst_system`)             | `Pi` (stab), `Pig` (diss) |
| DT Markov transition (multi-mode)         | `P`  *(ext, `mvst_system`)*      | `P`                     |
| switch / impulse reset (multi-mode)       | `J`  *(ext, `mvst_system`)*      | `J`                     |
| quadratic cost blocks                     | `Q`,`R`,`S`,`W`,`X0`,`ZT`        | (none)                  |
| disturbance input                         | (none)                           | `Bw`                    |
| disturbance-multiplicative noise          | (none)                           | `Dwc`                   |
| performance output / feedthrough (flow)   | (none)                           | `C`, `Dw`               |
| Poisson jump-disturbance channels         | (none)                           | `Bdf`                   |
| jump disturbance / output / feedthrough   | (none)                           | `Bd`, `CJ`, `DJ`        |

*(ext)* marks the optional extension fields added to `mvst_mode` /
`mvst_system`. Every extension defaults to **absent**: a struct built without
these arguments is byte-identical to the pre-extension output, and the
Riccati solvers (which read fields by name and ignore the rest) are
unaffected.

Multi-mode rich structs stack the per-mode data along the last dimension:
`A0` is `n x n x M`, `Ac` is `n x n x q x M`, `B0` is `n x m x M`, `Bc` is
`n x m x q x M`.

## Helper functions

### `mvst_torich(sys)` — thin (+extensions) -> rich

Idempotent and permissive. A struct that already carries `.A0` is returned
unchanged, so the SOS entry points call it unconditionally.

* thin single mode (`.A` present) -> single-mode rich: `A -> A0`, `B -> B0`,
  and `Ac`/`Bc`/`J`/`BJ`/`Jf`/`lamf`/`BJf` plus any dissipativity I/O fields
  already attached are copied verbatim;
* thin multi-mode system (`.modes` present) -> mode-indexed rich: modes are
  stacked, noise channels zero-padded to the common maximum `q`, the CT
  generator `Pi` exposed under both `Pi` (stability) and `Pig`
  (dissipativity), and optional `P` / `J` carried through.

Cost data has no SOS counterpart and is dropped.

### `mvst_tothin(rich, 'Q',Q, 'R',R, ...)` — rich -> thin

Inverse of `mvst_torich` on the shared fields. `A0 -> A`, `B0 -> B`, and the
noise / jump / reset fields are forwarded to `mvst_mode`. Cost data is
supplied through name/value pairs (`Q`,`R`,`S`,`W`,`X0`,`ZT`; defaults
`Q=0`, `R=I`, `S=0`, `W=0`, `X0=I`). A struct that is already thin is
returned unchanged.

Fields with no thin home make the conversion ill-posed and raise an error
listing them:

* mode-indexed data (`A0` with `M>1`, or a 4-D `Ac`): assemble each mode and
  wrap with `mvst_system` instead;
* dissipativity I/O blocks (`Bw`, `Dwc`, `C`, `Dw`, `Bdf`, `Bd`, `CJ`, `DJ`)
  and the generator/transition couplings (`Pig`, `P`).

## SOS entry-point routing

Each SOS entry point begins with

```matlab
if ~isfield(sys,'A0'); sys = mvst_torich(sys); end
```

so a thin `mvst_mode` / `mvst_system` struct is accepted transparently while
the rich fast path is untouched. The internals were not restructured.
