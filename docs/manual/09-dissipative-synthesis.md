# 9. Dissipative synthesis

[Manual index](index.md) | Prev: [Games](08-games.md) | Next: [Positive systems](10-positive-systems.md)

This chapter covers dissipative state-feedback synthesis (`mvst_dissyn_sos`,
method `S.dissyn`) and its output-feedback sibling (`mvst_dof_sos`). It
corresponds to the book chapter **Dissipativity-Based Synthesis** of the
companion monograph. The section below states the synthesis problem in full so
the chapter is self-contained.

## Mathematical setting

**The problem.** Given the ported system

$$
\mathrm{d}x = (A x + B_u u + B_w w)\,\mathrm{d}t
   + \sum_l (A_l x + B_l u)\,\mathrm{d}W_l ,
\qquad z = C x + D_u u + D_w w ,
$$

**find a state-feedback gain** `K` (and, when the system jumps, a jump gain
`K_J` acting on the pre-jump state) so that the closed loop under

$$
u = K x ,
\qquad
\dot x = (A + B_u K) x + B_w w ,
$$

is **dissipative** with respect to a prescribed supply, over a chosen dwell-time
regime. For the `L2`-gain supply this means: render the closed loop mean-square
stable and drive its `L2`-gain below a level `gamma`, then minimize `gamma`.

**Why it needs convexifying.** The dissipation certificate of the analysis
chapter, written for the closed loop, is the inequality (schematically)

$$
\operatorname{He}\big((A + B_u K)^{*} P\big) + (\text{supply and noise terms}) \preceq 0 ,
$$

which is **bilinear** in the storage matrix `S = -P` and the unknown gain `K`.
Bilinear matrix inequalities are not convex and cannot be handed to an SDP
solver directly.

**The convexifying change of variables.** Following the chapter, one congruences
the closed-loop inequality by the inverse storage matrix and substitutes new
variables

$$
W := S^{-1} , \qquad Y := K W , \qquad Y_J := K_J W^{-} ,
$$

after which the design inequality becomes **affine** in `(W, Y, Y_J)`, hence an
LMI (a "synthesis block"). With the shorthands `V_l = A_l W + B_l Y` and
`G_0 = C W + D_u Y`, the flow synthesis block has the schematic Schur form

$$
\Phi_f =
\begin{bmatrix}
  -\dot W + \operatorname{He}(A W + B_u Y) & B_w - G_0^{*}\Pi_2 & \operatorname{row}_l[V_l^{*}] & (N G_0)^{*} \\
  \ast & -\Pi_3 - \operatorname{He}(D_w^{*}\Pi_2) & \operatorname{row}_l[\cdot] & (N D_w)^{*} \\
  \ast & \ast & -I \otimes W & 0 \\
  \ast & \ast & \ast & -I_\nu
\end{bmatrix} \preceq 0 ,
$$

with a matching jump block `Phi_j` at resets, where `Pi = [Pi_1 Pi_2; Pi_2^* Pi_3]`
is the port weight pricing `(z; w)`. The **convexifiability** condition is
`Pi_1 <= 0`; factoring `-Pi_1 = N^{*} N` with `N` of full row rank (`nu` rows),
the fourth (performance) row and column above carry the output-strict channel
`(N G_0)^{*}`, `(N D_w)^{*}` closed by `-I_nu`. For the `L2`-gain supply
`Pi_1 = -I_pz` (so `N = I_pz`, `nu = pz`); for output-strict passivity
`N = sqrt(rho) I`, and for a non-output-strict supply `nu = 0` and the fourth row
and column drop out. Once feasible, the gains are recovered by

$$
K = Y W^{-1} , \qquad K_J = Y_J (W^{-})^{-1} .
$$

**The dwell-time scenarios.** The certificate `W(tau)` is allowed to depend on a
timer `tau` (the clock since the last jump), which is what makes causal
dwell-time synthesis possible. Four regimes are assembled:

- **arbitrary**: fixed jump sequence or CT LTI, constant `W`;
- **periodic**: inter-jump time `= Tp`, flow LMI on `[0, Tp]` plus a jump LMI;
- **minimum**: inter-jump time `>= Tmin`, flow LMI on `[0, Tmin]` plus a
  frozen-tail LMI (`Wdot = 0`) plus a jump LMI;
- **range**: inter-jump time in `[Tmin, Tmax]`, flow LMI on `[0, Tmax]` plus a
  jump LMI parametrized over the range.

**Standing assumptions:** `(A, B_u)` stabilizable; the port weight is
convexifiable (`Pi_1 <= 0`), which holds for the `L2`-gain and passivity
supplies; `W(tau) > 0` throughout.

## `S.dissyn`: dissipative state feedback

`S.dissyn(supply, dwell, mode)` searches (SOSTOOLS + SeDuMi) for a
state-feedback gain `K(tau)` (and a jump gain `K_J`) rendering the closed loop
dissipative with respect to `supply` over the dwell-time regime given by `dwell`, using exactly
the `W = S^{-1}`, `Y = K W` convexification above. `mode` is `'gamma'` (optimize
the supply level) or `'feas'` (certify feasibility at a fixed level).

Arguments default to the `L2`-gain supply and the arbitrary-dwell-time regime, so the
one-line form synthesizes the minimal-`L2`-gain state feedback:

```matlab
S = mvst.System.ct([0 1;-2 -1], [0;1], 'Bw', eye(2), 'C', [1 0;0 0]);
rs = S.dissyn();
fprintf('dissyn gammastar=%.6f  feas=%d\n', rs.gammastar, rs.feas);
```
```
dissyn gammastar=0.046472  feas=1
```

`res.gammastar` is the certified supply level, `res.feas` the feasibility flag,
`res.K` the gain function `@(tau) -> K`, `res.KJ` the jump gain (empty when
there is no scheduled jump), `res.W` the synthesis certificate `@(tau) -> W`, and
`res.S` the storage `@(tau) -> S = W^{-1} = -P`.

### The validation gate

This subsection is the one to read. An SDP solver reports "feasible" whenever it
finds a point satisfying the constraints **to its own numerical tolerance**, and
those tolerances can be loose enough to certify garbage: a certificate that is
barely indefinite, an LMI residual that is small only because the problem is
badly scaled, or a gain that does not in fact stabilize the true (unrelaxed)
closed loop. A raw `feas = true` from the SDP is therefore not, by itself,
trustworthy.

For this reason `feas = true` is **never** returned on the strength of the SOS
slack alone. Before returning, the routine recomputes the actual certificate and
LMI blocks at the returned `(W, Y, Y_J)` and runs four gates:

1. **Certificate positive definiteness.** The certificate `W(tau)` is checked to
   be PD on the timer grid: `min eig W(tau) >= epsW/2` at every
   sample. A near-singular `W` (which inverts into an enormous storage matrix)
   is rejected.
2. **LMI residual check.** The flow block (and the frozen-tail block for the
   minimum-dwell-time regime, and the jump block at resets) is recomputed at the
   returned variables and checked to be negative semidefinite within `restol`.
   This catches the case where the solver's internal residual was optimistic.
3. **Closed-loop stability.** The recovered gain `K = Y W^{-1}` is fed back and
   the closed loop's exact lifted spectral abscissa is computed; it must be
   negative (mean-square stable). A gain that satisfies the LMI numerically but
   destabilizes the true closed loop is rejected here.
4. **Achieved-vs-certified re-analysis.** For the exact continuous-time
   `L2`-gain case, the synthesized closed loop is handed back to the exact
   bounded-real backend `mvst_l2gain_ct`, and the **achieved** gain must be no
   worse than the certified level, within `gaptol`. This is the decisive check:
   the certificate is only accepted if an independent exact solver confirms the
   claimed performance.

A gate failure demotes the level to infeasible, so the gamma-bisection converges
to a certified level; if nothing passes, `res.feas = false` and
`res.info.gate` names the failing gate. In short, `feas = true` now means the
synthesized gain has been re-analyzed by the exact backend and confirmed, not
merely that an SDP reported a positive slack. The gate diagnostics live in
`res.info.gate`: `.ok`, `.reason`, `.mineigW`, `.condW`, `.abscissa`, and
`.achieved`. The standalone check `mvst_dissyn_sos('selftest')` re-runs this
re-analysis across a battery.

The `supply` struct accepts the same `'l2gain'` / `'passivity'` types as the
analysis routine (see [Dissipativity](07-dissipativity.md)), and the `dwell`
struct the same `'arbitrary'` / `'periodic'` / `'minimum'` / `'range'` regimes.

## Worked case study: synthesize, inspect the gate, close the loop, re-analyze

This end-to-end example synthesizes a minimal-`L2`-gain state feedback for an
oscillator whose open-loop gain is high, then opens up the validation gate to
see what was checked, reconstructs the gain, closes the loop, and re-analyzes it
with the independent exact backend. The point is to show that the certified
level is not taken on faith.

```matlab
A = [0 1; -2 -1]; Bu = [0;1]; Bw = eye(2); C = [1 0; 0 0]; Dw = [0 0; 0 1];
S  = mvst.System.ct(A, Bu, 'Bw', Bw, 'C', C, 'Dw', Dw);
rs = S.dissyn();
g  = rs.info.gate;
fprintf('certified gammastar = %.6f   feas=%d   slack=%.2e\n', ...
        rs.gammastar, rs.feas, rs.slack);
fprintf('gate.ok=%d (%s)\n', g.ok, g.reason);
fprintf('  (1) min eig W  = %.4f   cond(W) = %.1f\n', g.mineigW, g.condW);
fprintf('  (3) closed-loop spectral abscissa = %.4f\n', g.abscissa);
fprintf('  (4) achieved gain = %.6f  (<= certified %.6f)\n', g.achieved, rs.gammastar);
K   = rs.K(0);
Acl = A + Bu*K;
fprintf('gain K = %s\n', mat2str(K,5));
fprintf('eig(Acl) = %s\n', mat2str(sort(eig(Acl)),5));
reanalysis = mvst_l2gain_ct(Acl, Bw, C, Dw);
openloop   = mvst_l2gain_ct(A,   Bw, C, Dw);
fprintf('open-loop L2 gain    = %.6f\n', openloop);
fprintf('independent re-check = %.6f\n', reanalysis);
```
```
certified gammastar = 1.000100   feas=1   slack=3.87e-05
gate.ok=1 (pass)
  (1) min eig W  = 4.8554   cond(W) = 835.2
  (3) closed-loop spectral abscissa = -10.9409
  (4) achieved gain = 1.000062  (<= certified 1.000100)
gain K = [-796.51 -9.9409]
eig(Acl) = [-5.4705-27.723i;-5.4705+27.723i]
open-loop L2 gain    = 1.565858
independent re-check = 1.000062
```

Interpretation of the run: the synthesizer certifies `gamma* = 1.000100`, down from the
open-loop gain of 1.565858, so feedback improves the attenuation.
The gate then confirms the certificate is real: (1) `W` is comfortably positive
definite (smallest eigenvalue 4.86, condition number 835), (3) the closed loop
is deeply stable (spectral abscissa -10.94), and (4) the achieved gain
recomputed by the exact backend is 1.000062, no worse than the certified level.
Reconstructing `K = Y W^{-1}` and closing the loop independently gives a stable
`Acl` (a well-damped complex pair), and re-analyzing it from scratch reproduces
1.000062: the certified number and the true achieved number agree. Note the
large gain magnitude (about `-800` on the first state); pushing `gamma` to the
plant's floor demands high feedback authority, a real trade-off the certified
number does not hide.

## `mvst_dof_sos`: dynamic output feedback

`mvst_dof_sos(sys, dwell, opts)` is the output-feedback counterpart. It belongs
to the dynamic output-feedback development (treated separately from the
state-feedback chapter above), and synthesizes a timer-dependent dynamic
output-feedback controller reconstructed as `(A_c, B_c, C_c, D_c)`. Its signature
is deliberately shorter than the other SOS routines (no separate supply / cost
argument, no mode selector): the program is a fixed convex bracket in the timer,
solved without any inverse of the certificate. The example uses the scalar
timer-dependent plant of the output-feedback example under a periodic dwell-time:

```matlab
sys.a0 = 1.0; sys.a1 = 0.4; sys.sig = 0.5; sys.mu = 0.5;
sys.B0 = 1; sys.Cy = 1; sys.J0 = 0.5; sys.BJ = 1; sys.CJ = 1;
opts.alpha = 1e-3; opts.eps = 1e-3;
out = mvst_dof_sos(sys, struct('type','periodic','T',1.0,'deg',2), opts);
fprintf('periodic DOF solved=%d  max lam(M)=%+.3e  min lam(Om)=%+.3e\n', ...
        out.feas, out.maxM, out.minOm);
```
```
periodic DOF solved=1  max lam(M)=-3.416e-01  min lam(Om)=+5.925e-01
```

`out.feas` is the feasibility flag; `out.maxM` (the maximum eigenvalue of the
flow bracket, negative here) and `out.minOm` (the minimum eigenvalue of the
coupling block, positive here) are the two certificate margins verified over the
whole dwell-time interval, and `out.K` carries the reconstructed controller
`(A_c, B_c, C_c, D_c)` and jump gains. The `dwell` struct extends to `'minimum'`
and `'range'`.

Both routines require SOSTOOLS and an SDP backend.
