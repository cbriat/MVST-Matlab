# 7. Dissipativity

[Manual index](index.md) | Prev: [Steering](06-steering.md) | Next: [Games](08-games.md)

This chapter covers the L2-gain and passivity, and the general dissipativity SOS
program. It corresponds to the book chapter **Dissipativity** of the companion
monograph. The section below states the framework in full so the chapter is
self-contained.

## Mathematical setting

**The system with a port.** Add a disturbance input `w` and a performance
output `z` to the matrix-valued system:

$$
\mathrm{d}x = (A x + B_w w)\,\mathrm{d}t
   + \sum_l (A_l x + D_{w,l} w)\,\mathrm{d}W_l ,
\qquad
z = C x + D_w w .
$$

The toolbox carries `A` (`A0`), the disturbance input `Bw`, the performance
output `C`, the feedthrough `Dw`, and, when present, the multiplicative and jump
channels. All second-order information is again carried by the joint
second-moment matrix `X = X^* >= 0`, and every quantity below is **linear** in
`X`.

**The supply rate.** Dissipativity meters the power crossing the `(w,z)` port
through a quadratic **supply rate**, written compactly as a linear functional of
`X`:

$$
\sigma(w,z) = \langle M, X\rangle ,
\qquad
M = M^{*} =
\begin{bmatrix} Q_s & S_s \\ S_s^{*} & R_s \end{bmatrix}
  (\text{indefinite, in the state/disturbance partition}).
$$

The decisive modelling choice is that `M` carries **no sign**: the right
indefinite signature recovers each classical supply. In vector terms
`sigma(w,z) = [x;w]^* M [x;w]`, and `M` is built from the port weights and the
output map. Throughout, $\operatorname{He}(N) = N + N^{*}$ (the doubled
Hermitian-part operator of [Systems](02-systems.md)) and
$\widetilde C = \begin{bmatrix} C & D_w \end{bmatrix}$ is the **stacked output map**, so that
$z = \widetilde C [x;w]$.

**The storage and the dissipation inequality.** A **linear storage function**
is `V(X) = <-P, X>` for a Hermitian `P` (no sign imposed; the sign convention
`-P` makes the available storage equal to `<-P(0), X_0>`). The system is
**dissipative** with respect to `M` when stored energy never grows faster than
the supply delivers, i.e. the dissipation inequality

$$
V\big(X(t_1)\big) - V\big(X(t_0)\big)  \le  \int_{t_0}^{t_1} \sigma\,\mathrm{d}t
$$

holds along every admissible trajectory. Differentiating gives the pointwise
form `d/dt <-P,X> <= <M,X>`.

**The mean-square hybrid `L2 x ell2` gain.** For the disturbance-to-output map,
the performance measure is the **(mean-square) hybrid `L2 x ell2` gain**: the
smallest `gamma >= 0` such that

$$
\lVert z\rVert_{L_2\times\ell_2} \le \gamma \lVert w\rVert_{L_2\times\ell_2}
$$

along every trajectory started from rest (`X_1(0) = 0`), where the hybrid norm
splits into a flow (`L2`) energy and a jump (`ell2`) energy,

$$
\lVert z\rVert^2_{L_2\times\ell_2}
   = \mathbb{E}\left[\int_0^\infty  z_c^{*}z_c\,\mathrm{d}t
     + \sum_{k\ge 1} z_{d,k}^{*}z_{d,k}\right].
$$

For a purely continuous-time system the jump sum is empty and this is the usual
mean-square `L2` gain, the state-space face of the `H_infinity` norm.

**The two supplies.** Two indefinite signatures recur:

- **`L2`-gain (bounded real):** `sigma = gamma^2 |w|^2 - |z|^2`, whose supply
  matrix in the state/disturbance partition is

$$
M = \begin{bmatrix} 0 & 0 \\ 0 & \gamma^2 I \end{bmatrix} - \widetilde C^{*}\widetilde C
  = \begin{bmatrix} -C^{*}C & -C^{*}D_w \\ -D_w^{*}C & \gamma^2 I - D_w^{*}D_w \end{bmatrix}.
$$

  The system is dissipative with this `M` exactly when its `L2`-gain is at most
  `gamma`.

- **Passivity:** `sigma = 2\Re\langle z, w\rangle`, i.e.

$$
M = \begin{bmatrix} 0 & C^{*} \\ C & \operatorname{He}(D_w) \end{bmatrix}.
$$

  Strict refinements subtract an **input-strict** index `nu I` from the
  `w`-block and an **output-strict** index `rho` times the output quadratic.

**The KYP certificate.** Dissipativity is certified by a storage matrix `P`
(equivalently `S = -P >= 0`) for which the augmented flow operator is positive
semidefinite:

$$
\mathcal{C}_M(P) =
\begin{bmatrix}
  \operatorname{He}(A^{*}P) + \sum_l A_l^{*}P A_l + Q_s &
  P B_w + \sum_l A_l^{*}P D_{w,l} + S_s \\[2pt]
  \ast & \sum_l D_{w,l}^{*}P D_{w,l} + R_s
\end{bmatrix} \succeq 0 ,
$$

where the blocks of `M` are **added** according to the chosen supply, matching
$\mathcal{C}_M = \mathcal{F}^{*}(P) + M$ of the book and the `flow_cM` assembly
in `mvst_diss_sos`. As a scalar sanity check, take the `L2`-gain supply with no
feedthrough (`Dw = 0`) and no multiplicative channel: the `(2,2)` block reduces
to $R_s = \gamma^2$, which is positive, so the certificate can be feasible.
Writing a minus there would force $-\gamma^2 < 0$ and make the LMI infeasible
for every `gamma`. For the `L2`-gain supply this is the matrix bounded-real
lemma; for passivity it is the matrix positive-real (KYP) lemma. **Standing
assumptions:** `A` Hurwitz
(mean-square stable) for the gain to be finite; `Q >= 0`, `R > 0` where the
supply requires them.

## Supply-rate convention in the toolbox

`mvst_diss_sos` selects between the two supplies with
`supply.type = 'l2gain' | 'passivity'`. The `'l2gain'` supply builds the `M`
above with `gamma` either fixed (`mode = 'feas'`) or minimized
(`mode = 'gamma'`); the `'passivity'` supply builds the port-power `M`, with the
strict indices `nu` and `rho` available as fields.

## Exact L2-gain: `S.l2gain`

`S.l2gain` computes the (mean-square) `L2`-gain by routing to the exact
per-class backend, integrating the bounded-real Riccati rather than relaxing it.
Supply the disturbance input `Bw` and output `C` through the factory:

```matlab
S = mvst.System.ct([-1 0.5; 0 -2], zeros(2,0), 'Bw', [1;0], 'C', [1 0]);
fprintf('L2 gain = %.6f  method=%s\n', S.l2gain().gamma, S.l2gain().method);
```
```
L2 gain = 1.000001  method=ct
```

`res.gamma` is the gain and `res.method` the resolved class. The class is
inferred from the data (Poisson jump data, sampling period, generator, or
dwell-time), or forced with `OPTS.class`. The backends are
`mvst_l2gain_ct`, `_ct_mjls`, `_ct_dwell`, `_ct_poisson`, and `_sd_hybrid`.

## Passivity by frequency sweep: `mvst_passivity_ct`

For a square continuous-time LTI system, `mvst_passivity_ct` tests passivity by
a frequency sweep of the transfer matrix and returns a passivity index `nu`. It
forms `G(jw) = C(jw I - A)^{-1}Bw + Dw`, takes the smallest eigenvalue of the
Hermitian part over the sweep, and reports `nu = min_w lambda_min(G+G^*)`, the
book's convention (supply `He(w'z) - nu|w|^2`, shared with `mvst_passivity` and
`mvst_isp_margin_ct`). A positive `nu` certifies input-strict passivity with
that index:

```matlab
[passive,nu] = mvst_passivity_ct([-1 0;0 -2], [1;1], [1 1], 0.5);
fprintf('passive=%d  nu=%.6f\n', passive, nu);
```
```
passive=1  nu=1.000000
```

The signature is positional, `mvst_passivity_ct(A, Bw, C, Dw)`, and the system
must be square (`size(C,1) == size(Bw,2)`) with `A` Hurwitz. `passive` flags
(approximate) passivity; a negative `nu` quantifies the shortage of passivity.
The sweep is approximate; `mvst_passivity` (below) gives the exact Hamiltonian
test and is the recommended entry point.

## Dissipativity SOS: `mvst_diss_sos`

`mvst_diss_sos(sys, supply, dwell, mode)` solves the dual dissipativity DLMIs as
a matrix-SOS program in the timer, covering the `L2`-gain and passivity supplies
and every dwell-time type of the chapter. The storage certificate is the dual
variable of the chapter (storage matrix `S = -P >= 0`); the routine poses
`C_M(P) >= 0` on flow (and its jump counterpart at resets) as a matrix sum of
squares. `mode` is `'gamma'` (the supply level is the free variable, minimized)
or `'feas'` (only certify feasibility at a fixed level).

### L2-gain supply, compared to the exact gain

Under the arbitrary-dwell-time (constant certificate, CT LTI) regime the SOS dual
reproduces the exact bounded-real gain:

```matlab
A = [-3 0;0 -3]; Bw=[0;1]; C=[0.5 1]; Dw=0.5;
r = mvst_diss_sos(struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                  struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
gR = mvst_l2gain_ct_mjls({A},Bw,{C},Dw,0);
fprintf('SOS-dual gamma=%.6f  exact Riccati gamma=%.6f  |diff|=%.2e\n', ...
        r.gamma, gR, abs(r.gamma-gR));
```
```
SOS-dual gamma=0.833334  exact Riccati gamma=0.833986  |diff|=6.52e-04
```

The SOS dual (0.833334) matches the exact Riccati gain (0.833986) to about
6e-4, the residual being the SOS relaxation gap at this certificate degree.

### Passivity supply

```matlab
rp = mvst_diss_sos(struct('A0',[-2 1;0 -3],'Bw',[1;0],'C',[1 0],'Dw',1), ...
                   struct('type','passivity'), struct('type','arbitrary'), 'feas');
fprintf('passivity feasible=%d\n', rp.feas);
```
```
passivity feasible=1
```

The `dwell` struct extends to `'discrete'`, `'periodic'` (`Tbar`), `'minimum'`
(`Tbar`), and `'range'` (`[Tmin,Tmax]`), and the system may carry
state-multiplicative (`Ac`, `Dwc`), Poisson (`Jf`, `lamf`, `Bdf`), and jump
(`J`, `Bd`, `CJ`, `DJ`) data. The multi-mode variant is `mvst_diss_sos_modes`.
This routine requires SOSTOOLS and an SDP backend; `S.l2gain` and
`mvst_passivity_ct` need only base MATLAB.

## Worked case study: compute an L2 gain exactly, then certify it by SOS

This end-to-end example takes a single non-normal CT system, computes its
`L2`-gain three independent ways, and confirms they agree: the exact
bounded-real Riccati engine (`S.l2gain`), the SOS dissipativity dual
(`mvst_diss_sos`), and a direct frequency sweep of the transfer function
`|C(jw I - A)^{-1}Bw + Dw|`. Agreement across all three is the practical check
that the certificate is real.

```matlab
A = [-2 1; 0 -3]; Bw = [1;1]; C = [1 0]; Dw = 0;
S  = mvst.System.ct(A, zeros(2,0), 'Bw', Bw, 'C', C, 'Dw', Dw);
ex = S.l2gain();
sos = mvst_diss_sos(struct('A0',A,'Bw',Bw,'C',C,'Dw',Dw), ...
                    struct('type','l2gain'), struct('type','arbitrary'), 'gamma');
w = logspace(-2,2,20000); G = zeros(size(w));
for k = 1:numel(w); G(k) = abs(C*((1i*w(k)*eye(2)-A)\Bw) + Dw); end
fprintf('exact Riccati gamma = %.6f  (method=%s)\n', ex.gamma, ex.method);
fprintf('SOS-dual     gamma = %.6f  (|diff|=%.2e)\n', sos.gamma, abs(ex.gamma-sos.gamma));
fprintf('freq-sweep  Hinf   = %.6f\n', max(G));
```
```
exact Riccati gamma = 0.666687  (method=ct)
SOS-dual     gamma = 0.666667  (|diff|=2.05e-05)
freq-sweep  Hinf   = 0.666657
```

Interpretation of the run: the exact bounded-real Riccati gain is 0.666687, the SOS dual
returns 0.666667 (within 2e-5, the relaxation gap), and the brute-force peak of
the frequency response is 0.666657. Three methods that never touch each other's
code agree to five digits, so the reported `L2`-gain is trustworthy. In a design
loop one would use `S.l2gain` for the fast exact value and `mvst_diss_sos` when
a dwell-time-constrained or multiplicative-noise system takes the exact backend
out of reach.

For linear positive systems, dissipativity with a general linear supply is
handled by the copositive-LP routine `mvst_pos_supply` (see
[Positive systems](10-positive-systems.md)).

## Further dissipativity engines

### `mvst_kyp`: operator-domain KYP

`res = mvst_kyp(sys, opts)` solves the Kalman-Yakubovich-Popov problem for the
matrix-valued second-moment dynamics, in continuous time, discrete time, or the
periodic impulsive setting (`opts.domain` is `'ct'` (default), `'dt'`, or
`'hybrid'`). The supply is chosen by `opts.supply` (`'boundedreal'`,
`'passivity'`) or given as an arbitrary-sign matrix `opts.M`. The engine solves
both the LMI form and the zero-frequency (dual) form, reports which one
certifies, and returns a frequency-domain witness when the LMI is infeasible. Both
forms are SOSTOOLS programs (primal SDPs with the witness as the primal variable), so
every backend works, including `'nosdp'`; the legacy YALMIP formulation is available
as `opts.engine = 'yalmip'` (optional external dependency).

```matlab
A = [-1 0.5;0 -2]; Bw = [1;0.5]; C = [1 1]; Dw = 0.2;
r = mvst_kyp(struct('A',A,'Bw',Bw,'C',C,'Dw',Dw),struct('domain','ct','supply','boundedreal'));
fprintf('regular=%d gamma=%.4f certifies=%s\n',r.regular,r.gamma,r.certifies);
```
```
regular=1 gamma=1.5750 certifies=both
```

The bisected level `gamma = 1.5750` equals the exact `mvst_l2gain_ct` value.
When the noise-channel reachability hypothesis fails, the engine refuses to
conclude (`res.regular = false`).

### `mvst_isp_margin_ct`: input-strict-passivity index

`[nu, P] = mvst_isp_margin_ct(A, Bw, C, Dw, Cc, Dwc, opts)` returns the largest
`nu` such that the Ito system `dx = (A x + Bw w) dt + sum_l (Cc{l} x + Dwc{l} w)
dW_l`, `z = C x + Dw w` is dissipative for the supply `He(w'z) - nu |w|^2`, by an
SOS (constant-certificate) program; `Cc` and `Dwc` are cells (possibly empty).
Without noise, `nu` equals the minimum over frequency of the smallest
eigenvalue of `G(jw) + G(jw)'` (the convention shared with `mvst_passivity`).

```matlab
nu = mvst_isp_margin_ct([-1 0.5;0 -2],[1;0.5],[1 1],2,{},{});
fprintf('nu=%.4f\n',nu);
```
```
nu=4.0000
```

A positive `nu` certifies input strict passivity; `P` (second output) is the
constant storage certificate.

### `mvst_qsr_supply`: quadratic supply contract

`[Q,S,R] = mvst_qsr_supply(supply, p, m)` validates and unpacks the general
quadratic supply `s(z,w) = [z;w]'[Q S;S' R][z;w]` for an output of dimension `p`
and a disturbance of dimension `m`. Pass `'jump'` as a fourth argument to select
the jump weights (`supply.Qk, Sk, Rk`). It is the single validation point used
by the SOS dissipativity engines and `mvst.System.dissipativity`.

```matlab
sup = struct('type','qsr','Q',-1,'S',0,'R',4);   % L2-gain supply, gamma = 2
[Q,S,R] = mvst_qsr_supply(sup,1,1);              % -> Q = -1, S = 0, R = 4
```

The sign convention is that of the Dissipativity chapter: dissipativity means
`d/dt <P,x> <= s`; the L2-gain supply is `Q = -I`, `S = 0`, `R = gamma^2 I`, and
the passivity supply with indices `(nu,rho)` is `Q = -rho I`, `S = I`,
`R = -nu I`.

### `mvst_hinfnorm`: H-infinity norm with a uniform result

`res = mvst_hinfnorm(sys, opts)` returns the H-infinity norm (mean-square
L2-gain) through the class-routed exact backends of `mvst_l2gain`, in one
uniform struct: `res.gamma`, `res.certificate` (a storage struct `.S`, `.gamma`,
`.residual` where available, `[]` otherwise), `res.class`, `res.method` (the
backend that produced the number), and `res.diagnostics` (including `.exact`,
false for shooting lower bounds). The object method is `S.hinfnorm()`.

```matlab
r = mvst_hinfnorm(struct('A0',[-1 1;0 -2],'Bw',[1;0],'C',[1 0]));
fprintf('gamma=%.4f class=%s method=%s\n',r.gamma,r.class,r.method);
```
```
gamma=1.0000 class=ct method=mvst_l2gain_ct
```

### `mvst_passivity`: routed passivity analysis with indices

`res = mvst_passivity(sys, opts)` tests passivity of a square-port system and
returns the input-strict-passivity index `res.nu_isp` (the largest `nu` with
supply `He(w'z) - nu |w|^2` dissipative; for an LTI system the minimum over
frequency of the smallest eigenvalue of `G(jw) + G(jw)'`), the output-strict
index `res.rho_osp`, a flag `res.passive`, a storage certificate, and the
method used (`res.method`). The default for continuous-time LTI systems is the
exact imaginary-axis test of the positive-real Hamiltonian; noise, mean-square,
and dwell-time classes route to the SOS programs. The object method is
`S.passivity()`.

```matlab
r = mvst_passivity(struct('A0',[0 1;-1 -1],'Bw',[0;1],'C',[1 0],'Dw',0.5));
fprintf('passive=%d nu=%.4f rho=%.4f method=%s\n',r.passive,r.nu_isp,r.rho_osp,r.method);
```
```
passive=1 nu=0.3333 rho=0.7534 method=hamiltonian
```

### `mvst_l2gain_dt` and `mvst_l2gain_dt_mjls`: discrete-time gains

`gam = mvst_l2gain_dt(A, Bw, C, Dw)` returns the H-infinity norm over the unit
circle of `x+ = A x + Bw w`, `z = C x + Dw w` (`A` Schur), exactly, through the
bilinear map to continuous time and the imaginary-axis Hamiltonian bisection. No
Control System Toolbox is needed. `[gam, S] = mvst_l2gain_dt_mjls(Ac, Bw, Cc,
Dw, Ptr)` returns the coupled mean-square gain of a discrete-time Markov-jump
system (cells `Ac`, `Cc` of mode matrices, row-stochastic transition matrix
`Ptr`) by the coupled bounded-real Riccati recursion, together with the cell `S`
of mode storage matrices as certificate. Both are reached through `mvst_l2gain`
and `mvst_hinfnorm` when the system is flagged discrete-time.

```matlab
g  = mvst_l2gain_dt(0.5,1,1,0);                              % 1/(1-0.5)
gm = mvst_l2gain_dt_mjls({0.5,0.4},1,{1,1},0,[0.9 0.1;0.2 0.8]);
fprintf('dt=%.4f  dt-mjls=%.4f\n',g,gm);
```
```
dt=2.0000  dt-mjls=1.9417
```

The Markov-jump value is a lower bound that converges from below as the
recursion horizon grows (see `help mvst_l2gain_dt_mjls`).

