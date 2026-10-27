# 3. Stability

[Manual index](index.md) | Prev: [Systems](02-systems.md) | Next: [Stabilization](04-stabilization.md)

This chapter covers mean-square stability analysis: the exact spectral tests of
`mvst_stability` (method `S.stability`) and the dwell-time sum-of-squares tests
of `mvst_stab_sos`. It corresponds to the companion monograph chapter
**Stability**.

It opens with the mathematical setting (what mean-square stability means, why a
spectrum decides it, how dwell-time notions and monodromy operators work, and
the margin conventions the toolbox reports), documents the exact and SOS
engines, and closes with an end-to-end stability-versus-dwell-time case study whose
critical dwell-time is confirmed against a closed-form value.

## Mathematical setting

### What is being stabilized

Recall from [Systems](02-systems.md) that the state is a Hermitian matrix
$X(t)$, almost always the second moment $X=\mathbb{E}[x x^{*}]$ of an underlying
vector system. Because $\operatorname{tr}X=\mathbb{E}\|x\|^{2}$, decay of $X$ is
exactly **mean-square** decay of $x$.

### Mean-square exponential stability

The system is (uniformly) exponentially stable when the matrix state contracts
at a geometric rate. For an impulsive system with flow rate $\alpha>0$ and
per-jump contraction $\rho\in(0,1)$ the definition is the hybrid bound

$$\|X(t)\|\le\beta e^{-\alpha(t-s)} \rho^{ \kappa(t,s)} \|X(s)\|,\qquad t\ge s\ge t_0,$$

with $\beta\ge1$, $\kappa(t,s)$ the number of jumps in $(s,t]$, and
$\|X\|=\operatorname{tr}X$. Two degenerate cases specialize it: persistent
flowing ($\rho=1$) gives $\|X(t)\|\le\beta e^{-\alpha(t-s)}\|X(s)\|$, and
persistent jumping ($\alpha=0$) gives
$\|X(t)\|\le\beta\rho^{\kappa(t,s)}\|X(s)\|$.

Since $\|X\|=\operatorname{tr}X=\mathbb{E}\|x\|^{2}$, this is mean-square
exponential stability of the vector system. For the pure continuous-time
diffusion $\mathrm{d}x=A_0x\,\mathrm{d}t+\sum_l A_l x\,\mathrm{d}W_l$ it is
equivalent to the existence of $P\succ0$ with

$$A_0^{*}P+PA_0+\sum_{l=1}^{q}A_l^{*}P A_l\prec0,$$

and for the discrete recursion $x(k{+}1)=\sum_i\nu_i A_i x(k)$ (the $\nu_i$ the
i.i.d. unit-variance scalars of [Systems](02-systems.md)) to
$\sum_i A_i^{*}P A_i-P\prec0$. These are the Lyapunov identities the spectral
tests below decide exactly.

### The lift, and why a spectrum decides stability

The flow and jump operators are linear on the space of Hermitian matrices, so
they have a **matrix representation**. For the **real data** the engines handle
(see the data-scope note in [Systems](02-systems.md)), the identity
$\operatorname{vec}(A X B^{\top})=(B\otimes A)\operatorname{vec}(X)$ turns the
operators $\mathcal{F}(X)=A_0X+XA_0^{\top}+\sum_l A_l X A_l^{\top}$ and
$\mathcal{J}(X)=\sum_l J_l X J_l^{\top}$ into the $n^{2}\times n^{2}$ matrices

$$M_{\mathcal{F}}=I\otimes A_0+A_0\otimes I+\sum_{l}A_l\otimes A_l,\qquad M_{\mathcal{J}}=\sum_{l}J_l\otimes J_l,$$

with **no conjugation** (the general Hermitian form carries $\overline{A_0}$,
$\overline{A_l}$, $\overline{J_l}$, which coincide with these on real data). This
is exact for real data. Because the eigenmatrices are symmetric, one restricts to
the real subspace through an isometry $U$ (the scaled half-vectorization of
`mvst_svecbasis`, dimension $n(n+1)/2$), giving the **reduced** real
representations

$$\widetilde M_{\mathcal{F}}=U^{\top}M_{\mathcal{F}}U,\qquad \widetilde M_{\mathcal{J}}=U^{\top}M_{\mathcal{J}}U.$$

This is exactly what `mvst_liftsym` builds. The restriction is exact, so the
operator spectrum equals the ordinary matrix spectrum of the reduced matrix, and
stability reduces to a spectral condition:

- **Continuous-flow classes** (`ct`, `mjls`, `poisson`, `levy`): the solution of
  $\dot X=\mathcal{F}(X)$ is $e^{\widetilde M_{\mathcal{F}}t}\operatorname{svec}(X_0)$,
  so it decays iff the **spectral abscissa** is negative,

$$\mu(\mathcal{F})=\max\{\operatorname{Re}\lambda:\lambda\in\lambda(\widetilde M_{\mathcal{F}})\}<0.$$

- **Discrete / monodromy classes** (`dt`, `dtmjls`, `sampled`, `impulsive`,
  `switched`): the solution of $X(k{+}1)=\mathcal{J}(X(k))$ is
  $\widetilde M_{\mathcal{J}}^{ k}\operatorname{svec}(X_0)$, so it decays iff the
  **spectral radius** is below one, equivalently the log spectral radius is
  negative,

$$\varrho(\mathcal{J})=\max\{|\lambda|:\lambda\in\lambda(\widetilde M_{\mathcal{J}})\}<1\iff\log\varrho(\mathcal{J})<0.$$

### Monodromy operators for dwell-time classes

A dwell-time or periodic system alternates a jump-free flow of length $\theta$
with a reset. The pre-jump samples $\zeta_k=X(t_k^{-})$ obey the discrete
recursion $\zeta_{k+1}=\Psi_\theta(\zeta_k)$ with the **monodromy operator**

$$\Psi_\theta=\mathcal{U}(\theta,0)\circ\mathcal{J},$$

where $\mathcal{U}(\theta,0)$ is the flow over one interval. For constant data the
reduced monodromy matrix is the explicit product

$$\widetilde M_{\mathcal{J}} e^{\widetilde M_{\mathcal{F}}\theta},$$

and stability is decided by its spectral radius. The order of the two factors is
immaterial to the verdict: $\widetilde M_{\mathcal{J}}e^{\widetilde M_{\mathcal{F}}\theta}$
and $e^{\widetilde M_{\mathcal{F}}\theta}\widetilde M_{\mathcal{J}}$ share the
same spectrum (as $AB$ and $BA$ always do), so either product may be formed. This
is the matrix `mvst_stability` forms for the dwell-time classes; crossing $\varrho=1$
(log radius zero) marks the **critical dwell-time**.

### Dwell-time notions

The admissible jump sequences are constrained through the spacings
$\Delta_k=t_{k+1}-t_k$:

- **Arbitrary**: no constraint (certified by a single constant $P$).
- **Constant (periodic)**, $\Delta_k\equiv\bar T$: one periodic sequence.
- **Minimum**, $\Delta_k\ge\bar T$: jumps never closer than $\bar T$.
- **Range**, $\Delta_k\in[T_{\min},T_{\max}]$.

(The chapter also treats maximum $\Delta_k\le\bar T$, average
$\kappa(t,s)\le N_0+(t-s)/\tau^{*}$, and persistent dwell-time; the SOS engine
covers the first four.)

### Margin conventions

All margins the toolbox reports are **second-moment (mean-square) rates**,
because the analysis runs on the lifted second-moment operator. For the scalar
diffusion $\dot x=ax$ the covariance obeys $\dot X=2aX$: the second moment decays
at rate $2a$, twice the deterministic state exponent $a$. So a reported margin is
twice the corresponding first-moment exponent; halve it to recover the
deterministic rate. A **negative** margin always means mean-square stable:

- continuous-flow classes: margin $=\mu(\mathcal{F})$, the spectral abscissa;
- discrete / monodromy classes: margin $=\log\varrho$, the log spectral radius.

The `res.kind` field names the criterion used.

### Timer-dependent certificates

For dwell-time (not just periodic) regimes, stability is certified by a
**timer-dependent** Lyapunov matrix $P(\tau)$, $\tau=t-t_k$, satisfying a
differential LMI (DLMI). The representative periodic condition is

$$\dot P(\tau)+\mathcal{F}^{*}(\tau,P(\tau))\preceq-\alpha_3 I \text{ on } (0,T_p),\qquad \mathcal{J}^{*}(P(0))-P(T_p)\preceq-\alpha_4 I,$$

which is necessary and sufficient for exponential stability. The DLMI is made
finite by seeking $P(\tau)=\sum_{i=0}^{d}P_i\tau^{i}$ and enforcing the
inequalities as **matrix sum-of-squares** conditions in $\tau$, localized to the
interval by a multiplier $\tau(\bar T-\tau)\ge0$. This is what `mvst_stab_sos`
solves; the `.deg` option sets the polynomial degree $d$.

For deeper theory (Perron-Frobenius structure, looped functionals, the
delay case) see the companion monograph.

## Exact analysis: `S.stability`

`S.stability` decides mean-square (exponential) stability from an exact
spectral test on the relevant lifted second-moment operator (the reduced
representation on the symmetric subspace). No sum-of-squares and no
optimization are used. The system is treated as autonomous: the drift is the
flow or closed-loop generator, and any input matrix is ignored.

`res.stable` is the decision, `res.margin` the margin (abscissa or log-radius as
above), `res.spectrum` the lifted eigenvalues, and `res.kind` the criterion.

### Continuous-time with state-multiplicative noise

The lifted generator here is
$\widetilde M_{\mathcal{F}}$ of $\dot X=\operatorname{He}(AX)+A_c X A_c^{*}$, and
its abscissa is the margin.

```matlab
S = mvst.System.ct([0 1;-2 -1], zeros(2,0), 'Ac', 0.15*eye(2));
r = S.stability();
fprintf('stable=%d  margin=%.6f  kind=%s\n', r.stable, r.margin, r.kind);
```
```
stable=1  margin=-0.977500  kind=ct:spectral-abscissa
```

The margin $-0.9775$ is the second-moment rate; the deterministic state exponent
is about $-0.49$.

### Switched system on a dwell-time grid

For a switched system the decision is for the **single periodic schedule** that
visits the modes in index order (1, 2, ..., M) with the same dwell-time `T` in
each mode (the monodromy over one such cycle), not a certificate for arbitrary
or worst-case switching. Pass the dwell-time through `OPTS.T`; a vector produces
a stability-vs-dwell-time curve in `res.curve`.

```matlab
S = mvst.System.switched({[-1 1;0 -2],[-2 0;1 -1]}, {zeros(2,0),zeros(2,0)});
r = S.stability(struct('T',0.5));
fprintf('stable=%d  margin=%.6f  kind=%s\n', r.stable, r.margin, r.kind);
```
```
stable=1  margin=-2.000000  kind=switched:monodromy-log-radius
```

### Impulsive system

The margin is the log spectral radius of the one-cycle monodromy
$\widetilde M_{\mathcal{J}}e^{\widetilde M_{\mathcal{F}}T}$.

```matlab
S = mvst.System.impulsive([0 1;-1 0], zeros(2,0), 0.5*eye(2));
r = S.stability(struct('T',0.5));
fprintf('stable=%d  margin=%.6f  kind=%s\n', r.stable, r.margin, r.kind);
```
```
stable=1  margin=-1.386294  kind=impulsive:monodromy-log-radius
```

A pure discrete-time single-mode system has no dedicated factory; request it
with `OPTS.class = 'dt'`.

## Case study: stability versus dwell-time

This end-to-end study builds a system with a flow/jump trade-off, sweeps
the dwell-time to trace the stability curve, locates the critical dwell-time
where the monodromy radius crosses one, and confirms it against a closed-form
prediction.

Take an **unstable** oscillatory flow damped only by a **contractive periodic
reset**:

$$A=\begin{bmatrix}0&1\\2&0.1\end{bmatrix}\quad(\text{unstable}),\qquad J=0.4 I_2\quad(\text{contractive}).$$

Over one dwell-time interval of length $T$ the second-moment monodromy grows through
the flow (lifted abscissa $a=\mu(\mathcal{F})>0$) and contracts through the reset
(lifted radius $\rho=\varrho(\mathcal{J})<1$), so its radius is roughly
$\rho e^{aT}$ and the log radius crosses zero at the critical dwell-time
$T^{*}=\log(1/\rho)/a$: stable for $T<T^{*}$, unstable for $T>T^{*}$.

```matlab
A0 = [0 1; 2 0.1];  J = 0.4*eye(2);
S  = mvst.System.impulsive(A0, zeros(2,0), J);

a   = max(real(eig(mvst_liftsym(A0,{}))));       % lifted flow abscissa
rho = max(abs(eig(mvst_liftsym(zeros(2),{J}))));  % lifted reset radius
fprintf('lifted flow abscissa a   = %.6f\n', a);
fprintf('lifted reset radius  rho = %.6f\n', rho);

Tgrid = linspace(0.05, 1.2, 240);
r = S.stability(struct('T', Tgrid));       % stability-vs-dwell-time curve
idx = find(r.curve.margin < 0, 1, 'last');
Tstar = interp1(r.curve.margin(idx:idx+1), Tgrid(idx:idx+1), 0);
fprintf('critical dwell-time  T* (from curve)      = %.6f\n', Tstar);
fprintf('critical dwell-time  T* = log(1/rho)/a    = %.6f\n', log(1/rho)/a);
```
```
lifted flow abscissa a   = 2.930194
lifted reset radius  rho = 0.160000
critical dwell-time  T* (from curve)      = 0.625413
critical dwell-time  T* = log(1/rho)/a    = 0.625413
```

The critical dwell-time obtained from the swept curve matches the closed-form
$\log(1/\rho)/a$ to all printed digits. Sampling a few points across $T^{*}$
shows the sign change of the log radius:

```matlab
for T = [0.2 0.5 0.6244 0.8 1.0]
    rr = S.stability(struct('T', T));
    fprintf('  T=%.4f  logradius=%+.4f  stable=%d\n', T, rr.margin, rr.stable);
end
```
```
  T=0.2000  logradius=-1.2465  stable=1
  T=0.5000  logradius=-0.3675  stable=1
  T=0.6244  logradius=-0.0030  stable=1
  T=0.8000  logradius=+0.5116  stable=0
  T=1.0000  logradius=+1.0976  stable=0
```

Interpretation: frequent resets (small $T$) keep the second moment contracting
despite the unstable flow, but once the dwell-time exceeds $T^{*}\approx0.625$ the
flow growth over one interval outpaces the fixed per-jump contraction and the
mean square diverges. This is the mechanism behind every dwell-time theorem:
stability is a race between continuous growth and discrete contraction, decided
by the monodromy radius.

The same call structure handles a **sample-and-hold reset** that damps only the
velocity, $J=\operatorname{diag}(1,0.3)$ on a neutrally stable oscillator
$A=\begin{bmatrix}0&1\\-1&0\end{bmatrix}$. Sweeping the sampling period gives a non-monotone margin that
is stable across the whole grid, best near $T\approx0.87$ and degrading at both
short and long periods:

```matlab
Ss = mvst.System.impulsive([0 1;-1 0], zeros(2,0), [1 0;0 0.3]);
for T = [0.2 0.5 1.0 2.0 3.0]
    rr = Ss.stability(struct('T', T));
    fprintf('  T=%.3f  logradius=%+.4f  stable=%d\n', T, rr.margin, rr.stable);
end
```
```
  T=0.200  logradius=-0.0768  stable=1
  T=0.500  logradius=-0.6301  stable=1
  T=1.000  logradius=-1.2040  stable=1
  T=2.000  logradius=-1.2040  stable=1
  T=3.000  logradius=-0.0378  stable=1
```

## Dwell-time SOS analysis: `mvst_stab_sos`

The sum-of-squares stability test certifies stability of an impulsive system
under a dwell-time regime by a matrix-SOS program in the timer (SOSTOOLS +
SeDuMi). It is the computational form of the timer-dependent DLMI certificate
above: it searches for a polynomial $P(\tau)=\sum_i P_i\tau^{i}$ making the flow
and jump inequalities hold as SOS matrices on the dwell-time interval. The dwell-time
struct selects the regime: `'arbitrary'`, `'discrete'`, `'periodic'` (fixed
`Tbar`), `'minimum'` (`Tbar`), or `'range'` (`[Tmin,Tmax]`), and `.deg` sets the
timer polynomial degree of the certificate.

The following certifies a fixed periodic dwell-time for a marginally stable flow
with a contractive reset:

```matlab
sys = struct('A0',[0 1;-1 0], 'J', 0.5*eye(2));
res = mvst_stab_sos(sys, struct('type','periodic','Tbar',0.5,'deg',4));
fprintf('feas=%d  slack=%.3e\n', res.feas, res.slack);
```
```
feas=1  slack=5.000e-01
```

`res.feas` is the feasibility flag and `res.slack` the certificate slack (a
positive value here confirms a strictly feasible storage certificate). The
product-cone variant `mvst_stab_sos_modes` handles multi-mode switched systems.

This routine requires SOSTOOLS and an SDP backend (see
[Getting started](01-getting-started.md)). The exact `S.stability` needs only
base MATLAB, and the case study above uses only the exact engine.

## Further stability engines

These two engines are called directly as `mvst_*` functions (there is no
`mvst.System` method for them).

### `mvst_stab_critical`: critical parameter on the stability margin

`res = mvst_stab_critical(sys, param, range, opts)` finds, by bisection, the
value of a scalar parameter in `range = [lo hi]` at which the exact
mean-square stability margin of `sys` (from `mvst_stability`) crosses zero. By
default `param` is forwarded as an `opts` field to `mvst_stability`, so the
natural use is the dwell-time or sampling period `'T'`; for an arbitrary
parameter (a noise scale, a coupling rate) supply `opts.setter = @(sys,x)
sys_rebuilt_at_x`, and `param` is then only a label.

```matlab
sys = struct('A',[0.3 1;0 0.3],'J',0.5*eye(2));   % unstable flow, contracting reset
res = mvst_stab_critical(sys,'T',[0.01 5],struct('class','impulsive'));
fprintf('crit=%.4f margin_at=%.1e\n',res.crit,res.margin_at);
```
```
crit=2.3105 margin_at=2.1e-09
```

The system is stable for dwell-time periods below `res.crit = 2.3105` and unstable
above; `res.margin_at` is the margin at the returned value (zero up to the
bisection tolerance). If the margin does not change sign on `range`, the
function errors rather than returning an end point.

### `mvst_stab_delay`: delay-independent stability

`res = mvst_stab_delay(sys, opts)` certifies stability for all delays
`h_i >= 0` at once, by a Lyapunov-Krasovskii LMI and, equivalently, the
spectrum of the zero-delay operator. `opts.domain` is `'ct'` (default), `'dt'`,
or `'positive'`; `sys` carries the drift `A0`, the drift channels, and a cell
`Ad` of delayed blocks. It is a SOSTOOLS program by default (any backend, including `'nosdp'`); the legacy
YALMIP formulation is `opts.engine = 'yalmip'` (optional external dependency).

```matlab
C = struct('A0',diag([-2 -3]),'Ad',{{diag([0.5 0.5])}});
r = mvst_stab_delay(C,struct('domain','ct'));
fprintf('stable=%d margin=%.4f\n',r.stable,r.margin);
```
```
stable=1 margin=-3.7500
```

`res.margin` is the spectral abscissa of the zero-delay operator; it is
negative here, so the system is stable for every delay.
