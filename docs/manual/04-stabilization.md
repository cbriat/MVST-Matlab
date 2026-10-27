# 4. Stabilization

[Manual index](index.md) | Prev: [Stability](03-stability.md) | Next: [Optimal control](05-optimal-control.md)

This chapter covers state-feedback stabilization: the exact Riccati routes of
`mvst_stabilize` (method `S.stabilize`) and the sum-of-squares synthesis of
`mvst_stabz_sos`. It corresponds to the companion monograph chapter
**Stabilization**.

It opens with the mathematical setting (the feedback structure, the closed-loop
covariance generator, the Riccati equation the exact route solves, and the
convexified change of variables behind the dwell-time synthesis), documents both
engines, and closes with an end-to-end case study that designs, certifies, and
interprets a family of feedbacks up to the mean-square stabilizability limit.

## Mathematical setting

### State-feedback structure

The controlled diffusion of [Systems](02-systems.md),

$$\mathrm{d}x=(A x+B u)\,\mathrm{d}t+\sum_{l=1}^{q}(A_l x+B_l u)\,\mathrm{d}W_l,$$

is closed with a **linear state feedback** $u=K x$. Substituting collapses drift
and channels into closed-loop coefficients

$$A_{\mathrm{cl}}=A+BK,\qquad A_{l,\mathrm{cl}}=A_l+B_l K,$$

so the closed-loop second moment $X=\mathbb{E}[x x^{*}]$ obeys the covariance ODE

$$\dot X=\operatorname{He}\big(A_{\mathrm{cl}}X\big)+\sum_{l=1}^{q}A_{l,\mathrm{cl}} X A_{l,\mathrm{cl}}^{*}=\mathcal{F}_{\mathrm{cl}}(X),$$

with $\operatorname{He}(M)=M+M^{*}$. **Mean-square stabilization** is the design
of $K$ making $\mathcal{F}_{\mathrm{cl}}$ mean-square exponentially stable, that
is (from [Stability](03-stability.md)) making the lifted generator
$\widetilde M_{\mathcal{F}_{\mathrm{cl}}}$ have negative spectral abscissa. Note
that the input channels $B_l$ make the noise itself feedback-dependent: a gain
that shrinks the drift can amplify the multiplicative noise, so a plant need not
be mean-square stabilizable even when its deterministic part is.

### The Riccati route

For the classes with an exact route, the optimal stabilizing gain comes from the
stabilizing solution $P\succ0$ of the continuous-time **indefinite / stochastic
algebraic Riccati equation** (the I-GRE). This is exactly what the engines
`mvst_ric_ct` / `mvst_care_ct` compute on real data (see the data-scope note in
[Systems](02-systems.md)), so it is written here in the real (transpose) form the
code uses:

$$\operatorname{He}(A^{\top}P)+\sum_{l=1}^{q}A_l^{\top}P A_l+Q-S^{\top}\Theta^{-1}S=0,$$

with the noise-weighted quantities

$$\Theta=R+\sum_{l=1}^{q}B_l^{\top}P B_l,\qquad S=B^{\top}P+S_{\times}^{\top}+\sum_{l=1}^{q}B_l^{\top}P A_l,$$

where $Q=Q^{\top}\succeq0$ and $R=R^{\top}\succ0$ are the state and input weights
and $S_{\times}$ the cross weight. Here $\operatorname{He}(A^{\top}P)=A^{\top}P+PA$.
The stabilizing gain is the completion of squares

$$K=-\Theta^{-1}S,\qquad u=Kx.$$

Two features distinguish it from the deterministic ARE: the diffusion channels
add $\sum_l A_l^{\top}P A_l$ to the state term, and the input channels make the
control weight $\Theta=R+\sum_l B_l^{\top}P B_l$ state-solution-dependent. When
$R+\sum_l B_l^{\top}PB_l$ fails to stay positive definite the design has hit the
mean-square stabilizability limit and no stabilizing $P\succ0$ exists.

The toolbox obtains $P$ by the monotone-limit construction of the monograph:
integrate the differential I-GRE $\dot P=-\mathcal{R}(P)$ backward from
$P(T)=0$ over a growing horizon until the residual $\|\mathcal{R}(P)\|$ is
negligible; $P_T$ increases monotonically to $P_\infty$. This needs only ODE
integration.

Per-class variants keep the same shape: the `mjls` class solves the **coupled**
Riccati system (one $P_i$ per mode, coupled through $\sum_j\pi_{ji}P_j$ added to
the state term), the `poisson` class adds the jump term
$\sum_r\lambda_r(J_r^{*}P J_r-P)$, and the `sampled` class solves the periodic
sampled-data I-GRE over one sampling interval.

### Convexified dwell-time synthesis

The Riccati route gives a single constant gain. Under a dwell-time constraint the
gain is allowed to depend on the timer $\tau=t-t_k$, and the synthesis is made
convex by a linearizing change of variables. The closed-loop timer-dependent
DLMI in a Lyapunov matrix $P(\tau)$ is bilinear in $(P,K)$, so one works with the
**dual** variables

$$Y(\tau)=P(\tau)^{-1}\succ0,\qquad U(\tau)=K(\tau) Y(\tau),$$

which turns the closed-loop condition into an LMI/SOS in $(Y,U)$: schematically,
for the flow part,

$$-\dot Y(\tau)+\operatorname{He}\big(A Y(\tau)+B U(\tau)\big)+\Xi\big(Y,U\big)\preceq0,$$

where $\Xi$ collects the diffusion channels via a Schur complement, together with
a jump condition linking $Y(0)$ and $Y(\bar T)$ through the reset. The gain is
recovered as $K(\tau)=U(\tau)Y(\tau)^{-1}$. Seeking polynomial $Y(\tau)$,
$U(\tau)$ and enforcing the inequalities as matrix sum-of-squares in $\tau$ (as
in [Stability](03-stability.md)) makes the search finite. This is what
`mvst_stabz_sos` solves, returning the timer-dependent gain as a function handle.

For deeper theory (existence, the maximal solution, dwell-time necessity) see the
companion monograph.

## Exact synthesis: `S.stabilize`

`S.stabilize` designs a mean-square stabilizing feedback `u = K x` by the
closed-form / Riccati route of the resolved class, then **certifies** the closed
loop by an exact call to `mvst_stability`. The design is only accepted when
`res.check.stable` is true. Routing by class mirrors the book:

- `ct`: a single stochastic algebraic Riccati (`mvst_care_ct`);
- `mjls`: the coupled continuous-time MJLS Riccati (`mvst_solve_infinite`);
- `switched`: per-mode infinite-horizon gains (valid for a large-enough
  dwell-time);
- `sampled`: the periodic sampled-data impulsive Riccati
  (`mvst_igre_periodic_sd`);
- `poisson`: the GBM + Poissonian-jump algebraic Riccati (`mvst_care_ct_poisson`).

### Continuous-time with multiplicative noise

```matlab
S = mvst.System.ct([0 1;-1 0.2], [0;1], 'Q', eye(2), 'R', 1, ...
                   'Ac', 0.3*[0 0;1 0], 'Bc', [0;0.2]);
r = S.stabilize();
disp(r.K); fprintf('closed-loop stable=%d\n', r.check.stable);
```
```
   -0.5040   -1.5421
closed-loop stable=1
```

`res.K` is the gain (a matrix for single-mode classes, a cell of gains for the
multi-mode classes), `res.P` the Riccati solution $P$, `res.closedloop` the
autonomous closed-loop struct fed to the check, and `res.check` the stability
result on it (its `.margin` is the certified second-moment abscissa).

### Markov-jump system

```matlab
Pi = [-1 1; 2 -2];
S  = mvst.System.mjls({[0 1;-1 0.2],[0 1;-2 -0.1]}, {[0;1],[0;1.2]}, Pi, ...
                      'Q', eye(2), 'R', 1);
r = S.stabilize();
fprintf('closed-loop stable=%d\n', r.check.stable);
```
```
closed-loop stable=1
```

The discrete-time (`dt`, `dtmjls`) and the dwell-time-constrained impulsive
designs are handled by the SOS engine below, which `S.stabilize` deliberately
does not duplicate.

## Case study: the mean-square stabilizability limit

This end-to-end study designs a feedback for a family of increasingly noisy
plants, certifies each closed loop, and locates the noise level beyond which no
gain can stabilize the mean square. It illustrates the point made above: state-
and input-multiplicative noise scale with the design, so a deterministically
controllable plant can lose mean-square stabilizability.

Take the unstable double-integrator-like plant

$$A=\begin{bmatrix}0&1\\1&0.5\end{bmatrix},\qquad B=\begin{bmatrix}0\\1\end{bmatrix},$$

with open-loop eigenvalues $\{1.2808,-0.7808\}$ (deterministically controllable),
and attach a common noise level $\sigma$ to both a state channel and the input
channel, $A_c=\sigma\begin{bmatrix}0&0\\1&0\end{bmatrix}$ and $B_c=\sigma\begin{bmatrix}0\\1\end{bmatrix}$. We sweep $\sigma$,
solve the I-GRE each time, and read the certified closed-loop second-moment
margin.

```matlab
Abase = [0 1; 1 0.5];  B = [0;1];
fprintf('open-loop deterministic eig: %s\n', mat2str(round(eig(Abase),4)'));
for sigma = [0.0 0.3 0.6 0.9 1.2 1.5 1.8]
    S = mvst.System.ct(Abase, B, 'Q', eye(2), 'R', 1, ...
                       'Ac', sigma*[0 0;1 0], 'Bc', sigma*[0;1]);
    r = S.stabilize();
    fprintf('  sigma=%.2f  K=[% .3f % .3f]  cl-margin=%+.4f  stable=%d\n', ...
            sigma, r.K(1), r.K(2), r.check.margin, r.check.stable);
end
```
```
open-loop deterministic eig: [-0.7808 1.2808]
  sigma=0.00  K=[-2.414 -2.965]  cl-margin=-1.8162  stable=1
  sigma=0.30  K=[-2.228 -2.729]  cl-margin=-1.5483  stable=1
  sigma=0.60  K=[-1.716 -2.066]  cl-margin=-0.7905  stable=1
  sigma=0.90  K=[-1.126 -1.225]  cl-margin=-0.0985  stable=1
  sigma=1.20  K=[-1.000 -0.694]  cl-margin=+0.3056  stable=0
  sigma=1.50  K=[-1.000 -0.444]  cl-margin=+0.5556  stable=0
  sigma=1.80  K=[-1.000 -0.309]  cl-margin=+0.6914  stable=0
```

Interpretation of the sweep: for $\sigma\le0.9$ the I-GRE has a stabilizing solution and
the certified second-moment margin is negative (mean-square stable), but the
margin degrades toward zero as the noise grows, and the gain magnitude shrinks
because pushing harder would only inject more multiplicative noise. Between
$\sigma=0.9$ (margin $-0.0985$) and $\sigma=1.2$ (margin $+0.3056$) the closed
loop crosses into mean-square instability: the control weight
$R+\sum_l B_l^{\top}P B_l$ can no longer be balanced against the state term
enlarged by the noise, the Riccati residual blows up (the solver warns that the residual
exceeds tolerance), and the returned gain no longer certifies. The plant is
**deterministically controllable at every $\sigma$** yet loses **mean-square
stabilizability** past $\sigma\approx1$. This threshold is invisible to a
deterministic design and is exactly what the matrix-valued (second-moment)
formulation captures.

## SOS synthesis: `mvst_stabz_sos`

`mvst_stabz_sos` searches (SOSTOOLS + SeDuMi) for a stabilizing gain by the
bounded-real change of variables of the mathematical setting, covering every
dwell-time type of the chapter (`'arbitrary'`, `'discrete'`, `'periodic'`,
`'minimum'`, `'range'`). It returns `res.feas`, `res.slack`, `res.Kfun` (a
function handle `@(tau) -> K(tau)`, the recovered $U(\tau)Y(\tau)^{-1}$), and
`res.Kf` (the spontaneous-jump gain).

The example below stabilizes an unstable continuous-time drift with a
state-multiplicative Brownian channel under the arbitrary-dwell-time (constant
certificate) program, and certifies the design with the exact lifted-operator
abscissa:

```matlab
sysz = struct('A0',[0 1; 2 0.1], 'B0',[0;1], 'Ac', 0.3*eye(2));
rz = mvst_stabz_sos(sysz, struct('type','arbitrary'));
K  = rz.Kfun(0);
absOpen   = max(real(eig(mvst_liftsym([0 1;2 0.1], {0.3*eye(2)}))));
absClosed = max(real(eig(mvst_liftsym([0 1;2 0.1]+[0;1]*K, {0.3*eye(2)}))));
fprintf('feas=%d  K=[% .4f % .4f]\n', rz.feas, K(1), K(2));
fprintf('open-loop abscissa=%.4f  closed-loop abscissa=%.4f\n', absOpen, absClosed);
```
```
feas=1  K=[-3.5002 -1.4417]
open-loop abscissa=3.0202  closed-loop abscissa=-1.2517
```

The open-loop lifted abscissa is positive (mean-square unstable) and the
closed-loop one is negative (mean-square stable), so the synthesized gain is
certified. `mvst_liftsym` builds the reduced second-moment operator used for the
check.

This routine requires SOSTOOLS and an SDP backend; `S.stabilize` and the case
study above need only base MATLAB.

## Aperiodic sampled-data refresh gains

### `mvst_stabz_sd_range`: range sampling

`res = mvst_stabz_sd_range(plant, opts)` designs the refresh gain of a
held-input sampled-data plant whose sampling period is only known to lie in a
window `[Tmin, Tmax]`, as a matrix sum-of-squares program in the timer.
`opts.mode` is `'scheduled'` (a period-dependent gain `K_J(theta)`), `'constant'`
(a schedule-free gain), or `'slack'` (a constant gain with a period-dependent
certificate, less conservative than `'constant'`); `opts.deg` is the certificate
degree. The plant carries `A`, `Bu`, and optionally the Brownian channels `Ac`,
`Bc`. It requires SOSTOOLS and an SDP backend.

```matlab
plant = struct('A',[0 1;2 0.1],'Bu',[0;1],'Ac',0.2*eye(2),'Bc',0.1*[0;1]);
r = mvst_stabz_sd_range(plant,struct('Tmin',0.2,'Tmax',0.3,'mode','scheduled','deg',4));
fprintf('feas=%d slack=%.4g\n',r.feas,r.slack);
```
```
feas=1 slack=0.005661
```

`res.feas` and `res.slack` are the feasibility flag and certificate slack,
`res.KJfun` the gain as a function of the period, and `res.KJ` a representative
value of it (see `help mvst_stabz_sd_range` for the exact convention).
