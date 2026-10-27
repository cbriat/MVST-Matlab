# 11. Further performance measures

[Manual index](index.md) | Prev: [Positive systems](10-positive-systems.md) | Next: [Book reproduction](12-book-reproduction.md)

This chapter covers the further performance measures (`mvst_perf`, method
`S.perf`, and its dwell-time backend `mvst_perf_dwell`), which correspond to the
**Further Performance Measures** chapter of the companion monograph. Every
quantity is defined below as a display before the code that computes it, so the
chapter stands on its own.

## The performance measures

### Signal norms

The performance measures are induced gains and averages built on the standard
signal norms. For a signal $s$,

$$
\|s\|_{L_2}^2 = \int_0^\infty \|s(t)\|^2\, dt
\qquad(\text{energy}),
\qquad
\|s\|_{L_\infty} = \sup_{t \ge 0} \|s(t)\|
\qquad(\text{peak}),
$$

with $\ell_2$ and $\ell_\infty$ the discrete counterparts summed / maximized over
the jump instants. The hybrid norms combine the flow ($L$) and jump ($\ell$)
parts of a mixed continuous-discrete signal.

### The Gramian

The covariance measures are obtained from the **controllability (reachability)
Gramian** $\Pi$, the state second moment under a unit-intensity disturbance
$w$. For a Hurwitz continuous-time drift it is the steady solution of the
Lyapunov equation

$$
A\Pi + \Pi A^\top + B_w B_w^\top = 0,
$$

and for the constant dwell-time impulsive system it is the periodic solution
$\Pi(\tau) = e^{A\tau}\Pi_0 e^{A^\top\tau} + \int_0^\tau e^{As} B_w B_w^\top e^{A^\top s}\,ds$
with the reset $\Pi_0 = J\Pi(\bar T)J^\top$ closed over one period. The output
covariance is

$$
S(\tau) = C  \Pi(\tau)  C^\top,
$$

and the measures are its trace and its largest eigenvalue over a period.

### The measures, defined

`S.perf(measure)` computes a named performance measure, routing per system class
to the matching exact engine (Riccati / lifted-eigenvalue for the induced gains,
Lyapunov Gramian for the covariance measures). The measures are the following.

The **mean-square L2 gain** `'l2'` is the worst-case energy amplification (the
$H_\infty$-type gain), delegated to `mvst_l2gain`:

$$
\gamma_{L_2} = \sup_{w \ne 0} \frac{\|z\|_{L_2}}{\|w\|_{L_2}} .
$$

The **trace gain** `'tracegain'` is its square,

$$
\gamma_{\mathrm{tr}} = \gamma_{L_2}^2
= \sup_{w \ne 0} \frac{\|z\|_{L_2}^2}{\|w\|_{L_2}^2} .
$$

The **energy-to-peak (generalized H2) gain** `'e2p'` is the induced
$L_2 \to L_\infty$ gain,

$$
\gamma_{\mathrm{e2p}} = \sup_{w \ne 0} \frac{\|z\|_{L_\infty}}{\|w\|_{L_2}}
= \sqrt{\sup_\tau \lambda_{\max}\bigl(C \Pi(\tau) C^\top\bigr)} .
$$

The **peak gain** `'peak'` is its square,

$$
\gamma_{\mathrm{peak}} = \gamma_{\mathrm{e2p}}^2
= \sup_\tau \lambda_{\max}\bigl(C \Pi(\tau) C^\top\bigr) .
$$

The **H2 cost** `'h2'` is the average output power,

$$
\gamma_{H_2} = \sqrt{\frac{1}{\bar T}\int_0^{\bar T} \mathrm{tr}\bigl(C \Pi(\tau) C^\top\bigr)\,d\tau},
$$

which reduces to the classical $H_2$ norm in the jump-free case. The
**deterministic peak-to-peak (star-norm)** `'p2p'` is the induced
$L_\infty \to L_\infty$ gain; it is not covered by an exact engine (see below).

```matlab
S = mvst.System.ct([-1 0.5; 0 -2], zeros(2,0), 'Bw', [1;0], 'C', [1 0]);
fprintf('h2 = %.6f\n', S.perf('h2').gamma);
fprintf('e2p = %.6f\n', S.perf('e2p').gamma);
fprintf('peak = %.6f\n', S.perf('peak').gamma);
```
```
h2 = 0.707107
e2p = 0.707107
peak = 0.500000
```

Note `peak` (0.5) is the square of `e2p` (0.7071), as defined. For a scalar
output $C\Pi C^\top$ is a scalar, so its trace and largest eigenvalue coincide
and `h2 = e2p`; with a vector output they differ. `res.gamma` carries the measure
value and `res.measure` / `res.class` the resolved names.

### Coverage matrix

Not every measure is characterized for every class by the exact engines. The
dispatcher computes the pairs below and raises `mvst:perf:unsupported` (naming
the nearest supported pair) for the rest:

| measure | ct | mjls | switched | impulsive/sd | poisson |
|---|---|---|---|---|---|
| `l2` | yes | yes | yes | yes | yes |
| `tracegain` | yes | yes | yes | yes | yes |
| `e2p` | yes | no | no | yes | no |
| `peak` | yes | no | no | yes | no |
| `h2` | yes | no | no | yes | no |

The induced gains (`l2`, `tracegain`) inherit every class of `mvst_l2gain`. The
Gramian measures (`e2p`, `peak`, `h2`) are exact for the continuous-time system
and its constant dwell-time impulsive / sampled instance. Two further families
are routed to the dwell-time backend `mvst_perf_dwell` (next section): the
minimum and range dwell-time versions of `e2p`, `peak`, and `h2` (select with
`opts.dwell`), and the deterministic peak-to-peak (star-norm) measure `p2p`,
which needs an invariant-set LMI with a free decay rate. Both need an SOS or SDP
backend. Pairs that are covered by none of these (for example the multi-mode
stochastic covariance measures) are refused deliberately: rather than silently
returning an approximate or wrong number, the dispatcher raises
`mvst:perf:unsupported` and points at the closest supported measure.

## Dwell-time and peak-to-peak measures

### `mvst_perf_dwell`

`res = mvst_perf_dwell(sys, measure, opts)` is the backend of `mvst_perf` for
the timer-dependent measures; it is normally reached through `S.perf` or
`mvst_perf`. The measure is `'peak'`, `'e2p'`, or `'h2'` with `opts.dwell` in
`{'min','range'}` (and `opts.Tmin`, `opts.Tmax`, `opts.deg` for the degree of the
timer certificate), or `'p2p'` / `'peakfull'` for the deterministic peak-to-peak
and amplitude-bounded peak gains. The dwell-time tubes and, by default, the
decay-rate certificates (a two-parameter bisection) use SOSTOOLS; the legacy YALMIP
formulation of the decay-rate certificates is `opts.engine = 'yalmip'` (optional
external dependency, not bundled).

```matlab
s = struct('A',[-0.2 1;-1 -0.2],'Bw',[0.5;0.3],'C',[1 0],'J',0.6*eye(2));
r = mvst_perf(s,'h2',struct('dwell','min','Tmin',1,'deg',4));
fprintf('h2 = %.4f\n',r.gamma);
```
```
h2 = 0.7105
```

The same system under the deterministic peak-to-peak measure:

```matlab
r = mvst_perf(s,'p2p',struct('dwell','min','Tmin',1,'deg',4));
fprintf('p2p = %.4f (lower bound field gammaLB)\n',r.gamma);
```
```
p2p = 2.1523 (lower bound field gammaLB)
```

`res.gamma` is the certified value (an upper bound for `p2p`, with the DC-gain
companion lower bound in `res.gammaLB`), and `res.lambda`, `res.nu` the decay
parameters at the optimum. The `p2p` call takes noticeably longer than the
others because of the bisection.

## Worked case study: a performance sweep

We take a lightly damped second-order mode (a resonant flexible mode) with
natural frequency $\omega_0 = 1$ and damping ratio $\zeta = 0.15$,

$$
A = \begin{bmatrix} 0 & 1 \\ -\omega_0^2 & -2\zeta\omega_0 \end{bmatrix},
\qquad B_w = \begin{bmatrix} 0 \\ 1 \end{bmatrix},
\qquad C = \begin{bmatrix} 1 & 0 \end{bmatrix},
$$

with the disturbance entering the velocity channel, measure its performance, and
sweep the damping.

**Step 1: measure the performance.**

```matlab
w0 = 1.0; zeta = 0.15;
A = [0 1; -w0^2 -2*zeta*w0]; Bw = [0;1]; C = [1 0];
S = mvst.System.ct(A, zeros(2,0), 'Bw', Bw, 'C', C);
fprintf('h2   = %.6f\n', S.perf('h2').gamma);
fprintf('e2p  = %.6f\n', S.perf('e2p').gamma);
fprintf('peak = %.6f\n', S.perf('peak').gamma);
fprintf('l2   = %.6f\n', S.perf('l2').gamma);
```
```
h2   = 1.290994
e2p  = 1.290994
peak = 1.666667
l2   = 3.371521
```

**Step 2: sweep the damping ratio.** As $\zeta$ falls, the resonance sharpens
and every gain grows; the $L_2$ ($H_\infty$) gain, which sees the resonant peak
directly, grows fastest:

```matlab
fprintf('%6s %10s %10s %10s\n','zeta','h2','e2p','l2');
for z = [0.10 0.15 0.25 0.40 0.70]
    Az = [0 1; -w0^2 -2*z*w0];
    Sz = mvst.System.ct(Az, zeros(2,0), 'Bw', Bw, 'C', C);
    fprintf('%6.2f %10.6f %10.6f %10.6f\n', z, ...
            Sz.perf('h2').gamma, Sz.perf('e2p').gamma, Sz.perf('l2').gamma);
end
```
```
  zeta         h2        e2p         l2
  0.10   1.581139   1.581139   5.025208
  0.15   1.290994   1.290994   3.371521
  0.25   1.000000   1.000000   2.065613
  0.40   0.790569   0.790569   1.363892
  0.70   0.597614   0.597614   1.000245
```

The $H_2$ and energy-to-peak costs coincide because the output is scalar (see the note on scalar outputs above), while the $L_2$ gain climbs from $1.0$ at $\zeta = 0.7$ to $5.0$ at
$\zeta = 0.1$: a lightly damped mode is far more sensitive to worst-case
disturbance energy than the average-power measures suggest.

The exact engines used in this case study need only base MATLAB; the dwell-time
backend `mvst_perf_dwell` needs SOSTOOLS and an SDP backend.
