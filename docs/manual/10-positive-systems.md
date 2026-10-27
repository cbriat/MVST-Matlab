# 10. Positive systems

[Manual index](index.md) | Prev: [Dissipative synthesis](09-dissipative-synthesis.md) | Next: [Further performance measures](11-performance.md)

This chapter covers the linear positive-systems family: analysis and synthesis
by copositive linear programs. It corresponds to the **Linear Positive Systems**
chapter of the companion monograph. Unlike the SOS layer, this family uses
`linprog` and so needs the **Optimization Toolbox** rather than an SDP backend.

The sections below are written to stand on their own: every problem is stated as
a display before the code that solves it, so the chapter is usable without the
book.

## The mathematical setting

A **linear positive system** is a linear system whose state stays in the
nonnegative orthant: if the initial state is nonnegative, every future state is
nonnegative. Two entrywise matrix notions organize the whole family.

A square matrix $M$ is **nonnegative**, written $M \ge 0$, when all its entries
are nonnegative, and a matrix $A$ is **Metzler** when its off-diagonal entries
are nonnegative:

$$
M \ge 0 \iff M_{ij} \ge 0 \text{ for all } i,j,
\qquad
A \text{ Metzler} \iff A_{ij} \ge 0 \text{ for all } i \ne j .
$$

For the continuous-time system $\dot x = A x$ with resets $x(t_k^+) = J x(t_k^-)$,
positivity is characterized entrywise:

$$
e^{A t} \ge 0 \text{ for all } t \ge 0 \iff A \text{ is Metzler},
\qquad
J x \ge 0 \text{ for all } x \ge 0 \iff J \ge 0 .
$$

So a system is positive exactly when its drift is Metzler at every time and
every reset map is nonnegative. The discrete-time content is carried by the
nonnegative reset $J$ and by spectral-radius conditions such as
$\rho(J e^{A \bar T}) < 1$; there is no separate discrete state equation.

### Linear copositive Lyapunov functions

On the nonnegative orthant a **linear** functional already serves as a Lyapunov
function. For a Metzler $A$ the following are equivalent: $A$ is Hurwitz;
$-A^{-1} \ge 0$; there is a vector $d > 0$ with $A d < 0$; and there is a vector
$p > 0$ with

$$
A^\top p < 0,
\qquad\text{so that}\qquad
V(x) = p^\top x, \quad \dot V = (A^\top p)^\top x < 0 \text{ for } 0 \ne x \ge 0 .
$$

The function $V(x) = p^\top x$ is a **linear copositive Lyapunov function**: it is
positive on the orthant and decreasing along the dynamics. Its existence is a
**linear program** in the nonnegative weight $p$.

### Why linear programs replace semidefinite programs

On the cone of positive semidefinite matrices the certificates are quadratic
($V(X) = \langle P, X \rangle$) and their existence is a semidefinite program.
On the nonnegative orthant everything collapses one level: the quadratic
certificate becomes the **linear** functional $p^\top x$, the LMIs become linear
inequalities, and the trace norm becomes the $1$-norm. The positive-systems
theory is exactly the diagonal slice of the monograph's matrix-valued theory, and
on that slice a semidefinite program becomes a linear program. That is the reason
this family calls `linprog` and not an SDP solver.

The engines read the drift and I/O data at the top level of a flat struct
(`sys.A` a matrix for a single mode, a cell of drifts for a multi-mode system),
not through the `.modes` wrapper. The object method `S.pos(what, opts)` bridges a
`mvst.System` to that contract; the examples below call the `mvst_pos_*`
functions directly, matching the book's example scripts.

## Stability: `mvst_pos_stability`

The continuous-time stability certificate is the copositive linear program

$$
\text{find } p \ge \mathbf{1}
\quad\text{such that}\quad
A^\top p \le -\varepsilon \mathbf{1},
$$

with a small margin $\varepsilon > 0$ and $\mathbf 1$ the all-ones vector (the
lower bound $p \ge \mathbf 1$ normalizes the strictly positive $p$). A
Metzler-Hurwitz drift is certified:

```matlab
r = mvst_pos_stability(struct('A',[-1 0;3 -1]), struct('class','ct'));
fprintf('stable=%d\n', r.stable);
disp('p='); disp(r.p);
```
```
stable=1
p=
    3.0000
    1.0000
```

The solver returns $p = (3.000001, 1)^\top > 0$ (its printed values round to
$(3,1)$), and then $A^\top p = (-10^{-6}, -1)^\top$, strictly negative
entrywise by the LP margin $\varepsilon = 10^{-6}$, so $V(x) = 3 x_1 + x_2$
certifies stability. With the exactly rounded $p = (3,1)^\top$ the first entry is
$A^\top p|_1 = -3 + 3\cdot 1 = 0$, which meets the constraint only
non-strictly; the tiny shift the solver adds is exactly that margin. The
discrete-time analogue, selected with `opts.class = 'dt'`, tests a nonnegative
Schur transition matrix $M$ through

$$
\text{find } p \ge \mathbf{1}
\quad\text{such that}\quad
(M^\top - I) p \le -\varepsilon \mathbf{1}
\quad\Longleftrightarrow\quad M^\top p < p .
$$

## Dwell-time stability

When flow and reset are combined, stability holds only in a dwell-time window
and the certificate carries the **matrix exponential** of the flow over one
dwell-time period.

### Impulsive maximum dwell-time

Take an unstable Metzler flow caught by a contractive nonnegative reset. Over one
period the monodromy is $J e^{A \bar T}$; it is Schur only while the dwell-time
$\bar T$ stays below a maximum $\bar T^\star$. Below that threshold the
copositive LP

$$
p \ge \mathbf{1},
\qquad
A^\top p \ge \varepsilon \mathbf{1},
\qquad
\bigl(J e^{A \bar T}\bigr)^\top p \le p - \varepsilon \mathbf{1}
$$

is feasible: $V = p^\top x$ grows along the expanding flow but the reset more
than compensates over each cycle. With the book's running data
$A = \begin{psmallmatrix} 0.5 & 1 \\ 0 & 0.5 \end{psmallmatrix}$,
$J = \begin{psmallmatrix} 0.1 & 0.2 \\ 0 & 0.1 \end{psmallmatrix}$, the exact
threshold is $\bar T^\star = 2\ln 10 \approx 4.605$:

```matlab
A = [0.5 1; 0 0.5]; J = [0.1 0.2; 0 0.1];
r = mvst_pos_stability(struct('A',A,'J',J,'Tbar',1.5), ...
        struct('class','impulsive','threshold',true));
fprintf('at Tbar=1.5 stable=%d  Tstar=%.4f\n', r.stable, r.Tstar);
```
```
at Tbar=1.5 stable=1  Tstar=4.6052
```

`opts.threshold` also bisects the LP feasibility boundary and returns it in
`res.Tstar`, matching the analytic $2\ln 10$. (When the flow is Hurwitz the
routine solves instead the constant-dwell-time slice $A^\top p \le -\varepsilon\mathbf 1$,
$(J e^{A \bar T})^\top p \le p - \varepsilon\mathbf 1$.)

### Switched minimum dwell-time

For a positive switched system $\dot x = A_{\sigma} x$ with Metzler-Hurwitz modes
that share no common linear copositive Lyapunov function, stability is recovered
above a **minimum** dwell-time. The mode-dependent certificate uses one weight
per mode:

$$
p_i \ge \mathbf{1},
\qquad
A_i^\top p_i \le -\varepsilon \mathbf{1},
\qquad
e^{A_i^\top T_{\min}} p_j \le p_i - \varepsilon \mathbf{1}
\quad\text{for each admissible } i \to j .
$$

The two shear modes below are each Hurwitz, but $A_1 + A_2$ has an unstable
eigenvalue, so fast switching diverges:

```matlab
A1 = [-1 0; 3 -1]; A2 = [-1 3; 0 -1];
r = mvst_pos_stability(struct('A',{{A1,A2}}), ...
        struct('class','switched','dwell','min','Tmin',2.0,'threshold',true));
fprintf('stable=%d  Tstar=%.4f\n', r.stable, r.Tstar);
disp('p1='); disp(r.p{1}); disp('p2='); disp(r.p{2});
```
```
stable=1  Tstar=1.6221
p1=
    3.0000
    1.0000
p2=
    1.0000
    3.0000
```

At $T_{\min} = 2$ the program is feasible; the bisected threshold
$\bar T^\star \approx 1.62$ is the smallest dwell-time at which the copositive LP
stays feasible. Pass `opts.trans` (a $K \times 2$ list of ordered transitions
$[i, j]$) to restrict the admissible switches; the default is all $i \ne j$.

### Markov-jump systems

A positive Markov-jump system with generator $\Pi = (\pi_{ij})$ is certified by
the **coupled** copositive program

$$
A_i^\top p_i + \sum_{j} \pi_{ij}  p_j \le -\varepsilon \mathbf{1},
\qquad p_i \ge \mathbf{1},
\qquad i = 1,\dots,N .
$$

```matlab
A1 = [-2 0.5; 0.3 -2]; A2 = [-1.5 0.2; 0.6 -1.8]; Pi = [-1 1; 1 -1];
r = mvst_pos_stability(struct('A',{{A1,A2}},'Pi',Pi), struct('class','mjls'));
fprintf('stable=%d\n', r.stable);
```
```
stable=1
```

The generator couples the mode weights: unlike the switched case each mode gets
its own $p_i$, but the rows are linked through $\sum_j \pi_{ij} p_j$.

## Induced gains: `mvst_pos_gain`

For a stable positive LTI system

$$
\dot x = A x + E w, \qquad z = C x + F w,
$$

with $A$ Metzler-Hurwitz and $E, C, F \ge 0$, the induced gains are attained at
constant inputs and are computed by linear programs.

### The L1 gain

The induced $L_1 \times \ell_1$ gain is the worst-case amplification of the
integrated (total-mass) signals,

$$
\gamma_1
= \sup_{w \ne 0}
\frac{\displaystyle\int_0^\infty \mathbf{1}^\top z(t)\,dt}
     {\displaystyle\int_0^\infty \mathbf{1}^\top w(t)\,dt},
$$

and it is the value of the copositive linear program

$$
\min_{\gamma, \lambda \ge 0} \gamma
\quad\text{s.t.}\quad
A^\top \lambda + C^\top \mathbf{1} \le 0,
\qquad
E^\top \lambda + F^\top \mathbf{1} \le \gamma \mathbf{1} .
$$

The optimum is the vertex $\lambda^\star = -(A^\top)^{-1} C^\top \mathbf{1}$,
giving $\gamma_1 = \max$-channel of $E^\top \lambda^\star + F^\top \mathbf{1}$.
In the **scalar (SISO)** case this equals the DC gain $C(-A)^{-1} E$ (the area
under the nonnegative impulse response); for a MIMO system $C(-A)^{-1}E$ is a
matrix and $\gamma_1$ is its largest **column** channel, not a single scalar.

```matlab
r = mvst_pos_gain(struct('A',[-0.9 0;0.5 -0.6],'E',[1;0],'C',[1 1]), ...
                  struct('norm','L1'));
fprintf('L1 gain = %.6f  gammaDC=%.6f\n', r.gamma, r.gammaDC);
```
```
L1 gain = 2.037037  gammaDC=2.037037
```

The LP value and the closed-form DC gain agree, as they must for a positive
system. `res.weights` returns the dual certificate $\lambda^\star$.

### The Linf (peak-to-peak) gain

The induced $L_\infty \times \ell_\infty$ (peak-to-peak) gain amplifies the
peak amplitudes,

$$
\gamma_\infty
= \sup_{w \ne 0}
\frac{\sup_t \|z(t)\|_\infty}{\sup_t \|w(t)\|_\infty},
$$

and is the value of the dual (primal-tube) copositive program

$$
\min_{\gamma, \xi \ge 0} \gamma
\quad\text{s.t.}\quad
A \xi + E \mathbf{1} \le 0,
\qquad
C \xi + F \mathbf{1} \le \gamma \mathbf{1},
$$

with optimum $\xi^\star = (-A)^{-1} E \mathbf{1}$ and
$\gamma_\infty = \max(C \xi^\star + F \mathbf{1})$. In the **scalar (SISO)** case
this is again the DC gain, so $\gamma_1 = \gamma_\infty = C(-A)^{-1}E$; for a
MIMO system $\gamma_\infty$ is the largest **row** channel of the DC-gain matrix,
so in general $\gamma_1$ (worst column) and $\gamma_\infty$ (worst row) differ
and coincide only in the scalar case shown here:

```matlab
r = mvst_pos_gain(struct('A',[-0.9 0;0.5 -0.6],'E',[1;0],'C',[1 1]), ...
                  struct('norm','Linf'));
fprintf('Linf gain = %.6f\n', r.gamma);
```
```
Linf gain = 2.037037
```

Select the norm with `opts.norm = 'L1'` or `'Linf'`.

## Dissipativity with a general linear supply: `mvst_pos_supply`

Dissipativity of a positive system is tested against a **general linear** supply
rate (rather than a quadratic one),

$$
\sigma(w,z) = q^\top z + r^\top w,
$$

with weights $q, r$ of arbitrary sign. A linear storage $V = \lambda^\top x$,
$\lambda \ge 0$, certifies dissipativity through the copositive program

$$
A^\top \lambda - C^\top q \le 0,
\qquad
E^\top \lambda - F^\top q - r \le 0
$$

(with two further rows $(J^\top - I)\lambda - C_d^\top q_d \le 0$ and
$E_d^\top \lambda - F_d^\top q_d - r_d \le 0$ when a jump channel is present).
Choosing $q = -\mathbf 1$, $r = \mathbf 1$ recovers the mass-balance
(passivity) supply $\sigma = \mathbf 1^\top w - \mathbf 1^\top z$; choosing
$q = -\mathbf 1$, $r = \gamma \mathbf 1$ recovers the $L_1$-gain program.

The supply is given through its coefficients `qc` (on the output) and `rc` (on
the disturbance). For the compartmental system the supply is satisfiable only
once the disturbance weight `rc` clears the induced-gain threshold:

```matlab
sys = struct('A',[-0.9 0;0.5 -0.6],'Ec',[1;0],'Cc',[1 1]);
r1 = mvst_pos_supply(sys, struct('supply',struct('qc',-1,'rc',1)));
r3 = mvst_pos_supply(sys, struct('supply',struct('qc',-1,'rc',3)));
fprintf('rc=1 feas=%d   rc=3 feas=%d\n', r1.feas, r3.feas);
```
```
rc=1 feas=0   rc=3 feas=1
```

At `rc = 1` the supply is not satisfiable; at `rc = 3` it is, consistent with the
induced gain (2.037) lying between those two levels. When feasible, `res.lambda`
is the storage weight.

## Positivity-preserving stabilization: `mvst_pos_stabilize`

The goal is a state feedback $u = K x$ that makes the closed loop **both**
positive and stable. Positivity requires the closed-loop drift $A + B_c K$ to be
Metzler, which is bilinear in $(K, \text{certificate})$. The copositive change of
variables $X = \mathrm{diag}(\lambda)$ with $\lambda > 0$ and $U_c = K X$
linearizes it: $A + B_c K$ is Metzler iff $A X + B_c U_c$ has nonnegative
off-diagonal entries. The design is then the linear program

$$
(A X + B_c U_c)_{ij} \ge 0 \quad (i \ne j),
\qquad
(A X + B_c U_c)\mathbf{1} \le -\varepsilon \mathbf{1},
$$

in $(\lambda, U_c)$, with $X = \mathrm{diag}(\lambda)$. The first rows force the
closed loop Metzler (positivity preserved); the last is the copositive Lyapunov
decrease $(A + B_c K)\lambda < 0$. The gain is recovered as $K = U_c X^{-1}$:

```matlab
r = mvst_pos_stabilize(struct('A',[1 2; 0 1],'Bc',eye(2)), struct('class','ct'));
fprintf('feas=%d\n', r.feas); fprintf('K =\n'); disp(r.K);
```
```
feas=1
K =
   -1.0010   -2.0000
         0   -1.0010
```

The closed-loop drift $A + B_c K$ is Metzler with negative diagonal, hence
positive and Hurwitz. When `opts.certify` is true (the default) the design is
re-checked with `mvst_pos_stability` and returned in `res.cert`. A single shared
input column cannot in general make an arbitrary $A$ Metzler AND Hurwitz, so
adequate actuation is required for feasibility. The impulsive class additionally
designs a jump gain $K_d$ so that $J + B_d K_d \ge 0$ under a constant
dwell-time.

## Worked case study: an epidemic-like positive system

We put the pieces together on a small compartmental model. Three coupled regions
carry an infection; near the disease-free equilibrium the linearized prevalences
$x = (x_1, x_2, x_3)$ obey $\dot x = A x + E w$ with

$$
A = \begin{bmatrix} -1.6 & 0.5 & 0.2 \\ 0.4 & -1.8 & 0.6 \\ 0.3 & 0.4 & -2.0 \end{bmatrix},
\qquad E = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix},
\qquad C = \begin{bmatrix} 1 & 1 & 1 \end{bmatrix} .
$$

The off-diagonal entries (cross-region transmission and mobility) are
nonnegative, so $A$ is Metzler and the model is positive. The disturbance $w$ is
an imported-infection rate into region 1, and the output $z = x_1 + x_2 + x_3$ is
the total prevalence.

**Step 1: stability below the epidemic threshold.** The copositive stability LP
certifies the disease-free equilibrium:

```matlab
Ae = [-1.6 0.5 0.2; 0.4 -1.8 0.6; 0.3 0.4 -2.0];
rs = mvst_pos_stability(struct('A',Ae), struct('class','ct'));
fprintf('stable=%d, margin=%.4f\n', rs.stable, rs.margin);
disp('copositive weight p='); disp(rs.p);
```
```
stable=1, margin=0.9000
copositive weight p=
     1
     1
     1
```

The uniform weight $p = \mathbf 1$ works: $V(x) = x_1 + x_2 + x_3$ (the total
prevalence itself) is a Lyapunov function, and $A^\top \mathbf 1 = (-0.9, -0.9,
-1.2)^\top < 0$ gives the decrease margin $0.9$.

**Step 2: the L1 gain (total infection burden).** The induced $L_1$ gain is the
total infection-days produced per unit of imported infection:

```matlab
rg = mvst_pos_gain(struct('A',Ae,'E',[1;0;0],'C',[1 1 1]), struct('norm','L1'));
fprintf('L1 gain = %.6f\n', rg.gamma);
```
```
L1 gain = 1.061947
```

**Step 3: an above-threshold outbreak, stabilized with positivity preserved.**
Weaken recovery so the epidemic is now unstable, and design a per-region
intervention $u = K x$ (added removal, actuated by $B_c = I$) that returns the
system to a positive, stable regime. A stronger decrease margin `eps` is asked
for so the controlled outbreak decays visibly:

```matlab
Au = [-0.3 0.5 0.2; 0.4 -0.4 0.6; 0.3 0.4 -0.5];
fprintf('open-loop eig = %s\n', mat2str(round(eig(Au)',4)));
rz = mvst_pos_stabilize(struct('A',Au,'Bc',eye(3)), struct('class','ct','eps',0.5));
disp('gain K='); disp(rz.K);
Acl = Au + rz.K;
fprintf('closed-loop Metzler offdiag min = %.4f\n', min(Acl(~eye(3))));
fprintf('closed-loop eig = %s\n', mat2str(round(eig(Acl)',4)));
fprintf('re-certified stable = %d\n', rz.cert.stable);
```
```
open-loop eig = [0.4098 -0.6699 -0.9398]
gain K=
   -0.2000   -0.5000   -0.2000
   -0.4000   -0.1000   -0.6000
   -0.3000         0   -0.4000
closed-loop eig = [-0.5 -0.9 -0.5]
re-certified stable = 1
```

The open loop has the unstable mode $+0.41$ (an ongoing outbreak). The designed
gain yields the closed-loop drift

$$
A_u + K = \begin{bmatrix} -0.5 & 0 & 0 \\ 0 & -0.5 & 0 \\ 0 & 0.4 & -0.9 \end{bmatrix},
$$

which is Metzler (its off-diagonal minimum is $0$) and Hurwitz (eigenvalues
$-0.5, -0.9, -0.5$). Positivity is preserved and the outbreak is driven to
extinction; `mvst_pos_stability` re-certifies the result.

## Multi-mode systems and the object bridge

For a switched or Markov-jump positive system, pass a cell of per-mode drifts;
the coupled Markov-jump program additionally reads the generator `Pi`. A
multi-mode system carrying a common reset has no home in this family and is
refused with `mvst:System:posUnsupported` rather than being silently misrouted.
The internal LP wrapper is `mvst_pos_lp`, and the analysis front end that
mirrors the book's example scripts is the `S.pos(what, opts)` sugar method on a
`mvst.System`, which routes to the `mvst_pos_*` family.

## Further positive-systems engines

The following engines complete the family. They are called directly, or reached
through the routers `mvst_pos_gain`, `mvst_pos_stabilize`, and `mvst_pos_ms`,
which dispatch to the timer-dependent (SOS) backends when `opts.dwell` is
`'min'` or `'range'`.

### `mvst_pos_oc`: linear-cost optimal control

`res = mvst_pos_oc(sys, opts)` solves the linear-cost optimal control of a
positive system by the copositive costate programs: the value is linear,
`V(x) = p'x`, and the optimal costate `p` solves a copositive program.
`opts.mode` selects `'finite'` (horizon `opts.T`, jump times `opts.tk`, adjoint
costate recursion), `'inf'` (LTI, bounded costate `p = -(A')^{-1} q`), or
`'dwell'` (periodic dwell-time); `opts.q` is the nonnegative running cost
weight and `opts.x0` the initial state.

```matlab
r = mvst_pos_oc(struct('A',[-1 0;1 -2]),struct('mode','inf','q',[1;1],'x0',[1;1]));
fprintf('feas=%d Jopt=%.4f p0=%s\n',r.feas,r.Jopt,mat2str(r.p0.',4));
```
```
feas=1 Jopt=2.0000 p0=[1.5 0.5]
```

### `mvst_pos_ms`: second-moment certificates for stochastic positive systems

`res = mvst_pos_ms(sys, opts)` certifies mean-square properties of positive
systems driven by Brownian and Poisson noise, which the mean-level LP layer
cannot see. `opts.what` is `'stability'` (default), `'tracegain'`, or
`'peakgain'`, and `opts.cert` is `'spn'` (default, a split into a
semidefinite and an entrywise-nonnegative part) or `'diag'` (diagonal quadratic
certificate). The constant certificates are SOSTOOLS programs by default (any backend,
including `'nosdp'`); the legacy YALMIP formulation is `opts.engine = 'yalmip'`
(optional external dependency).

```matlab
sys = struct('A0',[-1 0;1 -2],'Ac',sqrt(1.5)*eye(2));
r = mvst_pos_ms(sys,struct('what','stability','cert','diag'));
fprintf('feas=%d\n',r.feas);
```
```
feas=1
```

`res.feas` is the certificate flag, `res.D` the diagonal certificate (or `res.P`
for `'spn'`). For this diagonal system the certificate exists for noise
intensities up to a threshold computed in closed form (the test suite checks it
at `s* = 1.6667` for `Ac = sqrt(s) I`).

### `mvst_pos_ms_dwell`: timer-dependent second-moment certificates

`mvst_pos_ms_dwell(sys, what, cert, opts)` is the minimum / range dwell-time
backend of `mvst_pos_ms`, with `opts.dwell`, `opts.Tmin`, `opts.Tmax`, and the
polynomial degree `opts.deg` of the timer certificate. It is normally reached
through `mvst_pos_ms`:

```matlab
s = struct('A0',[-1 0.5;0.2 -1.5],'Ac',{{0.2*eye(2)}},'J',{{[0.6 0;0.1 0.5]}});
r = mvst_pos_ms(s,struct('what','stability','cert','diag','dwell','min','Tmin',0.5,'deg',2));
fprintf('feas=%d\n',r.feas);
```
```
feas=1
```

With `deg = 0` the timer certificate collapses to the constant (arbitrary
dwell-time) certificate, so the verdict then equals that of the constant
program; a higher degree can only enlarge the feasible set. It requires
SOSTOOLS and an SDP backend.

### `mvst_pos_gain_sos`: timer-dependent induced gains

`mvst_pos_gain_sos(sys, opts)` computes the `L_inf` x `ell_inf` (peak) and `L_1`
x `ell_1` (trace) induced gains of a positive impulsive or switched system under
a minimum or range dwell-time, by the copositive-SOS forms of the book's sliced
conditions (forward reachability tube for the peak gain, backward adjoint for the
trace gain). It is reached through `mvst_pos_gain`:

```matlab
s = struct('A',[-1 0;1 -2],'Ec',[0.1;1.1],'Cc',[0 1],'Fc',0.3, ...
           'J',[0.1 1;0.1 0.1],'Ed',[0.3;0.3],'Cd',[0 1],'Fd',0.03);
r = mvst_pos_gain(s,struct('norm','Linf','class','impulsive','dwell','min','Tmin',1,'deg',4));
fprintf('feas=%d gamma=%.4f\n',r.feas,r.gamma);
```
```
feas=1 gamma=1.1356
```

The gain is non-increasing in `Tmin` (a longer guaranteed dwell-time can only help).

### `mvst_pos_stabilize_sos`: timer-dependent copositive synthesis

`mvst_pos_stabilize_sos(sys, opts)` designs a flow gain `Kc(tau)` and a jump
gain `Kd` that make the closed loop both positive and exponentially stable under
a constant, minimum, or range dwell-time, through the copositive change of
variables of the book (so the conditions are linear at each timer value). With
`opts.norm = 'Linf'` and disturbance/output data present it additionally
minimizes the closed-loop induced `L_inf` gain. It is reached through
`mvst_pos_stabilize`:

```matlab
s = struct('A',[3 -1;2 -1],'Ec',[0.5;0.1],'Cc',[0 1],'Fc',0.1,'J',[2 1;0 1.2], ...
           'Ed',[0.8;0.3],'Cd',[0 1],'Fd',0.1,'Bc',[1;0],'Bd',[1;1]);
r = mvst_pos_stabilize(s,struct('class','impulsive','dwell','min','Tmin',0.2,'deg',2,'norm','Linf'));
fprintf('feas=%d gamma=%.4f\n',r.feas,r.gamma);
```
```
feas=1 gamma=0.5981
```

### `mvst_sosmatrvar` and `mvst_sosposmatrvar`: matrix variables for SOSTOOLS programs

These two helpers declare polynomial matrix variables in a SOSTOOLS program and
are used by the copositive-SOS engines above. `[prog,P] = mvst_sosmatrvar(prog,
Z, [sp sq])` declares a free (sign-indefinite) `sp`-by-`sq` polynomial matrix in
the monomial vector `Z`. `[prog,Q] = mvst_sosposmatrvar(prog, n, d, vars)`
declares an `n`-by-`n` matrix of degree `d` that is positive semidefinite for
every value of `vars` by construction (a scalar sum-of-squares polynomial when
`n = 1`). Both are independent reimplementations, on the modern SOSTOOLS API, of
helpers of M. Peet's MatrixTools; the attribution is in their help text.

```matlab
syms t; prog = sosprogram(t);
[prog,Q] = mvst_sosposmatrvar(prog,2,2,t);              % 2x2 SOS matrix in t
[prog,w] = mvst_sosmatrvar(prog,monomials(t,0:2),[2 1]); % free 2x1 polynomial vector
```
