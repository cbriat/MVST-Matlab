# 6. Steering

[Manual index](index.md) | Prev: [Optimal control](05-optimal-control.md) | Next: [Dissipativity](07-dissipativity.md)

This chapter covers finite-time optimal covariance steering: driving the state
second-moment covariance to a prescribed target over a horizon. It corresponds
to the book chapter **Steering** of the companion monograph. The section below
states the problem in full so the chapter is self-contained.

## Mathematical setting

**The problem.** Take the controlled diffusion of the previous chapter,

$$
\mathrm{d}x = (A x + B u)\,\mathrm{d}t
  + \sum_l (A_l x + B_l u)\,\mathrm{d}W_l ,
\qquad x(0) \sim ( \cdot , \Sigma_0) ,
$$

and track only its state covariance

$$
\Sigma(t) = \mathbb{E}\big[x(t)x(t)^{*}\big] \in \mathbb{H}^{n}_{\succeq 0} .
$$

Covariance steering asks for a control policy that moves this Hermitian
positive-semidefinite matrix from a given start to a given end over a fixed
horizon:

$$
\text{find } u(\cdot) \text{ such that }
\Sigma(0) = \Sigma_0  \longrightarrow  \Sigma(T) = \Sigma_f ,
$$

with `Sigma_f = Sigma_f^* >= 0` the prescribed terminal covariance. Under a
linear state feedback `u = K(t) x` the closed-loop covariance obeys the matrix
Lyapunov (Riccati-type) flow

$$
\dot\Sigma = (A + B K)\Sigma + \Sigma (A + B K)^{*}
   + \sum_l (A_l + B_l K)\Sigma (A_l + B_l K)^{*} ,
\qquad \Sigma(0) = \Sigma_0 .
$$

**Reachability.** The classical controllability Gramian below is only the
noiseless, vector-state **intuition**: the true reachability object here is
operator-valued (the reachable set of the matrix covariance flow), not a single
Gramian matrix. As a special-case guide, for the noiseless controllable model
the reachable set is governed by the controllability Gramian

$$
W_c(T) = \int_0^T e^{A t} B B^{*} e^{A^{*} t}\,\mathrm{d}t \succ 0 .
$$

When `(A,B)` is controllable, `W_c(T) > 0` for every `T > 0` and a large family
of targets `Sigma_f` is reachable; targets that are "too small" relative to the
process noise, or that ask for contraction faster than the input authority
allows, sit outside the feasible set and no exact steering law exists. Because
the underlying object is operator-valued, reachability is in practice confirmed by
**convergence of the solver** rather than by forming `W_c(T)` explicitly.

**The feedback construction.** Steering is solved as a two-point boundary-value
problem. It minimizes the control energy $\int_0^T u^{*} R u\,\mathrm{d}t$
subject to the boundary covariances, so `R = R^* > 0` is the **control weight**
of that underlying minimum-energy cost. A symmetric terminal multiplier `M` seeds
the **backward** generalized Riccati equation

$$
-\dot P = \operatorname{He}(A^{*}P) + \sum_l A_l^{*} P A_l
   - \big(P B + \textstyle\sum_l A_l^{*}P B_l\big) \Theta^{-1} 
     \big(P B + \textstyle\sum_l A_l^{*}P B_l\big)^{*} ,
\qquad P(T) = M ,
$$

with the **noise-weighted** control term
$\Theta = R + \sum_l B_l^{*} P B_l$: the engine `mvst_ric_ct` inverts exactly
this $\Theta$, not the bare `R`, so that the input-multiplicative channels raise
the effective control weight. The equation produces the time-varying gain

$$
K(t) = -\Theta^{-1}\big(B^{*}P(t) + \textstyle\sum_l B_l^{*}P(t)A_l\big) ,
\qquad
\Theta = R + \textstyle\sum_l B_l^{*}P(t)B_l .
$$

The state term is displayed here in its **minimum-energy specialization** (no
state penalty); the engine reads the mode's cost `md.Q` and adds it to the state
block when present.

The gain is then propagated **forward** through the closed-loop Lyapunov flow
above to obtain `Sigma(T; M)`. The steering law is the multiplier `M` at which
the terminal covariance hits the target:

$$
\Phi(M) := \Sigma(T; M) - \Sigma_f = 0 .
$$

This root-finding map `M -> Sigma(T; M)` is the forward-backward (FB) structure;
its Jacobian is dense and the map is nonlinear, so the toolbox solves
`Phi(M) = 0` by a **damped Newton iteration**. **Standing assumptions:** `(A,B)`
controllable (or stabilizable with `Sigma_f` in the reachable set), `Sigma_0`
and `Sigma_f` Hermitian positive semidefinite, `R > 0`.

## `S.steer`

`S.steer(Xf, T)` steers the state covariance from its initial value to the
target `Xf` over `[0,T]` by exactly this forward-backward algorithm: for a
symmetric terminal multiplier `M` it runs the backward I-GRE from `P(T) = M` to
obtain the time-varying gains, propagates the closed-loop covariance forward,
and solves `Sigma(T; M) = Xf` by the damped Newton core `mvst_steer_newton`.

The initial covariance defaults to the mode's `X0`. The example steers a
two-state diffusion from `Sigma(0) = I` to `Sigma(T) = 0.5 I` over one time
unit, with a full-rank input so the target is reachable:

```matlab
S = mvst.System.ct([-1 0;0 -1], eye(2), 'X0', eye(2));
r = S.steer(0.5*eye(2), 1);
fprintf('iters=%d  final gap=%.3e  class=%s\n', r.iters, r.res, r.class);
fprintf('Sigma(0) =\n'); disp(r.Sigma(:,:,1));
fprintf('Sigma(T) =\n'); disp(r.Sigma(:,:,end));
```
```
iters=8  final gap=2.218e-08  class=ct
Sigma(0) =
     1     0
     0     1
Sigma(T) =
    0.5000         0
         0    0.5000
```

Newton converges in 8 iterations to a terminal gap of about 2e-8, and the final
covariance is the requested `0.5 I`.

Result fields: `res.M` (optimal terminal multiplier), `res.t` and `res.Sigma`
(the covariance trajectory, an `n x n x nt` array), `res.K` (the optimal gain
object, shape depending on class), `res.res` (final gap norm), `res.iters`, and
`res.hist` (gap-norm history).

## Worked case study: an anisotropic target and the covariance trajectory

The isotropic example above hides the interesting behaviour, because the two
axes move together. This case study steers a non-normal, weakly damped diffusion
to an **anisotropic, correlated** target and inspects the whole covariance path.
The target has off-diagonal correlation and eigenvalues `0.19` and `0.51`, so
the two principal axes must contract by different amounts.

The default Newton seed `M0 = 0` is isotropic and can stall on a hard,
anisotropic target (the trial step drives the backward Riccati into finite
escape and the line search stalls); a seed in the solution basin restores
convergence, so `OPTS.M0` is set here.

```matlab
S  = mvst.System.ct([-0.5 1; 0 -0.5], eye(2), 'X0', eye(2));
Xf = [0.4 0.15; 0.15 0.3];
r  = S.steer(Xf, 2, struct('M0', -0.5*eye(2)));
fprintf('iters=%d  final gap=%.3e  class=%s\n', r.iters, r.res, r.class);
fprintf('Sigma(0) eig = %s\n', mat2str(sort(eig(r.Sigma(:,:,1)))',4));
fprintf('Sigma(T) eig = %s\n', mat2str(sort(eig(r.Sigma(:,:,end)))',4));
fprintf('Sigma(T) =\n'); disp(r.Sigma(:,:,end));
fprintf('||Sigma(T) - Xf|| = %.3e\n', norm(r.Sigma(:,:,end)-Xf,'fro'));
```
```
iters=13  final gap=9.182e-08  class=ct
Sigma(0) eig = [1 1]
Sigma(T) eig = [0.1919 0.5081]
Sigma(T) =
    0.4000    0.1500
    0.1500    0.3000
||Sigma(T) - Xf|| = 1.179e-07
```

Interpretation of the run: the covariance starts isotropic (both eigenvalues `1`), and by
`t = T` it has been shaped into the requested correlated target, with the two
principal axes now at `0.19` and `0.51`. The terminal covariance reproduces
`Xf` to a Frobenius error of 1.2e-7, and the damped Newton needs 13 iterations
(more than the isotropic case, as expected for the harder off-diagonal target).
The full trajectory `r.Sigma` is available for plotting; the companion helpers
below draw it.

## Other classes and hard targets

`S.steer` dispatches on class: `ct`, `mjls`, `switched`, `sd` (sampled-data),
and `poisson` are supported, with time-varying variants detected automatically
when a drift is a function handle. Pass options through the fourth argument:

- `OPTS.X0`: initial covariance (required for `switched`, `sd`, `poisson`);
- `OPTS.M0`: initial Newton multiplier. As the worked example shows, for hard
  targets the default `M0 = 0` can leave the damped Newton stranded at a large
  terminal gap; an `M0` in the solution basin restores quadratic convergence.
  The sampled-data and Poissonian-jump steering examples set it for this reason;
- `OPTS.tol` (default 1e-7) and `OPTS.maxit` (default 60).

The companion figure helpers `mvst_steer_ellipse` and `mvst_steer_components`
draw the book's steering figures from a `res` struct: the first plots the
covariance confidence ellipses at start, middle, and end, the second plots the
individual covariance entries against time. Steering uses only base MATLAB (no
SDP backend).
