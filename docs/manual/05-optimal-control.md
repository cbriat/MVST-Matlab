# 5. Optimal control

[Manual index](index.md) | Prev: [Stabilization](04-stabilization.md) | Next: [Steering](06-steering.md)

This chapter covers linear-cost optimal control over the three horizons and the
sum-of-squares (SOS) layer that reproduces it by a dual storage-function
program. It corresponds to the book chapter **Linear-Cost Optimal Control** of
the companion monograph. The section below states the problem in full so that
the chapter is self-contained; no access to the book is assumed.

## Mathematical setting

**The system.** The object of study is a linear Hermitian matrix-valued system.
The concrete case is a controlled diffusion whose state `x(t)` in `C^n` and
control `u(t)` in `C^m` evolve as

$$
\mathrm{d}x = (A x + B u)\,\mathrm{d}t
  + \sum_{l} (A_l x + B_l u)\,\mathrm{d}W_l ,
$$

with independent scalar Brownian increments `dW_l`; spontaneous (Poisson) jumps
`x^+ = J_r x + B_{J,r} u` firing at rate `lambda_r` may be added. The toolbox
carries `A`, `B` (drift and control input), the multiplicative-noise channels
`A_l`, `B_l` (`Ac`, `Bc`), and the jump data (`Jf`, `lamf`, `BJf`).

**Why matrices.** All second-order information is carried by the joint
second-moment matrix

$$
X(t) = \mathbb{E}\begin{bmatrix} x \\ u \end{bmatrix}
                    \begin{bmatrix} x \\ u \end{bmatrix}^{*}
     = \begin{bmatrix} X_1 & X_2 \\ X_2^{*} & X_3 \end{bmatrix}
     \succeq 0 ,
\qquad
X_1 = \mathbb{E}[x x^{*}],  X_3 = \mathbb{E}[u u^{*}] .
$$

`X(t)` is a Hermitian positive-semidefinite matrix trajectory, and every cost
and constraint of the theory is **linear** in `X`.

**The linear cost.** With the augmented cost weight

$$
Z = \begin{bmatrix} Q & S \\ S^{*} & R \end{bmatrix} \succeq 0,
\qquad Q = Q^{*} \succeq 0,  R = R^{*} \succ 0,
$$

the running cost is the inner product `<Z,X> = tr(ZX)`, which unfolds to the
familiar LQ integrand

$$
\langle Z, X\rangle = \operatorname{tr}(Q X_1) + 2 \Re\operatorname{tr}(S^{*} X_2)
   + \operatorname{tr}(R X_3)
   = \mathbb{E}\left[ x^{*}Q x + 2 \Re(x^{*}S u) + u^{*}R u \right].
$$

The **finite-horizon** cost over `[0,T]`, with terminal weight `Z_T` on the
state block `X_1(T)`, is

$$
J_T = \int_{0}^{T} \langle Z, X(t)\rangle\,\mathrm{d}t
      + \langle Z_T, X_1(T)\rangle .
$$

The **infinite-horizon** cost and the **long-term average** cost rate are

$$
J_\infty = \int_{0}^{\infty} \langle Z, X(t)\rangle\,\mathrm{d}t ,
\qquad
\rho = \lim_{T\to\infty} \frac{1}{T} J_T
     = \sum_{i} \pi_i^{\mathrm{ss}} \langle P_i, W_i\rangle ,
$$

where `W_i` is the additive process-noise intensity of mode `i` and
`pi^ss` is the stationary distribution of the switching Markov chain (for a
single mode the sum is one term `<P,W>`). The average-cost problem is the right
one under persistent excitation: `J_infty` is then infinite, but `rho` is
finite.

**The impulsive generalized Riccati equation (I-GRE).** Optimality is certified
by a dual variable `P = P^*` in `C^{n x n}`, using the Hermitian-part operator
`He(M) = M + M^*`, through the three flow blocks

$$
\begin{aligned}
(\mathcal{C}_Z)_1 &= \operatorname{He}(A^{*}P) + Q
   + \textstyle\sum_l A_l^{*} P A_l
   + \textstyle\sum_r \lambda_r (J_r^{*} P J_r - P) , \\
(\mathcal{C}_Z)_2 &= P B + S
   + \textstyle\sum_l A_l^{*} P B_l
   + \textstyle\sum_r \lambda_r J_r^{*} P B_{J,r} , \\
(\mathcal{C}_Z)_3 &= R
   + \textstyle\sum_l B_l^{*} P B_l
   + \textstyle\sum_r \lambda_r B_{J,r}^{*} P B_{J,r} ,
\end{aligned}
$$

In the nonsingular case `(C_Z)_3 > 0` the finite- and
infinite-horizon I-GRE are the backward differential and the algebraic Riccati
equations

$$
-\dot P = (\mathcal{C}_Z)_1 - (\mathcal{C}_Z)_2 (\mathcal{C}_Z)_3^{-1} (\mathcal{C}_Z)_2^{*},
\quad P(T) = Z_T
\qquad\text{and}\qquad
(\mathcal{C}_Z)_1 - (\mathcal{C}_Z)_2 (\mathcal{C}_Z)_3^{-1} (\mathcal{C}_Z)_2^{*} = 0 .
$$

The optimal state-feedback gain and the optimal cost are

$$
K^{\star} = -(\mathcal{C}_Z)_3^{-1} (\mathcal{C}_Z)_2^{*} ,
\qquad
J^{\star} = \langle P, X_0\rangle = \operatorname{tr}(P X_0) ,
$$

with `X_0 = X_1(0)` the initial second moment. The singular case
(`(C_Z)_3` rank-deficient, free or singular inputs) is handled by the
pseudoinverse form; the book develops it, and the SOS layer certifies it
without inverting anything. **Standing assumptions:** `(A,B)` stabilizable and
`(A, Q^{1/2})` detectable for the infinite-horizon problem, `Q >= 0`, `R > 0`.

**What the solver returns.** Each engine returns the dual `P` (the I-GRE
solution), the gain `K^star`, and the cost: `.cost = <P,X0>` for the finite and
infinite horizons, `.rho` for the average horizon.

## The three horizons

A system carries the quadratic-cost blocks `Q`, `R` (and optionally `S`, the
terminal weight `ZT`, the initial covariance `X0`, and the noise intensity `W`).
Three methods price the three horizons of the chapter:

- `S.ocfinite(T)`: the finite-horizon differential I-GRE on `[0,T]` integrated
  backward from `P(T) = Z_T`, via `mvst_solve_finite`;
- `S.ocinfinite()`: the infinite-horizon (algebraic, coupled) I-GRE, via
  `mvst_solve_infinite`;
- `S.ocaverage()`: the long-term average cost rate `rho`, via
  `mvst_solve_average`.

### Finite and infinite horizon

The finite-horizon cost is `<P(0), X0>` for the built initial covariance; here
`X0 = I`, so the reported cost is `trace(P(0))`. As the horizon grows it
approaches the infinite-horizon cost:

```matlab
S = mvst.System.ct([-1 0;0 -2], [1;1], 'Q', eye(2), 'R', 1, 'X0', eye(2));
fprintf('finite cost(T=3) = %.6f\n', S.ocfinite(3).cost);
fprintf('infinite cost = %.6f\n', S.ocinfinite().cost);
```
```
finite cost(T=3) = 0.660876
infinite cost = 0.661064
```

The finite-horizon value (0.660876) sits just below the infinite-horizon value
(0.661064), as expected for a converging horizon: the finite horizon accrues
strictly less cost than the infinite one, and closes the gap as `T` grows.

### Average cost

The long-term average cost rate `rho` weights the stationary running cost by the
noise intensity `W`:

```matlab
S = mvst.System.ct([-1 0;0 -2], [1;1], 'Q', eye(2), 'R', 1, 'W', eye(2));
fprintf('average rho = %.6f\n', S.ocaverage().rho);
```
```
average rho = 0.661064
```

`S.ocaverage` returns the infinite-horizon solution plus the extra field
`rho`. Result structs additionally carry the optimal gain `K` and the Riccati
solution `P` (both `1 x M` cells, one entry per mode).

## Worked case study: horizon convergence and the closed loop

This end-to-end example prices a single non-normal mode over a sweep of
horizons, confirms the monotone convergence to the infinite-horizon value,
extracts the optimal gain, verifies that it stabilizes the closed loop, and
cross-checks the whole solution against the exact algebraic Riccati engine
`mvst_care_ct`. Every number below is produced by the run.

```matlab
S = mvst.System.ct([-1 0.2; 0 -2], [1;0.5], ...
                   'Q', eye(2), 'R', 2, 'X0', eye(2), 'W', 0.3*eye(2));
for T = [0.5 1 2 5 10]
    fprintf('T=%5.1f  cost=%.6f\n', T, S.ocfinite(T).cost);
end
ri = S.ocinfinite();
fprintf('infinite cost = %.6f  (iters=%d, residual=%.2e)\n', ...
        ri.cost, ri.iters, ri.residual);
fprintf('average rho    = %.6f\n', S.ocaverage().rho);
K = ri.K{1};
fprintf('optimal gain K = %s\n', mat2str(K,5));
Acl = [-1 0.2;0 -2] + [1;0.5]*K;
fprintf('eig(A)   = %s\n', mat2str(sort(eig([-1 0.2;0 -2])),5));
fprintf('eig(Acl) = %s\n', mat2str(sort(eig(Acl)),5));
[Pr,Kr] = mvst_care_ct(mvst_mode([-1 0.2;0 -2],[1;0.5],eye(2),2,'X0',eye(2)));
fprintf('||K-Kcare|| = %.1e   ||P-Pcare|| = %.1e\n', ...
        norm(K-Kr,'fro'), norm(ri.P{1}-Pr,'fro'));
```
```
T=  0.5  cost=0.523734
T=  1.0  cost=0.650952
T=  2.0  cost=0.693441
T=  5.0  cost=0.697055
T= 10.0  cost=0.697057
infinite cost = 0.697057  (iters=2, residual=1.53e-11)
average rho    = 0.209117
optimal gain K = [-0.2286 -0.071784]
eig(A)   = [-2;-1]
eig(Acl) = [-2.0173;-1.2472]
||K-Kcare|| = 0.0e+00   ||P-Pcare|| = 0.0e+00
```

Interpretation of the run: the finite-horizon cost climbs monotonically from 0.523734 at
`T = 0.5` and settles at 0.697057 by `T = 10`, matching the infinite-horizon
value to six digits (the algebraic solve itself needs only two Gauss-Seidel
iterations to a residual of 1.5e-11). The average-cost rate `rho = 0.209117` is
the noise-weighted stationary running cost, a different quantity from the
transient cost. The optimal gain shifts both eigenvalues left (from `-1, -2` to
`-1.25, -2.02`), so the closed loop is stable. Finally the Riccati solution `P`
and the gain `K` recovered from the infinite-horizon engine agree with the
exact algebraic Riccati solver to machine precision, so the horizon result and
the closed-form result are the same object.

## SOS layer: `mvst_oc_sos`

`mvst_oc_sos` solves the same infinite-horizon problem by the dual
(storage-function) SOS / LMI program of the chapter. Instead of integrating the
Riccati equation it maximizes `<P,X0>` over all `P` for which the augmented KYP
matrix

$$
\mathcal{C}_Z(P) =
\begin{bmatrix} (\mathcal{C}_Z)_1 & (\mathcal{C}_Z)_2 \\
                (\mathcal{C}_Z)_2^{*} & (\mathcal{C}_Z)_3 \end{bmatrix} \succeq 0
$$

is a matrix sum of squares in the indeterminates `[x;u]`. This is the
pre-Schur (bounded-real) form of the I-GRE, and the maximizer is exactly the
stabilizing solution `P^star`, from which the gain is recovered as
`K^star = -(C_Z)_3^{-1}(C_Z)_2^*`. The routine returns a primal cost, a dual
(moment) cost, and the recovered `P` and `K`. The multiplicative-noise and
Poisson terms enter `P` linearly, so the program stays affine. The example
below builds a thin mode, solves the exact algebraic Riccati with
`mvst_care_ct`, then the SOS dual, and compares:

```matlab
md = mvst_mode([-1 0;0 -2], [1;1], eye(2), 1, 'X0', eye(2));
[Pr,Kr] = mvst_care_ct(md);
r1 = mvst_oc_sos(md);
fprintf('SOS cost primal=%.6f  dual=%.6f\n', r1.cost.primal, r1.cost.dual);
fprintf('||P_sos-P_ric||=%.2e  ||K_sos-K_ric||=%.2e\n', ...
        norm(r1.P-Pr,'fro'), norm(r1.K-Kr,'fro'));
```
```
SOS cost primal=0.661064  dual=0.661064
||P_sos-P_ric||=4.79e-10  ||K_sos-K_ric||=5.32e-10
```

The SOS primal and dual costs coincide (0.661064) and match the
infinite-horizon value above, and the recovered `P` and `K` agree with the
Riccati solution to about 1e-10. The gap between the two is the SOS relaxation
gap at this certificate degree, negligible here. The multi-mode / dwell-time
variant is `mvst_oc_sos_modes`, which handles the periodic, minimum, and range
dwell-time regimes by a timer-dependent certificate. The SOS layer requires
SOSTOOLS and an SDP backend; the Riccati horizons need only base MATLAB.

## Dwell-time optimal control

### `mvst_oc_dwell`: causal minimum and range dwell-time control

`res = mvst_oc_dwell(md, jump, opts)` computes the causal optimal control of a
single-cone linear impulsive system whose inter-jump times obey a minimum or a
range dwell-time constraint. `md` is a mode struct from `mvst_mode` (flow data
`A,B,Q,R`), `jump` carries the reset `J`, the jump input `Bd`, and the jump cost
`Qd`, `Rd`; `opts.dwell` is `'min'` or `'range'` with `opts.Tmin` (and
`opts.Tmax`), and `opts.X10` the initial second moment. The certificate is the
fixed point of the backward flow integration of the differential Riccati
equation composed with the jump Riccati map.

```matlab
md   = mvst_mode([-1 0.5;0 -0.8],[0;1],eye(2),1);
jump = struct('J',1.2*eye(2),'Bd',[0;1],'Qd',0.05*eye(2),'Rd',1);
r = mvst_oc_dwell(md,jump,struct('dwell','min','Tmin',1,'X10',eye(2)));
fprintf('Jbound=%.4f feas=%d Kd=%s\n',r.Jbound,r.feas,mat2str(r.Kd,4));
```
```
Jbound=1.0655 feas=1 Kd=[-0.0974 -0.4166]
```

`res.Jbound` is the certified cost bound realized by the recovered policy for
every admissible dwell-time sequence, `res.P0` the certificate at the start of a
dwell-time interval, `res.Kc` and `res.Kd` the flow and jump gains. As `Tmin` grows,
`res.P0` converges to the infinite-horizon solution of `mvst_care_ct`.
