# 8. Games

[Manual index](index.md) | Prev: [Dissipativity](07-dissipativity.md) | Next: [Dissipative synthesis](09-dissipative-synthesis.md)

This chapter covers the state-feedback zero-sum `H_infinity` game: the exact
critical level of `mvst_game_gamma` (method `S.gamegamma`) and the SOS
relaxation of `mvst_game_sos`. It corresponds to the book chapter **Differential
Games** of the companion monograph. The section below states the game in full so
the chapter is self-contained.

## Mathematical setting

**The two-player system.** The state `x` is driven by a **minimizing** control
`u` (mu channels) and a **maximizing** disturbance `w` (mw channels):

$$
\mathrm{d}x = (A x + B_u u + B_w w)\,\mathrm{d}t
   + \sum_l (A_l x + B_l u)\,\mathrm{d}W_l .
$$

The toolbox carries `A` (`A0`), the control input `Bu`, the disturbance input
`Bw`, the state weight `Q`, and the control weight `Ru`.

**The minimax payoff.** With `Q >= 0`, `Ru > 0`, and a disturbance penalty
`gamma^2`, the players share the quadratic cost, linear in the joint second
moment `X`,

$$
J(u,w) = \int_0^\infty \Big(\langle Q, X_1\rangle + \langle R_u, X_{22}\rangle
    - \gamma^2 \langle I, X_{33}\rangle\Big)\,\mathrm{d}t
   = \mathbb{E}\int_0^\infty \big(\lVert z\rVert^2 - \gamma^2 \lVert w\rVert^2\big)\,\mathrm{d}t ,
$$

where `z = [Q^{1/2} x;  R_u^{1/2} u]` is the performance output and
`X_{22}, X_{33}` are the control and disturbance second-moment blocks. The
minimizer picks `u` to reduce `J`, the maximizer picks `w` to increase it. The
negative sign on the `w`-block encodes the opposing objective. The game value is
the saddle

$$
J^{\circ} = \min_{u} \max_{w}  J(u,w) .
$$

**Isaacs condition.** The value is well defined (min-max equals max-min, a
saddle point exists) when **Isaacs' condition** holds: the bottom-right Isaacs
block of the augmented weight,

$$
\mathcal{G} = \begin{bmatrix} R_u & 0 \\ 0 & -\gamma^2 I \end{bmatrix}
\qquad (\mathcal{G}_{uu} = R_u \succ 0,  \mathcal{G}_{ww} = -\gamma^2 I \prec 0),
$$

has the saddle inertia that lets the order of min and max be exchanged.
`G` is **indefinite** (positive on the control block, negative on the
disturbance block); this indefiniteness is what separates the game from ordinary
optimal control.

**The game generalized Riccati equation (game GRE).** The saddle is
characterized by a stabilizing Hermitian solution `P >= 0` of the indefinite
(Hamilton-Jacobi-Isaacs) Riccati equation

$$
\operatorname{He}(A^{*}P) + Q
   - P B_u R_u^{-1} B_u^{*} P
   + \tfrac{1}{\gamma^2}  P B_w B_w^{*} P = 0 ,
$$

whose defining feature is the `+ (1/gamma^2)` disturbance term (opposite in sign
to the control term). The saddle policies are the two central gains

$$
K_u = -R_u^{-1} B_u^{*} P \quad (\text{optimal control}),
\qquad
K_w = \tfrac{1}{\gamma^2} B_w^{*} P \quad (\text{worst-case disturbance}).
$$

**The bounded-real interpretation of `H_infinity`.** As `gamma` decreases toward a
critical level `gamma*`, the `+ (1/gamma^2)` term grows and the stabilizing
solution of the game GRE ceases to exist. That critical level is exactly the
closed-loop bounded-real (`H_infinity`) level: the smallest disturbance
attenuation the controller can guarantee. For `gamma > gamma*` the game GRE has
a stabilizing `P >= 0` (the controller wins); for `gamma < gamma*` it does not.
The toolbox finds `gamma*` by **bisecting** on `gamma`, at each trial testing
whether the game GRE (or its SOS certificate) is solvable. **Standing
assumptions:** `(A,B_u)` stabilizable, `Q >= 0`, `R_u > 0`.

## The bounded-real interpretation, in the toolbox

The critical level `gamma*` is what `S.gamegamma` returns. A feasibility solve
at a fixed `gamma` (`mode = 'feas'`) certifies the bounded-real property at that
level: feasible above `gamma*`, infeasible below it. The worked example below
shows exactly this crossover.

## Exact level: `S.gamegamma`

`S.gamegamma` computes `gamma*` by gamma-bisection over the game Riccati of the
resolved class. Supply the state weight `Q`, control weight `R`, control input
`B`, and disturbance input `Bw` through the factory:

```matlab
S = mvst.System.ct([0 1;-1 -0.5], [0;1], 'Q', eye(2), 'R', 1, 'Bw', [1;0]);
r = S.gamegamma();
fprintf('gamma* = %.6f  method=%s\n', r.gammastar, r.method);
```
```
gamma* = 1.293481  method=ct
```

`res.gammastar` is the critical level, `res.P` the settled dual(s), and
`res.method` the resolved class. The backends are `mvst_game_gamma_ct`,
`_ct_mjls`, `_ct_dwell`, and `_ct_poisson`; the sampled-data game shares the
`ct` backend (use `OPTS.class = 'ct'`).

## SOS relaxation: `mvst_game_sos`

`mvst_game_sos(sys, game, dwell, mode)` returns the critical level and central
gains over every dwell-time type of the chapter.

**Caveat.** The book's game Riccati inequality is not jointly convex in
the dual variable (the Isaacs block `G` is indefinite, so the generalized Schur
complement is rational in `P` and paired with the saddle-inertia condition), so
it cannot be handed to an SDP solver as stated. This routine instead solves the
sufficient convex relaxation from the classical bounded-real change of
variables. For the separated continuous-time LTI game (zero cross term) the
relaxation is **exact** and reproduces the true `gamma*`; for the
timer-dependent / dwell-time certificates it returns a certified **upper bound**
that tightens as the certificate degree grows. The exact per-class solvers above
always return the exact `gamma*`.

The example below solves the same system by the arbitrary-dwell-time SOS program and
compares to the exact level:

```matlab
Ag=[0 1;-1 -0.5]; Bu=[0;1]; Bw=[1;0]; Q=eye(2); Ru=1;
rg = mvst_game_sos(struct('A0',Ag,'Bu',Bu,'Bw',Bw), ...
                   struct('Q',Q,'Ru',Ru), struct('type','arbitrary'), 'gamma');
ex = mvst_game_gamma(struct('A',Ag,'Bu',Bu,'Bw',Bw,'Q',Q,'Ru',Ru));
fprintf('SOS gamma=%.6f  exact gamma*=%.6f\n', rg.gamma, ex.gammastar);
```
```
SOS gamma=1.293300  exact gamma*=1.293481
```

The SOS value (1.293300) matches the exact level (1.293481) to about 2e-4, as
expected for this separated game. The `dwell` struct extends to `'discrete'`,
`'periodic'`, `'minimum'`, and `'range'`, and the system may carry
state-multiplicative and Poisson jump data. This routine requires SOSTOOLS and
an SDP backend; `S.gamegamma` needs only base MATLAB.

## Worked case study: the critical level and the bounded-real crossover

This end-to-end example takes an undamped-plus-feedback oscillator, computes the
critical level exactly, extracts both saddle gains from the settled game GRE
solution, and then demonstrates the bounded-real interpretation directly: a feasibility
solve **below** `gamma*` fails and one **above** `gamma*` succeeds.

```matlab
A = [0 1; -2 -1]; Bu = [0;1]; Bw = [1;0]; Q = eye(2); Ru = 1;
S  = mvst.System.ct(A, Bu, 'Q', Q, 'R', Ru, 'Bw', Bw);
ex = S.gamegamma();
sos = mvst_game_sos(struct('A0',A,'Bu',Bu,'Bw',Bw), ...
                    struct('Q',Q,'Ru',Ru), struct('type','arbitrary'), 'gamma');
P  = ex.P;
Ku = -Ru\(Bu'*P);           % optimal control gain
Kw = (1/ex.gammastar^2)*(Bw'*P);  % worst-case disturbance gain
fprintf('exact gamma* = %.6f   SOS gamma = %.6f   |diff| = %.2e\n', ...
        ex.gammastar, sos.gamma, abs(ex.gammastar-sos.gamma));
fprintf('Ku = %s\nKw = %s\n', mat2str(Ku,5), mat2str(Kw,5));
for g = [0.9 1.1]*ex.gammastar
    rf = mvst_game_sos(struct('A0',A,'Bu',Bu,'Bw',Bw), ...
             struct('Q',Q,'Ru',Ru,'gamma',g), struct('type','arbitrary'), 'feas');
    fprintf('  gamma=%.4f (=%.2f*gamma*)  feas=%d\n', g, g/ex.gammastar, rf.feas);
end
```
```
exact gamma* = 1.331323   SOS gamma = 1.336055   |diff| = 4.73e-03
Ku = [-2.7763 -2.4252]
Kw = [3.2336 1.5664]
  gamma=1.1982 (=0.90*gamma*)  feas=0
  gamma=1.4645 (=1.10*gamma*)  feas=1
```

Interpretation of the run: the exact critical level is `gamma* = 1.331323`, and the SOS
relaxation gives 1.336055, a certified upper bound within 5e-3. The two central
gains follow from the settled game GRE solution `P`: `Ku` is the minimizing
controller, `Kw` the direction the disturbance would push. The last two lines
are the bounded-real interpretation made concrete: at `gamma = 0.90 gamma*` (below the
critical level) the game is infeasible (`feas=0`, no controller can attenuate
that hard), while at `gamma = 1.10 gamma*` (above it) the game is feasible
(`feas=1`). This is exactly what the bisection in `S.gamegamma` certifies: it
squeezes `gamma` down to the boundary between these two regimes.

## Further game engines

These engines are called directly as `mvst_*` functions; the object method
`S.gamegamma` covers the critical level only.

### `mvst_game_average`: ergodic zero-sum value

`res = mvst_game_average(sys, opts)` returns the long-run average value `rho` of
the zero-sum game under a persistent forcing `sys.W` (a constant covariance, or
a handle `@(t)` with period `sys.Tp`), at a level `sys.gamma` above the critical
level. The homogeneous part is the H-infinity game dual; the value is linear in
the forcing.

```matlab
sys = struct('A',[-1 0.5;0 -2],'Bu',[1;0.5],'Bw',[0.3;0.6],'Q',eye(2),'Ru',1);
g = mvst_game_gamma(sys); sys.gamma = 1.3*g.gammastar; sys.W = eye(2);
r = mvst_game_average(sys);
fprintf('gammastar=%.4f rho=%.4f residual=%.1e\n',g.gammastar,r.rho,r.residual);
```
```
gammastar=0.3840 rho=0.7495 residual=3.7e-10
```

`res.Pinf` is the homogeneous game dual and `res.Ku`, `res.Kw` the saddle gains.

### `mvst_game_mixed`: mixed H2 / H-infinity design

`res = mvst_game_mixed(sys, opts)` solves the mixed design as a two-player
nonzero-sum Nash game: the designer carries the H2 control cost, the fictitious
disturbance carries the bounded-real cost at level `sys.gamma` (required). The
two coupled Riccati legs are solved by a Gauss-Seidel iteration. As `gamma`
grows, the solution tends to the pure H2 solution and `Kw` to zero.

```matlab
sys = struct('A',[-1 0.5;0 -2],'Bu',[1;0.5],'Bw',[0.3;0.6],'Q',eye(2),'Ru',1,'Su',zeros(2,1));
sys.gamma = 2*0.3840;                       % twice the critical level
r = mvst_game_mixed(sys);
fprintf('iters=%d residual=%.1e h2cost=%.4f\n',r.iters,r.residual,r.h2cost);
```
```
iters=8 residual=1.1e-10 h2cost=0.7082
```

The gains are `res.Ku = [-0.4606 -0.2072]` and `res.Kw = [0.1978 0.2655]`. Below
the critical level the iteration does not converge.

### `mvst_game_nash`: N-player Nash equilibrium

`res = mvst_game_nash(sys, opts)` solves the coupled generalized Riccati system
of an N-player game by a Gauss-Seidel sweep. `sys.Bu`, `sys.Q` are 1-by-N cells
and `sys.R{i}{j}` is player `i`'s weight on control `j`. A two-player zero-sum
instance reproduces the H-infinity saddle.

```matlab
sys = struct('A',[-1 0.5;0 -2],'Bu',{{[1;0.5],[0.3;0.6]}},'Q',{{eye(2),-eye(2)}},'R',{{{1,-4},{-1,4}}});
r = mvst_game_nash(sys);
fprintf('iters=%d residual=%.1e cost=%s\n',r.iters,r.residual,mat2str(r.cost,4));
```
```
iters=3 residual=6.9e-11 cost=[0.665 -0.665]
```

`res.K{i}` and `res.P{i}` are the equilibrium gains and duals; the zero-sum
costs are opposite, as expected.

### `mvst_game_sd_held` and `mvst_game_gamma_sd_held`: held-input sampled-data game

`out = mvst_game_sd_held(A,Bu,Bw,Q,Ru,gam,Tp,Dwc,Nmax,tol)` solves by backward
iteration the periodic HJI Riccati equation of the zero-sum game on a
sampled-data plant whose control is held between samples of period `Tp`
(`Dwc` is the optional disturbance-multiplicative noise matrix, `[]` for none).
`[gstar,lo] = mvst_game_gamma_sd_held(A,Bu,Bw,Q,Ru,Tp,Dwc)` bisects the
critical level, returning the feasible upper end and the infeasible lower end of
the final bracket.

```matlab
A=[-1 0.5;0 -2]; Bu=[1;0.5]; Bw=[0.3;0.6]; Q=eye(2); Ru=1; Tp=0.2;
[gstar,lo] = mvst_game_gamma_sd_held(A,Bu,Bw,Q,Ru,Tp,[]);
o = mvst_game_sd_held(A,Bu,Bw,Q,Ru,1.2*gstar,Tp,[]);
fprintf('gstar=%.4f ok=%d reason=%s Ku=%s\n',gstar,o.ok,o.reason,mat2str(o.Ku,4));
```
```
gstar=0.3839 ok=1 reason=converged Ku=[-0.4623 -0.2479]
```

`out.ok` is true when the iteration converged with the Isaacs conditions in
force; `out.reason` otherwise names the failure (`'escape'`, `'isaacs_flow'`,
`'isaacs_jump'`, `'blowup'`, `'maxiter'`). `out.Ku` is the sample-time gain.
