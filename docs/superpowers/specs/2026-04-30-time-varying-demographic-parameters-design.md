# Time-Varying Demographic Parameters in discoal

**Status:** Draft for review
**Issue:** [#82 — Demes importer emits unfaithful migration events](https://github.com/kr-colab/discoal/issues/82)
**Branch base:** `origin/nsp-yaml-revamp`
**Working branch:** `feature/issue-82-time-varying-demography`
**Date:** 2026-04-30

## 1. Motivation

### 1.1 The proximate trigger: issue #82

The demes importer in `src/core/demesInterface.c:402-443` translates each
`migrations:` window into a pair of `'M'` events (start with rate $r$, end
with rate 0). These events are not handled by the runtime event-dispatch
switch in `src/core/discoal_multipop.c:228`. A back-derivation hack at
`src/core/discoalFunctions.c:222-261` scans the entire event list to
reconstruct what `migMat` should be at $t = 0$, then sets it. That hack
is the only path by which demes-emitted migration events influence the
simulation.

The hack works for single-window migration. It produces incorrect
results whenever the demes graph has multi-window migration (e.g., a
deme born partway through the simulation): the rate-0 events still fire
at their original times and zero out migration the demes spec says
should be active. The detailed reproducer is in issue #82.

### 1.2 The deeper problem

The current discoal engine assumes **rates are constant between
events**. It supports piecewise-constant population sizes via `'n'`
events (mutating `currentSize[popID]`) but has no mechanism to vary
migration rates over time and no mechanism for continuously-varying
sizes (exponential or linear growth). The `population_structure.rst`
docs explicitly note: "Time-varying migration rates are not currently
implemented." The demes importer rejects exponential and linear epochs
outright (`demesInterface.c:323-345`).

This blocks msprime parity: any demes graph with growth or
multi-window migration cannot be simulated faithfully by discoal.

### 1.3 Goal

Add first-class support for time-varying demographic parameters with
the shape vocabulary that demes uses (constant, exponential, linear),
applied to both population sizes and pairwise migration rates. The
result should:

- Faithfully simulate any demes graph that demes-c can parse.
- Achieve statistical parity with msprime for neutral simulations
  under demes models.
- Remain backward-compatible with all existing CLI flags.
- Preserve sweep-simulation correctness for piecewise-constant
  demography (statistical parity with current discoal).
- Remove the back-derivation hack as a side effect.

## 2. Goals and Non-Goals

### In scope

- Constant, exponential, and linear shape primitives for population
  sizes and migration rates.
- Continuous integration of these shapes in both the neutral phase
  and the sweep phase.
- New event types `'g'` (size-shape change) and `'em'` (migration-shape
  change) in the runtime.
- New CLI flags `-eg`, `-eG`, `-em`, `-eM` mirroring ms.
- Importer changes: emit faithful shape events; remove rejection of
  exponential/linear epochs.
- Removal of the back-derivation hack at `discoalFunctions.c:222-261`.
- Test-driven development for every component.
- Statistical parity tests against existing discoal (sweeps under
  piecewise-constant demography) and msprime (neutral under demes).

### Out of scope

- Arbitrary user-supplied $\lambda(t)$ shapes (thinning sampler).
- Time-varying recombination rate, mutation rate, gene-conversion rate,
  selection coefficient.
- Sweep behavior under time-varying selection.
- Migration during sweep phases (currently suspended; remains suspended).
- A linear-shape CLI primitive (linear shapes are accepted via demes/YAML
  only at first cut; no precedent in ms vocabulary).
- Changes to existing semantics of `-en`, `-em` (in the sense of single
  pairwise migration set), `-M`, `-m`.

## 3. Mathematical Foundation

### 3.1 Non-homogeneous Poisson process draws

The neutral phase is a competing-exponentials sampler today. With each
rate component $\lambda_k$ constant, the next inter-event time is
$T \sim \text{Exp}(\sum_k \lambda_k)$. With time-varying rates each
component becomes its own non-homogeneous Poisson process. To draw a
waiting time $T_k$ for component $k$, sample $\xi_k = -\log U_k$ with
$U_k \sim \text{Uniform}(0, 1)$ and solve

$$\xi_k = \int_0^{T_k} \lambda_k(s)\,ds.$$

Take $T^* = \min_k T_k$ as the next event time, capped at the next
shape-change boundary. By memorylessness, components that did not fire
get fresh draws at the next iteration; or equivalently, all $T_k$ are
redrawn after each event because the lineage state may have changed.

### 3.2 Closed-form integrals for the three shapes

#### Sign convention

discoal follows the **msprime forward-time per-generation rate
convention** for `Shape.rate_param` across all shape types. Each
population has an initial absolute size $N_e$ and a per-generation
forward-time rate of change ($\alpha$ for exponential, $\gamma$ for
linear). The size at backward-time $s$ in the simulator's natural
direction is $N(s) = N_e\,e^{-\alpha s}$ for exponential and
$N(s) = N_e - \gamma s$ for linear. Migration rates use the same
convention. A *positive* rate parameter means forward-time growth,
which means the past held a *smaller* value.

#### Setup

Let $\binom{k}{2}$ denote the within-population pair count and $k_i$
the lineage count in population $i$. Coalescent rate within
population $i$ at backward-time $s$ is

$$\lambda_C^{(i)}(s) = \binom{k_i}{2} \big/ N_i(s),$$

migration rate from $i$ to $j$ for the lineages in $i$ is

$$\lambda_M^{(i \to j)}(s) = k_i \cdot m_{ij}(s),$$

and recombination/gene-conversion are constant per lineage. The three
shapes give:

**Constant.** $N(s) = N_0$ or $m(s) = m_0$. $T = \xi / \lambda$.
This is the existing case.

**Exponential.** $N(s) = N_0 e^{-\alpha s}$ where $\alpha$ is the
forward-time growth rate. $\alpha > 0$ ⇒ backward decline (past was
smaller).

For the coalescent rate, $\lambda_C(s) = \binom{k}{2} e^{\alpha s}/N_0$
and

$$\int_0^T \lambda_C(s)\,ds = \frac{\binom{k}{2}}{N_0\alpha}\big(e^{\alpha T} - 1\big) = \xi$$

inverts to

$$T = \frac{1}{\alpha}\,\log\!\left(1 + \frac{N_0\alpha}{\binom{k}{2}}\,\xi\right).$$

For exponential migration $m(s) = m_0 e^{-\beta s}$ where $\beta$ is
the forward-time growth rate of migration, $\lambda_M(s) = k\,m_0\,
e^{-\beta s}$ and

$$\int_0^T \lambda_M(s)\,ds = \frac{k\,m_0}{\beta}\,(1 - e^{-\beta T}) = \xi$$

inverts to

$$T = -\frac{1}{\beta}\,\log\!\left(1 - \frac{\beta\,\xi}{k\,m_0}\right).$$

If $\beta\xi/(k\,m_0) \ge 1$ the integrated hazard over $[0,\infty)$
is finite and $\xi$ exceeds it — the draw is unreachable; return
$-1$ at the call site.

**Linear.** $N(s) = N_0 - \gamma s$ where $\gamma$ is the forward-time
growth rate. $\gamma > 0$ ⇒ backward decline (past was smaller); 
$\gamma < 0$ ⇒ backward growth.

$$\int_0^T \frac{\binom{k}{2}}{N_0 - \gamma s}\,ds = \frac{\binom{k}{2}}{\gamma}\,\log\!\left(\frac{N_0}{N_0 - \gamma T}\right) = \xi$$

inverts to

$$T = \frac{N_0}{\gamma}\left(1 - \exp\!\left(-\frac{\gamma\xi}{\binom{k}{2}}\right)\right).$$

The $\gamma = 0$ degenerate case reduces to the constant formula
$T = \xi N_0/\binom{k}{2}$.

(Care: if $\gamma > 0$ — the population was growing forward in time
— $N(s) = N_0 - \gamma s$ goes non-positive at the zero-crossing
$T^* = N_0/\gamma$. The integrated hazard from 0 to $T^*$ diverges
to $+\infty$, so any finite $\xi$ produces $T < T^*$ from the closed
form — but the inversion must clamp at $T^*$ via $\texttt{nextafter}$
to guard against floating-point ties. If $\gamma < 0$ — past was
larger — $N$ never reaches zero and the closed form is well-behaved
for all $T > 0$.)

The migration linear case has integrand $k\,(m_0 - \delta s)$ where
$\delta$ is the forward-time rate of migration change.
$\int_0^T k\,(m_0 - \delta s)\,ds = k\,m_0\,T - \tfrac{1}{2}\,k\,\delta\,T^2 = \xi$
is a quadratic in $T$ with positive root

$$T = \frac{m_0 - \sqrt{m_0^2 - 2\,\delta\,\xi/k}}{\delta}.$$

If the discriminant is negative the draw is unreachable; return
$-1$ at the call site. (When $\delta > 0$ the migration rate
declines backward to zero at $T^* = m_0/\delta$ and the integrated
hazard caps at $k\,m_0\,T^*/2$.) The $\delta = 0$ degenerate case
reduces to $T = \xi/(k\,m_0)$.

### 3.3 Sweep-phase math under continuous $N(t)$

The sweep phase is already an Euler-style small-dt forward integration
(`proposeTrajectory` in `discoalFunctions.c:1764-1875`,
`sweepPhaseEventsConditionalTrajectory` in 2143+, and
`sweepPhaseEventsGeneralPopNumber` in 1881+). At every grid step it
re-evaluates $N$, $\alpha_{\text{eff}} = \alpha \cdot \text{sizeRatio}$,
and $\text{tInc} = 1/(\text{deltaTMod} \cdot N)$. Today these are
pulled from a scalar `currentSizeRatio` mutated by `'n'` events; under
B2 they will be pulled from `sizeAt(popID, t)` evaluating the current
shape.

Three sweep modes:

- **Stochastic forward (`'s'`)**: WF SDE Euler step. Algorithm unchanged;
  just reads $N$ from `sizeAt` per step.
- **Neutral stochastic (`'N'`)**: drift-diffusion Euler step. Same.
- **Deterministic (`'d'`)**: today uses `detSweepFreq(\tau, \alpha_{\text{eff}})`,
  the Stephan et al. 1992 closed-form logistic. This formula assumes
  *constant* $\alpha_{\text{eff}}$. Under continuous $N(t)$ the
  closed form **still applies** — the deterministic ODE
  $dx/d\tau = -\alpha_{\text{eff}}(\tau)\,x(1-x)$ is separable, so

  $$x(\tau) = \frac{x_0\,e^{-A(\tau)}}{1 - x_0 + x_0\,e^{-A(\tau)}}$$

  where $A(\tau) := \int_0^{\tau}\alpha_{\text{eff}}(s)\,ds$. For
  constant $\alpha_{\text{eff}}$, $A(\tau)=\alpha\tau$ recovers the
  Stephan form bit-exactly. For our shape vocabulary $A(\tau)$ has
  closed forms (just $\alpha_0$ times the integral of `sizeAt`):

  - CONSTANT: $A = \alpha_0\,N_0\,\tau$ (with $N_0 = N_{\text{ref}}$)
  - EXPONENTIAL ($N(s)=N_0 e^{-\alpha_g s}$): $A = (\alpha_0 N_0/\alpha_g)(1-e^{-\alpha_g\tau})$, $\alpha_g\to 0$ degenerate
  - LINEAR ($N(s)=N_0 - \gamma_g s$): $A = \alpha_0(N_0\tau - \tfrac{1}{2}\gamma_g\tau^2)$

  Phase 5 introduces `detSweepFreqGeneral(x_0, A_\tau, A_{\tau_s})` and
  a sister `integratedSizeRatio(popID, t_0, T)` that returns
  $\int_{t_0}^{t_0+T}\text{sizeAt}(s)\,ds$. The dispatch on shape type
  for the deterministic path collapses to a single closed-form call;
  no Euler step, no truncation error.

Migration is suspended during sweep phases (existing behavior). Remains
suspended; time-varying migration is silently ignored across sweep
durations. This is documented and out of scope to change.

The acceptance probability returned by `proposeTrajectory` is
`currentSizeRatio / Nmax` with `Nmax = max(sizeRatio over walk)`. For
continuous shapes, `Nmax` is the running max as the trajectory walk
progresses. For monotone shapes within an epoch this is just the
endpoint of larger value; across multi-epoch walks it is the running
max as today. Generalization is straightforward.

## 4. Architecture

### 4.1 Shape state

Add per-population and per-pair shape state holding the *currently
active* shape for that pop or pair. Updated by `'n'`, `'g'`, `'em'`
events when they fire.

```c
typedef enum { SHAPE_CONSTANT, SHAPE_EXPONENTIAL, SHAPE_LINEAR } ShapeType;

typedef struct {
    ShapeType type;
    double anchor_value;   // size or rate at anchor_time
    double rate_param;     // alpha for EXP, gamma for LIN, unused for CONST
    double anchor_time;    // time at which anchor_value applies
} Shape;

extern Shape popShape[MAXPOPS];
extern Shape migShape[MAXPOPS][MAXPOPS];
```

Both arrays initialize from the parser/importer to the t=0 shape (see
section 4.6 on initialization).

### 4.2 Accessors

```c
double sizeAt(int popID, double t);
double migAt(int srcPopID, int dstPopID, double t);
double integratedHazardSize(int popID, double t0, double T, int k);
double integratedHazardMig(int srcPopID, int dstPopID, double t0, double T, int k);
double drawWaitingTimeSize(int popID, double t0, double xi, int k);
double drawWaitingTimeMig(int srcPopID, int dstPopID, double t0, double xi, int k);
```

Each is a small switch on `Shape.type` calling the appropriate
closed-form. `sizeAt` and `migAt` are the only accessors the sweep
phase needs. The neutral phase additionally calls
`drawWaitingTime{Size,Mig}` for NHPP draws.

`integratedHazard*` is exported for tests (verification of the closed
forms against numerical quadrature) but the sampler uses
`drawWaitingTime*` which inverts directly.

### 4.3 Event vocabulary

Today's runtime events: `'n'`, `'s'`, `'p'`, `'a'`, `'A'`. Adding two:

- `'g'`: at time $t$, set `popShape[popID]` to a new shape.
  - Carries: `popID`, new `Shape` (type, anchor_value, rate_param)
  - `anchor_time` of the new shape is set to $t$.
  - Convention 1: when `'g'` fires, the previous shape's value at $t$
    becomes the new shape's `anchor_value` if the importer doesn't
    override it (i.e., shape changes are continuous unless an `'n'`
    event explicitly snaps the size).
- `'em'`: at time $t$, set `migShape[i][j]` to a new shape.
  - Carries: `popID` (dst), `popID2` (src), new `Shape`.
  - Same continuity convention.

The existing `'n'` event remains: it sets `popShape[popID]` to
`(CONSTANT, value, 0, t)`, anchored at $t$ — the existing semantic.

Existing CLI flag `-en` continues to emit `'n'` events. New flags emit
`'g'`, `'em'`.

### 4.4 Event struct extension

```c
typedef struct event {
    double time, popnSize;    // popnSize used by 'n' for the snap-to value
    char type;
    int popID, popID2, popID3;
    int lineageNumber;
    double admixProp;
    Shape newShape;           // used by 'g' and 'em' only
} event;
```

`Shape` is fixed-size (3 doubles + an int), 32 bytes. Per-event memory
goes up; harmless.

### 4.4.1 Relation to nspope's suggestions in issue #82

The issue lists two proposed fix directions: (1) importer-side
synthesis of the active migration matrix per interval, removing the
runtime back-derivation; (2) adding a `-em` CLI flag for time-varying
migration. This design implements both. The importer produces faithful
shape events (1), and the CLI exposes the same primitive directly (2).
The runtime back-derivation is removed unconditionally.

### 4.5 Inner-loop changes

`neutralPhaseGeneralPopNumber` (`discoalFunctions.c:1589+`):

```c
// today: total constant rate, draw Exp(total)
// new: NHPP competing-events draw

double T_min = LOCAL_NEXT_TIME - currentTime;  // upper bound = next event
int winnerKind = NONE;  // {COAL_i, MIG_ij, RECOMB, GC, ...}
int winnerArg = -1;

for each population i {
    if (popnSizes[i] >= 2) {
        double xi = -log(ranf());
        double T = drawWaitingTimeSize(i, currentTime, xi, popnSizes[i]);
        // T may be infinity if xi > total integrated hazard over [t, infty)
        if (T < T_min) { T_min = T; winnerKind = COAL_i; winnerArg = i; }
    }
    for each j != i {
        if (popnSizes[i] >= 1 && migAt(i, j, currentTime) > 0) {
            double xi = -log(ranf());
            double T = drawWaitingTimeMig(i, j, currentTime, xi, popnSizes[i]);
            if (T < T_min) { T_min = T; winnerKind = MIG_ij; winnerArg = pack(i,j); }
        }
    }
}

// recomb, gene conv: rate constant in time per lineage; draw Exp directly
// rRate = rho * sum_i popnSizes[i] / 2
// double T_recomb = -log(ranf()) / rRate;   if (T_recomb < T_min) ...

if (winnerKind == NONE) {
    // no event before the next epoch boundary; advance to it
    currentTime = LOCAL_NEXT_TIME;
} else {
    currentTime += T_min;
    fire(winnerKind, winnerArg, currentTime);
}
```

Recombination and gene conversion remain constant-in-time per lineage,
so they stay homogeneous Poisson and are drawn with a single
$\text{Exp}(\rho/2 \cdot \sum_i k_i)$ call.

`recurrentSweepPhaseGeneralPopNumber` and the conditional-trajectory
sweep functions are *not* changed at the sampler level — they remain
small-dt Euler — but their inner per-step rate computations switch to
read from `sizeAt(i, t)` instead of `sizeRatio[i]` and from `migAt(i, j, t)`
where applicable (which is currently nowhere, since migration is
suspended during sweeps).

`proposeTrajectory` likewise reads `sizeAt(0, t)` per step (and
`sizeAt(i, t)` if it tracks other pops; it does not in current code).

Deterministic-mode sweep: replace the `detSweepFreq` call with the
shape-aware closed form using the integrated selection coefficient
$A(\tau) = \alpha_0 \cdot S(0,\tau)$:

```c
case 'd':
    /* x_0 was set at the start of the sweep walk to a value just below 1
     * (e.g., 1 - 1/(2*N(0))). A_now = alpha_0 * integratedSizeRatio(0, 0, ttau)
     * accumulates incrementally from A_prev and the per-step closed-form
     * size integral. */
    A_now = A_prev + alpha * integratedSizeRatio(0, currentTime + ttau - tIncOrig, tIncOrig);
    x = detSweepFreqGeneral(x_0, A_now);
    A_prev = A_now;
    break;
```

For SHAPE_CONSTANT, `integratedSizeRatio` returns `anchor_value * dt` and
$A(\tau) = \alpha_0 \cdot N_0 \cdot \tau$, recovering `detSweepFreq`
bit-exactly when expressed in the boundary-condition form. For
SHAPE_EXPONENTIAL and SHAPE_LINEAR, the per-step `integratedSizeRatio`
is the closed-form integral above. No Euler step appears; no truncation
error. The Q1 verification (section 6.3) ruled out plain Euler, but the
separable-ODE closed form sidesteps the issue entirely.

### 4.6 Initialization

The simulation starts at $t = 0$ with `migMatConst` and `currentSize[]`
populated by command-line flags or the importer. Under B2:

- `popShape[i]` initializes to `(CONSTANT, currentSize[i], 0, 0)` if no
  growth shape applies at $t=0$, else to the active shape with anchor
  values at $t=0$ for that population.
- `migShape[i][j]` initializes to `(CONSTANT, migMatConst[i][j], 0, 0)`
  if no migration window applies at $t=0$, else to the active shape.

The back-derivation hack at `discoalFunctions.c:222-261` is **deleted**.
The importer is responsible for writing the t=0 active matrix into
`migMatConst` directly (just as it does for sizes today via
`currentSize[]`), and emitting `'em'` events only for *transitions*.

### 4.7 Importer changes (`demesInterface.c`)

Replace the current paired-`'M'` emission with shape-aware events:

For each demes deme:
- For each epoch: if `size_function == EXPONENTIAL`, compute
  $\alpha = -\log(\text{end\_size}/\text{start\_size}) / (\text{end\_time} - \text{start\_time})$
  in coalescent time units, and emit a `'g'` event at the epoch's
  more-recent boundary anchoring the new shape.
- For `LINEAR`, similarly compute $\gamma$ and emit `'g'` with shape
  type `LINEAR`.
- For `CONSTANT` epochs across boundaries, emit `'n'` as today.
- Remove the rejection at `demesInterface.c:323-345`.

For each demes migration:
- Compute the windows-active migration matrix for each piecewise
  interval implied by the union of all migration windows. Within each
  interval, every pair has a constant rate (that is what demes
  guarantees: rates are constant within their own window, and windows
  do not overlap on the same pair).
- Emit one `'em'` event per pair per interval boundary, anchoring the
  per-pair shape.
- The t=0 active matrix (interval ending at $t=0$) gets written to
  `migMatConst[i][j]` directly, not as an event.

### 4.8 CLI surface

New flags (mirroring ms):

- `-eg <time> <popID> <alpha>`: emit `'g'` for population `popID` at
  `time` with `(EXPONENTIAL, sizeAt(popID, time), alpha, time)`. Anchor
  value is computed by evaluating the prior shape at `time` (Convention
  1 continuity).
- `-eG <time> <alpha>`: same as `-eg` for every population.
- `-em <time> <i> <j> <rate>`: emit `'em'` for pair (i,j) at `time`
  with `(CONSTANT, rate, 0, time)`. (Constant within the new window;
  to make a window with start and end, supply two `-em` events.)
- `-eM <time> <rate>`: same as `-em` for all off-diagonal pairs.

No CLI surface for linear shapes at first cut. Linear is reachable via
demes/YAML.

Existing flags unchanged: `-en`, `-m`, `-M`, `-ed`, `-ej`, `-ea`, `-A`,
`-w`, `-l`, all sweep flags.

### 4.9 What gets removed

- The back-derivation block at `discoalFunctions.c:222-261` (the
  fprintf-laced scan that infers the $t=0$ matrix). Replaced by direct
  initialization from `migMatConst`.
- The importer rejection at `demesInterface.c:323-345`.
- The current importer's paired-`'M'`-event emission. Replaced by
  per-interval `'em'` events plus direct write to `migMatConst` for
  the $t=0$ interval.

## 5. Implementation Phases (TDD-driven)

Every phase is **test-first**. For each component below, the order is:
write tests describing intended behavior; verify they fail; implement
the minimum to pass; refactor.

### Phase 0: Branch hygiene and scaffolding

- Confirm branch is `feature/issue-82-time-varying-demography` based on
  `origin/nsp-yaml-revamp`.
- Add `test/unit/test_shapes.c` (math primitives) and
  `test/parity/` directory (for parity test harnesses).
- Wire `make test-shapes` and `make test-parity` into the build.

### Phase 1: Shape primitives and accessors

Tests first:
- `test_shape_constant`: `sizeAt`/`migAt` returns anchor value for CONST.
- `test_shape_exponential`: `sizeAt(t)` matches $N_0 e^{-\alpha (t-t_0)}$
  to $10^{-12}$ relative tolerance for a battery of inputs.
- `test_shape_linear`: similarly for $N_0 - \gamma (t - t_0)$.
- `test_integratedHazard_quadrature`: closed-form integrated hazards
  agree with high-resolution numerical quadrature (Simpson's rule, 1024
  steps) to $10^{-9}$ for all three shapes, sizes and migrations,
  random parameters.
- `test_drawWaitingTime_distribution`: empirical CDF of $T_k$ samples
  matches theoretical CDF (Kolmogorov-Smirnov, $n = 10^5$, $p > 0.05$)
  for all three shapes.
- `test_drawWaitingTime_inverse`: round-trip $T \to \xi \to T$ recovers
  to $10^{-12}$ relative tolerance.

Implementation:
- Add `Shape`, `popShape`, `migShape` to `discoal.h`.
- Add `sizeAt`, `migAt`, `integratedHazardSize`, `integratedHazardMig`,
  `drawWaitingTimeSize`, `drawWaitingTimeMig` to a new
  `src/core/shapes.c` / `shapes.h`.

### Phase 2: Q1 verification — Euler vs detSweepFreq under constant $N$

Tests first:
- `test_sweep_deterministic_constant_N_parity`: under
  $\alpha \in \{50, 200, 1000\}$ and constant $N$ (no shape change),
  run $10^4$ deterministic sweep replicates with the closed-form path
  and $10^4$ with the Euler path. Compare distributions of:
  - Number of segregating sites
  - $\pi$ (nucleotide diversity)
  - Tajima's D
  - SFS bin frequencies (per-bin chi-squared)
- Statistical thresholds: KS test on continuous statistics,
  $p > 0.01$ Bonferroni-corrected; chi-squared on SFS bins, same
  threshold.

If the test passes (expected): collapse the branch in 4.5 to
always-Euler. If it fails (drift visible at large $\alpha$): keep the
branch on shape type so constant-$N$ runs use the closed form.

### Phase 3: Inner-loop NHPP sampler with shape=CONSTANT only

Tests first:
- `test_neutral_phase_constant_regression`: the rewritten neutral phase
  with shape state initialized to all-CONSTANT must produce **byte-for-byte
  identical output** to the current `neutralPhaseGeneralPopNumber`
  given identical RNG seeds, for a battery of:
  - Single-pop neutral
  - Two-pop with `-en` size changes
  - Two-pop with `-m` constant migration
  - Two-pop with `-ed` split
- This regression is the safety net for the algorithmic refactor before
  any new shape support is added.

Implementation:
- Refactor `neutralPhaseGeneralPopNumber` to use the NHPP sampler with
  shape state. With CONST-only shapes the NHPP draws collapse to
  exponential draws and bit-equality should hold.

If bit-equality cannot be preserved (e.g., due to RNG-call-order
differences), fall back to statistical parity at $p > 0.01$ on the same
summary statistics as Phase 2, plus exact agreement of `eventNumber`,
`tDiv`, segregating-sites count, lineage-count trajectory at fixed
checkpoints. Document the divergence in the design doc and CHANGELOG.

### Phase 4: Add EXPONENTIAL shape

Tests first:
- `test_shape_exp_single_pop_msprime_parity`: single population with
  exponential growth (matched parameters in discoal and msprime). Run
  $10^4$ replicates each, compare:
  - SFS (per-bin chi-squared)
  - $\pi$, Tajima's D, segregating sites distributions (KS)
  - Pairwise coalescent-time distribution (KS)
- Threshold: $p > 0.01$ Bonferroni-corrected across statistics.
- `test_shape_exp_neutral_no_migration_msprime_parity`: 2-pop with
  split + exp growth in one branch, no migration. Same statistics.

Implementation:
- Wire `SHAPE_EXPONENTIAL` through the closed-form drawWaitingTime.
- Add CLI `-eg`, `-eG` parsing.
- Update YAML config to accept growth-rate fields.

### Phase 5: Sweep accessor wiring

Tests first:
- `test_sweep_constant_N_regression`: piecewise-constant `-en` sweep
  configurations produce byte-equal or statistically-indistinguishable
  output (per Phase 2 thresholds) before and after the `sizeAt` swap.
- `test_sweep_exp_growth_internal_consistency`: under exponential
  growth, sweep replicates produce sensible SFS shapes (sanity, not a
  parity test — there is no msprime gold standard for sweeps).
- `test_sweep_recurrent_constant_N_regression`: same regression for
  the recurrent-sweep code path.

Implementation:
- Replace `currentSizeRatio` and `sizeRatio[i]` reads in
  `proposeTrajectory`, `sweepPhaseEventsConditionalTrajectory`,
  `sweepPhaseEventsGeneralPopNumber` with `sizeAt` calls.
- Apply Q1 outcome from Phase 2 to the deterministic-mode dispatch.
- Audit `currentSize[]` usage. If `sizeAt` becomes the sole accessor
  in the inner loop, remove `currentSize[]` and any `'n'` dispatch
  that only mutates it. Single source of truth in `popShape[]`.

### Phase 6: Add LINEAR shape

Tests first:
- `test_shape_linear_quadrature` (already in Phase 1; revisit
  edge cases including $\gamma$ pushing $N \to 0$).
- `test_shape_linear_neutral_msprime_parity`: 2-pop with linear-growth
  branch, parity vs msprime. Same statistics as Phase 4.
- `test_shape_linear_zero_crossing_robustness`: linear shape with
  $\gamma > 0$ such that $N_0 - \gamma T \to 0^+$ at finite $T$ within
  the epoch. Verify the integrator forces a coalescence and does not
  segfault or produce NaN times.

Implementation:
- Wire `SHAPE_LINEAR` through the closed-form drawWaitingTime.
- Linear shape exposed via demes/YAML only at this stage.

### Phase 7: Migration shapes (`'em'` events)

Tests first:
- `test_migration_constant_window_msprime_parity`: 3-pop with single
  constant-rate migration window covering full simulation, parity vs
  msprime. (Sanity baseline — should match the Phase 4 results.)
- `test_migration_multi_window_msprime_parity`: the issue #82 fixture
  (`config_examples/demes_example.demes.yaml`). Compare against msprime
  running the same demes graph. **This is the closing of issue #82.**
- `test_migration_back_derivation_hack_removed`: the existing
  `migMat`-debug fprintf output is gone; importer writes `migMatConst`
  directly for $t=0$.

Implementation:
- Add `'em'` event parsing (CLI `-em`, `-eM`).
- Add `'em'` dispatch in the inner loop.
- Importer rewrite: piecewise migration matrices per interval, emit
  one `'em'` per pair per boundary.
- Delete `discoalFunctions.c:222-261`.

### Phase 8: Documentation and full-vocabulary parity

- Update `docs/population_structure.rst`: remove "Time-varying migration
  rates are not currently implemented" note. Document `-eg`, `-eG`,
  `-em`, `-eM`.
- Update `docs/demes_integration.md`: exponential and linear epochs are
  now supported; document any caveats.
- Update `docs/yaml_configuration.md` with growth-rate fields.
- Re-run the full demes-vs-msprime parity suite on a battery of demes
  examples (single pop, splits, migration, growth combinations).

## 6. Parity Testing Strategy

Two separate parity test suites, both runnable as Make targets and CI
jobs.

### 6.1 Sweep parity vs current discoal (regression)

**Purpose**: ensure the new code does not regress sweep behavior under
the demography current discoal supports (piecewise-constant).

**Reference**: the tip of `origin/nsp-yaml-revamp` at branch creation
time, built and committed as a reference binary into `test/parity/bin/`.

**Configurations** (full grid):
- Sweep mode: `'d'` (deterministic), `'s'` (stochastic forward),
  `'N'` (neutral stochastic). Recurrent sweep variant for each.
- Selection strength $\alpha \in \{50, 200, 1000\}$.
- Sweep timing $\tau \in \{0.01, 0.1, 0.5\}$.
- Demography: single-pop constant; single-pop with `-en` size changes
  at 2 epochs; 2-pop with split + migration.

**Replicates**: $10^4$ per configuration. Same RNG seeds for both
binaries.

**Statistics compared** per replicate:
- Segregating sites count
- $\pi$ (nucleotide diversity)
- Tajima's D
- Per-bin SFS frequencies
- Number of recombination events (if available in trees output)

**Tests**:
- For each summary statistic distribution: KS test, $p > 0.01$
  Bonferroni-corrected across the configuration grid.
- For per-bin SFS: chi-squared, same threshold.
- Phase 3 attempts byte-equality where the RNG sequence is preserved;
  Phase 5 falls back to statistical-equality once the sweep code
  changes.

### 6.2 Neutral parity vs msprime under demes models

**Purpose**: validate that discoal under the new shape vocabulary
faithfully simulates the demes models that msprime is the gold standard
for.

**Reference**: msprime ≥ 1.3 with `msprime.sim_ancestry` and
`msprime.Demography.from_demes`.

**Configurations**:
- Single pop, constant size — sanity baseline.
- Single pop, exponential growth (one epoch with `size_function: exponential`).
- Single pop, linear growth.
- Two pops with split + constant migration.
- Two pops with split + exp growth + constant migration.
- Three pops with two splits + multi-window migration (the issue #82
  fixture, `config_examples/demes_example.demes.yaml`).
- Three pops with growth + multi-window migration (full B2 stress test).

**Replicates**: $10^4$ per configuration, matched random seeds when
possible (independent runs otherwise; rely on distributional tests).

**Statistics compared**:
- SFS (per-bin)
- $\pi$ overall and per-pop
- $F_{\text{ST}}$ between pop pairs
- Tajima's D per-pop
- Pairwise coalescent-time distribution (extracted from trees)
- Number of segregating sites distribution

**Tests**:
- KS for continuous statistics, chi-squared for SFS, $p > 0.01$
  Bonferroni-corrected.
- For each demes model, document the test in
  `test/parity/test_msprime_<model_name>.py` as a runnable
  `pytest`-style harness with msprime as a dev dependency.

### 6.3 Q1 verification: detSweepFreq vs Euler under constant $N$

Specific to Phase 2. Goal: decide whether the deterministic-sweep
closed form is empirically distinguishable from the Euler step under
the regimes discoal currently uses.

**Setup**: single population, constant $N$, `'d'` mode, range of
$\alpha$, $\tau$, recombination rates. Run $10^4$ replicates with each
method.

**Statistics**: same as 6.1 plus the sweep frequency trajectory itself
(record the $x(\tau)$ at fixed grid points and compare distributions
pointwise via KS).

**Outcome**:
- If $p > 0.01$ Bonferroni across all configurations: the dispatch on
  shape type in section 4.5 collapses to always-Euler.
- If any configuration rejects: keep the dispatch; document where
  closed-form retention is required.

This is a one-time test that informs the design; document the result
in this design doc as an addendum after Phase 2 completes.

## 7. TDD Discipline

This project is large and touches the simulation engine. **Every commit
on this branch follows red-green-refactor.**

- New behavior gets a failing test first. The test names what the
  behavior is.
- Implementation does the minimum to pass. No speculative generality.
- Refactor with the test suite green.
- A commit may not introduce code without an accompanying test, except
  for non-functional changes (formatting, comments, doc-only) and
  test-harness scaffolding.
- The two parity suites (6.1, 6.2) and the unit test suite all run in
  CI. PRs may not merge with any of them red.

The `superpowers:test-driven-development` skill provides the discipline.
This design is the spec; the writing-plans skill produces the
implementation plan; subagent-driven-development executes phase-by-phase
with each phase ending green.

## 8. Risks and Open Questions

### 8.1 NHPP sampler performance

The neutral-phase inner loop today draws one $\text{Exp}$ per event;
NHPP draws one log per rate component per event. For $P$ populations
the migration draw count is $O(P^2)$, up from $O(P)$ today (we sum
total migration). For typical $P \le 10$ this is fine; for large $P$
(deep demes graphs) it may matter.

**Mitigation**: benchmark Phase 3 against the current code on a
stress-test config; if regression > 50%, profile and consider summing
total-coalescent-rate and total-migration-rate via a tighter NHPP
formulation.

### 8.2 Numerical stability

Exponential growth with very large $\alpha$ in coalescent units, or
linear growth with $\gamma$ near zero, can produce numerically tricky
$\log$ and $\exp$ arguments. The integrators must:
- Clip $\xi$-driven $T$ to the next epoch boundary explicitly.
- Detect $N(t) \to 0$ in linear shape and force coalescent.
- Not produce NaN times under any input the importer can emit.

Tested in Phase 6 (`test_shape_linear_zero_crossing_robustness`) and
Phase 4 (large-$\alpha$ growth).

### 8.3 Sweep shape branch decision

Phase 2's outcome (Q1 verification) determines the sweep dispatch
structure. If closed-form retention is required for constant $N$ to
preserve current discoal's exact distribution, the branch in 4.5
stays. Document either way.

### 8.4 Importer shape extraction for migrations

demes' migration windows are constant-rate within their window. But
multiple windows on overlapping pairs require the importer to compute
the *piecewise constant* matrix per disjoint interval and emit
boundary events. The interval algorithm:
1. Collect all window boundaries $\{t_b\}$ for all pairs.
2. Sort, deduplicate.
3. For each interval $(t_b, t_{b+1})$ and each pair $(i, j)$, set rate
   to the demes rate of any window covering this interval (windows for
   the same pair are disjoint by demes spec).
4. Emit `'em'` events at each $t_b$ for pairs whose rate changes.

Tested in Phase 7.

### 8.5 Convention 1 anchoring corner case

When a user mixes `-eg t p alpha` with a subsequent `-en t' p N`, the
size at $t'$ might be discontinuous (the exponential extrapolation does
not necessarily pass through `N`). This is intentional: `-en` is an
instantaneous snap. Document in `docs/population_structure.rst`.

When the importer composes shape transitions, it must compute the
anchor value of the new shape from the prior shape evaluated at the
transition time. Tested via the multi-epoch demes parity tests.

## 9. Out-of-Scope for This Design

- Time-varying selection coefficient $s(t)$ during sweeps. Sweep math
  with $s(t)$ is materially harder; the deterministic logistic does
  not have a closed form for non-constant $s$ even before adding $N(t)$.
- Time-varying recombination rate. Demes does not specify it; users
  can approximate with multiple discoal runs.
- Migration during sweep phase. Suspended today, suspended after this
  change.
- Shape-aware tree-sequence output annotations. Trees output already
  records discrete events; shape state is engine-internal.
- Removing CLI duplicates. `-en` stays exactly as is.

## 10. Concrete Deliverables

- `src/core/shapes.h`, `src/core/shapes.c` — shape state, accessors,
  drawWaitingTime functions. Unit-tested.
- `src/core/discoal.h` — extended `event` struct with `Shape newShape`
  and shape-state globals.
- `src/core/discoalFunctions.c` — refactored `neutralPhaseGeneralPopNumber`,
  `proposeTrajectory`, `sweepPhaseEventsConditionalTrajectory`,
  `sweepPhaseEventsGeneralPopNumber`, `recurrentSweepPhaseGeneralPopNumber`
  to use shape accessors. Back-derivation hack at lines 222-261 removed.
- `src/core/discoal_multipop.c` — CLI parsing of `-eg`, `-eG`, `-em`,
  `-eM`. Inner-loop dispatch case for `'g'` and `'em'`.
- `src/core/demesInterface.c` — emit `'g'`/`'em'` events; remove
  exp/linear rejection; piecewise migration matrix construction.
- `src/core/configInterface.c` — YAML acceptance of growth-rate and
  multi-window migration.
- `test/unit/test_shapes.c` — Phase 1 tests.
- `test/parity/test_sweep_regression.{c,sh}` — Phase 2 + Phase 5.
- `test/parity/test_msprime_*.py` — Phase 4, 6, 7 parity vs msprime.
- `docs/population_structure.rst`, `docs/demes_integration.md`,
  `docs/yaml_configuration.md` — documentation updates.
- This design document, committed at start of work.

## 11. Acceptance

- All parity tests in 6.1 and 6.2 green.
- Issue #82 fixture (`config_examples/demes_example.demes.yaml`)
  produces statistical parity with msprime.
- Existing YAML validation suite (`testing/yaml_validation_suite.sh`)
  remains 100% green.
- Existing unit test suite green.
- Documentation updated to reflect new flags and demes coverage.
- Back-derivation hack at `discoalFunctions.c:222-261` is gone.

## Addendum (2026-05-01): Q1 Verification Result

Phase 2 of the foundations implementation plan ran the Q1 verification
harness at `test/parity/q1_detsweep_verification.sh` across 9
(alpha, tau) configurations with $10^3$ replicates each (n=1000;
the plan called for $10^4$, scaled down for harness runtime —
$10^3$ still gives K-S statistical power well past the
Bonferroni-corrected threshold for 108 comparisons). The harness
compared the closed-form `detSweepFreq` and the Euler-step
`detSweepFreqEuler` paths under constant N via the
`--det-sweep-mode {closed,euler}` runtime flag.

**Configuration grid:** alpha in {50, 200, 1000}, tau in {0.01, 0.1, 0.5},
n=10, theta=10, rho=10, nsites=10000.

**Result:** **FAIL.**

29 of 108 comparisons rejected distributional equality at
Bonferroni-corrected $p < 9.26 \times 10^{-5}$
($\alpha = 0.01 / 108$). The pattern of failures is informative:

| alpha | tau=0.01 | tau=0.1 | tau=0.5 |
|---|---|---|---|
| 50 | 11/12 reject | 4/12 reject | 0/12 reject |
| 200 | 11/12 reject | 4/12 reject | 0/12 reject |
| 1000 | 0/12 reject | 0/12 reject | 0/12 reject |

- Recent sweeps (tau=0.01) with smaller alpha (50, 200) reject
  heavily — the Euler integrator's $O(\text{dt}^2)$ error per step
  accumulates over a relatively short trajectory whose tail-shape
  has strong influence on observable summary statistics when the
  sweep finished recently.
- Larger alpha (1000) does not reject at any tau, because the per-step
  error for an Euler step on $dx/d\tau = \alpha x(1-x)$ scales as
  $\alpha^2 x^2 (1-x)^2 \cdot \text{dt}^2$, but the trajectory duration
  also shrinks like $1/\alpha$, so total error scales like $\alpha$
  for fixed total integration time. With our fixed `tIncOrig`, the
  *number* of steps grows with alpha, mitigating the per-step error.
  The empirical effect is that closed-form and Euler track each
  other closely at large alpha despite the formal $O(\text{dt}^2)$
  truncation.
- Older sweeps (tau=0.5) wash out trajectory differences in
  post-sweep neutral coalescent.

**Conclusion (initial reading): keep dispatch on shape type, use
Euler only when N is non-constant.**

**Conclusion (revised, post-Q1): don't use Euler at all.** The
deterministic-sweep ODE $dx/d\tau = -\alpha_{\text{eff}}(\tau)\,x(1-x)$
is separable, so under arbitrary time-varying $N$ it has a closed-form
solution: $x(\tau) = x_0 e^{-A(\tau)}/(1-x_0+x_0 e^{-A(\tau)})$ with
$A(\tau) = \int_0^\tau \alpha_{\text{eff}}(s)\,ds = \alpha_0\,S(0,\tau)$,
where $S$ is the integral of `sizeAt`. For our shape vocabulary $S$ has
closed forms for all three types (see section 4.5). Phase 5 implements
`detSweepFreqGeneral` and a sister `integratedSizeRatio` and replaces
the deterministic-mode sweep call with a single closed-form invocation
that handles constant, exponential, and linear shapes uniformly. No
Euler step, no truncation error, no Q1-style discrepancy possible.

The Q1 verification, the harness in `test/parity/`, the
`--det-sweep-mode` runtime flag, and the `detSweepFreqEuler` function
are kept on the branch as a record of how the design pivoted. They
will be deleted in Phase 5 when the closed-form general path lands.

**Raw results:** `test/parity/q1_results/analysis.txt` on
`feature/issue-82-time-varying-demography`. Re-running the harness
(`./test/parity/q1_detsweep_verification.sh`) regenerates the .ms
and .stats files (gitignored due to size; ~30 MB total).

## Addendum (2026-05-01): Phase 4b msprime Parity Result

The Phase 4 SHAPE_EXPONENTIAL implementation was validated against
msprime 1.4.1 via the parity harness at `test/parity/phase4b_msprime/`.

**Configuration discovery:** The conversion between discoal's CLI/internal
units and msprime's parameters has three plausible (ploidy, population_size)
combinations. A trial sweep at the no-demography baseline (theta=5, rho=5,
n=6, reps=200) decisively picked **Convention B: ploidy=2, popsize=Ne** with
relative error 0.025 on mean pi. The other conventions failed:
A (ploidy=1, popsize=Ne) had rel_err 0.503 (half the coalescent timescale);
C (ploidy=1, popsize=2*Ne) had rel_err 0.023 (mathematically equivalent to B).
The discoal `EFFECTIVE_POPN_SIZE` parameter is therefore best interpreted as
diploid Ne — consistent with the documented "4N convention" since 4*Ne_diploid
generations per coalescent unit.

**Two implementation gotchas surfaced during conversion discovery, both
caught by the test harness:**

1. msprime's default mutation model (JC69 with `discrete_genome=True`) produces
   multi-allelic sites; discoal uses infinite-sites binary. The harness uses
   `model=msprime.BinaryMutationModel(), discrete_genome=False`.
2. msprime's `samples=n` under `ploidy=2` returns `2n` haploid samples; discoal's
   `n` is already haploid. The harness passes `n // ploidy` individuals.

**EXP parity tests at $10^3$ replicates each:**

- Single-pop EXP (`-eg 0.5 0 50` in discoal; matched
  `add_population_parameters_change` in msprime): 4 comparisons (segregating
  sites, pi, Tajima's D, folded SFS chi-squared) at Bonferroni-corrected
  $\alpha = 0.01/4 = 2.5 \times 10^{-3}$. **Result: PASS**, smallest $p = 6.7 \times 10^{-3}$
  (SFS chi-squared — the most discriminating statistic), all others $p > 10^{-2}$.

- Two-pop split + EXP in pop0 (`-p 2 6 6 -eg 0.3 0 50 -ed 1.0 0 1`; matched
  msprime `Demography` with `population_split` + `parameters_change`):
  3 comparisons (ss, pi, tajD) at $\alpha = 3.3 \times 10^{-3}$. **Result: PASS**,
  smallest $p = 0.31$.

**Overall: PASS.** The discoal SHAPE_EXPONENTIAL implementation is statistically
indistinguishable from msprime under the tested conditions, validating the
NHPP per-component sampler and the `'g'` event handler.

Detailed per-comparison output: `test/parity/phase4b_msprime/results/analysis.txt`.

### 2026-05-01 update: broadened single-pop sweep + alpha-conversion fix

After the initial single-point single-pop test passed, a configuration
sweep was run with alpha in {50, 200, 1000} x eg_time_cli in {0.1, 0.5, 1.0}
and 11 statistics per cell (segregating sites, pi, Tajima's D, Watterson's
theta, haplotype diversity, number of distinct haplotypes, plus per-bin
folded SFS). Initial run at REPS=1000 rejected 4 statistics in the
`a50_t0.1` cell (small alpha, recent growth onset) at Bonferroni p < 1e-4,
with msprime systematically more diverse than discoal — a 2x parameter
mis-mapping somewhere in the conversions.

Diagnosis: the bug was in `discoal_alpha_to_msp_growth`, which was dividing
by `4*Ne` but should divide by `2*Ne`. Discoal stores alpha as the
`rate_param` of an exponential shape evaluated in *internal* time units
(`size(t) = anchor * exp(-alpha * (t - t0))`, see `shapes.c`). The pair
coalescent rate in `neutralPhase` is `n*(n-1)/2` per internal unit, so 1
internal unit equals `2*Ne_diploid` generations. Therefore the per-generation
growth rate `g` that produces the same exponential is
`g = alpha / (2*Ne)`, not `alpha / (4*Ne)`. The CLI-time conversion is
unaffected: discoal multiplies user-supplied event times by 2.0 to enter
internal units, so `t_gen = t_cli * 4 * Ne` is correct (one factor of 2
from CLI->internal, one factor of 2 from internal->generations).

The bug was masked at the original single-point test (alpha=50, t_cli=0.5)
because the deeper growth event left less integrated growth contribution to
diversity, and was only revealed by the broadened sweep — small alpha plus
recent growth makes the absolute exponential rate strongly determining of
the post-event coalescent time distribution.

Fix: divide alpha by `2*Ne` instead of `4*Ne` in
`discoal_alpha_to_msp_growth`. Re-running the sweep at REPS=5000 across 99
comparisons (Bonferroni alpha = 1.01e-4): PASS, no rejections. Two-pop
test (REPS=5000) also PASSes.

This corroborates Convention B (ploidy=2, popsize=Ne) for the rate
mapping, the existing time conversion (`t_cli * 4 * Ne`), and the corrected
growth-rate conversion (`alpha / (2 * Ne)`). The broadened sweep is a
stronger validation than the original single-point test.
