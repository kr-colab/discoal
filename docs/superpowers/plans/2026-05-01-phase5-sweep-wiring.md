# Phase 5: Sweep Accessor Wiring with Closed-Form detSweepFreqGeneral

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire `sizeAt(popID, t)` into the sweep code paths (`proposeTrajectory`, `sweepPhaseEventsConditionalTrajectory`, `sweepPhaseEventsGeneralPopNumber`) so sweep simulations work correctly under time-varying $N$. For deterministic mode, replace the per-step `detSweepFreq(\tau, \alpha_{\text{eff}})` calls with `detSweepFreqGeneral(\alpha, A(\tau))` based on the integrated selection coefficient $A(\tau) = \alpha \cdot \int_0^\tau \text{sizeAt}(s)\,ds$ (per the spec revision: deterministic ODE is separable, no Euler step needed). Bit-equality preserved for `SHAPE_CONSTANT`-only configs by dispatching on `allShapesConstant()`. The Euler-step variant (`detSweepFreqEuler`) and the `--det-sweep-mode` flag added in Phase 2 are deleted as part of this phase.

**Architecture:** Add `integratedSizeRatio(popID, t0, T)` to `src/core/shapes.{h,c}` returning $\int_{t_0}^{t_0+T} \text{sizeAt}(s)\,ds$ (closed form per shape). Add `detSweepFreqGeneral(alpha, A)` to `src/core/alleleTraj.{h,c}` that reduces to `detSweepFreq(τ, alpha)` when `A = alpha * τ` (i.e., constant shape). In `proposeTrajectory`, save the global `popShape[]` state at entry, walk events forward updating `popShape[]` per `'n'`/`'g'` events as it iterates, restore at exit. In each step inside the trajectory loop, dispatch on `allShapesConstant()`: constant-shape → existing `detSweepFreq` path (preserves bit-equality); else → `detSweepFreqGeneral(alpha, A_prev + alpha*integratedSizeRatio(0, t_now, dt))`. Same dispatch in `sweepPhaseEventsGeneralPopNumber`. The conditional-trajectory consumer just needs `sizeAt(i, t)` reads for non-sweep-pop coalescent rates.

**Tech Stack:** C99, Unity, GNU make, bash for parity tests, niceStats for sweep-statistic comparisons.

**Spec:** `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md` §3.3 (sweep math), §4.5 (deterministic-mode dispatch revision), §6.1 (sweep parity tests).

**Deliverables:**
- `integratedSizeRatio` accessor in shapes module with unit tests.
- `detSweepFreqGeneral` in alleleTraj with unit test verifying it reduces to `detSweepFreq` exactly under constant shape.
- Sweep code paths read `sizeAt(i, t)` instead of `sizeRatio[i]` and `currentSizeRatio`.
- Deterministic mode: closed-form general path under non-constant shapes; bit-equal under constant shapes.
- `discoal_pre_phase5` reference binary + sweep bit-equality regression harness with 6 sweep configurations.
- Smoke test demonstrating sweep + EXP growth runs to completion.
- Deletion: `detSweepFreqEuler` function from `alleleTraj.{h,c}`, the `--det-sweep-mode` CLI flag, the Q1 verification harness in `test/parity/q1_*`, the Q1 result files (gitignored anyway).
- All Phase 1-4 regressions still PASS.
- Tag `phase5-sweep-wiring-complete`.

**Out of scope:**
- `SHAPE_LINEAR` for sweeps beyond a sanity check (Phase 6 broadens).
- `'em'` migration shapes (Phase 7).
- Importer changes (Phase 7).
- Full msprime parity for sweeps — there's no msprime gold standard for sweeps; sanity checks only.

**Convention reminder:** Never mention Claude or AI in commits/code/docs. Never use emojis. Stay on `feature/issue-82-time-varying-demography`. Do not push.

---

## Conversion math reminder

For the deterministic-sweep ODE $dx/d\tau = -\alpha_{\text{eff}}(\tau)\,x(1-x)$ with $\alpha_{\text{eff}}(\tau) = \alpha \cdot \text{sizeAt}(0, \tau)$:

$$x(\tau) = \frac{x_0\,e^{-A(\tau)}}{1 - x_0 + x_0\,e^{-A(\tau)}} \qquad A(\tau) = \int_0^\tau \alpha_{\text{eff}}(s)\,ds = \alpha \cdot \int_0^\tau \text{sizeAt}(0, s)\,ds$$

In the Stephan boundary-condition form (matching `detSweepFreq`):

$$x(\tau) = \frac{\varepsilon}{\varepsilon + (1-\varepsilon)\,e^{A(\tau) - A_{\tau_s}}}$$

with $\varepsilon = 0.05/\alpha$ and $A_{\tau_s} = -2\log\varepsilon$. For constant $\alpha$, $A(\tau) = \alpha\tau$ and $A_{\tau_s} = \alpha\tau_s$ — recovers `detSweepFreq(τ, α)` exactly.

**`integratedSizeRatio` closed forms** (sister to `integratedHazardSize` but integrating `sizeAt` directly instead of `1/sizeAt`):

- `SHAPE_CONSTANT`: $\int_{t_0}^{t_0+T} N_0 \,ds = N_0 \cdot T = \text{sizeAt}(0, t_0) \cdot T$
- `SHAPE_EXPONENTIAL`: $\int_{t_0}^{t_0+T} N_0 e^{-\alpha (s-t_a)}\,ds = \text{sizeAt}(0,t_0) \cdot (1 - e^{-\alpha T})/\alpha$ (degenerate to $N_0 T$ for $\alpha=0$)
- `SHAPE_LINEAR`: $\int_{t_0}^{t_0+T} (N_0 - \gamma(s-t_a))\,ds = \text{sizeAt}(0, t_0) \cdot T - \gamma T^2/2$

(Each uses `sizeAt(0, t_0)` so the formula handles `t_0 ≠ anchor_time` correctly.)

---

## Task 1: `integratedSizeRatio` accessor

**Files:**
- Modify: `src/core/shapes.h`
- Modify: `src/core/shapes.c`
- Modify: `test/unit/test_shapes.c`

- [ ] **Step 1: Failing tests**

Add to `test/unit/test_shapes.c` (before `main`):

```c
void test_integratedSizeRatio_constant(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 2.5;
    /* T=4 -> 2.5 * 4 = 10 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 10.0, integratedSizeRatio(0, 0.0, 4.0));
}

void test_integratedSizeRatio_exponential(void) {
    /* sizeAt(s) = 1.0 * exp(-0.5 * s); integral over [0, 4] = (1 - exp(-2))/0.5 */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    double expected = (1.0 - exp(-2.0)) / 0.5;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected, integratedSizeRatio(0, 0.0, 4.0));
}

void test_integratedSizeRatio_exponential_offset_t0(void) {
    /* If t0 != anchor_time, the relevant N_0 is sizeAt(t0). */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.2;
    popShape[0].anchor_time = 0.0;
    /* At t0=5, N(t0) = exp(-1). Integral over T=2: N(t0) * (1 - exp(-0.4))/0.2 */
    double N_t0 = exp(-1.0);
    double expected = N_t0 * (1.0 - exp(-0.4)) / 0.2;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected, integratedSizeRatio(0, 5.0, 2.0));
}

void test_integratedSizeRatio_exponential_alpha_zero(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 3.0;
    popShape[0].rate_param = 0.0;
    popShape[0].anchor_time = 0.0;
    /* alpha=0 -> matches constant: 3 * 5 = 15 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 15.0, integratedSizeRatio(0, 0.0, 5.0));
}

void test_integratedSizeRatio_linear(void) {
    /* sizeAt(s) = 5 - 0.5 * s. Integral over [0, 4] = 5*4 - 0.5*16/2 = 20 - 4 = 16 */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 5.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 16.0, integratedSizeRatio(0, 0.0, 4.0));
}

void test_integratedSizeRatio_quadrature_match(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.5;
    popShape[0].rate_param = 0.7;
    popShape[0].anchor_time = 2.0;
    double t0 = 3.5;
    double T = 4.0;
    int N = 16384;
    double dt = T / N;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double s_mid = (i + 0.5) * dt;
        sum += sizeAt(0, t0 + s_mid) * dt;
    }
    double closed = integratedSizeRatio(0, t0, T);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, sum, closed);
}
```

Register all 6 in `main`.

- [ ] **Step 2: Verify failure** — 6 link errors.

- [ ] **Step 3: Add declaration to `shapes.h`**

Before `#endif`:

```c
/* Integrated size ratio: returns int_{t0}^{t0+T} sizeAt(popID, s) ds.
 * Sister to integratedHazardSize but integrates sizeAt directly rather than
 * 1/sizeAt. Used by detSweepFreqGeneral to compute the integrated selection
 * coefficient A(tau) = alpha * integratedSizeRatio(0, sweep_start, tau). */
double integratedSizeRatio(int popID, double t0, double T);
```

- [ ] **Step 4: Implement in `shapes.c`**

After the existing functions:

```c
double integratedSizeRatio(int popID, double t0, double T) {
    Shape *s = &popShape[popID];
    double N0 = sizeAt(popID, t0);
    switch (s->type) {
        case SHAPE_CONSTANT:
            return N0 * T;
        case SHAPE_EXPONENTIAL: {
            double a = s->rate_param;
            if (a == 0.0) return N0 * T;
            return N0 * (1.0 - exp(-a * T)) / a;
        }
        case SHAPE_LINEAR: {
            double g = s->rate_param;
            return N0 * T - 0.5 * g * T * T;
        }
        default:
            return 0.0;
    }
}
```

- [ ] **Step 5: Verify pass** — all 6 new tests PASS, total tests now 57+.

- [ ] **Step 6: Commit**

```bash
git add src/core/shapes.h src/core/shapes.c test/unit/test_shapes.c
git commit -m "Add integratedSizeRatio for sweep-deterministic-mode integration

Returns int_{t0}^{t0+T} sizeAt(popID, s) ds in closed form for all
three shape types. Used by Phase 5's detSweepFreqGeneral to compute
the integrated selection coefficient A(tau) under time-varying N
without numerical integration."
```

---

## Task 2: `detSweepFreqGeneral` in alleleTraj

**Files:**
- Modify: `src/core/alleleTraj.h`
- Modify: `src/core/alleleTraj.c`
- Modify: `test/unit/test_alleleTraj.c`

- [ ] **Step 1: Add declaration to `alleleTraj.h`**

After the existing declarations:

```c
/* Generalized deterministic-sweep frequency under (possibly time-varying)
 * alpha_eff. A is the integrated selection coefficient int_0^tau alpha_eff(s) ds.
 * For constant alpha_eff = alpha, A = alpha * tau and this reduces to
 * detSweepFreq(tau, alpha) bit-exactly. */
double detSweepFreqGeneral(double alpha, double A);
```

- [ ] **Step 2: Implement in `alleleTraj.c`**

```c
double detSweepFreqGeneral(double alpha, double A) {
    double epsilon = 0.05 / alpha;
    /* For constant alpha_eff, sweep duration tau_s satisfies alpha*tau_s = -2*log(eps),
     * so A_ts = -2*log(eps). For time-varying alpha_eff, this anchors the boundary
     * to the user-specified alpha (interpretation: alpha is the "intended" sweep
     * strength; the trajectory follows the time-varying alpha_eff). */
    double A_ts = -2.0 * log(epsilon);
    double denom = epsilon + ((1.0 - epsilon) * exp(A - A_ts));
    return epsilon / denom;
}
```

- [ ] **Step 3: Add unit test**

Replace the body of `test_alleleTraj.c` with (preserving the existing test plus a new equivalence test):

```c
#include "unity.h"
#include "alleleTraj.h"
#include <math.h>

#ifndef TEST_RUNNER_MODE
void setUp(void) { }
void tearDown(void) { }
#endif

void test_detSweepFreqEuler_against_closed_form(void) {
    /* (existing test from Phase 2) — keeping until detSweepFreqEuler is removed in Task 11 */
    double alpha = 200.0;
    double epsilon = 0.05 / alpha;
    double ts = -2.0 * log(epsilon) / alpha;
    int N_steps = 100000;
    double dt = ts / N_steps;
    double x_euler = detSweepFreq(0.0, alpha);
    for (int i = 0; i < N_steps; i++) {
        x_euler = detSweepFreqEuler(x_euler, dt, alpha);
    }
    double x_closed = detSweepFreq(ts, alpha);
    TEST_ASSERT_DOUBLE_WITHIN(1e-3 * x_closed, x_closed, x_euler);
}

void test_detSweepFreqGeneral_reduces_to_detSweepFreq_constant(void) {
    /* For constant alpha, detSweepFreqGeneral(alpha, alpha*tau) == detSweepFreq(tau, alpha). */
    double alpha = 200.0;
    for (double tau = 0.0; tau < 0.1; tau += 0.005) {
        double x_old = detSweepFreq(tau, alpha);
        double x_new = detSweepFreqGeneral(alpha, alpha * tau);
        TEST_ASSERT_DOUBLE_WITHIN(1e-12, x_old, x_new);
    }
}

void test_detSweepFreqGeneral_handles_zero_A(void) {
    /* At A=0, x should be near 1 (start of sweep). */
    double x = detSweepFreqGeneral(200.0, 0.0);
    TEST_ASSERT_TRUE(x > 0.99);
    TEST_ASSERT_TRUE(x < 1.0);
}

void test_detSweepFreqGeneral_handles_large_A(void) {
    /* At A = -2*log(epsilon) = sweep duration, x ≈ epsilon. */
    double alpha = 200.0;
    double epsilon = 0.05 / alpha;
    double A_ts = -2.0 * log(epsilon);
    double x = detSweepFreqGeneral(alpha, A_ts);
    TEST_ASSERT_DOUBLE_WITHIN(epsilon * 0.01, epsilon, x);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_detSweepFreqEuler_against_closed_form);
    RUN_TEST(test_detSweepFreqGeneral_reduces_to_detSweepFreq_constant);
    RUN_TEST(test_detSweepFreqGeneral_handles_zero_A);
    RUN_TEST(test_detSweepFreqGeneral_handles_large_A);
    return UNITY_END();
}
#endif
```

- [ ] **Step 4: Build and verify**

```bash
make test_alleleTraj && ./build/test_alleleTraj
```

Expected: 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/alleleTraj.h src/core/alleleTraj.c test/unit/test_alleleTraj.c
git commit -m "Add detSweepFreqGeneral for time-varying alpha_eff sweeps

Generalizes detSweepFreq to take A = int alpha_eff(s) ds rather
than alpha and tau separately. Reduces to detSweepFreq exactly when
A = alpha*tau (constant). Boundary condition anchored to the
user-specified alpha (epsilon = 0.05/alpha, A_ts = -2*log(epsilon)).

Unit test verifies bit-equivalence to detSweepFreq for constant
alpha across a range of tau values."
```

---

## Task 3: Pre-Phase-5 reference binary

**Files:**
- Modify: `Makefile`

- [ ] **Step 1: Identify SHA**

```bash
git log -1 --pretty=format:'%h'
```

Use this short SHA as `PHASE5_REF_SHA`.

- [ ] **Step 2: Add Makefile recipe**

After the existing `discoal_pre_phase3:` recipe, add:

```makefile
# Phase 5 regression reference: discoal at the SHA where Phase 5 began.
PHASE5_REF_SHA = <fill in from Step 1>
.PHONY: discoal_pre_phase5
discoal_pre_phase5: build/discoal_pre_phase5

build/discoal_pre_phase5:
	@mkdir -p build
	@WT=$$(mktemp -d) && \
	  git worktree add --detach "$$WT" $(PHASE5_REF_SHA) && \
	  $(MAKE) -C "$$WT" discoal && \
	  cp "$$WT/build/discoal" build/discoal_pre_phase5 && \
	  git worktree remove "$$WT"
	@echo "Built pre-Phase-5 reference: build/discoal_pre_phase5"
```

- [ ] **Step 3: Build the reference**

```bash
make discoal_pre_phase5
ls -la build/discoal_pre_phase5
./build/discoal_pre_phase5 4 1 100 -t 1.0 -d 12345 67890 | head -5
```

- [ ] **Step 4: Commit**

```bash
git add Makefile
git commit -m "Add discoal_pre_phase5 Makefile recipe for sweep regression"
```

---

## Task 4: Sweep bit-equality regression harness

**Files:**
- Create: `test/parity/phase5_sweep_bit_equality.sh`
- Create: `test/parity/phase5_sweep_bit_equality_results/.gitignore`

- [ ] **Step 1: Create the harness**

```bash
#!/usr/bin/env bash
# Phase 5 regression: verify build/discoal sweep output is byte-identical to
# build/discoal_pre_phase5 for SHAPE_CONSTANT-only configurations.

set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

PRE="$ROOT/build/discoal_pre_phase5"
CUR="$ROOT/build/discoal"
[[ -x "$PRE" ]] || { echo "FAIL: build/discoal_pre_phase5 missing (make discoal_pre_phase5)"; exit 1; }
[[ -x "$CUR" ]] || { echo "FAIL: build/discoal missing (make discoal)"; exit 1; }

# Sweep configurations covering deterministic / stochastic-forward / neutral-stochastic
# modes, varying alpha and tau, single-pop and two-pop.
CONFIGS=(
  "single_pop_neutral_sweep|6 5 1000 -t 5 -r 5 -wn 0.05 -a 200 -x 0.5 -d 12345 67890"
  "single_pop_det_sweep|6 5 1000 -t 5 -r 5 -wd 0.05 -a 200 -x 0.5 -d 12345 67890"
  "single_pop_stoch_sweep|6 5 1000 -t 5 -r 5 -ws 0.05 -a 200 -x 0.5 -d 12345 67890"
  "single_pop_det_sweep_alpha50|6 5 1000 -t 5 -r 5 -wd 0.5 -a 50 -x 0.5 -d 12345 67890"
  "single_pop_det_sweep_size_change|6 5 1000 -t 5 -r 5 -wd 0.05 -a 200 -x 0.5 -en 0.1 0 0.5 -d 12345 67890"
  "two_pop_det_sweep|8 5 1000 -t 5 -r 5 -p 2 4 4 -wd 0.05 -a 200 -x 0.5 -ed 0.5 0 1 -d 12345 67890"
)

OUT="$HERE/phase5_sweep_bit_equality_results"
mkdir -p "$OUT"

fails=0
for cfg in "${CONFIGS[@]}"; do
  tag="${cfg%%|*}"
  args="${cfg##*|}"
  echo "=== $tag ==="
  # shellcheck disable=SC2086
  "$PRE" $args > "$OUT/${tag}_pre.ms" 2>"$OUT/${tag}_pre.err"; pre_rc=$?
  # shellcheck disable=SC2086
  "$CUR" $args > "$OUT/${tag}_cur.ms" 2>"$OUT/${tag}_cur.err"; cur_rc=$?
  if [[ $pre_rc -ne $cur_rc ]]; then
    echo "  FAIL exit code mismatch: pre=$pre_rc cur=$cur_rc"
    fails=$((fails+1))
    continue
  fi
  if diff <(sed '1d' "$OUT/${tag}_pre.ms") <(sed '1d' "$OUT/${tag}_cur.ms") > "$OUT/${tag}.diff"; then
    echo "  PASS"
  else
    echo "  FAIL  (see $OUT/${tag}.diff)"
    head -10 "$OUT/${tag}.diff" || true
    fails=$((fails+1))
  fi
done

if [[ $fails -eq 0 ]]; then
  echo
  echo "PASS: all ${#CONFIGS[@]} sweep configurations are byte-equal."
  exit 0
else
  echo
  echo "FAIL: $fails of ${#CONFIGS[@]} sweep configurations differ."
  exit 1
fi
```

Make executable.

- [ ] **Step 2: Create gitignore**

```
*.ms
*.err
*.diff
```

at `test/parity/phase5_sweep_bit_equality_results/.gitignore`.

- [ ] **Step 3: Smoke run** — at this point both binaries are from the same SHA so all 6 should PASS.

```bash
./test/parity/phase5_sweep_bit_equality.sh
```

Expected: 6/6 PASS.

- [ ] **Step 4: Commit**

```bash
git add test/parity/phase5_sweep_bit_equality.sh test/parity/phase5_sweep_bit_equality_results/.gitignore
git commit -m "Add Phase 5 sweep bit-equality regression harness

6 sweep configurations covering deterministic / stochastic-forward /
neutral-stochastic modes, single-pop / two-pop with size change,
varying alpha and tau."
```

---

## Task 5: Wire `sizeAt` reads into `proposeTrajectory`

This is the substantive change. `proposeTrajectory` walks events forward, currently tracking `currentSizeRatio` for `'n'` events. We replace this with `popShape[]` mutation (saved/restored), so future tasks can use `sizeAt(0, t)` directly.

**Files:**
- Modify: `src/core/discoalFunctions.c`

- [ ] **Step 1: Read the function**

```bash
sed -n '1764,1880p' src/core/discoalFunctions.c
```

Identify:
- The events-walk loop at the top (`for(i=currentEventNumber...`).
- The `if(events[i].type == 'n')` block updating `currentSizeRatio`.
- The inner trajectory loop and the `case 'd':`/'s'/'N' switch.

- [ ] **Step 2: Save popShape state at function entry**

After the existing local declarations near the top of the function body, add:

```c
    Shape saved_popShape[MAXPOPS];
    memcpy(saved_popShape, popShape, sizeof(saved_popShape));
```

(`memcpy` requires `<string.h>`; verify it's already included or add it.)

- [ ] **Step 3: Update events walk to mutate popShape**

In the events-walk loop, find:

```c
		if(events[i].type == 'n'){
			currentSizeRatio = events[i].popnSize;
			N = floor(N_0 *events[i].popnSize);
			if(currentSizeRatio > Nmax) Nmax = currentSizeRatio;
		}
```

Replace with:

```c
		if(events[i].type == 'n'){
			popShape[events[i].popID].type = SHAPE_CONSTANT;
			popShape[events[i].popID].anchor_value = events[i].popnSize;
			popShape[events[i].popID].rate_param = 0.0;
			popShape[events[i].popID].anchor_time = events[i].time;
			currentSizeRatio = events[i].popnSize;  /* keep for legacy code paths in this function */
			N = floor(N_0 * events[i].popnSize);
			if(currentSizeRatio > Nmax) Nmax = currentSizeRatio;
		}
		if(events[i].type == 'g'){
			popShape[events[i].popID].type = SHAPE_EXPONENTIAL;
			popShape[events[i].popID].anchor_value = sizeAt(events[i].popID, events[i].time);
			popShape[events[i].popID].rate_param = events[i].popnSize;
			popShape[events[i].popID].anchor_time = events[i].time;
			/* Nmax under EXP needs the running max sizeRatio; for now,
			 * approximate by sampling at the event time. Phase 6 may refine. */
			double sr_now = sizeAt(events[i].popID, events[i].time);
			if(sr_now > Nmax) Nmax = sr_now;
		}
```

- [ ] **Step 4: Restore popShape at function exit**

Find every `return` statement in `proposeTrajectory` (there's likely one main return at the end and possibly an `exit(1)` for the trajectory-too-bigly case). Before each, add:

```c
        memcpy(popShape, saved_popShape, sizeof(saved_popShape));
```

Be careful: the trajectory-too-bigly path calls `unlink(tempFilename)` and `exit(1)`, which terminates the process — restoration is moot there. Only restore before the normal `return` at the bottom of the function.

- [ ] **Step 5: Build and test**

```bash
make discoal
make discoal_pre_phase5
./test/parity/phase5_sweep_bit_equality.sh
```

Expected: all 6 sweep configs PASS. The `proposeTrajectory` events walk now updates `popShape` in addition to `currentSizeRatio`, but doesn't yet read from `popShape` — so output should be unchanged.

- [ ] **Step 6: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Save/restore popShape across proposeTrajectory; track 'g' events

proposeTrajectory now mutates popShape as it walks events forward
(in addition to the existing currentSizeRatio mutation), saving
and restoring the global state around the function call so the
caller's view of popShape is unchanged. Future tasks read sizeAt
inside the trajectory loop, which depends on this state being
correct during the walk.

The existing currentSizeRatio updates remain as a backward-
compatible cache for the trajectory loop's inner step computations
(line ~1820), to be removed in Task 6 once those reads are
migrated to sizeAt."
```

---

## Task 6: Wire `sizeAt` into the trajectory inner loop

Replace `currentSizeRatio` reads inside the inner trajectory loop with `sizeAt(0, currentTime + ttau)`. For `case 'd':` (deterministic mode), keep existing `detSweepFreq(ttau, alpha * currentSizeRatio)` for now; replacement to `detSweepFreqGeneral` lands in Task 8.

**Files:**
- Modify: `src/core/discoalFunctions.c`

- [ ] **Step 1: Inspect the inner loop**

The relevant block is around line 1820 (after the events `for` walk):

```c
		if(minF < 1.0/(2.*N)) minF = 1.0/(2.*N);
		tInc = 1.0 / (deltaTMod * N);
		while( x > 1.0/(2.*N) && (currentTime+ttau) < localNextTime){
			ttau += tIncOrig;
			if(x > minF && insweepphase){
				switch(sweepMode){
					case 'd':
					x = detSweepFreq(ttau, alpha * currentSizeRatio);
					break;
					case 's':
					x = 1.0 - genicSelectionStochasticForwardsOptimized(tInc, (1.0 - x), alpha * currentSizeRatio);
					break;
					case 'N':
					x = neutralStochasticOptimized(tInc, x);
					break;
				}
			}
			else{
				...
			}
```

`currentSizeRatio` is used as the size at the current step. We replace with `sizeAt(0, currentTime + ttau)`. For `SHAPE_CONSTANT` (the only shape currently produced by the events walk in this function), `sizeAt` returns `popShape[0].anchor_value` which equals the just-set `currentSizeRatio`. Bit-equal.

For `case 's'` and `case 'N'`, also replace any use of `N` inside that depends on `currentSizeRatio`. Specifically `tInc = 1.0 / (deltaTMod * N)` — `N` was set from `currentSizeRatio` earlier. We can read it directly: `tInc = 1.0 / (deltaTMod * N_0 * sizeAt(0, currentTime + ttau))`.

To preserve bit-equality, the order of operations matters. Probably easiest:

```c
		while( x > 1.0/(2.*N) && (currentTime+ttau) < localNextTime){
			ttau += tIncOrig;
			double sr_now = sizeAt(0, currentTime + ttau);
			N = floor(N_0 * sr_now);
			tInc = 1.0 / (deltaTMod * N);
			if(x > minF && insweepphase){
				switch(sweepMode){
					case 'd':
					x = detSweepFreq(ttau, alpha * sr_now);
					break;
					case 's':
					x = 1.0 - genicSelectionStochasticForwardsOptimized(tInc, (1.0 - x), alpha * sr_now);
					break;
					case 'N':
					x = neutralStochasticOptimized(tInc, x);
					break;
				}
			}
			else{
				insweepphase = 0;
				tInc = 1.0 / (deltaTMod * N );
				x = neutralStochasticOptimized(tInc, x);
			}
			...
		}
```

This reads `sizeAt` once per iteration and uses it consistently. Under `SHAPE_CONSTANT` (set by the events walk), `sr_now = currentSizeRatio` (the same value).

But wait — `N` and `tInc` were set OUTSIDE the inner loop in the original code (only updated when 'n' fired). Inside the loop, they're constant for the epoch. With the new code, we recompute them every step.

For `SHAPE_CONSTANT`, `sr_now` is constant within an epoch, so `N` and `tInc` are also constant. Bit-equal. (Floating-point recomputation may produce identical bits if the operations are deterministic, which they are.)

For `SHAPE_EXPONENTIAL`, `sr_now` varies, so `N` and `tInc` vary. New behavior; that's the point.

- [ ] **Step 2: Apply the change**

Replace the inner loop carefully per the structure above.

- [ ] **Step 3: Build and test**

```bash
make discoal
./test/parity/phase5_sweep_bit_equality.sh
```

Expected: all 6 PASS. Bit-equality should hold under the constant-shape configs.

If FAIL: investigate where the floating-point order diverged. The most likely culprit is the moment `N = floor(N_0 * sr_now)` or `tInc = ...` happens once per loop instead of once per epoch.

- [ ] **Step 4: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Read size via sizeAt in proposeTrajectory inner loop

Inner trajectory loop now reads sizeAt(0, currentTime + ttau)
each step, replacing the per-epoch currentSizeRatio scalar.
Under SHAPE_CONSTANT (the only shape currently produced by the
events walk), sizeAt returns the anchor_value which equals the
just-set currentSizeRatio, preserving bit-equality.

Under SHAPE_EXPONENTIAL the size now varies per-step within
an epoch — the intended behavior. Phase 5 Task 8 will replace
the deterministic-mode detSweepFreq call with the closed-form
detSweepFreqGeneral for the variable-size case."
```

---

## Task 7: Wire `sizeAt` into `sweepPhaseEventsConditionalTrajectory`

The trajectory consumer uses `sizeRatio[i]` to scale per-pop coalescent rates. Replace with `sizeAt(i, cTime + ttau)`. Function is around line 2143.

**Files:**
- Modify: `src/core/discoalFunctions.c`

- [ ] **Step 1: Locate the rate computations**

```bash
grep -n 'sizeRatio\[' src/core/discoalFunctions.c | grep -v '^15\|^17\|^18'
```

This skips proposeTrajectory and the neutralPhase functions. Look at lines ~2247 (in `sweepPhaseEventsConditionalTrajectory`):

```c
				cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 * tIncOrig / sizeRatio[i];
```

And similar at lines ~2225-2230 for `pCoalB`, `pCoalb` etc. (per-pop in the sweep population).

- [ ] **Step 2: Replace reads**

Replace each `sizeRatio[i]` (or `sizeRatio[0]` for the sweep pop) with `sizeAt(i, cTime + ttau)` (or `sizeAt(0, cTime + ttau)`).

Be systematic — every `sizeRatio[...]` reference in the function body becomes a `sizeAt(...)` call.

- [ ] **Step 3: Build and test**

```bash
make discoal
./test/parity/phase5_sweep_bit_equality.sh
make run_tests
```

Expected: all PASS.

- [ ] **Step 4: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Read sizes via sizeAt in sweepPhaseEventsConditionalTrajectory

Replaces all sizeRatio[i] reads in the function body with
sizeAt(i, cTime + ttau). Under SHAPE_CONSTANT, sizeAt returns
popShape[i].anchor_value which equals the prior sizeRatio[i]
seed; output is bit-equal."
```

---

## Task 8: Replace deterministic-mode `detSweepFreq` with `detSweepFreqGeneral` in `proposeTrajectory`

**Files:**
- Modify: `src/core/discoalFunctions.c`

The deterministic-mode line in `proposeTrajectory` is:

```c
case 'd':
    x = detSweepFreq(ttau, alpha * sr_now);
    break;
```

(After Task 6's edit. `sr_now = sizeAt(0, currentTime + ttau)`.)

Under SHAPE_CONSTANT, `detSweepFreq(ttau, alpha * sr_now) == detSweepFreq(ttau, alpha * sizeRatio_current)`, exactly the prior behavior.

For non-constant shapes, we need the closed-form general path. The dispatch:

```c
case 'd':
    if (allShapesConstant()) {
        /* Bit-equal with pre-Phase-5: existing per-step detSweepFreq with current alpha_eff. */
        x = detSweepFreq(ttau, alpha * sr_now);
    } else {
        /* Time-varying alpha_eff: use the closed-form general formula with the
         * incrementally-tracked integrated alpha. */
        A_now = A_prev + alpha * integratedSizeRatio(0, currentTime + ttau - tIncOrig, tIncOrig);
        x = detSweepFreqGeneral(alpha, A_now);
        A_prev = A_now;
    }
    break;
```

Add `double A_now = 0.0, A_prev = 0.0;` to the local variable declarations near the top of `proposeTrajectory`.

(Reset `A_prev = 0.0` whenever a new sweep starts. Since proposeTrajectory is called once per sweep, the local initialization at function entry is sufficient.)

- [ ] **Step 1: Add A tracking variables**

Add `double A_now = 0.0, A_prev = 0.0;` to the local declarations in `proposeTrajectory`.

- [ ] **Step 2: Replace the deterministic-mode case**

Apply the dispatch from above.

- [ ] **Step 3: Build and test**

```bash
make discoal
./test/parity/phase5_sweep_bit_equality.sh
make run_tests
```

Expected: all 6 sweep configs PASS (constant shapes route through `detSweepFreq`, identical to pre-Phase-5).

- [ ] **Step 4: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Dispatch deterministic-mode sweep on shape constancy

For SHAPE_CONSTANT-only configs, keep the existing
detSweepFreq(ttau, alpha*sr_now) path bit-equal with pre-Phase-5.
For non-constant shapes, use detSweepFreqGeneral with
A_now = A_prev + alpha*integratedSizeRatio(0, t-dt, dt)
incrementally tracked.

Mirrors the Phase 4 NHPP dispatch pattern in the neutral phase:
common case stays exact; non-constant case uses the time-varying
closed form with no truncation error."
```

---

## Task 9: Same dispatch in `sweepPhaseEventsGeneralPopNumber` (recurrent sweep)

The recurrent sweep path (function at line 1885) has its own deterministic-mode line at ~1958 and ~1965:

```c
case 'd':
    if (detSweepMode == 0) {
        x = detSweepFreq(ttau, alpha * sizeRatio[0]);
    } else {
        x = detSweepFreqEuler(x, tIncOrig, alpha * sizeRatio[0]);
    }
    break;
```

(This was the Phase 2 dispatch on `--det-sweep-mode`.)

Replace with the same structure as Task 8 (allShapesConstant dispatch using detSweepFreq vs detSweepFreqGeneral). The Euler branch is removed entirely (deletion lands in Task 11).

**Files:**
- Modify: `src/core/discoalFunctions.c`

- [ ] **Step 1: Find and update**

The function body needs:
- A_now, A_prev locals (initialized to 0 at function entry).
- A loop-internal `sr_now = sizeAt(0, cTime + ttau)` similar to Task 6.
- Dispatch: constant → existing `detSweepFreq`; non-constant → `detSweepFreqGeneral`.

Apply the changes.

- [ ] **Step 2: Build and test**

```bash
make discoal
./test/parity/phase5_sweep_bit_equality.sh
make run_tests
```

Expected: all PASS.

- [ ] **Step 3: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Same constant/general dispatch in sweepPhaseEventsGeneralPopNumber

Deterministic-mode sweep in the recurrent-sweep code path now
routes through the same allShapesConstant dispatch as
proposeTrajectory: constant -> detSweepFreq (bit-equal),
non-constant -> detSweepFreqGeneral with incremental A tracking.

The previous --det-sweep-mode dispatch (Phase 2 verification
artifact) is gone from this site; the runtime flag itself is
removed in a follow-up task."
```

---

## Task 10: Smoke test sweep + EXP growth

A simple test that runs `discoal` with both `-w*` (sweep) and `-eg` (EXP growth) flags, confirming the simulator runs to completion and produces ms-format output.

**Files:**
- Create: `test/parity/phase5_sweep_exp_smoke.sh`

- [ ] **Step 1: Create the script**

```bash
#!/usr/bin/env bash
# Phase 5 smoke: confirm sweep + -eg runs to completion.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"
[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

CONFIGS=(
  "det_sweep_with_eg|6 5 1000 -t 5 -r 5 -wd 0.05 -a 200 -x 0.5 -eg 0.5 0 50 -d 12345 67890"
  "stoch_sweep_with_eg|6 5 1000 -t 5 -r 5 -ws 0.05 -a 200 -x 0.5 -eg 0.5 0 50 -d 12345 67890"
)

for cfg in "${CONFIGS[@]}"; do
  tag="${cfg%%|*}"
  args="${cfg##*|}"
  echo "=== $tag ==="
  if "$DISCOAL" $args > /tmp/$tag.ms 2>/tmp/$tag.err; then
    n_segsites=$(grep -c '^segsites:' /tmp/$tag.ms || echo 0)
    if [[ $n_segsites -ge 1 ]]; then
      echo "  PASS ($n_segsites segsites entries)"
    else
      echo "  FAIL: no segsites lines"
      exit 1
    fi
  else
    echo "  FAIL: discoal exited nonzero"
    cat /tmp/$tag.err >&2
    exit 1
  fi
done

echo
echo "PASS: all sweep+eg smoke configurations ran to completion."
```

Make executable.

- [ ] **Step 2: Run**

```bash
./test/parity/phase5_sweep_exp_smoke.sh
```

Expected: PASS for both configs.

- [ ] **Step 3: Commit**

```bash
git add test/parity/phase5_sweep_exp_smoke.sh
git commit -m "Add sweep + EXP growth smoke test

Confirms discoal runs to completion when -wd/-ws and -eg are
combined. Sanity check that the time-varying alpha_eff path
in proposeTrajectory and sweepPhaseEvents* doesn't crash on
typical parameters."
```

---

## Task 11: Remove `detSweepFreqEuler` and `--det-sweep-mode`

Per the spec revision, the Euler step and runtime flag are no longer needed. Remove them.

**Files:**
- Modify: `src/core/alleleTraj.h` (remove declaration)
- Modify: `src/core/alleleTraj.c` (remove implementation)
- Modify: `src/core/discoal.h` (remove `detSweepMode` global)
- Modify: `src/core/discoal_multipop.c` (remove `detSweepMode` definition + `--det-sweep-mode` parser block)
- Modify: `src/core/discoalFunctions.c` (remove Euler dispatch leftovers, if any)
- Modify: `test/unit/test_alleleTraj.c` (remove the Euler test)
- Modify: `Makefile` (no change expected; just verify still builds)
- Delete: `test/parity/q1_detsweep_verification.sh`
- Delete: `test/parity/q1_detsweep_analyze.py`

The Q1 verification harness and result files are no longer relevant once the Euler path is gone. Keep the design-spec addendum about the Q1 result (it's historical context); just remove the harness scripts.

- [ ] **Step 1: Remove `detSweepFreqEuler` declaration and implementation**

```bash
grep -n 'detSweepFreqEuler' src/core/alleleTraj.{h,c}
```

Remove the declaration block from `alleleTraj.h` and the implementation block from `alleleTraj.c`. Also remove the related comment block above the implementation.

- [ ] **Step 2: Remove `detSweepMode` global**

```bash
grep -n 'detSweepMode' src/core/{discoal.h,discoal_multipop.c,discoalFunctions.c}
```

Remove:
- Declaration in `discoal.h`.
- Definition (`int detSweepMode = 0;`) in `discoal_multipop.c`.
- The `--det-sweep-mode` long-option parser block in `getParameters` (in `discoal_multipop.c`).
- Any stale `if (detSweepMode == 0)` conditions in `discoalFunctions.c` (Tasks 8 and 9 should have already removed them).

- [ ] **Step 3: Remove the Euler unit test**

In `test/unit/test_alleleTraj.c`, delete `test_detSweepFreqEuler_against_closed_form` and its `RUN_TEST` line. The remaining tests (`test_detSweepFreqGeneral_*`) stay.

- [ ] **Step 4: Delete Q1 harness files**

```bash
rm test/parity/q1_detsweep_verification.sh
rm test/parity/q1_detsweep_analyze.py
```

The `test/parity/q1_results/` directory and its `.gitignore` + `analysis.txt` stay as a historical record.

- [ ] **Step 5: Build and test**

```bash
make discoal
make test_alleleTraj && ./build/test_alleleTraj
make run_tests
./test/parity/phase5_sweep_bit_equality.sh
./test/parity/phase5_sweep_exp_smoke.sh
```

Expected: all PASS. The flag is gone, but observable behavior is unchanged for the regression configs.

- [ ] **Step 6: Commit**

```bash
git add -A src/core/ test/unit/ test/parity/q1_*
git commit -m "Remove detSweepFreqEuler and --det-sweep-mode

Per the spec revision, the deterministic-sweep ODE under time-
varying N has a closed-form solution via detSweepFreqGeneral; no
Euler step is needed. The Phase 2 Q1 verification was conclusive
(detSweepFreq != Euler under constant N), and the Euler path was
never used at runtime after Phase 5 wired in the closed-form
general dispatch.

Deleted:
- detSweepFreqEuler function and unit test
- detSweepMode global and the --det-sweep-mode CLI flag
- Q1 verification harness scripts (test/parity/q1_*.sh, *.py)

Kept:
- The Q1 spec addendum (historical context for why we chose the
  closed-form approach over Euler)
- test/parity/q1_results/analysis.txt (the recorded result)"
```

---

## Task 12: Final regression sweep + tag

- [ ] **Step 1: Full regression**

```bash
make discoal discoal_pre_phase3 discoal_pre_phase5
./test/parity/phase3_bit_equality.sh
./test/parity/phase4_eg_smoke.sh
./test/parity/phase5_sweep_bit_equality.sh
./test/parity/phase5_sweep_exp_smoke.sh
make run_tests
```

Expected: all PASS.

- [ ] **Step 2: Tag**

```bash
git tag -a phase5-sweep-wiring-complete -m "$(cat <<'TAGEOF'
Phase 5 of issue-82 design complete

Sweep code paths now read sizes via sizeAt(i, t) instead of
sizeRatio[i] / currentSizeRatio. Deterministic mode dispatches
on allShapesConstant():
- Constant: existing detSweepFreq(tau, alpha*sr_now). Bit-equal
  with pre-Phase-5 across the 6-config sweep regression grid.
- Non-constant: closed-form detSweepFreqGeneral(alpha, A_now)
  with A_now = A_prev + alpha*integratedSizeRatio(0, t-dt, dt)
  incrementally tracked. No Euler step; the deterministic ODE's
  separability gives an exact solution under time-varying N.

Sweep + EXP growth runs to completion (smoke test).

Removed: detSweepFreqEuler, --det-sweep-mode flag, Q1 verification
harness scripts. The closed-form path supersedes the Euler-step
approach validated as inferior in Phase 2.

Out of scope for Phase 5: SHAPE_LINEAR sweep parity (Phase 6),
'em' migration shapes (Phase 7), importer changes (Phase 7).
TAGEOF
)"
```

- [ ] **Step 3: List tags**

```bash
git tag --list 'phase*'
```

Expected: 6 tags including `phase5-sweep-wiring-complete`.

(Push not done yet; user pushes when ready.)

---

## Self-Review Checklist

- [ ] All tasks have explicit file paths.
- [ ] Every code step contains the actual code an engineer needs.
- [ ] Bit-equality preserved for constant-shape sweep configs at every step.
- [ ] `detSweepFreqGeneral` reduces to `detSweepFreq` exactly under constant alpha.
- [ ] `integratedSizeRatio` handles all three shapes including alpha=0 degenerate.
- [ ] popShape state is saved/restored across `proposeTrajectory`.
- [ ] `detSweepFreqEuler` and `--det-sweep-mode` are gone after Task 11.
- [ ] No emojis, no Claude/AI references.

## What's Next After This Plan

- **Plan: Phase 6** — Add `SHAPE_LINEAR` support (already in the math primitives; need engine wiring + msprime parity if applicable). Reuses the Phase 4b harness.
- **Plan: Phase 7** — `'em'` migration shape events, importer rewrite, removal of back-derivation hack. **Closes issue #82.**
- **Plan: Phase 8** — Documentation, full-vocabulary parity sweep, CHANGELOG, PR prep.
