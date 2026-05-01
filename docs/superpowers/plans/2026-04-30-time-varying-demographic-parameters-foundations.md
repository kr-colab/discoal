# Time-Varying Demographic Parameters — Foundations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the `Shape` state, accessors, and NHPP draw primitives, plus run the Q1 verification deciding whether deterministic-sweep mode under constant $N$ can collapse to always-Euler.

**Architecture:** Self-contained `shapes` module (`src/core/shapes.{h,c}`) holding per-population/per-pair shape state and closed-form integrated-hazard inversions for `CONSTANT`, `EXPONENTIAL`, `LINEAR`. Unity-based unit tests against numerical quadrature and KS distribution checks. Q1 is a runtime toggle on `proposeTrajectory` plus a parity harness comparing summary-statistic distributions.

**Tech Stack:** C99 (existing), Unity test framework (`extern/Unity`), Python+SciPy for distribution comparisons (`scripts/validation/`), bash for orchestration.

**Spec:** `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md`

**Deliverables of this plan:**
- A tested `shapes` module (no integration into the simulation engine yet — the runtime still uses `currentSize[]` and `migMatConst[]`).
- A Q1 verification outcome: addendum to the design spec stating whether `detSweepFreq` and Euler are statistically indistinguishable under constant $N$.
- All new tests wired into `make run_tests`.
- Foundation for follow-up plans (NHPP sampler refactor, sweep accessor wiring, importer changes).

**Out of scope for this plan:** any change to `neutralPhaseGeneralPopNumber`, `proposeTrajectory` behavior beyond a temporary CLI toggle, importer changes, removal of the back-derivation hack, new event types, new CLI flags beyond the temporary `--det-sweep-mode` toggle.

---

## Convention for every task

Every task is one TDD cycle: write the failing test, run it, implement, run again, commit. Steps are explicit. **Commit messages do not mention Claude or AI.**

If a step's "Run" command fails for an unexpected reason (compiler error, missing dependency), stop and investigate. Do not skip ahead.

---

## Task 1: Scaffolding — `Shape` type and module skeleton

**Files:**
- Modify: `src/core/discoal.h` (add `Shape` typedef + extern globals)
- Create: `src/core/shapes.h`
- Create: `src/core/shapes.c`
- Modify: `src/core/discoal_multipop.c` (define globals once, no behavior change)
- Modify: `Makefile` (add `shapes.c` to discoal build dependencies)

- [ ] **Step 1: Add `Shape` typedef to `discoal.h`**

In `src/core/discoal.h`, immediately after the `event` typedef (around line 90), add:

```c
/******************************************************************************/
/* Shape state for time-varying demographic parameters                         */
/******************************************************************************/

typedef enum {
    SHAPE_CONSTANT = 0,
    SHAPE_EXPONENTIAL = 1,
    SHAPE_LINEAR = 2
} ShapeType;

typedef struct {
    ShapeType type;
    double anchor_value;   /* size or rate at anchor_time */
    double rate_param;     /* alpha for EXP, gamma for LIN, unused for CONST */
    double anchor_time;    /* time at which anchor_value applies */
} Shape;
```

Then, in the globals section (around line 158, near `migMat`), add:

```c
Shape popShape[MAXPOPS];
Shape migShape[MAXPOPS][MAXPOPS];
```

(Plain definitions are fine — discoal already uses `-fcommon` so duplicate definitions in TUs are merged. Match existing style.)

- [ ] **Step 2: Create `src/core/shapes.h`**

```c
#ifndef SHAPES_H
#define SHAPES_H

#include "discoal.h"

/* Evaluate the size/rate of a shape at absolute time t. */
double sizeAt(int popID, double t);
double migAt(int srcPopID, int dstPopID, double t);

/* Evaluate the integrated hazard for k lineages from t0 over duration T.
 * Used in tests to verify drawWaitingTime against numerical quadrature. */
double integratedHazardSize(int popID, double t0, double T, int k);
double integratedHazardMig(int srcPopID, int dstPopID, double t0, double T, int k);

/* Draw a waiting time T such that integratedHazard(t0, T, k) == xi.
 * Returns -1.0 if xi exceeds the integrated hazard over [t0, infinity)
 * or if the linear shape would drive N(t) to zero before the draw resolves. */
double drawWaitingTimeSize(int popID, double t0, double xi, int k);
double drawWaitingTimeMig(int srcPopID, int dstPopID, double t0, double xi, int k);

#endif /* SHAPES_H */
```

- [ ] **Step 3: Create `src/core/shapes.c` with stubs**

```c
#include <math.h>
#include <stdlib.h>
#include "shapes.h"

double sizeAt(int popID, double t) {
    (void)popID; (void)t;
    return 0.0;  /* TBD: implemented in Tasks 3-5 */
}

double migAt(int srcPopID, int dstPopID, double t) {
    (void)srcPopID; (void)dstPopID; (void)t;
    return 0.0;
}

double integratedHazardSize(int popID, double t0, double T, int k) {
    (void)popID; (void)t0; (void)T; (void)k;
    return 0.0;
}

double integratedHazardMig(int srcPopID, int dstPopID, double t0, double T, int k) {
    (void)srcPopID; (void)dstPopID; (void)t0; (void)T; (void)k;
    return 0.0;
}

double drawWaitingTimeSize(int popID, double t0, double xi, int k) {
    (void)popID; (void)t0; (void)xi; (void)k;
    return -1.0;
}

double drawWaitingTimeMig(int srcPopID, int dstPopID, double t0, double xi, int k) {
    (void)srcPopID; (void)dstPopID; (void)t0; (void)xi; (void)k;
    return -1.0;
}
```

The "TBD" comment is fine here because it scopes to one stub function being filled in by the next named tasks. Remove every "TBD" by Task 16.

- [ ] **Step 4: Wire `shapes.c` into the discoal build**

In `Makefile`, find the `discoal:` recipe (line 55-57). Add `$(SRC_CORE)/shapes.c $(SRC_CORE)/shapes.h` to the dependency list and `$(SRC_CORE)/shapes.c` to the compiler invocation.

The recipe currently reads:

```makefile
discoal: libyaml demes-c libcyaml $(SRC_CORE)/discoal_multipop.c $(SRC_CORE)/discoalFunctions.c $(SRC_CORE)/discoal.h $(SRC_CORE)/discoalFunctions.h $(SRC_CORE)/ancestrySegment.c ...
```

After modification, append `$(SRC_CORE)/shapes.h $(SRC_CORE)/shapes.c` to the dependencies and `$(SRC_CORE)/shapes.c` to the gcc command line. Apply the same change to the `discoal_legacy_rng:` and `discoal_edited:` recipes (lines 61-68) — these are alternate builds used by parity tests.

- [ ] **Step 5: Verify build succeeds**

Run: `make discoal`
Expected: build succeeds with no warnings about `shapes.c`.

If GCC complains about unused `Shape` definitions in `discoal.h`, that's expected and harmless because `Shape popShape[MAXPOPS]` etc. will be referenced in Task 3.

- [ ] **Step 6: Smoke-test that nothing regresses**

Run: `./build/discoal 4 1 100 -t 1.0 > /tmp/discoal_smoke.txt && head -5 /tmp/discoal_smoke.txt`
Expected: typical ms-format output (`//`, segsites, positions, haplotypes). Confirms the shapes module did not break the main build.

- [ ] **Step 7: Commit scaffolding**

```bash
git add src/core/discoal.h src/core/shapes.h src/core/shapes.c Makefile
git commit -m "Add Shape type and shapes module scaffold

Stubs for sizeAt, migAt, integratedHazard*, drawWaitingTime*
implementations to be filled in subsequent tasks. No runtime
behavior change yet."
```

---

## Task 2: Test runner scaffolding for `test_shapes`

**Files:**
- Create: `test/unit/test_shapes.c`
- Modify: `Makefile` (add `test_shapes` target, add to `run_tests`)

- [ ] **Step 1: Create `test/unit/test_shapes.c` with empty Unity runner**

```c
#include "unity.h"
#include "discoal.h"
#include "shapes.h"
#include <math.h>
#include <stdlib.h>

extern long seed1, seed2;
extern void setall(long iseed1, long iseed2);

#ifndef TEST_RUNNER_MODE
void setUp(void) {
    /* Initialize the RNG so ranf() is usable from KS tests later */
    seed1 = 12345;
    seed2 = 67890;
    setall(seed1, seed2);

    /* Reset shape state for every test */
    for (int i = 0; i < MAXPOPS; i++) {
        popShape[i].type = SHAPE_CONSTANT;
        popShape[i].anchor_value = 1.0;
        popShape[i].rate_param = 0.0;
        popShape[i].anchor_time = 0.0;
        for (int j = 0; j < MAXPOPS; j++) {
            migShape[i][j].type = SHAPE_CONSTANT;
            migShape[i][j].anchor_value = 0.0;
            migShape[i][j].rate_param = 0.0;
            migShape[i][j].anchor_time = 0.0;
        }
    }
}

void tearDown(void) { }
#endif

void test_placeholder(void) {
    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, 0);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_placeholder);
    return UNITY_END();
}
#endif
```

- [ ] **Step 2: Add `test_shapes` build target to `Makefile`**

After the `test_event:` target (line 248-250), add:

```makefile
test_shapes: $(TEST_DIR)/test_shapes.c $(SRC_CORE)/shapes.c $(SRC_CORE)/shapes.h $(SRC_CORE)/discoal.h $(TEST_DIR)/test_globals.c
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUSE_XOSHIRO256PP -DUNITY_INCLUDE_DOUBLE -o build/test_shapes $(TEST_DIR)/test_shapes.c $(SRC_CORE)/shapes.c $(TEST_DIR)/test_globals.c $(SRC_RNG)/xoshiro256pp_compat.c $(UNITY_SOURCES) -I$(UNITY_DIR) -lm -fcommon
```

- [ ] **Step 3: Add `test_shapes` to `run_tests` aggregator**

Find `run_tests:` recipe (around line 290). Add `test_shapes` to the dependency list and `./build/test_shapes` to the run sequence.

- [ ] **Step 4: Build and run the empty runner**

Run: `make test_shapes && ./build/test_shapes`
Expected: `1 Tests 0 Failures 0 Ignored OK` (or equivalent Unity success output).

- [ ] **Step 5: Commit**

```bash
git add test/unit/test_shapes.c Makefile
git commit -m "Add test_shapes runner scaffold

Empty Unity test target wired into run_tests. Resets popShape
and migShape arrays in setUp so each test starts from a clean
slate."
```

---

## Task 3: `sizeAt` for `SHAPE_CONSTANT`

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

- [ ] **Step 1: Write the failing test**

Replace `test_placeholder` and its `RUN_TEST` line with:

```c
void test_sizeAt_constant_returns_anchor(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 2.5;
    popShape[0].anchor_time = 10.0;
    /* For CONSTANT, t should be irrelevant */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.5, sizeAt(0, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.5, sizeAt(0, 50.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.5, sizeAt(0, 1e6));
}

void test_sizeAt_constant_per_population(void) {
    popShape[1].type = SHAPE_CONSTANT;
    popShape[1].anchor_value = 0.5;
    popShape[2].type = SHAPE_CONSTANT;
    popShape[2].anchor_value = 4.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, sizeAt(1, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 4.0, sizeAt(2, 0.0));
}
```

In `main`:

```c
RUN_TEST(test_sizeAt_constant_returns_anchor);
RUN_TEST(test_sizeAt_constant_per_population);
```

- [ ] **Step 2: Run to verify failure**

Run: `make test_shapes && ./build/test_shapes`
Expected: both tests FAIL — `sizeAt` returns 0.0 from the stub.

- [ ] **Step 3: Implement `sizeAt` for `SHAPE_CONSTANT`**

In `src/core/shapes.c`, replace the `sizeAt` stub:

```c
double sizeAt(int popID, double t) {
    Shape *s = &popShape[popID];
    switch (s->type) {
        case SHAPE_CONSTANT:
            return s->anchor_value;
        default:
            return 0.0;  /* other shapes implemented in subsequent tasks */
    }
}
```

- [ ] **Step 4: Verify pass**

Run: `make test_shapes && ./build/test_shapes`
Expected: both new tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement sizeAt for SHAPE_CONSTANT

Returns anchor_value regardless of t. Tests verify per-population
isolation."
```

---

## Task 4: `sizeAt` for `SHAPE_EXPONENTIAL`

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

The exponential shape: $N(t) = N_0 \cdot e^{-\alpha (t - t_0)}$ where $N_0$ = `anchor_value`, $\alpha$ = `rate_param`, $t_0$ = `anchor_time`. (Forward-time growth rate $\alpha > 0$ ⇒ backward-time decline.)

- [ ] **Step 1: Write the failing test**

Add to `test/unit/test_shapes.c`:

```c
void test_sizeAt_exponential_at_anchor_time(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 2.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 10.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0, sizeAt(0, 10.0));
}

void test_sizeAt_exponential_decays_backward(void) {
    /* alpha > 0 means N decreases as t grows past anchor_time */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    /* t=2: N = 1 * exp(-0.5*2) = exp(-1) */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, exp(-1.0), sizeAt(0, 2.0));
    /* t=4: N = exp(-2) */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, exp(-2.0), sizeAt(0, 4.0));
}

void test_sizeAt_exponential_grows_forward(void) {
    /* For t < anchor_time, N is larger */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 5.0;
    /* t=3: N = 1 * exp(-0.5 * (3-5)) = exp(1) */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, exp(1.0), sizeAt(0, 3.0));
}

void test_sizeAt_exponential_zero_alpha_equals_constant(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 3.0;
    popShape[0].rate_param = 0.0;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 3.0, sizeAt(0, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 3.0, sizeAt(0, 1e6));
}
```

In `main`:

```c
RUN_TEST(test_sizeAt_exponential_at_anchor_time);
RUN_TEST(test_sizeAt_exponential_decays_backward);
RUN_TEST(test_sizeAt_exponential_grows_forward);
RUN_TEST(test_sizeAt_exponential_zero_alpha_equals_constant);
```

- [ ] **Step 2: Verify failure**

Run: `make test_shapes && ./build/test_shapes`
Expected: 4 new tests FAIL.

- [ ] **Step 3: Implement**

Add a `case` to `sizeAt` in `src/core/shapes.c`:

```c
case SHAPE_EXPONENTIAL:
    return s->anchor_value * exp(-s->rate_param * (t - s->anchor_time));
```

- [ ] **Step 4: Verify pass**

Run: `make test_shapes && ./build/test_shapes`
Expected: all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement sizeAt for SHAPE_EXPONENTIAL

N(t) = anchor_value * exp(-rate_param * (t - anchor_time)).
Forward growth rate alpha > 0 maps to backward decline in
coalescent time. Tests cover anchor, forward, backward, and
the alpha=0 degenerate case."
```

---

## Task 5: `sizeAt` for `SHAPE_LINEAR`

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

Linear shape (msprime forward-time convention): $N(t) = N_0 - \gamma (t - t_0)$. $\gamma$ is the forward-time growth rate. $\gamma > 0$ means $N$ *declines* backward-in-time (past was smaller); $\gamma < 0$ means $N$ grows backward.

- [ ] **Step 1: Failing test**

```c
void test_sizeAt_linear_at_anchor(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.5;
    popShape[0].rate_param = 0.25;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.5, sizeAt(0, 0.0));
}

void test_sizeAt_linear_declines_backward_with_positive_gamma(void) {
    /* gamma > 0 = forward growth = backward decline (past smaller) */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 5.0;
    popShape[0].rate_param = 2.0;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, sizeAt(0, 2.0));   /* 5 - 2*2 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 3.0, sizeAt(0, 1.0));   /* 5 - 2*1 */
}

void test_sizeAt_linear_grows_backward_with_negative_gamma(void) {
    /* gamma < 0 = forward decline = backward growth (past larger) */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = -1.0;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0, sizeAt(0, 1.0));   /* 1 - (-1)*1 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 5.0, sizeAt(0, 4.0));   /* 1 - (-1)*4 */
}
```

Register all three in `main`.

- [ ] **Step 2: Verify failure**: `make test_shapes && ./build/test_shapes` — 3 tests FAIL.

- [ ] **Step 3: Implement**

Add to `sizeAt` BEFORE `default:`:

```c
case SHAPE_LINEAR:
    return s->anchor_value - s->rate_param * (t - s->anchor_time);
```

- [ ] **Step 4: Verify pass**: tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement sizeAt for SHAPE_LINEAR

N(t) = anchor_value - rate_param * (t - anchor_time). rate_param is
the forward-time growth rate (msprime convention); positive value
means past was smaller."
```

---

## Task 6: `migAt` for all three shapes

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

`migAt` is structurally identical to `sizeAt` but indexes into `migShape[srcPopID][dstPopID]`. Implement and test all three shapes in one task since the logic is parallel.

- [ ] **Step 1: Failing tests**

```c
void test_migAt_constant(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.5;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, migAt(0, 1, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, migAt(0, 1, 100.0));
}

void test_migAt_exponential(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 1.0;
    migShape[0][1].rate_param = 0.5;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, exp(-1.0), migAt(0, 1, 2.0));
}

void test_migAt_linear_declines_backward_with_positive_delta(void) {
    /* delta > 0 = forward growth = backward decline. m(t) = 0.5 - 0.05*t. */
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.5;
    migShape[0][1].rate_param = 0.05;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.3, migAt(0, 1, 4.0));
}

void test_migAt_linear_grows_backward_with_negative_delta(void) {
    /* delta < 0 = forward decline = backward growth. m(t) = 0.1 + 0.05*t. */
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.1;
    migShape[0][1].rate_param = -0.05;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.6, migAt(0, 1, 10.0));
}

void test_migAt_pair_isolation(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.3;
    migShape[1][0].type = SHAPE_CONSTANT;
    migShape[1][0].anchor_value = 0.7;
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.3, migAt(0, 1, 0.0));
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.7, migAt(1, 0, 0.0));
}
```

Register in `main`.

- [ ] **Step 2: Verify failure** — 5 tests FAIL.

- [ ] **Step 3: Implement `migAt` in `src/core/shapes.c`**

```c
double migAt(int srcPopID, int dstPopID, double t) {
    Shape *s = &migShape[srcPopID][dstPopID];
    switch (s->type) {
        case SHAPE_CONSTANT:
            return s->anchor_value;
        case SHAPE_EXPONENTIAL:
            return s->anchor_value * exp(-s->rate_param * (t - s->anchor_time));
        case SHAPE_LINEAR:
            return s->anchor_value - s->rate_param * (t - s->anchor_time);
        default:
            return 0.0;
    }
}
```

- [ ] **Step 4: Verify pass**: tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement migAt for all shapes

Mirrors sizeAt structure but indexes migShape[src][dst]. Same
forward-time msprime convention: positive rate_param means past
held a smaller migration rate. Tested per shape and for src/dst
pair isolation."
```

---

## Task 7: `integratedHazardSize` for `SHAPE_CONSTANT`

The integrated hazard for $k$ lineages is

$$H(T) = \int_0^T \frac{\binom{k}{2}}{N(t_0 + s)}\,ds.$$

For `CONSTANT`: $H(T) = \binom{k}{2} T / N_0$.

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

- [ ] **Step 1: Failing test**

```c
void test_integratedHazardSize_constant(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 2.0;
    /* k=2, T=4: H = 1*4/2 = 2 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0, integratedHazardSize(0, 0.0, 4.0, 2));
    /* k=4, T=10: H = 6*10/2 = 30 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 30.0, integratedHazardSize(0, 0.0, 10.0, 4));
}
```

- [ ] **Step 2: Verify failure** — FAIL.

- [ ] **Step 3: Implement**

Replace the `integratedHazardSize` stub in `src/core/shapes.c`:

```c
double integratedHazardSize(int popID, double t0, double T, int k) {
    Shape *s = &popShape[popID];
    double pairs = (double)k * (k - 1) / 2.0;
    if (pairs == 0.0) return 0.0;
    switch (s->type) {
        case SHAPE_CONSTANT:
            return pairs * T / s->anchor_value;
        default:
            return 0.0;
    }
}
```

- [ ] **Step 4: Verify pass** — PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement integratedHazardSize for SHAPE_CONSTANT

H(T) = (k choose 2) * T / N_0. Returns 0 when k < 2."
```

---

## Task 8: `integratedHazardSize` for `SHAPE_EXPONENTIAL`

For exponential $N(t) = N_0 e^{-\alpha s}$ (with $s = t - t_0$, $N_0$ = `anchor_value` *evaluated at* $t_0$):

$$H(T) = \int_0^T \binom{k}{2}\,\frac{e^{\alpha s}}{N_0}\,ds = \frac{\binom{k}{2}}{N_0\,\alpha}\,(e^{\alpha T} - 1).$$

Special case: $\alpha = 0$ gives $H(T) = \binom{k}{2} T / N_0$ (matches CONSTANT).

**Critical detail**: the integral must use $N_0$ = `sizeAt(popID, t0)`, not the raw `anchor_value`, because `t0` may not be the shape's anchor time.

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

- [ ] **Step 1: Failing tests**

```c
void test_integratedHazardSize_exponential_at_anchor(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 2.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    /* k=2, T=4, alpha=0.5, N_0 at t=0 = 2.0 */
    /* H = 1/(2*0.5) * (exp(2) - 1) = (exp(2)-1) */
    double expected = (exp(2.0) - 1.0);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected, integratedHazardSize(0, 0.0, 4.0, 2));
}

void test_integratedHazardSize_exponential_offset_t0(void) {
    /* If t0 != anchor_time, the relevant N_0 is sizeAt(t0) */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.2;
    popShape[0].anchor_time = 0.0;
    /* At t0=5, N(t0) = exp(-1). Then H over T=2 with that as anchor:
     * H = 1/(N(t0)*0.2) * (exp(0.2*2) - 1) = 1/(exp(-1)*0.2) * (exp(0.4)-1) */
    double N_t0 = exp(-1.0);
    double expected = 1.0 * (exp(0.4) - 1.0) / (N_t0 * 0.2);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected, integratedHazardSize(0, 5.0, 2.0, 2));
}

void test_integratedHazardSize_exponential_alpha_zero(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 3.0;
    popShape[0].rate_param = 0.0;
    popShape[0].anchor_time = 0.0;
    /* alpha=0 should match constant: H = 1*5/3 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 5.0/3.0, integratedHazardSize(0, 0.0, 5.0, 2));
}

void test_integratedHazardSize_exponential_quadrature_match(void) {
    /* High-resolution numerical quadrature should match the closed form */
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.5;
    popShape[0].rate_param = 0.7;
    popShape[0].anchor_time = 2.0;
    double t0 = 3.5;
    double T = 4.0;
    int k = 5;
    double pairs = k*(k-1)/2.0;
    int N = 16384;
    double dt = T / N;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double s_lo = i * dt;
        double s_hi = (i+1) * dt;
        double s_mid = 0.5 * (s_lo + s_hi);
        sum += pairs / sizeAt(0, t0 + s_mid) * dt;  /* midpoint rule */
    }
    double closed = integratedHazardSize(0, t0, T, k);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, sum, closed);
}
```

- [ ] **Step 2: Verify failure** — 4 tests FAIL.

- [ ] **Step 3: Implement**

In `src/core/shapes.c`, expand `integratedHazardSize`:

```c
double integratedHazardSize(int popID, double t0, double T, int k) {
    Shape *s = &popShape[popID];
    double pairs = (double)k * (k - 1) / 2.0;
    if (pairs == 0.0) return 0.0;
    double N0 = sizeAt(popID, t0);
    switch (s->type) {
        case SHAPE_CONSTANT:
            return pairs * T / N0;
        case SHAPE_EXPONENTIAL: {
            double a = s->rate_param;
            if (a == 0.0) return pairs * T / N0;
            return pairs * (exp(a * T) - 1.0) / (N0 * a);
        }
        default:
            return 0.0;
    }
}
```

- [ ] **Step 4: Verify pass** — all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement integratedHazardSize for SHAPE_EXPONENTIAL

Closed form (k choose 2)*(exp(alpha*T) - 1) / (N(t0)*alpha)
verified to 1e-6 against midpoint-rule quadrature with 16k
subintervals. Includes alpha=0 degenerate case and offset t0."
```

---

## Task 9: `integratedHazardSize` for `SHAPE_LINEAR`

For linear $N(s) = N_0 - \gamma s$ (forward-time convention; $\gamma > 0$ means past was smaller):

$$H(T) = \int_0^T \frac{\binom{k}{2}}{N_0 - \gamma s}\,ds = \frac{\binom{k}{2}}{\gamma}\,\log\!\left(\frac{N_0}{N_0 - \gamma T}\right).$$

Special cases:
- $\gamma = 0$ ⇒ matches CONSTANT.
- $N_0 - \gamma T \le 0$ (size hits zero before $T$ when $\gamma > 0$) ⇒ hazard diverges; return $+\infty$ (the caller in `drawWaitingTime` clamps $T$ at the crossing).

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

- [ ] **Step 1: Failing tests**

```c
void test_integratedHazardSize_linear_basic(void) {
    /* gamma = -0.5 (forward decline => backward growth):
     * N(s) = 1 - (-0.5)*s = 1 + 0.5 s, k=2, T=4
     * H = 1/(-0.5) * log(1/(1+0.5*4)) = -2 * log(1/3) = 2*log(3) */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = -0.5;
    popShape[0].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 2.0 * log(3.0), integratedHazardSize(0, 0.0, 4.0, 2));
}

void test_integratedHazardSize_linear_zero_gamma(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 2.0;
    popShape[0].rate_param = 0.0;
    popShape[0].anchor_time = 0.0;
    /* matches constant */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 1.0 * 5.0 / 2.0, integratedHazardSize(0, 0.0, 5.0, 2));
}

void test_integratedHazardSize_linear_diverges_at_zero_crossing(void) {
    /* gamma = 0.5 (forward growth => backward decline):
     * N(t) = 1 - 0.5 t hits zero at t=2. Integration to T=2.5 should diverge. */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    double H = integratedHazardSize(0, 0.0, 2.5, 2);
    TEST_ASSERT_TRUE(isinf(H) || H > 1e15);
}

void test_integratedHazardSize_linear_quadrature_match(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 2.0;
    popShape[0].rate_param = 0.3;
    popShape[0].anchor_time = 1.0;
    double t0 = 1.5;
    double T = 6.0;
    int k = 4;
    double pairs = k*(k-1)/2.0;
    int N = 16384;
    double dt = T / N;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double s_mid = (i + 0.5) * dt;
        sum += pairs / sizeAt(0, t0 + s_mid) * dt;
    }
    double closed = integratedHazardSize(0, t0, T, k);
    /* 1e-5 not 1e-6 -- midpoint-rule error is ~4e-6 here because the
     * integrand 1/N(s) ranges from 0.54 to 20 (N near the far end is 0.05).
     * The closed form is exact; this test still catches any algebraic
     * mistake in the implementation. */
    TEST_ASSERT_DOUBLE_WITHIN(1e-5, sum, closed);
}
```

- [ ] **Step 2: Verify failure** — 4 tests FAIL.

- [ ] **Step 3: Implement**

Add to `integratedHazardSize`. Use the `log1p` form for numerical
stability — `log(N0/end) = -log1p(-frac)` where `frac = g*T/N0`,
mathematically identical and more accurate when `frac` is small.

```c
case SHAPE_LINEAR: {
    double g = s->rate_param;
    if (g == 0.0) return pairs * T / N0;
    double frac = g * T / N0;
    if (frac >= 1.0) return INFINITY;  /* end <= 0 case */
    return -pairs * log1p(-frac) / g;
}
```

- [ ] **Step 4: Verify pass** — all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement integratedHazardSize for SHAPE_LINEAR

H = (k choose 2) * log(N0 / (N0 - gamma*T)) / gamma. Returns
+infinity when N(t) crosses zero within T (gamma > 0) to
signal forced coalescence to the waiting-time draw."
```

---

## Task 10: `integratedHazardMig` for all three shapes

For migration, the integrand is $\lambda_M(s) = k_i\,m_{ij}(s)$, *not* $\binom{k}{2}/N$. All shapes follow the msprime forward-time convention ($\beta$, $\delta$ are forward-time growth rates):

- CONST: $m(s) = m_0$, $H(T) = k\,m_0\,T$.
- EXP: $m(s) = m_0\,e^{-\beta s}$, $H(T) = k\,m_0\,(1 - e^{-\beta T}) / \beta$.
- LIN: $m(s) = m_0 - \delta s$, $H(T) = k\,(m_0\,T - \tfrac{1}{2}\,\delta\,T^2)$.

(The migration linear case integrates as a simple polynomial because the integrand is affine in $s$, not the reciprocal-of-affine form size has.)

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

- [ ] **Step 1: Failing tests**

```c
void test_integratedHazardMig_constant(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.5;
    /* k=3, T=4: H = 3*0.5*4 = 6 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 6.0, integratedHazardMig(0, 1, 0.0, 4.0, 3));
}

void test_integratedHazardMig_exponential(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 1.0;
    migShape[0][1].rate_param = 0.5;
    migShape[0][1].anchor_time = 0.0;
    /* k=2, T=4: m(s) = exp(-0.5 s), H = 2 * (1 - exp(-2))/0.5 */
    double expected = 2.0 * (1.0 - exp(-2.0)) / 0.5;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, expected, integratedHazardMig(0, 1, 0.0, 4.0, 2));
}

void test_integratedHazardMig_linear(void) {
    /* delta = -0.05 (forward decline => backward growth):
     * m(s) = 0.1 - (-0.05)*s = 0.1 + 0.05 s, k=2, T=4
     * H = 2*(0.1*4 - (-0.05)*16/2) = 2*(0.4 + 0.4) = 1.6 */
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.1;
    migShape[0][1].rate_param = -0.05;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 1.6, integratedHazardMig(0, 1, 0.0, 4.0, 2));
}

void test_integratedHazardMig_quadrature_match(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 0.3;
    migShape[0][1].rate_param = 0.7;
    migShape[0][1].anchor_time = 1.0;
    double t0 = 1.5;
    double T = 3.0;
    int k = 4;
    int N = 16384;
    double dt = T / N;
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        double s_mid = (i + 0.5) * dt;
        sum += k * migAt(0, 1, t0 + s_mid) * dt;
    }
    double closed = integratedHazardMig(0, 1, t0, T, k);
    TEST_ASSERT_DOUBLE_WITHIN(1e-6, sum, closed);
}
```

- [ ] **Step 2: Verify failure** — 4 tests FAIL.

- [ ] **Step 3: Implement**

In `src/core/shapes.c`, replace the `integratedHazardMig` stub:

```c
double integratedHazardMig(int srcPopID, int dstPopID, double t0, double T, int k) {
    Shape *s = &migShape[srcPopID][dstPopID];
    double m0 = migAt(srcPopID, dstPopID, t0);
    if (k <= 0 || m0 == 0.0) return 0.0;
    switch (s->type) {
        case SHAPE_CONSTANT:
            return k * m0 * T;
        case SHAPE_EXPONENTIAL: {
            double b = s->rate_param;
            if (b == 0.0) return k * m0 * T;
            return k * m0 * (1.0 - exp(-b * T)) / b;
        }
        case SHAPE_LINEAR: {
            double d = s->rate_param;
            return k * (m0 * T - 0.5 * d * T * T);
        }
        default:
            return 0.0;
    }
}
```

- [ ] **Step 4: Verify pass** — all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement integratedHazardMig for all shapes

Constant: k*m0*T. Exponential: k*m0*(1-exp(-beta*T))/beta.
Linear: k*(m0*T - delta*T^2/2). Forward-time convention
(msprime); quadrature match to 1e-6."
```

---

## Task 11: `drawWaitingTimeSize` — invert the closed forms

Given $\xi$, solve $H(T) = \xi$ for $T$:

- CONST: $T = \xi N_0 / \binom{k}{2}$
- EXP: $T = \log\!\big(1 + N_0 \alpha \xi / \binom{k}{2}\big) / \alpha$ if $\alpha \ne 0$, else CONST formula
- LIN: $T = (N_0/\gamma)(1 - \exp(-\gamma \xi / \binom{k}{2}))$ if $\gamma \ne 0$, else CONST formula. **Plus zero-crossing protection**: if $\gamma > 0$ (forward growth ⇒ backward decline) the crossing $T^* = N_0/\gamma$ is finite. The closed form returns $T < T^*$ for any finite $\xi$ in exact arithmetic, but we clamp via `nextafter` to defend against floating-point ties.

Returns $-1$ if $k < 2$ (no possible coalescence).

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

- [ ] **Step 1: Failing tests**

```c
void test_drawWaitingTimeSize_constant_inverts(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 2.0;
    /* xi = 3, k=2: T = 3*2/1 = 6 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 6.0, drawWaitingTimeSize(0, 0.0, 3.0, 2));
}

void test_drawWaitingTimeSize_exponential_inverts(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    double xi = 2.0;
    int k = 2;
    /* T such that H(T)=xi: T = log(1 + N0*alpha*xi/pairs)/alpha = log(1+1)/0.5 = 2*log(2) */
    double T = drawWaitingTimeSize(0, 0.0, xi, k);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, 2.0 * log(2.0), T);
    /* Verify round trip: H(T) ≈ xi */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, xi, integratedHazardSize(0, 0.0, T, k));
}

void test_drawWaitingTimeSize_linear_inverts(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    double xi = 1.0;
    int k = 2;
    double T = drawWaitingTimeSize(0, 0.0, xi, k);
    /* Round trip */
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, xi, integratedHazardSize(0, 0.0, T, k));
}

void test_drawWaitingTimeSize_returns_negative_for_k_lt_2(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 1.0;
    TEST_ASSERT_TRUE(drawWaitingTimeSize(0, 0.0, 1.0, 1) < 0.0);
    TEST_ASSERT_TRUE(drawWaitingTimeSize(0, 0.0, 1.0, 0) < 0.0);
}

void test_drawWaitingTimeSize_linear_zero_crossing(void) {
    /* gamma = 0.5 (forward growth => backward decline):
     * N(t) = 1 - 0.5 t hits zero at t=2. H(2) = +inf, so any finite xi maps to T<2. */
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    int k = 2;
    /* Try a battery of xi values; T must always be < 2. */
    for (double xi = 0.1; xi < 100.0; xi *= 1.5) {
        double T = drawWaitingTimeSize(0, 0.0, xi, k);
        TEST_ASSERT_TRUE(T > 0.0);
        TEST_ASSERT_TRUE(T < 2.0);
    }
}
```

- [ ] **Step 2: Verify failure** — 5 tests FAIL.

- [ ] **Step 3: Implement**

```c
double drawWaitingTimeSize(int popID, double t0, double xi, int k) {
    if (k < 2) return -1.0;
    Shape *s = &popShape[popID];
    double pairs = (double)k * (k - 1) / 2.0;
    double N0 = sizeAt(popID, t0);
    if (N0 <= 0.0) return -1.0;
    switch (s->type) {
        case SHAPE_CONSTANT:
            return xi * N0 / pairs;
        case SHAPE_EXPONENTIAL: {
            double a = s->rate_param;
            if (a == 0.0) return xi * N0 / pairs;
            return log1p(N0 * a * xi / pairs) / a;
        }
        case SHAPE_LINEAR: {
            double g = s->rate_param;
            if (g == 0.0) return xi * N0 / pairs;
            double T = (N0 / g) * (1.0 - exp(-g * xi / pairs));
            if (g > 0.0) {
                double T_cross = N0 / g;
                if (T >= T_cross) T = nextafter(T_cross, 0.0);
            }
            return T;
        }
        default:
            return -1.0;
    }
}
```

- [ ] **Step 4: Verify pass** — all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement drawWaitingTimeSize for all shapes

Closed-form inversion of integrated hazard. log1p used for
numerical stability near alpha*xi -> 0. Linear shape (msprime
forward-time convention: gamma > 0 means past was smaller)
clips T at the zero-crossing boundary using nextafter to
ensure strict less-than."
```

---

## Task 12: `drawWaitingTimeMig` — invert migration closed forms

Same idea for migration (msprime forward-time convention: $\beta$, $\delta$ are forward-time growth rates of migration):
- CONST: $T = \xi / (k\,m_0)$
- EXP: $m(s) = m_0 e^{-\beta s}$, $T = -\log(1 - \beta \xi / (k m_0))/\beta$ when $\beta \ne 0$. (If $\beta\,\xi/(k m_0) \ge 1$, the integrated hazard over $[0,\infty)$ is finite and $\xi$ exceeds it — return $-1$.)
- LIN: integrand $k(m_0 - \delta s)$; integral $k(m_0 T - \delta T^2/2) = \xi$. Quadratic root: $T = (m_0 - \sqrt{m_0^2 - 2 \delta \xi/k}) / \delta$ when $\delta \ne 0$. Negative discriminant (when $\delta > 0$ and $\xi > k m_0^2/(2\delta)$, i.e. xi exceeds the cap reached at the migration zero-crossing) ⇒ unreachable, return $-1$.

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c`

- [ ] **Step 1: Failing tests**

```c
void test_drawWaitingTimeMig_constant_inverts(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.5;
    /* xi=2, k=4: T = 2/(4*0.5) = 1 */
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, drawWaitingTimeMig(0, 1, 0.0, 2.0, 4));
}

void test_drawWaitingTimeMig_exponential_round_trip(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 0.5;
    migShape[0][1].rate_param = 0.3;
    migShape[0][1].anchor_time = 0.0;
    double xi = 1.0;
    int k = 3;
    double T = drawWaitingTimeMig(0, 1, 0.0, xi, k);
    TEST_ASSERT_TRUE(T > 0.0);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, xi, integratedHazardMig(0, 1, 0.0, T, k));
}

void test_drawWaitingTimeMig_exponential_unreachable_xi(void) {
    /* m(s) = 0.5 exp(-0.3 s). Total integrated hazard over [0,inf) for k=3:
     * 3 * 0.5 / 0.3 = 5. Any xi >= 5 should be unreachable. */
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 0.5;
    migShape[0][1].rate_param = 0.3;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_TRUE(drawWaitingTimeMig(0, 1, 0.0, 5.5, 3) < 0.0);
    TEST_ASSERT_TRUE(drawWaitingTimeMig(0, 1, 0.0, 100.0, 3) < 0.0);
}

void test_drawWaitingTimeMig_linear_round_trip(void) {
    /* delta = -0.05 (forward decline => backward growth) so m(s) = 0.1+0.05s
     * grows monotonically and the round-trip is well-defined for any xi. */
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.1;
    migShape[0][1].rate_param = -0.05;
    migShape[0][1].anchor_time = 0.0;
    double xi = 0.8;
    int k = 4;
    double T = drawWaitingTimeMig(0, 1, 0.0, xi, k);
    TEST_ASSERT_TRUE(T > 0.0);
    TEST_ASSERT_DOUBLE_WITHIN(1e-10, xi, integratedHazardMig(0, 1, 0.0, T, k));
}

void test_drawWaitingTimeMig_zero_rate_returns_negative(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.0;
    TEST_ASSERT_TRUE(drawWaitingTimeMig(0, 1, 0.0, 1.0, 4) < 0.0);
}
```

- [ ] **Step 2: Verify failure** — 5 tests FAIL.

- [ ] **Step 3: Implement**

```c
double drawWaitingTimeMig(int srcPopID, int dstPopID, double t0, double xi, int k) {
    if (k <= 0) return -1.0;
    Shape *s = &migShape[srcPopID][dstPopID];
    double m0 = migAt(srcPopID, dstPopID, t0);
    if (m0 <= 0.0) return -1.0;
    double km0 = k * m0;
    switch (s->type) {
        case SHAPE_CONSTANT:
            return xi / km0;
        case SHAPE_EXPONENTIAL: {
            double b = s->rate_param;
            if (b == 0.0) return xi / km0;
            double arg = b * xi / km0;
            if (arg >= 1.0) return -1.0;  /* unreachable */
            return -log1p(-arg) / b;
        }
        case SHAPE_LINEAR: {
            double d = s->rate_param;
            if (d == 0.0) return xi / km0;
            double disc = m0 * m0 - 2.0 * d * xi / k;
            if (disc < 0.0) return -1.0;
            double T = (m0 - sqrt(disc)) / d;
            if (T < 0.0) return -1.0;
            return T;
        }
        default:
            return -1.0;
    }
}
```

- [ ] **Step 4: Verify pass** — all PASS.

- [ ] **Step 5: Commit**

```bash
git add src/core/shapes.c test/unit/test_shapes.c
git commit -m "Implement drawWaitingTimeMig for all shapes

Constant: T = xi/(k*m0). Exponential: T = -log1p(-b*xi/(k*m0))/b
with unreachable-xi detection. Linear: smaller positive root
of the quadratic k*m0*T - k*delta*T^2/2 = xi. Forward-time
convention (msprime)."
```

---

## Task 13: KS distribution test for `drawWaitingTimeSize`

Sample $\xi \sim \text{Exp}(1)$ many times, compute $T_i = \texttt{drawWaitingTimeSize}(\xi_i)$, and verify the empirical CDF of $T_i$ matches the theoretical $F_T(t) = 1 - \exp(-H(t))$ via a Kolmogorov-Smirnov statistic.

For ease of testing, compute the K-S statistic in C against a discretized theoretical CDF, with the threshold set generously ($D < 0.05$ for $n = 10^4$ — KS critical at 5% is $\approx 1.36/\sqrt{n} = 0.0136$, but we want a forgiving threshold that catches gross errors without false-positiving on rare samples).

**Files:**
- Modify: `test/unit/test_shapes.c`
- Modify: `src/core/shapes.c` (no changes; uses existing functions)

- [ ] **Step 1: Failing test**

Add to `test/unit/test_shapes.c`. The `setUp` in Task 2 already calls `setall(12345, 67890)` so `ranf()` is ready to use.

```c
extern double ranf(void);  /* from xoshiro256pp_compat */

static int compare_doubles(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

static double ks_statistic_size(int popID, int k, int n) {
    double *T = malloc(sizeof(double) * n);
    for (int i = 0; i < n; i++) {
        double xi = -log(ranf());
        T[i] = drawWaitingTimeSize(popID, 0.0, xi, k);
    }
    qsort(T, n, sizeof(double), compare_doubles);
    double max_d = 0.0;
    for (int i = 0; i < n; i++) {
        double F_emp = (double)(i + 1) / n;
        double F_theory = 1.0 - exp(-integratedHazardSize(popID, 0.0, T[i], k));
        double d = fabs(F_emp - F_theory);
        if (d > max_d) max_d = d;
    }
    free(T);
    return max_d;
}

void test_drawWaitingTimeSize_distribution_constant(void) {
    popShape[0].type = SHAPE_CONSTANT;
    popShape[0].anchor_value = 1.0;
    int n = 10000;
    double D = ks_statistic_size(0, 2, n);
    /* For n=10000, KS critical at 5% is ~0.0136. Use 0.05 as forgiving threshold. */
    TEST_ASSERT_TRUE(D < 0.05);
}

void test_drawWaitingTimeSize_distribution_exponential(void) {
    popShape[0].type = SHAPE_EXPONENTIAL;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.5;
    popShape[0].anchor_time = 0.0;
    int n = 10000;
    double D = ks_statistic_size(0, 4, n);
    TEST_ASSERT_TRUE(D < 0.05);
}

void test_drawWaitingTimeSize_distribution_linear(void) {
    popShape[0].type = SHAPE_LINEAR;
    popShape[0].anchor_value = 1.0;
    popShape[0].rate_param = 0.3;
    popShape[0].anchor_time = 0.0;
    int n = 10000;
    double D = ks_statistic_size(0, 3, n);
    TEST_ASSERT_TRUE(D < 0.05);
}
```

Note: the `ranf` symbol is exported by `xoshiro256pp_compat.c`; the test_shapes Makefile recipe already links it.

- [ ] **Step 2: Verify pass** (these tests should already pass given correct implementations)

Run: `make test_shapes && ./build/test_shapes`
Expected: 3 new tests PASS. If any fails, investigate the closed-form math — the round-trip tests in earlier tasks pass, so failure here likely indicates an RNG seeding issue or a transcription error.

- [ ] **Step 3: If pass, commit**

```bash
git add test/unit/test_shapes.c
git commit -m "Add KS distribution tests for drawWaitingTimeSize

10k samples compared to theoretical CDF 1 - exp(-H(t)) for
each shape. Threshold D < 0.05 is loose enough to avoid
false positives but tight enough to catch any inversion error."
```

---

## Task 14: KS distribution test for `drawWaitingTimeMig`

Mirror Task 13 for migration.

**Files:**
- Modify: `test/unit/test_shapes.c`

- [ ] **Step 1: Add test (same shape as Task 13)**

```c
static double ks_statistic_mig(int src, int dst, int k, int n) {
    double *T = malloc(sizeof(double) * n);
    int valid = 0;
    for (int i = 0; i < n; i++) {
        double xi = -log(ranf());
        double t = drawWaitingTimeMig(src, dst, 0.0, xi, k);
        if (t > 0.0) T[valid++] = t;
    }
    if (valid < n / 2) { free(T); return 1.0; }  /* too many unreachable; fail */
    qsort(T, valid, sizeof(double), compare_doubles);
    double max_d = 0.0;
    for (int i = 0; i < valid; i++) {
        double F_emp = (double)(i + 1) / valid;
        double F_theory = 1.0 - exp(-integratedHazardMig(src, dst, 0.0, T[i], k));
        /* Note: F_theory here is conditional-on-finite, which matches valid sample */
        double d = fabs(F_emp - F_theory);
        if (d > max_d) max_d = d;
    }
    free(T);
    return max_d;
}

void test_drawWaitingTimeMig_distribution_constant(void) {
    migShape[0][1].type = SHAPE_CONSTANT;
    migShape[0][1].anchor_value = 0.3;
    TEST_ASSERT_TRUE(ks_statistic_mig(0, 1, 4, 10000) < 0.05);
}

void test_drawWaitingTimeMig_distribution_exponential(void) {
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    migShape[0][1].anchor_value = 1.0;
    migShape[0][1].rate_param = 0.5;
    migShape[0][1].anchor_time = 0.0;
    /* H(infty) = k*m0/beta = 4*1/0.5 = 8. Most xi will resolve. */
    TEST_ASSERT_TRUE(ks_statistic_mig(0, 1, 4, 10000) < 0.05);
}

void test_drawWaitingTimeMig_distribution_linear(void) {
    migShape[0][1].type = SHAPE_LINEAR;
    migShape[0][1].anchor_value = 0.5;
    migShape[0][1].rate_param = 0.1;
    migShape[0][1].anchor_time = 0.0;
    TEST_ASSERT_TRUE(ks_statistic_mig(0, 1, 3, 10000) < 0.05);
}
```

- [ ] **Step 2: Run and verify pass**

Run: `make test_shapes && ./build/test_shapes`

- [ ] **Step 3: Commit**

```bash
git add test/unit/test_shapes.c
git commit -m "Add KS distribution tests for drawWaitingTimeMig"
```

---

## Task 15: Linear-shape edge cases — explicit zero-crossing test

The earlier `test_drawWaitingTimeSize_linear_zero_crossing` ensures $T < T^*$ for valid $\xi$. This task adds a stress test exercising the integrator on the boundary.

**Files:**
- Modify: `test/unit/test_shapes.c`

- [ ] **Step 1: Add test**

```c
void test_drawWaitingTimeSize_linear_zero_crossing_stress(void) {
    /* Battery of (gamma, N0) configurations that drive N -> 0 going backward.
     * Under msprime forward-time convention, gamma > 0 (forward growth) means
     * past was smaller, hitting zero at T_cross = N0/gamma. */
    struct { double N0; double gamma; } configs[] = {
        {1.0, 0.1}, {1.0, 1.0}, {1.0, 10.0},
        {0.01, 0.001}, {100.0, 50.0}
    };
    for (int c = 0; c < (int)(sizeof(configs)/sizeof(configs[0])); c++) {
        popShape[0].type = SHAPE_LINEAR;
        popShape[0].anchor_value = configs[c].N0;
        popShape[0].rate_param = configs[c].gamma;
        popShape[0].anchor_time = 0.0;
        double T_cross = configs[c].N0 / configs[c].gamma;
        for (int trial = 0; trial < 1000; trial++) {
            double xi = -log(ranf());
            double T = drawWaitingTimeSize(0, 0.0, xi, 2);
            TEST_ASSERT_TRUE(T > 0.0);
            TEST_ASSERT_TRUE(T < T_cross);
            TEST_ASSERT_FALSE(isnan(T));
            TEST_ASSERT_FALSE(isinf(T));
        }
    }
}
```

- [ ] **Step 2: Run and verify pass**

If this test fails for some configuration, investigate the `nextafter` clamp in `drawWaitingTimeSize` linear case.

- [ ] **Step 3: Commit**

```bash
git add test/unit/test_shapes.c
git commit -m "Add linear-shape zero-crossing stress test

5 (N0, gamma) configurations x 1000 random xi each. T must
always be in (0, T_cross), no NaN, no inf."
```

---

## Task 16: Replace `TBD` placeholder with confirmation comment

**Files:**
- Modify: `src/core/shapes.c`

- [ ] **Step 1: Audit `shapes.c` for any remaining `TBD` strings**

Run: `grep -n TBD src/core/shapes.c`
Expected: no matches (the stub TBD was replaced by the implementations of Tasks 3-12).

If any TBD remains, finish the corresponding stub before proceeding.

- [ ] **Step 2: Verify all tests still pass**

Run: `make test_shapes && ./build/test_shapes`
Expected: all tests pass.

- [ ] **Step 3: Verify discoal still builds and runs**

```bash
make discoal
./build/discoal 4 1 100 -t 1.0 > /tmp/discoal_smoke.txt && head -5 /tmp/discoal_smoke.txt
```

Expected: typical ms-format output. Confirms shapes module integrates with the main build.

No commit — this is a sanity audit.

---

## Task 17: Commit Phase 1 milestone tag

**Files:** none

- [ ] **Step 1: Tag the foundations milestone**

```bash
git tag -a phase1-shapes-complete -m "Phase 1 of issue-82 design: shape primitives complete

src/core/shapes.{h,c} provides sizeAt, migAt, integratedHazardSize,
integratedHazardMig, drawWaitingTimeSize, drawWaitingTimeMig for
SHAPE_CONSTANT, SHAPE_EXPONENTIAL, SHAPE_LINEAR. Unit-tested
against numerical quadrature and KS distribution checks. No
integration into the simulation engine yet."
```

- [ ] **Step 2: Verify tag is on local branch only (not pushed)**

Run: `git tag --list 'phase1-*'`
Expected: shows the tag.

(We do not push the tag yet; Phase 2 lands first.)

---

## Task 18: Q1 verification — add `detSweepFreqEuler` alongside `detSweepFreq`

Phase 2 of the design spec is the empirical decision: does the Euler step on the deterministic logistic ODE produce a statistically indistinguishable trajectory from the closed-form `detSweepFreq` under constant $N$?

We add a parallel function `detSweepFreqEuler(currentX, ttau, dt, alpha)` that takes the *current* $x$ and advances it by one Euler step, returning the new $x$. This is structurally different from `detSweepFreq` (which is closed-form, returns $x$ at absolute time $\tau$); the verification test will call them in two different binaries with the same RNG seeds and configurations.

**Files:**
- Modify: `src/core/alleleTraj.h`
- Modify: `src/core/alleleTraj.c`
- Create: `test/unit/test_alleleTraj.c`
- Modify: `Makefile`

- [ ] **Step 1: Add declaration to `alleleTraj.h`**

After the existing `detSweepFreq` declaration:

```c
/* Euler-step variant of detSweepFreq for time-varying N support.
 * Given current frequency x and per-step time increment dt,
 * advance by one logistic-ODE Euler step:
 *   x_{t+dt} = x_t + alpha_eff * x * (1 - x) * dt
 * Caller is responsible for clamping x in [0, 1]. */
double detSweepFreqEuler(double x, double dt, double alpha_eff);
```

- [ ] **Step 2: Implement in `alleleTraj.c`**

```c
double detSweepFreqEuler(double x, double dt, double alpha_eff) {
    double dx = alpha_eff * x * (1.0 - x) * dt;
    double x_new = x + dx;
    if (x_new < 0.0) x_new = 0.0;
    if (x_new > 1.0) x_new = 1.0;
    return x_new;
}
```

- [ ] **Step 3: Failing unit test — ODE matches closed-form to leading order**

Create `test/unit/test_alleleTraj.c`:

```c
#include "unity.h"
#include "alleleTraj.h"
#include <math.h>

#ifndef TEST_RUNNER_MODE
void setUp(void) { }
void tearDown(void) { }
#endif

void test_detSweepFreqEuler_against_closed_form(void) {
    /* Run both methods over the entire sweep and compare endpoints */
    double alpha = 200.0;
    /* closed-form sweep duration */
    double epsilon = 0.05 / alpha;
    double ts = -2.0 * log(epsilon) / alpha;
    /* fine Euler grid */
    int N_steps = 100000;
    double dt = ts / N_steps;
    double x_euler = epsilon / (epsilon + (1.0 - epsilon));  /* matches detSweepFreq(0, alpha) */
    for (int i = 0; i < N_steps; i++) {
        x_euler = detSweepFreqEuler(x_euler, dt, alpha);
    }
    double x_closed = detSweepFreq(ts, alpha);
    /* Euler with 100k steps should match closed form to 1e-3 relative */
    TEST_ASSERT_DOUBLE_WITHIN(1e-3 * x_closed, x_closed, x_euler);
}

#ifndef TEST_RUNNER_MODE
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_detSweepFreqEuler_against_closed_form);
    return UNITY_END();
}
#endif
```

- [ ] **Step 4: Wire into Makefile**

After the `test_shapes:` target, add:

```makefile
test_alleleTraj: $(TEST_DIR)/test_alleleTraj.c $(SRC_CORE)/alleleTraj.c $(SRC_CORE)/alleleTraj.h
	@mkdir -p build
	$(CC) $(TEST_CFLAGS) -DUNITY_INCLUDE_DOUBLE -o build/test_alleleTraj $(TEST_DIR)/test_alleleTraj.c $(SRC_CORE)/alleleTraj.c $(SRC_RNG)/xoshiro256pp_compat.c $(UNITY_SOURCES) -I$(UNITY_DIR) -lm -fcommon
```

Add `test_alleleTraj` to `run_tests` dependency list and `./build/test_alleleTraj` to its run sequence.

- [ ] **Step 5: Build and verify**

Run: `make test_alleleTraj && ./build/test_alleleTraj`
Expected: PASS. If FAIL, investigate dt or step-count tuning.

- [ ] **Step 6: Commit**

```bash
git add src/core/alleleTraj.h src/core/alleleTraj.c test/unit/test_alleleTraj.c Makefile
git commit -m "Add detSweepFreqEuler for Q1 verification

Parallel implementation that advances sweep frequency by one
Euler step on the deterministic logistic ODE. Unit test
confirms agreement with closed-form detSweepFreq to 1e-3
over a full sweep with 100k steps."
```

---

## Task 19: Add CLI toggle `--det-sweep-mode` for runtime selection

The Q1 verification harness needs to run the same discoal binary in two modes — closed-form vs Euler — without rebuilding. Add a runtime flag.

**Files:**
- Modify: `src/core/discoal.h` (add global)
- Modify: `src/core/discoal_multipop.c` (parse flag, define global)
- Modify: `src/core/discoalFunctions.c` (dispatch in `proposeTrajectory` and `sweepPhaseEventsGeneralPopNumber`)

- [ ] **Step 1: Add global to `discoal.h`**

In the globals section:

```c
int detSweepMode;  /* 0 = closed-form (default), 1 = Euler */
```

- [ ] **Step 2: Define and parse in `discoal_multipop.c`**

In the global definitions section (where `int npops` etc. are defined):

```c
int detSweepMode = 0;
```

In `getParameters`, after the existing `case 'F'` etc., add a long-option-style parse. Discoal currently uses single-letter flags; we add a special-case prefix check at the top of the option loop. Around line 832 where `case 'F'` lives, the dispatching `switch (argv[args][1])` reads the second char of the flag. For the new flag `--det-sweep-mode` we need a string compare. Add:

```c
/* Long-option pre-check */
if (strcmp(argv[args], "--det-sweep-mode") == 0) {
    args++;
    if (strcmp(argv[args], "closed") == 0) detSweepMode = 0;
    else if (strcmp(argv[args], "euler") == 0) detSweepMode = 1;
    else {
        fprintf(stderr, "--det-sweep-mode: expected 'closed' or 'euler', got '%s'\n", argv[args]);
        exit(1);
    }
    args++;
    continue;  /* skip the switch */
}
```

Place this immediately inside the `while` loop, before the `switch (argv[args][1])`.

- [ ] **Step 3: Dispatch in `proposeTrajectory`**

In `src/core/discoalFunctions.c`, find the `case 'd':` inside the trajectory-walk loop (line 1822):

```c
case 'd':
    x = detSweepFreq(ttau, alpha * currentSizeRatio);
    break;
```

Replace with:

```c
case 'd':
    if (detSweepMode == 0) {
        x = detSweepFreq(ttau, alpha * currentSizeRatio);
    } else {
        x = detSweepFreqEuler(x, tIncOrig, alpha * currentSizeRatio);
    }
    break;
```

Apply the same change at the second occurrence (line 1958, in `sweepPhaseEventsGeneralPopNumber`):

```c
case 'd':
    if (detSweepMode == 0) {
        x = detSweepFreq(ttau, alpha * sizeRatio[0]);
    } else {
        x = detSweepFreqEuler(x, tIncOrig, alpha * sizeRatio[0]);
    }
    break;
```

- [ ] **Step 4: Rebuild and smoke-test both modes**

```bash
make discoal
./build/discoal 4 1 1000 -t 5.0 -r 5.0 -p 1 4 -ws 0.05 -a 200 -x 0.5 -d 12345 67890 > /tmp/q1_closed.txt
./build/discoal 4 1 1000 -t 5.0 -r 5.0 -p 1 4 -ws 0.05 -a 200 -x 0.5 -d 12345 67890 --det-sweep-mode closed > /tmp/q1_closed2.txt
./build/discoal 4 1 1000 -t 5.0 -r 5.0 -p 1 4 -ws 0.05 -a 200 -x 0.5 -d 12345 67890 --det-sweep-mode euler > /tmp/q1_euler.txt
diff /tmp/q1_closed.txt /tmp/q1_closed2.txt
```

Expected: `closed` and the default produce **byte-identical** output (the flag defaulting to `0` must not change behavior). `euler` may differ — that's the question Q1 is asking.

Note: `-ws` invokes a stochastic sweep, which doesn't go through the deterministic-mode path. To exercise it, use `-wd` (deterministic). Adjust the smoke test:

```bash
./build/discoal 4 1 1000 -t 5.0 -r 5.0 -p 1 4 -wd 0.05 -a 200 -x 0.5 -d 12345 67890 > /tmp/q1_closed.txt
./build/discoal 4 1 1000 -t 5.0 -r 5.0 -p 1 4 -wd 0.05 -a 200 -x 0.5 -d 12345 67890 --det-sweep-mode euler > /tmp/q1_euler.txt
```

- [ ] **Step 5: Commit**

```bash
git add src/core/discoal.h src/core/discoal_multipop.c src/core/discoalFunctions.c
git commit -m "Add --det-sweep-mode runtime flag for Q1 verification

Default 'closed' preserves existing detSweepFreq behavior
exactly. 'euler' switches deterministic-mode sweeps to the
detSweepFreqEuler integrator. Used by the Q1 parity harness
to compare distributions empirically."
```

---

## Task 20: Q1 verification harness script

Bash + python to run discoal under both modes across a configuration grid, compute summary statistics with niceStats, compare distributions with KS via the existing `compare_nicestats_distributions.py`.

**Files:**
- Create: `test/parity/q1_detsweep_verification.sh`
- Create: `test/parity/q1_detsweep_analyze.py`

- [ ] **Step 1: Create the directory and harness**

```bash
mkdir -p test/parity
```

Write `test/parity/q1_detsweep_verification.sh`:

```bash
#!/usr/bin/env bash
# Q1 verification: detSweepFreq closed-form vs Euler step under constant N.
#
# Runs discoal in two modes across a grid of (alpha, tau, rho, sample size,
# nsites) configs, generates niceStats summary statistics, and compares
# distributions with KS tests.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
OUT="$HERE/q1_results"
mkdir -p "$OUT"

DISCOAL="$ROOT/build/discoal"
NICESTATS="$ROOT/build/niceStats"
[[ -x "$DISCOAL" ]] || { echo "build discoal first"; exit 1; }
[[ -x "$NICESTATS" ]] || { echo "build niceStats first"; exit 1; }

REPS=10000
N=10
NSITES=10000
THETA=10
RHO=10

ALPHAS=(50 200 1000)
TAUS=(0.01 0.1 0.5)

for alpha in "${ALPHAS[@]}"; do
  for tau in "${TAUS[@]}"; do
    cfg="alpha${alpha}_tau${tau}"
    for mode in closed euler; do
      out="$OUT/${cfg}_${mode}.ms"
      stats="$OUT/${cfg}_${mode}.stats"
      echo "=== $cfg $mode ==="
      "$DISCOAL" "$N" "$REPS" "$NSITES" \
        -t "$THETA" -r "$RHO" \
        -wd "$tau" -a "$alpha" -x 0.5 \
        -d 12345 67890 \
        --det-sweep-mode "$mode" > "$out"
      "$NICESTATS" "$N" "$REPS" < "$out" > "$stats"
    done
  done
done

echo
echo "All configurations run. Pass result directory to q1_detsweep_analyze.py:"
echo "  python3 $HERE/q1_detsweep_analyze.py $OUT"
```

Make executable:

```bash
chmod +x test/parity/q1_detsweep_verification.sh
```

- [ ] **Step 2: Create the analysis script**

Write `test/parity/q1_detsweep_analyze.py`:

```python
#!/usr/bin/env python3
"""Compare niceStats output distributions between closed-form and Euler
deterministic-sweep modes. Bonferroni-corrected KS at p > 0.01.

Usage: q1_detsweep_analyze.py <results_dir>
"""
import sys
import os
import re
from pathlib import Path

import numpy as np
from scipy import stats

def parse_stats(path):
    cols = {}
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for h in header:
            cols[h] = []
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != len(header):
                continue
            for h, p in zip(header, parts):
                try:
                    cols[h].append(float(p))
                except ValueError:
                    pass
    return cols

def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    results = Path(sys.argv[1])
    configs = set()
    for p in results.glob('*_closed.stats'):
        configs.add(p.name.replace('_closed.stats', ''))
    configs = sorted(configs)

    # collect all (config, statistic) pairs and run KS
    comparisons = []
    for cfg in configs:
        a = parse_stats(results / f'{cfg}_closed.stats')
        b = parse_stats(results / f'{cfg}_euler.stats')
        common = sorted(set(a.keys()) & set(b.keys()))
        for stat in common:
            arr_a = np.asarray(a[stat], dtype=float)
            arr_b = np.asarray(b[stat], dtype=float)
            if len(arr_a) < 100 or len(arr_b) < 100:
                continue
            if np.all(arr_a == arr_a[0]) and np.all(arr_b == arr_b[0]):
                continue
            D, p = stats.ks_2samp(arr_a, arr_b)
            comparisons.append((cfg, stat, D, p, len(arr_a), len(arr_b)))

    if not comparisons:
        print('no comparisons made')
        sys.exit(1)

    n = len(comparisons)
    alpha = 0.01 / n  # Bonferroni
    print(f'{n} comparisons, Bonferroni alpha = {alpha:.2e}')
    print(f"{'config':<24} {'stat':<16} {'D':>10} {'p':>10} {'reject?'}")
    rejected = 0
    for cfg, stat, D, p, na, nb in comparisons:
        flag = '  **' if p < alpha else ''
        if p < alpha:
            rejected += 1
        print(f'{cfg:<24} {stat:<16} {D:>10.4f} {p:>10.2e}{flag}')

    print()
    if rejected == 0:
        print(f'PASS: all {n} comparisons within Bonferroni-corrected p > {alpha:.2e}')
        print('CONCLUSION: detSweepFreq and Euler are statistically indistinguishable')
        print('            under constant N for the tested grid.')
        sys.exit(0)
    else:
        print(f'FAIL: {rejected}/{n} comparisons reject equality at Bonferroni p < {alpha:.2e}')
        sys.exit(1)

if __name__ == '__main__':
    main()
```

Make executable:

```bash
chmod +x test/parity/q1_detsweep_analyze.py
```

- [ ] **Step 3: Smoke-test the harness with a tiny config**

Create a one-config dry-run:

```bash
cd /tmp
python3 -c "
import subprocess
DISCOAL='${HOME}/discoal/build/discoal'
NICESTATS='${HOME}/discoal/build/niceStats'
for mode in ('closed', 'euler'):
    p = subprocess.run([DISCOAL, '6', '50', '1000', '-t', '5', '-r', '5',
        '-wd', '0.1', '-a', '200', '-x', '0.5', '-d', '12345', '67890',
        '--det-sweep-mode', mode], capture_output=True, text=True, check=True)
    print(f'{mode}: {len(p.stdout)} bytes, first line: {p.stdout.splitlines()[0]}')
"
```

Expected: both modes run; output sizes may differ slightly. If discoal fails on either mode, debug before proceeding.

- [ ] **Step 4: Commit harness**

```bash
git add test/parity/q1_detsweep_verification.sh test/parity/q1_detsweep_analyze.py
git commit -m "Add Q1 verification harness for detSweepFreq vs Euler

Runs discoal under closed-form and Euler deterministic-sweep
modes across alpha x tau grid, computes niceStats summary
statistics, runs KS comparison with Bonferroni-corrected
threshold. Output PASS or FAIL with per-statistic detail."
```

---

## Task 21: Run Q1 verification end-to-end

**Files:** none modified; produces results in `test/parity/q1_results/`.

- [ ] **Step 1: Build prerequisites**

```bash
make discoal niceStats
```

Expected: both binaries built into `build/`.

- [ ] **Step 2: Run the harness**

```bash
./test/parity/q1_detsweep_verification.sh
```

Expected: 9 configurations × 2 modes = 18 discoal runs, ~10000 reps each. Total runtime depends on alpha and tau but estimate 10-30 minutes on an A100 box. **GPU not used by discoal — runs on CPU.** Confirm a free CPU before kicking off.

- [ ] **Step 3: Analyze**

```bash
python3 test/parity/q1_detsweep_analyze.py test/parity/q1_results
```

Possible outcomes:

- **PASS**: all comparisons p > Bonferroni-corrected threshold. Q1 result is "indistinguishable", future plans may collapse the dispatch in section 4.5 of the design.
- **FAIL**: one or more configurations reject. Q1 result is "distinguishable in some regimes". Future plans must keep the closed-form path for constant $N$.

- [ ] **Step 4: Save analysis output**

```bash
python3 test/parity/q1_detsweep_analyze.py test/parity/q1_results > test/parity/q1_results/analysis.txt 2>&1 || true
```

- [ ] **Step 5: Commit results**

```bash
git add test/parity/q1_results/analysis.txt
git commit -m "Q1 verification result: detSweepFreq vs Euler under constant N

Full results in test/parity/q1_results/analysis.txt. Outcome
documented as addendum to the design spec."
```

(The `.ms` files and per-config `.stats` are not committed — they are large and reproducible.)

---

## Task 22: Document Q1 outcome as addendum to the design spec

**Files:**
- Modify: `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md`

- [ ] **Step 1: Append addendum**

Append a new section at the end of the design spec:

```markdown
## Addendum (2026-MM-DD): Q1 Verification Result

Phase 2 of the implementation plan ran the Q1 verification harness at
`test/parity/q1_detsweep_verification.sh` across 9 (alpha, tau)
configurations with $10^4$ replicates each, comparing the closed-form
`detSweepFreq` and the Euler-step `detSweepFreqEuler` under constant $N$.

**Configuration grid:** alpha in {50, 200, 1000}, tau in {0.01, 0.1, 0.5},
n=10, theta=10, rho=10, nsites=10000.

**Result:** [PASS | FAIL — fill in based on Task 21 output]

[If PASS:]
All [N] comparisons across all (alpha, tau) configurations and all
niceStats summary statistics yielded p > [Bonferroni-corrected threshold].
Conclusion: the two methods are statistically indistinguishable under
constant N. Section 4.5's dispatch on shape type can collapse to
always-Euler in subsequent phases. The design's main body still
specifies the conservative dispatch; subsequent plans may remove it.

[If FAIL:]
[K] of [N] comparisons rejected equality at Bonferroni-corrected
p < [threshold]. The configurations and statistics that rejected:
[list]. Conclusion: the closed-form and Euler step are not statistically
indistinguishable in some regimes. Section 4.5's dispatch on shape type
must be kept; constant-N runs use the closed form, only continuous-N
runs use the Euler step. This is the conservative behavior already
specified in the design.

Raw results: `test/parity/q1_results/analysis.txt` on
`feature/issue-82-time-varying-demography`.
```

Fill in the brackets based on the actual Task 21 output.

- [ ] **Step 2: Commit addendum**

```bash
git add docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md
git commit -m "Spec addendum: Q1 verification result"
```

---

## Task 23: Final phase-1+2 milestone

**Files:** none

- [ ] **Step 1: Run all unit tests**

```bash
make run_tests
```

Expected: all unit-test binaries pass, including `test_shapes`, `test_alleleTraj`, and the existing tests.

- [ ] **Step 2: Smoke-test discoal**

```bash
make discoal
./build/discoal 6 5 1000 -t 5 -r 5 -p 1 6 -d 12345 67890 | head -20
```

Expected: typical ms-format output, no crashes.

- [ ] **Step 3: Tag milestone**

```bash
git tag -a phase2-q1-complete -m "Phases 1+2 of issue-82 design complete

Phase 1: shapes module unit-tested. Phase 2: Q1 verification
result documented as addendum to the design spec. Foundation
ready for the NHPP sampler refactor (next plan)."
```

- [ ] **Step 4: Push branch and tags**

Wait for user confirmation before pushing.

```bash
echo "Ready to push. Run:"
echo "  git push -u origin feature/issue-82-time-varying-demography"
echo "  git push origin phase1-shapes-complete phase2-q1-complete"
```

---

## Self-Review Checklist (run before declaring this plan complete)

- [ ] All tasks have explicit file paths.
- [ ] Every code step contains the actual code an engineer needs.
- [ ] Every test step includes the expected output of `make` / `./build/test_*`.
- [ ] No "TBD" markers except the one in Task 1, which is resolved by Task 16.
- [ ] Function names used in later tasks match the names defined in earlier tasks (`sizeAt`, `migAt`, `integratedHazardSize`, `integratedHazardMig`, `drawWaitingTimeSize`, `drawWaitingTimeMig`, `detSweepFreqEuler`).
- [ ] `Shape` struct fields used everywhere are `type`, `anchor_value`, `rate_param`, `anchor_time`.
- [ ] Commit messages contain no Claude/AI references and no emojis.
- [ ] All tasks live on branch `feature/issue-82-time-varying-demography`.

## Subsequent Plans

Phases 3-8 of the design get their own plans, written after the deliverables of this plan are reviewed:

- **Plan 2 — Phase 3**: NHPP sampler refactor in `neutralPhaseGeneralPopNumber` with shape=CONST only. Bit-equality regression vs current discoal.
- **Plan 3 — Phase 4**: Add EXPONENTIAL shape support, msprime parity for single-pop and 2-pop exp growth.
- **Plan 4 — Phase 5**: Sweep accessor wiring (`sizeAt`/`migAt` in proposeTrajectory, sweepPhase*). Apply Q1 outcome.
- **Plan 5 — Phase 6**: Add LINEAR shape, msprime parity.
- **Plan 6 — Phase 7**: Migration shapes (`'em'` events), importer rewrite, removal of back-derivation hack at `discoalFunctions.c:222-261`. **Closes issue #82.**
- **Plan 7 — Phase 8**: Documentation, full-vocabulary parity sweep, CHANGELOG, PR prep.
