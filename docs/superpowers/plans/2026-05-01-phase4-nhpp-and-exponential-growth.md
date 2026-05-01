# Phase 4: NHPP Sampler + Exponential Growth Shape

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `SHAPE_EXPONENTIAL` support to the simulation engine via a non-homogeneous Poisson process (NHPP) per-component sampler, a new `'g'` event type, and CLI flags `-eg`/`-eG` mirroring ms. Bit-equality for `SHAPE_CONSTANT`-only configs is preserved by dispatching on a per-iteration `all_constant` flag.

**Architecture:** In `neutralPhaseGeneralPopNumber` the inner-loop sampler computes rates as before, then checks whether every active `popShape`/`migShape` is `SHAPE_CONSTANT`. If yes, fall through to the existing `Exp(total_rate)` sampler — bit-identical to Phase 3. If no, draw $T_k$ per rate component via `drawWaitingTimeSize` / `drawWaitingTimeMig` (or `Exp(\lambda)` for time-constant components like recombination), pick the minimum, and fire that event. A new `'g'` event handler updates `popShape[popID]` to `(SHAPE_EXPONENTIAL, anchor_value, alpha, t_event)` where `anchor_value = sizeAt(popID, t_event)` (continuity convention). CLI parses `-eg t pop alpha` and `-eG t alpha` to emit `'g'` events.

**Tech Stack:** C99, Unity, GNU make, bash for parity tests. msprime parity tests are deferred to Phase 4b.

**Spec:** `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md` §3.1, §4.3, §4.5, §4.8 (CLI), and Phase 4 in §5. Note: Phase 3 already wired `sizeAt`/`migAt` into reads; this phase adds the NHPP sampler that those reads enable.

**Deliverables:**
- A flag function `allShapesConstant()` in `src/core/shapes.{h,c}`.
- NHPP per-component sampler in `neutralPhaseGeneralPopNumber` (non-constant path).
- `'g'` event type with dispatch handler in `discoal_multipop.c`.
- CLI flags `-eg time popID alpha` and `-eG time alpha`.
- Bit-equality regression for CONST-only configs (8 configs from Phase 3) still PASS.
- Smoke tests demonstrating that `-eg` produces distributionally different output from `-en` (sanity check that EXP is actually being used).
- All 75 existing unit tests still pass, plus new tests for `'g'` event handling and the NHPP sampler.
- Tag `phase4-exp-complete` on the local branch.

**Out of scope:**
- msprime parity validation (Phase 4b).
- `SHAPE_LINEAR` support (Phase 6).
- Sweep code refactor (Phase 5).
- `'em'` event type (Phase 7).
- Importer changes (Phase 7).
- Removing the `discoalFunctions.c:222-261` back-derivation hack (Phase 7).

**Convention reminder:** Never mention Claude or AI in commits/code/docs. Never use emojis. Stay on `feature/issue-82-time-varying-demography`. Do not push (push is a separate manual step at end).

---

## Task 1: `allShapesConstant()` accessor

A small predicate function returning 1 iff every entry of `popShape[]` and `migShape[][]` is `SHAPE_CONSTANT`. The inner-loop sampler will branch on this.

**Files:**
- Modify: `src/core/shapes.h`
- Modify: `src/core/shapes.c`
- Modify: `test/unit/test_shapes.c`

- [ ] **Step 1: Failing tests**

Add to `test/unit/test_shapes.c` (before `main`):

```c
void test_allShapesConstant_default_state_is_true(void) {
    /* setUp resets all shapes to SHAPE_CONSTANT */
    npops = 3;
    TEST_ASSERT_EQUAL(1, allShapesConstant());
}

void test_allShapesConstant_false_when_a_pop_is_exponential(void) {
    npops = 3;
    popShape[1].type = SHAPE_EXPONENTIAL;
    TEST_ASSERT_EQUAL(0, allShapesConstant());
}

void test_allShapesConstant_false_when_a_pop_is_linear(void) {
    npops = 3;
    popShape[2].type = SHAPE_LINEAR;
    TEST_ASSERT_EQUAL(0, allShapesConstant());
}

void test_allShapesConstant_false_when_a_pair_is_exponential(void) {
    npops = 2;
    migShape[0][1].type = SHAPE_EXPONENTIAL;
    TEST_ASSERT_EQUAL(0, allShapesConstant());
}

void test_allShapesConstant_ignores_pops_outside_npops(void) {
    /* Shapes for pops 5..MAXPOPS that are not in use should not affect the answer */
    npops = 2;
    popShape[5].type = SHAPE_EXPONENTIAL;
    TEST_ASSERT_EQUAL(1, allShapesConstant());
}
```

Register all 5 in `main`.

- [ ] **Step 2: Verify failure**

Run: `make test_shapes && ./build/test_shapes`
Expected: 5 link errors or test failures (`undefined reference to allShapesConstant`).

- [ ] **Step 3: Add declaration to `shapes.h`**

Before `#endif`, add:

```c
/* Returns 1 if every popShape[i] for i in [0, npops) and every migShape[i][j]
 * for i,j in [0, npops) has type SHAPE_CONSTANT. Returns 0 otherwise.
 * Used by the inner-loop sampler to dispatch between the constant-rate
 * Exp(total) path (bit-equal with pre-Phase-4 behavior) and the
 * non-homogeneous Poisson process per-component path. */
int allShapesConstant(void);
```

- [ ] **Step 4: Implement in `shapes.c`**

After the existing functions, add:

```c
int allShapesConstant(void) {
    extern int npops;
    for (int i = 0; i < npops; i++) {
        if (popShape[i].type != SHAPE_CONSTANT) return 0;
        for (int j = 0; j < npops; j++) {
            if (migShape[i][j].type != SHAPE_CONSTANT) return 0;
        }
    }
    return 1;
}
```

- [ ] **Step 5: Verify pass**

Run: `make test_shapes && ./build/test_shapes`
Expected: all tests pass (51 total: 46 prior + 5 new).

- [ ] **Step 6: Commit**

```bash
git add src/core/shapes.h src/core/shapes.c test/unit/test_shapes.c
git commit -m "Add allShapesConstant predicate for sampler dispatch

Returns 1 iff every active shape (popShape[i] for i<npops and
migShape[i][j] for i,j<npops) has type SHAPE_CONSTANT. Used by
the inner-loop sampler in Phase 4 to keep the Exp(total) fast
path for the common all-constant case while routing non-constant
shapes through the NHPP per-component sampler."
```

---

## Task 2: NHPP per-component sampler

Refactor `neutralPhaseGeneralPopNumber` to add a non-constant code path. The constant path is unchanged; the non-constant path draws $T_k$ per rate component using the closed-form `drawWaitingTime*` accessors and selects the minimum.

**Files:**
- Modify: `src/core/discoalFunctions.c`

The existing inner loop (around lines 1607-1700) computes `cRate[i]`, `rRate[i]`, `gcRate[i]`, `mRate[i]`, sums to `totRate`, draws `waitTime = genexp(1.0)/totRate`, then samples a category and fires the corresponding event. Under `SHAPE_CONSTANT`, drawWaitingTimeSize and drawWaitingTimeMig are mathematically identical to `Exp(rate)` draws — but with different RNG sequences. To preserve bit-equality for the all-constant case we keep the existing path and only add a new path when at least one shape varies.

- [ ] **Step 1: Read the existing function**

Run:
```bash
sed -n '1589,1764p' src/core/discoalFunctions.c
```

Identify:
- The outer `while` condition (drives the loop).
- The rate-computation block (lines 1607-1632 area).
- The waitTime draw (line 1632).
- The event-firing logic (lines 1641-1700+) — recomb, GC, migration, coalescence cases.

- [ ] **Step 2: Add include for shapes.h**

If not already present, add `#include "shapes.h"` to the includes block at the top of `discoalFunctions.c`. (It should already be there from Phase 3 Task 4.)

- [ ] **Step 3: Add the all_constant dispatch**

After the rate-computation block (after the line that sets `totRate`), before the existing `waitTime = genexp(1.0) * (1.0/ totRate);` line, add a branch:

```c
		if (allShapesConstant()) {
			/* CONSTANT-only path: existing Exp(total_rate) sampler.
			 * Bit-equal with pre-Phase-4 behavior. */
			waitTime = genexp(1.0) * (1.0/ totRate);
			/* ...existing event-selection code unchanged... */
		} else {
			/* NHPP per-component path: draw T_k per rate component,
			 * take the minimum. Used for SHAPE_EXPONENTIAL / SHAPE_LINEAR. */
			waitTime = drawNHPPWaitingTime(cRate, rRate, gcRate, mRate, totRRate, totGCRate, &winnerKind, &winnerArg);
			/* ... event-selection uses winnerKind/winnerArg ... */
		}
```

The non-constant path needs a helper function `drawNHPPWaitingTime` that returns the minimum waiting time across components AND tells us which component fired (so the caller can dispatch the corresponding event without re-sampling).

- [ ] **Step 4: Define `drawNHPPWaitingTime` helper**

Add a new helper function in `discoalFunctions.c` BEFORE `neutralPhaseGeneralPopNumber`. The function signature:

```c
/* Draw the next event time across all rate components using NHPP per-component
 * draws. Returns the minimum waiting time and writes the firing event class
 * + index into *winnerKind / *winnerArg.
 *
 * winnerKind values: 0 = recomb, 1 = gc, 2 = migration, 3 = coalescence
 * winnerArg meaning depends on winnerKind:
 *   recomb/gc: winnerArg is the population index whose rate fired
 *   migration: winnerArg encodes (src_pop * MAXPOPS + dst_pop)
 *   coalescence: winnerArg is the population index */
static double drawNHPPWaitingTime(double *cRate, double *rRate, double *gcRate,
                                  double *mRate, double totRRate, double totGCRate,
                                  int *winnerKind, int *winnerArg) {
    double T_min = INFINITY;
    *winnerKind = -1;
    *winnerArg = -1;

    /* Recombination: total rate is constant in time within an interval (depends only
     * on total lineage count and rho), so an Exp(totRRate) draw is exact under
     * any shape mix. Random pop is selected proportional to rRate[i]. */
    if (totRRate > 0.0) {
        double T = -log(ranf()) / totRRate;
        if (T < T_min) {
            T_min = T;
            *winnerKind = 0;
            /* winnerArg is decided at firing time by sampling rRate[i]/totRRate. */
            *winnerArg = -1;
        }
    }

    /* Gene conversion: same structure. */
    if (totGCRate > 0.0) {
        double T = -log(ranf()) / totGCRate;
        if (T < T_min) {
            T_min = T;
            *winnerKind = 1;
            *winnerArg = -1;
        }
    }

    /* Coalescence: per-population, time-varying under SHAPE_EXPONENTIAL/LINEAR. */
    for (int i = 0; i < npops; i++) {
        if (popnSizes[i] < 2) continue;
        double xi = -log(ranf());
        double T = drawWaitingTimeSize(i, currentTime, xi, popnSizes[i]);
        if (T > 0.0 && T < T_min) {
            T_min = T;
            *winnerKind = 3;
            *winnerArg = i;
        }
    }

    /* Migration: per-pair, time-varying. */
    for (int i = 0; i < npops; i++) {
        if (popnSizes[i] < 1) continue;
        for (int j = 0; j < npops; j++) {
            if (i == j) continue;
            double m = migAt(i, j, currentTime);
            if (m <= 0.0) continue;
            double xi = -log(ranf());
            double T = drawWaitingTimeMig(i, j, currentTime, xi, popnSizes[i]);
            if (T > 0.0 && T < T_min) {
                T_min = T;
                *winnerKind = 2;
                *winnerArg = i * MAXPOPS + j;
            }
        }
    }

    return T_min;
}
```

(Note: the recomb and GC branches use total rates because per-population recomb/GC rates are mathematically equivalent to drawing one Exp at the total. For coalescence and migration we draw per component because shapes vary. This is a hybrid that minimizes RNG consumption while remaining correct.)

- [ ] **Step 5: Wire dispatch in `neutralPhaseGeneralPopNumber`**

Replace the existing `waitTime = genexp(1.0) * (1.0/ totRate);` and the subsequent event-selection block. The event-selection block currently has nested `if (r < threshold/totRate)` chains — under the constant path we keep that. Under the non-constant path, `winnerKind`/`winnerArg` from `drawNHPPWaitingTime` directly identify the event.

This is the substantial refactor. The cleanest approach: extract the event-firing logic into a small set of helper functions (one per event class), then call the right helper based on either the existing-r-threshold logic (CONST path) or the winnerKind (NHPP path).

Concretely, study lines ~1641-1720 of `discoalFunctions.c`. The block does:
1. `r = ranf()` — uniform
2. If `r < totRRate/totRate`: recomb event (then sample which pop proportional to rRate).
3. Elif `r < (totRRate+totGCRate)/totRate`: GC event.
4. Elif `r < (totRRate+totGCRate+totMRate)/totRate`: migration event (sample pop, then pair).
5. Else: coalescence event (sample pop proportional to cRate).

For the NHPP path, we already know which class won. We just need to fire the right event; sampling within-class (which lineage, which pair) follows the existing logic.

Implementation strategy: keep the existing big switch unchanged for the constant path. For the non-constant path, after `drawNHPPWaitingTime` returns, dispatch on `winnerKind`:

- `winnerKind == 0` (recomb): draw a uniform $r$ and sample pop proportional to `rRate[i]/totRRate`, then call the existing recomb-firing code.
- `winnerKind == 1` (GC): same structure.
- `winnerKind == 2` (migration): pop-pair already known from `winnerArg`; call the existing migration-firing code with that pair.
- `winnerKind == 3` (coalescence): pop already known from `winnerArg`; call the existing coalescence-firing code in that pop.

To minimize duplication, factor the event-firing bodies (which call functions like `coalesceAtTime`, `recombineAtTime`, `migrateAtTime`, etc.) into helper functions. A reasonable factoring:

```c
static void fireRecombEvent(double *bpArray, double currentTime, double *rRate, double totRRate);
static void fireGCEvent(double *bpArray, double currentTime, double *gcRate, double totGCRate);
static void fireCoalEvent(double currentTime, int popID);
static void fireMigEvent(double currentTime, int srcPop, int dstPop);
```

If extracting helpers is too risky (the existing code uses many closures over local variables), an alternative is to inline the dispatch with `if (winnerKind == 0) { ... } else if ... { ... }` inside the same block, leveraging the fact that the original code has clearly-delimited per-event blocks already.

**Recommendation**: do the inline dispatch first (less code churn), get bit-equality regression passing, commit. Refactor to helpers in a follow-up if the inline version is too tangled to maintain.

- [ ] **Step 6: Build**

Run: `make discoal`
Expected: builds cleanly.

- [ ] **Step 7: Bit-equality regression**

```bash
./test/parity/phase3_bit_equality.sh
```
Expected: PASS for all 8 configurations. The constant path is unchanged in this task; under all-constant configs `allShapesConstant()` returns 1 and the existing Exp(total) sampler runs.

If any config FAILs, the dispatch is misrouting CONST configs through the NHPP path. Investigate `allShapesConstant()` and the dispatch placement.

- [ ] **Step 8: Run unit tests**

Run: `make run_tests`
Expected: all 75-80 tests pass (no test changes; this is a refactor).

- [ ] **Step 9: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Add NHPP per-component sampler dispatch for non-constant shapes

neutralPhaseGeneralPopNumber now branches on allShapesConstant():
- All-constant: existing Exp(total_rate) sampler. Bit-equal preserved.
- Otherwise: drawNHPPWaitingTime selects minimum across per-component
  closed-form draws (drawWaitingTimeSize / drawWaitingTimeMig for
  coalescence and migration; Exp(total) for recomb/gc which are
  time-constant). The winning event class and target are returned
  to the caller for direct dispatch.

No SHAPE_EXPONENTIAL / SHAPE_LINEAR shapes are wired into the
runtime yet, so this commit exercises only the constant path.
The non-constant path becomes active in Task 3 once the 'g'
event handler starts setting popShape.type = SHAPE_EXPONENTIAL."
```

---

## Task 3: `'g'` event type and dispatch handler

Define the `'g'` event semantics: at time $t$, set `popShape[popID]` to `(SHAPE_EXPONENTIAL, sizeAt(popID, t), alpha, t)`. The `anchor_value` is set to the population's *current* size at the event time (Convention 1 continuity from the spec §4.5), so the population size is continuous across the shape-change boundary.

**Files:**
- Modify: `src/core/discoal.h` (no changes expected; `Shape` and event types already there)
- Modify: `src/core/discoal_multipop.c` (add `case 'g':` in main event dispatch)

- [ ] **Step 1: Find the main event dispatch**

In `src/core/discoal_multipop.c`, locate the main event-loop switch (around line 228). It has cases for `'n'`, `'s'`, `'p'`, `'a'`, `'A'`. We add `'g'` here.

- [ ] **Step 2: Add `'g'` case after `'n'`**

After the closing `break;` of the `case 'n':` block (and any associated logic; preserve all existing code), insert a new case:

```c
			case 'g':
				currentTime = events[j].time;
				popShape[events[j].popID].type = SHAPE_EXPONENTIAL;
				popShape[events[j].popID].anchor_value = sizeAt(events[j].popID, currentTime);
				popShape[events[j].popID].rate_param = events[j].popnSize;  /* alpha is stored in popnSize field */
				popShape[events[j].popID].anchor_time = currentTime;
				/* Now run the inter-event interval the same way 'n' does. */
				if(activeSweepFlag == 0){
					if(recurSweepMode == 0){
						currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				else{
					if(recurSweepMode == 0){
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
								currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
							currentTime = neutralPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, currentSize);
					}
					else{
						currentTime = sweepPhaseEventsConditionalTrajectory(breakPoints, currentTime, nextTime, sweepSite, \
								currentFreq, &currentFreq, &activeSweepFlag, alpha, currentSize, sweepMode, f0, uA);
						if (currentTime < nextTime)
							currentTime = recurrentSweepPhaseGeneralPopNumber(breakPoints, currentTime, nextTime, &currentFreq, alpha, sweepMode, currentSize);
					}
				}
				break;
```

The body is the same as `case 'n':` MINUS the `currentSize[events[j].popID] = events[j].popnSize;` mutation. The shape change replaces what 'n' would do for size; the population's runtime size at any future time is computed by `sizeAt(popID, t)` from the new EXP shape.

(Note: `events[j].popnSize` is conventionally where 'n' stored the new size value. For 'g' events we reuse this field to store `alpha` (the rate_param). The CLI parser in Task 4 writes alpha there.)

- [ ] **Step 3: Build**

Run: `make discoal`
Expected: builds cleanly.

- [ ] **Step 4: Bit-equality regression**

```bash
./test/parity/phase3_bit_equality.sh
```
Expected: PASS for all 8 configurations. None of the regression configs use 'g' events, so nothing should change.

- [ ] **Step 5: Run unit tests**

Run: `make run_tests`
Expected: all 75-80 tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/core/discoal_multipop.c
git commit -m "Add 'g' event handler for SHAPE_EXPONENTIAL transitions

When a 'g' event fires, popShape[popID] is set to (EXPONENTIAL,
sizeAt(popID, t), alpha, t). The anchor_value uses sizeAt at the
event time so the size is continuous across the shape-change
boundary (msprime forward-time convention; spec section 4.5
Convention 1). The alpha rate_param comes from events[j].popnSize.

The case body otherwise mirrors 'n' — runs the inter-event
interval through neutralPhase or sweepPhase as appropriate.
No CLI / importer flag emits 'g' events yet; that comes in Task 4."
```

---

## Task 4: CLI flags `-eg` and `-eG`

Add CLI parsing for `-eg time popID alpha` (single pop) and `-eG time alpha` (matrix-wide / all pops). Both emit one or more `'g'` events with the time and rate_param set appropriately.

**Files:**
- Modify: `src/core/discoal_multipop.c` (add to the `case 'e':` parser dispatch around line 946)

- [ ] **Step 1: Find the existing `-e<x>` parser**

In `src/core/discoal_multipop.c`, search for the parser block that handles `case 'e':` which dispatches on `argv[args][2]` to handle `-en`, `-ed`, `-ej`, `-ea`. It's around line 946-980.

The existing structure (per Task 1's earlier exploration):

```c
			case 'e' :
				switch(argv[args][2]){
					case 'n':
						/* parse -en time popID size */
						ensureEventsCapacity();
						events[eventNumber].time = atof(argv[++args]) * 2.0;
						events[eventNumber].popID = atoi(argv[++args]);
						events[eventNumber].popnSize = atof(argv[++args]);
						events[eventNumber].type = 'n';
						eventNumber++;
						break;
					case 'd':
					case 'j':
						/* ... */
					case 'a':
						/* ... */
				}
				break;
```

- [ ] **Step 2: Add `case 'g':` and `case 'G':` parsers**

Inside the inner `switch(argv[args][2])`, add:

```c
					case 'g':
						/* -eg time popID alpha */
						ensureEventsCapacity();
						events[eventNumber].time = atof(argv[++args]) * 2.0;
						events[eventNumber].popID = atoi(argv[++args]);
						events[eventNumber].popnSize = atof(argv[++args]);  /* alpha stored in popnSize field */
						events[eventNumber].type = 'g';
						eventNumber++;
						break;
					case 'G':
						/* -eG time alpha — applies to all populations */
						{
							double t = atof(argv[++args]) * 2.0;
							double alpha_val = atof(argv[++args]);
							for (int p = 0; p < npops; p++) {
								ensureEventsCapacity();
								events[eventNumber].time = t;
								events[eventNumber].popID = p;
								events[eventNumber].popnSize = alpha_val;
								events[eventNumber].type = 'g';
								eventNumber++;
							}
						}
						break;
```

The time is multiplied by 2.0 for backwards compatibility with the existing `-en`/`-ed` semantics (CLI time is in 2N units, internal time is in 4N units).

The `-eG` form requires `npops` to be set BEFORE `-eG` appears in argv. discoal already requires `-p npops samp1 samp2 ...` to come before `-en` etc. for the same reason — the parser is sequential.

- [ ] **Step 3: Build**

Run: `make discoal`
Expected: builds cleanly.

- [ ] **Step 4: Smoke-test the new flag**

```bash
./build/discoal 6 1 1000 -t 5 -r 5 -eg 0.5 0 50 -d 12345 67890 | head -10
```
Expected: ms-format output. Compare against the no-`-eg` case to confirm behavior changes:

```bash
./build/discoal 6 1 1000 -t 5 -r 5 -d 12345 67890 > /tmp/no_eg.txt
./build/discoal 6 1 1000 -t 5 -r 5 -eg 0.5 0 50 -d 12345 67890 > /tmp/with_eg.txt
diff /tmp/no_eg.txt /tmp/with_eg.txt | head -20
```

Expected: substantial diff — the simulations differ because the `-eg` triggers a 'g' event setting popShape to EXPONENTIAL with alpha=50 at t=1.0 (after the 2x conversion).

- [ ] **Step 5: Bit-equality regression for CONST configs**

```bash
./test/parity/phase3_bit_equality.sh
```
Expected: PASS for all 8 configs. None use `-eg`.

- [ ] **Step 6: Run unit tests**

Run: `make run_tests`
Expected: all 75-80 tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/core/discoal_multipop.c
git commit -m "Add -eg and -eG CLI flags for exponential-growth events

-eg time popID alpha emits one 'g' event setting popShape[popID]
to SHAPE_EXPONENTIAL with the given growth rate at the given time
(in 2N units, multiplied by 2.0 for internal 4N representation
matching -en convention).

-eG time alpha emits one 'g' event per active population, all
with the same alpha. -p must come before -eG for npops to be set.

alpha is the per-generation forward-time growth rate per the
msprime convention. The CLI matches ms's -eg / -eG flags."
```

---

## Task 5: Sanity test that `-eg` produces different output than `-en` in equivalent configs

A quick test demonstrating that `-eg` actually exercises the EXPONENTIAL shape rather than silently no-op'ing.

**Files:**
- Create: `test/parity/phase4_eg_smoke.sh`

- [ ] **Step 1: Create the test**

Write `test/parity/phase4_eg_smoke.sh`:

```bash
#!/usr/bin/env bash
# Phase 4 smoke test: confirm -eg produces distributionally different output
# from comparable -en configurations. This is a sanity check that the EXP
# shape is actually being used by the simulation engine.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"

[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

# Three pairs of configs; each pair (constant vs eg) should yield different output
# at the same RNG seeds.

NREPS=20
configs=(
  "no_demography|6 $NREPS 1000 -t 5 -r 5 -d 12345 67890"
  "with_eg_alpha50|6 $NREPS 1000 -t 5 -r 5 -eg 0.5 0 50 -d 12345 67890"
  "with_en_size_change|6 $NREPS 1000 -t 5 -r 5 -en 0.5 0 0.5 -d 12345 67890"
)

OUT="$HERE/phase4_eg_smoke_results"
mkdir -p "$OUT"

for cfg in "${configs[@]}"; do
  tag="${cfg%%|*}"
  args="${cfg##*|}"
  echo "=== $tag ==="
  # shellcheck disable=SC2086
  "$DISCOAL" $args > "$OUT/${tag}.ms" 2>"$OUT/${tag}.err" || {
    echo "  FAIL: discoal exited non-zero (see $OUT/${tag}.err)"
    exit 1
  }
done

# Confirm with_eg_alpha50 differs from no_demography
if diff <(sed '1d' "$OUT/no_demography.ms") <(sed '1d' "$OUT/with_eg_alpha50.ms") > /dev/null; then
  echo "FAIL: -eg config produced byte-identical output to no-demography config (eg appears to be a no-op)"
  exit 1
fi

# Confirm with_eg_alpha50 differs from with_en_size_change
if diff <(sed '1d' "$OUT/with_eg_alpha50.ms") <(sed '1d' "$OUT/with_en_size_change.ms") > /dev/null; then
  echo "FAIL: -eg and -en configs produced byte-identical output (eg may be silently treated as en)"
  exit 1
fi

echo
echo "PASS: -eg produces distributionally different output from constant-N and from -en configurations."
```

Make executable: `chmod +x test/parity/phase4_eg_smoke.sh`.

- [ ] **Step 2: Add `.gitignore`**

Create `test/parity/phase4_eg_smoke_results/.gitignore`:

```
*.ms
*.err
```

- [ ] **Step 3: Run the smoke test**

```bash
./test/parity/phase4_eg_smoke.sh
```
Expected: PASS.

If the smoke test fails (`-eg` produces same output as constant-N or as `-en`), the EXP shape isn't being exercised. Investigate:
- Is the 'g' event actually being parsed and added to events[]? (Add `fprintf(stderr, ...)` to the parser to confirm.)
- Is the 'g' case in the main switch being hit? (Add `fprintf(stderr, ...)` to confirm.)
- Is `popShape[popID].type` actually SHAPE_EXPONENTIAL after the event fires?
- Is `allShapesConstant()` returning 0 after the event fires?

- [ ] **Step 4: Commit**

```bash
git add test/parity/phase4_eg_smoke.sh test/parity/phase4_eg_smoke_results/.gitignore
git commit -m "Add -eg smoke test to confirm EXP shape is exercised

Three configs at the same RNG seed: no demography, -eg with
alpha=50 at t=0.5, and -en with size change at t=0.5. The smoke
test confirms -eg differs from both no-demography and -en —
sanity check that the new event type is actually wired into the
simulation engine and not silently no-op'ing."
```

---

## Task 6: Final regression sweep + tag

- [ ] **Step 1: Full regression**

```bash
make discoal discoal_pre_phase3
./test/parity/phase3_bit_equality.sh
./test/parity/phase4_eg_smoke.sh
```

Expected: both PASS. The Phase 3 regression confirms CONST configs are still bit-equal; the Phase 4 smoke confirms `-eg` is exercised.

- [ ] **Step 2: Full unit tests**

```bash
make run_tests
```

Expected: all tests pass.

- [ ] **Step 3: Tag the milestone**

```bash
git tag -a phase4-exp-complete -m "$(cat <<'TAGEOF'
Phase 4 of issue-82 design complete

Adds SHAPE_EXPONENTIAL support to the simulation engine:

- allShapesConstant() predicate gates the inner-loop sampler.
- neutralPhaseGeneralPopNumber dispatches on the predicate:
  - All-constant: existing Exp(total_rate) sampler. Bit-equal
    with Phase 3 (8/8 configs in the regression grid pass).
  - Otherwise: drawNHPPWaitingTime selects minimum across
    per-component closed-form draws using drawWaitingTimeSize
    and drawWaitingTimeMig.
- 'g' event type sets popShape[popID] = (EXPONENTIAL, sizeAt(t),
  alpha, t). Continuity preserved across the shape boundary
  (msprime forward-time convention).
- CLI flags -eg time popID alpha and -eG time alpha emit 'g'
  events.

Out of scope for Phase 4: msprime parity validation (Phase 4b),
SHAPE_LINEAR (Phase 6), sweep-phase EXP support (Phase 5),
'em' migration shapes (Phase 7), importer changes (Phase 7).
TAGEOF
)"
```

- [ ] **Step 4: List tags**

Run: `git tag --list 'phase*'`
Expected: `phase1-shapes-complete`, `phase2-q1-complete`, `phase3-wiring-complete`, `phase4-exp-complete`.

(Push not done yet; user typically pushes when ready.)

---

## Self-Review Checklist (run before declaring this plan complete)

- [ ] All tasks have explicit file paths.
- [ ] Every code step contains the actual code or precise instructions an engineer needs.
- [ ] Tests are included for each new function/behavior.
- [ ] `allShapesConstant()`, `drawNHPPWaitingTime()` names consistent across tasks.
- [ ] `popShape`/`migShape`/`Shape` field names used uniformly.
- [ ] Bit-equality preserved for CONST-only configs at every step.
- [ ] Commit messages contain no Claude/AI references and no emojis.
- [ ] All work lands on `feature/issue-82-time-varying-demography`.

## What's Next After This Plan

- **Plan: Phase 4b** — msprime parity validation. Build a python harness using msprime that simulates the same demography. Run discoal under `-eg` and msprime in matched configurations; compare SFS, π, Tajima's D, segregating sites, pairwise coalescent-time distributions via KS / chi-squared at Bonferroni-corrected p > 0.01. This validates that SHAPE_EXPONENTIAL produces statistically correct results.

- **Plan: Phase 5** — Sweep accessor wiring using closed-form `detSweepFreqGeneral` per the spec revision. Removes `--det-sweep-mode` flag and `detSweepFreqEuler` function as side effect.

- **Plan: Phase 6** — Add `SHAPE_LINEAR` support, msprime parity for linear-growth configs.

- **Plan: Phase 7** — `'em'` migration shape events, importer rewrite, removal of back-derivation hack. **Closes issue #82.**

- **Plan: Phase 8** — Documentation, full-vocabulary parity sweep, CHANGELOG, PR prep.
