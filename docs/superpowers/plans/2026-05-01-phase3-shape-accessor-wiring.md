# Phase 3: Shape Accessor Wiring (Bit-Equal Foundation)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire `sizeAt(popID, t)` and `migAt(srcPopID, dstPopID, t)` accessors into `neutralPhaseGeneralPopNumber` and the `'n'` event handler. With `SHAPE_CONSTANT` only, the simulation output is **bit-equal** to the pre-Phase-3 binary across the regression configuration grid.

**Architecture:** Add `initializeShapesFromGlobals()` to `src/core/shapes.c` that copies the `currentSize[]` and `migMatConst[][]` globals into `popShape[]` and `migShape[][]` as `SHAPE_CONSTANT` shapes anchored at $t=0$. Call it from `initialize()` in `src/core/discoalFunctions.c` (the per-replicate init). Refactor the rate computation in `neutralPhaseGeneralPopNumber` to read sizes via `sizeAt(i, currentTime)` and migration rates via `migAt(i, j, currentTime)`. Update the `'n'` event handler in `src/core/discoal_multipop.c` to also write `popShape[popID]` whenever it writes `currentSize[popID]`. The `Exp(total_rate)` sampler is unchanged — under `SHAPE_CONSTANT` the rate computations are mathematically identical, so RNG-call order is preserved and bit-equality holds.

**Tech Stack:** C99, Unity test framework, bash for parity tests, GNU make.

**Spec:** `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md` §4.5 (narrowed Phase 3 scope: accessor wiring only; the per-component NHPP draw is deferred to Phase 4 when EXPONENTIAL is added and the rate-varies-within-interval problem actually requires it).

**Deliverables:**
- `initializeShapesFromGlobals()` in `src/core/shapes.{h,c}`, unit-tested.
- `neutralPhaseGeneralPopNumber` reads through `sizeAt` / `migAt`.
- `'n'` event handler writes `popShape`.
- A regression harness at `test/parity/phase3_bit_equality.sh` that builds two binaries (pre-Phase-3 reference and current) and `diff`s their output across 8 configurations.
- Bit-equality verified across the full configuration grid.
- All 75 existing unit tests still pass.
- Tag `phase3-wiring-complete` on the local branch.

**Out of scope:**
- The NHPP per-component draw (deferred to Phase 4).
- Sweep code refactor (Phase 5).
- New event types `'g'` and `'em'` (Phase 4 / Phase 7).
- Importer changes (Phase 7).
- Removing the `discoalFunctions.c:222-261` back-derivation hack (Phase 7).
- Any user-facing behavior change.

**Convention reminder:** Never mention Claude or AI in commits/code/docs. Never use emojis. Stay on `feature/issue-82-time-varying-demography`. Do not push.

---

## Task 1: Build pre-Phase-3 reference binary

The regression target is the discoal binary as it stands at the start of Phase 3 (HEAD = `18ddde1` "Spec revision: deterministic sweep uses closed-form, not Euler"). We build it once, save it at `build/discoal_pre_phase3`, and diff against it after each refactor task.

**Files:**
- Modify: `Makefile` (add `discoal_pre_phase3` target)

- [ ] **Step 1: Identify the pre-Phase-3 SHA**

Run: `git log -1 --pretty=format:'%H' HEAD`
Expected: a single commit SHA. Record this SHA and the short-SHA in the next step.

- [ ] **Step 2: Add `discoal_pre_phase3` Makefile target**

In `Makefile`, after the `discoal_edited:` recipe (around line 66-68), add a new recipe that builds discoal from the recorded SHA via a temp git worktree:

```makefile
# Phase 3 regression reference: discoal binary at the SHA where Phase 3 began.
# This recipe checks out that SHA in a temp worktree, builds discoal there,
# copies the binary to build/discoal_pre_phase3, and removes the worktree.
PHASE3_REF_SHA = 18ddde1
discoal_pre_phase3:
	@mkdir -p build
	@if [ -x build/discoal_pre_phase3 ]; then \
		echo "build/discoal_pre_phase3 already present; remove it to rebuild"; \
		exit 0; \
	fi
	@WT=$$(mktemp -d) && \
	  git worktree add --detach "$$WT" $(PHASE3_REF_SHA) && \
	  $(MAKE) -C "$$WT" discoal && \
	  cp "$$WT/build/discoal" build/discoal_pre_phase3 && \
	  git worktree remove "$$WT"
	@echo "Built pre-Phase-3 reference: build/discoal_pre_phase3"
```

(Replace `18ddde1` with the actual current HEAD short-SHA from Step 1 if different.)

- [ ] **Step 3: Build the reference binary**

Run: `make discoal_pre_phase3`
Expected: build succeeds, `ls -la build/discoal_pre_phase3` shows an executable file.

- [ ] **Step 4: Smoke-test the reference**

Run:
```bash
./build/discoal_pre_phase3 4 1 100 -t 1.0 -d 12345 67890 | head -5
```
Expected: typical ms-format output identical to what current `./build/discoal` produces (since they're built from the same SHA).

- [ ] **Step 5: Commit Makefile change**

```bash
git add Makefile
git commit -m "Add discoal_pre_phase3 Makefile recipe for Phase 3 regression

Builds discoal from the Phase 3 base SHA in a temp git worktree
and copies the binary to build/discoal_pre_phase3. Used as the
bit-equality regression target during Phase 3 wiring tasks."
```

---

## Task 2: Bit-equality regression harness

A small bash script that runs both binaries with identical seeds and configs, then `diff`s their outputs ignoring the line-3 command-line echo (which differs because the path prefix differs).

**Files:**
- Create: `test/parity/phase3_bit_equality.sh`

- [ ] **Step 1: Create the harness**

Write `test/parity/phase3_bit_equality.sh`:

```bash
#!/usr/bin/env bash
# Phase 3 regression: verify build/discoal output is byte-identical to
# build/discoal_pre_phase3 across a configuration grid. Excludes the
# line-3 command-line echo (which differs because the binary paths differ).

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

PRE="$ROOT/build/discoal_pre_phase3"
CUR="$ROOT/build/discoal"
[[ -x "$PRE" ]] || { echo "FAIL: build/discoal_pre_phase3 missing (make discoal_pre_phase3)"; exit 1; }
[[ -x "$CUR" ]] || { echo "FAIL: build/discoal missing (make discoal)"; exit 1; }

# Configuration grid. Each line is: tag <discoal args>
CONFIGS=(
  "single_pop_neutral|6 5 1000 -t 5 -r 5 -d 12345 67890"
  "single_pop_size_change|6 5 1000 -t 5 -r 5 -en 0.5 0 0.1 -en 1.0 0 1.0 -d 12345 67890"
  "two_pop_constant_mig|8 5 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -m 1 0 0.3 -d 12345 67890"
  "two_pop_split|8 5 1000 -t 5 -r 5 -p 2 4 4 -ed 0.5 0 1 -d 12345 67890"
  "two_pop_split_mig|8 5 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -m 1 0 0.3 -ed 1.0 0 1 -d 12345 67890"
  "three_pop_combo|9 5 1000 -t 5 -r 5 -p 3 3 3 3 -m 0 1 0.5 -ed 0.5 0 1 -ed 1.5 1 2 -d 12345 67890"
  "single_pop_with_admixture|8 5 1000 -t 5 -r 5 -p 2 4 4 -ea 0.2 0 0 1 0.7 -ed 1.0 0 1 -d 12345 67890"
  "ancient_sample|6 5 1000 -t 5 -r 5 -p 2 4 0 -A 2 1 0.3 -ed 0.5 0 1 -d 12345 67890"
)

OUT="$HERE/phase3_bit_equality_results"
mkdir -p "$OUT"

fails=0
for cfg in "${CONFIGS[@]}"; do
  tag="${cfg%%|*}"
  args="${cfg##*|}"
  echo "=== $tag ==="
  # shellcheck disable=SC2086
  "$PRE" $args > "$OUT/${tag}_pre.ms" 2>"$OUT/${tag}_pre.err"
  # shellcheck disable=SC2086
  "$CUR" $args > "$OUT/${tag}_cur.ms" 2>"$OUT/${tag}_cur.err"
  # diff excluding line 3 (the command-line echo)
  if diff <(sed '3d' "$OUT/${tag}_pre.ms") <(sed '3d' "$OUT/${tag}_cur.ms") > "$OUT/${tag}.diff"; then
    echo "  PASS"
  else
    echo "  FAIL  (see $OUT/${tag}.diff)"
    head -20 "$OUT/${tag}.diff" || true
    fails=$((fails+1))
  fi
done

if [[ $fails -eq 0 ]]; then
  echo
  echo "PASS: all ${#CONFIGS[@]} configurations are byte-equal."
  exit 0
else
  echo
  echo "FAIL: $fails of ${#CONFIGS[@]} configurations differ."
  exit 1
fi
```

Make executable: `chmod +x test/parity/phase3_bit_equality.sh`.

- [ ] **Step 2: Smoke-run the harness**

Build both binaries, then run:
```bash
make discoal discoal_pre_phase3
./test/parity/phase3_bit_equality.sh
```

Expected: at this point both binaries are built from the same SHA, so all 8 configurations PASS. This validates the harness mechanics (runs the configs, parses output, ignores line 3) before any Phase 3 changes land.

- [ ] **Step 3: Add `.gitignore` for the result files**

Create `test/parity/phase3_bit_equality_results/.gitignore`:

```
*.ms
*.err
*.diff
```

(The result files are reproducible from the harness; we don't commit them.)

- [ ] **Step 4: Commit harness**

```bash
git add test/parity/phase3_bit_equality.sh test/parity/phase3_bit_equality_results/.gitignore
git commit -m "Add Phase 3 bit-equality regression harness

8 configurations covering single-pop neutral, single-pop with size
changes, two-pop with constant migration, two-pop with split,
two-pop with split + migration, three-pop combo, admixture, and
ancient samples. Each config runs both binaries with identical
seeds and diffs outputs ignoring the line-3 command-line echo."
```

---

## Task 3: `initializeShapesFromGlobals` accessor

Add a function that copies `currentSize[]` and `migMatConst[][]` into `popShape[]` and `migShape[][]` respectively, all as `SHAPE_CONSTANT` shapes anchored at $t=0$.

**Files:**
- Modify: `src/core/shapes.h`
- Modify: `src/core/shapes.c`
- Modify: `test/unit/test_shapes.c`

- [ ] **Step 1: Write the failing test**

Add to `test/unit/test_shapes.c` (after the existing test functions, before `main`):

```c
extern double currentSize[];
extern double migMatConst[MAXPOPS][MAXPOPS];
extern int npops;

void test_initializeShapesFromGlobals_copies_currentSize(void) {
    /* Set globals as if parameters were just parsed */
    npops = 3;
    currentSize[0] = 1.0;
    currentSize[1] = 0.5;
    currentSize[2] = 2.0;

    /* Pre-condition: popShape entries set to defaults by setUp */
    /* Pre-condition: migMatConst entries are zero by default */

    initializeShapesFromGlobals();

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, popShape[0].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 1.0, popShape[0].anchor_value);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, popShape[0].rate_param);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, popShape[0].anchor_time);

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, popShape[1].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.5, popShape[1].anchor_value);

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, popShape[2].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 2.0, popShape[2].anchor_value);
}

void test_initializeShapesFromGlobals_copies_migMatConst(void) {
    npops = 2;
    currentSize[0] = 1.0;
    currentSize[1] = 1.0;
    migMatConst[0][1] = 0.7;
    migMatConst[1][0] = 0.3;
    migMatConst[0][0] = 0.0;  /* diagonals are zero */
    migMatConst[1][1] = 0.0;

    initializeShapesFromGlobals();

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, migShape[0][1].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.7, migShape[0][1].anchor_value);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, migShape[0][1].rate_param);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, migShape[0][1].anchor_time);

    TEST_ASSERT_EQUAL(SHAPE_CONSTANT, migShape[1][0].type);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.3, migShape[1][0].anchor_value);

    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, migShape[0][0].anchor_value);
    TEST_ASSERT_DOUBLE_WITHIN(1e-12, 0.0, migShape[1][1].anchor_value);
}
```

Register both in `main`:
```c
RUN_TEST(test_initializeShapesFromGlobals_copies_currentSize);
RUN_TEST(test_initializeShapesFromGlobals_copies_migMatConst);
```

- [ ] **Step 2: Verify failure**

Run: `make test_shapes && ./build/test_shapes`
Expected: 2 tests FAIL (function `initializeShapesFromGlobals` is not declared yet — link error or compile error).

If you get a compile error rather than a runtime FAIL, that's fine — TDD with C handles missing-symbol errors as the failing-test signal.

- [ ] **Step 3: Add declaration to `shapes.h`**

In `src/core/shapes.h`, before the `#endif`, add:

```c
/* Initialize popShape[] and migShape[][] from the current values in
 * currentSize[] and migMatConst[][]. Sets all shapes to SHAPE_CONSTANT
 * anchored at t=0. Called once per replicate from initialize(). */
void initializeShapesFromGlobals(void);
```

- [ ] **Step 4: Implement in `shapes.c`**

In `src/core/shapes.c`, after the existing functions, add:

```c
void initializeShapesFromGlobals(void) {
    extern double currentSize[];
    extern double migMatConst[MAXPOPS][MAXPOPS];
    extern int npops;

    for (int i = 0; i < npops; i++) {
        popShape[i].type = SHAPE_CONSTANT;
        popShape[i].anchor_value = currentSize[i];
        popShape[i].rate_param = 0.0;
        popShape[i].anchor_time = 0.0;
        for (int j = 0; j < npops; j++) {
            migShape[i][j].type = SHAPE_CONSTANT;
            migShape[i][j].anchor_value = migMatConst[i][j];
            migShape[i][j].rate_param = 0.0;
            migShape[i][j].anchor_time = 0.0;
        }
    }
}
```

- [ ] **Step 5: Verify pass**

Run: `make test_shapes && ./build/test_shapes`
Expected: all tests PASS, including the two new ones (77 tests total).

- [ ] **Step 6: Commit**

```bash
git add src/core/shapes.h src/core/shapes.c test/unit/test_shapes.c
git commit -m "Add initializeShapesFromGlobals to seed shape state from globals

Copies currentSize[] -> popShape[] and migMatConst[][] -> migShape[][]
as SHAPE_CONSTANT shapes anchored at t=0. Will be called from
initialize() in the next task to seed per-replicate shape state."
```

---

## Task 4: Wire `initializeShapesFromGlobals` into `initialize()`

The per-replicate `initialize()` function in `src/core/discoalFunctions.c` is the right call site — it runs after parameters are parsed and before the main event loop. Currently it sets `migMat[][]` from `migMatConst[][]` at line 220-225; we add the shape-state init right after.

**Files:**
- Modify: `src/core/discoalFunctions.c`

- [ ] **Step 1: Add include**

In `src/core/discoalFunctions.c`, find the existing `#include` block at the top of the file. Add:

```c
#include "shapes.h"
```

- [ ] **Step 2: Call the initializer**

In `initialize()`, find this block (around line 218-225):

```c
	if (npops>1){
		if(tDiv==666 && migFlag == 0){
			fprintf(stderr,"tDiv or migration not set in population split model\n");
			exit(1);
		}
		//initialize migration matrix
		for(i=0;i<npops;i++){
			for(j=0;j<npops;j++){
				migMat[i][j]=migMatConst[i][j];
			}
		}
```

Add immediately after the `migMat` initialization loop (after the second `}` closing `for(j=...)`):

```c
		//initialize shape state for time-varying parameter framework
		initializeShapesFromGlobals();
```

This places the shape init inside the `npops > 1` branch, which is fine for multi-pop runs. Single-pop runs don't need migration shapes but still need population-size shapes — handle that next:

- [ ] **Step 3: Also init shapes for single-pop case**

After the closing `}` of the `if (npops > 1)` block (around line 280-290 — search for the matching brace), add a single-pop fallback. Look for code that runs after the migration init / before any further per-rep setup.

Actually a cleaner approach: move the `initializeShapesFromGlobals()` call to BEFORE the `if (npops > 1)` block so it always runs. Replace the change from Step 2 with this: find the block

```c
	activeSites = nSites;
	if (npops>1){
```

and insert before the `if`:

```c
	//initialize shape state for time-varying parameter framework
	initializeShapesFromGlobals();
	activeSites = nSites;
	if (npops>1){
```

Revert the placement from Step 2 if you used it (the Step 3 location is the correct one).

- [ ] **Step 4: Verify build**

Run: `make discoal`
Expected: build succeeds, no warnings about `initializeShapesFromGlobals`.

- [ ] **Step 5: Smoke-test discoal still runs**

Run:
```bash
./build/discoal 4 1 100 -t 1.0 -d 12345 67890 | head -5
```
Expected: typical ms-format output.

- [ ] **Step 6: Bit-equality regression**

```bash
./test/parity/phase3_bit_equality.sh
```
Expected: PASS for all 8 configurations. The shape state is computed but not yet read by any code path, so output is unchanged.

- [ ] **Step 7: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Call initializeShapesFromGlobals from per-replicate initialize()

Seeds popShape[] and migShape[][] from the current parameter values
at the start of each replicate. Shape state is computed but not
yet read; bit-equality preserved across the regression grid."
```

---

## Task 5: Refactor coalescent rate read in `neutralPhaseGeneralPopNumber`

Replace `sizeRatio[i]` reads in the coalescent-rate computation with `sizeAt(i, currentTime)` calls. Under `SHAPE_CONSTANT`, `sizeAt` returns the anchor value, which equals `currentSize[i]` (which is what was passed into `sizeRatio`); the rate computation is bit-identical.

**Files:**
- Modify: `src/core/discoalFunctions.c`

- [ ] **Step 1: Inspect the call site**

Look at `src/core/discoalFunctions.c:1614`:

```c
		cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 / sizeRatio[i];
```

`sizeRatio` is the parameter passed in by the caller (`currentSize` array). `currentTime` is also accessible (it's the `startTime` parameter or a local variable; verify by reading the function signature and locals).

Actually `currentTime` is a global — confirm by `grep -n 'double currentTime' src/core/discoal.h`. (It is.)

- [ ] **Step 2: Replace the sizeRatio read**

Change line 1614 from:

```c
		cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 / sizeRatio[i];
```

to:

```c
		cRate[i] = popnSizes[i] * (popnSizes[i] - 1) * 0.5 / sizeAt(i, currentTime);
```

(Note: the parameter `sizeRatio` is no longer read here. It's still used elsewhere in the function — leave those for the next task.)

- [ ] **Step 3: Build**

Run: `make discoal`
Expected: builds cleanly.

- [ ] **Step 4: Bit-equality regression**

```bash
./test/parity/phase3_bit_equality.sh
```
Expected: PASS for all 8 configurations. `sizeAt(i, currentTime)` for `SHAPE_CONSTANT` equals `popShape[i].anchor_value`, which equals `currentSize[i]` (the seed value). Bit-equality holds.

If any config FAILs, investigate whether the shape state was initialized correctly (Task 4) and whether `currentTime` is the right value at the call site.

- [ ] **Step 5: Run unit tests**

Run: `make run_tests`
Expected: all 75 tests pass (no test changes; this is a refactor).

- [ ] **Step 6: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Read coalescent rate denominator via sizeAt in neutralPhaseGeneralPopNumber

cRate[i] now uses sizeAt(i, currentTime) instead of sizeRatio[i].
Under SHAPE_CONSTANT (the only shape type currently set by
initializeShapesFromGlobals), sizeAt returns the anchor value
which equals the prior sizeRatio seed; output is bit-equal."
```

---

## Task 6: Refactor migration rate read in `neutralPhaseGeneralPopNumber`

The function reads `migMat[i][j]` in two places (line 1618 for the `mRate[i]` row sum, and around line 1668-1673 for the per-pair migration sampling). Replace both with `migAt(i, j, currentTime)`.

**Files:**
- Modify: `src/core/discoalFunctions.c`

- [ ] **Step 1: Replace `mRate[i]` row sum**

Find line 1618:

```c
		for(j=0;j<npops;j++) mRate[i]+=migMat[i][j];
```

Note: the sum includes `migMat[i][i]` which is conventionally zero (diagonal). That's preserved by `migAt` because `migShape[i][i]` is initialized to anchor_value=0.

Change to:

```c
		for(j=0;j<npops;j++) mRate[i]+=migAt(i, j, currentTime);
```

- [ ] **Step 2: Replace per-pair migration sampling**

Find the block around line 1668-1674. The exact code looks like:

```c
						eSum = migMat[i][0]* popnSizes[i] * 0.5;
						j=0;
						while(eSum < eventProb){
							eSum += migMat[i][++j] * popnSizes[i] * 0.5;
```

Change both `migMat[i][...]` reads to `migAt(i, ..., currentTime)`:

```c
						eSum = migAt(i, 0, currentTime) * popnSizes[i] * 0.5;
						j=0;
						while(eSum < eventProb){
							eSum += migAt(i, ++j, currentTime) * popnSizes[i] * 0.5;
```

(There may be additional `migMat[i][...]` references in this same block — search the function for any remaining `migMat\[` occurrences and replace each with `migAt(...)`.)

- [ ] **Step 3: Build**

Run: `make discoal`
Expected: builds cleanly.

- [ ] **Step 4: Bit-equality regression**

```bash
./test/parity/phase3_bit_equality.sh
```
Expected: PASS for all 8 configurations. `migAt(i, j, t)` returns `migShape[i][j].anchor_value` for `SHAPE_CONSTANT`, which equals `migMatConst[i][j]` (which is what `migMat[i][j]` was set to by `initialize()`). Bit-equality holds.

- [ ] **Step 5: Run unit tests**

Run: `make run_tests`
Expected: all 75 tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/core/discoalFunctions.c
git commit -m "Read migration rate via migAt in neutralPhaseGeneralPopNumber

Replaces all migMat[i][j] reads in the function (row-sum mRate[i]
and per-pair event sampling) with migAt(i, j, currentTime). Under
SHAPE_CONSTANT, migAt returns the anchor value which equals the
prior migMat seed; output is bit-equal."
```

---

## Task 7: Update `'n'` event handler to write `popShape`

When an `'n'` event fires (size change for a population), the dispatch in `discoal_multipop.c` updates `currentSize[popID]`. Now we also need to update `popShape[popID]` so the next call to `sizeAt` returns the new value.

**Files:**
- Modify: `src/core/discoal_multipop.c`

- [ ] **Step 1: Find the `'n'` dispatch**

In `src/core/discoal_multipop.c`, search for `case 'n':` in the main event-dispatch switch (around line 229). The relevant code:

```c
			case 'n':
				currentTime = events[j].time;
				currentSize[events[j].popID] = events[j].popnSize;
```

(There may be additional logic in this case for sweep handling — leave that alone.)

- [ ] **Step 2: Add the popShape update**

Right after the `currentSize[events[j].popID] = events[j].popnSize;` line, add:

```c
				popShape[events[j].popID].type = SHAPE_CONSTANT;
				popShape[events[j].popID].anchor_value = events[j].popnSize;
				popShape[events[j].popID].rate_param = 0.0;
				popShape[events[j].popID].anchor_time = events[j].time;
```

(For an `'n'` event the shape becomes `(CONSTANT, new_value, 0, t_event)` — the rate_param is unused for CONST so it's set to 0, and anchor_time is the event time.)

- [ ] **Step 3: Add include if missing**

Make sure `src/core/discoal_multipop.c` has access to `popShape` and `SHAPE_CONSTANT`. Both are declared in `discoal.h`, which is already included. No new include needed.

- [ ] **Step 4: Build**

Run: `make discoal`
Expected: builds cleanly.

- [ ] **Step 5: Bit-equality regression**

```bash
./test/parity/phase3_bit_equality.sh
```
Expected: PASS for all 8 configurations, including the `single_pop_size_change` config which exercises `-en` events. After the `'n'` event fires, `currentSize[popID]` and `popShape[popID].anchor_value` agree, so subsequent `sizeAt` reads match the prior `currentSize` reads bit-for-bit.

- [ ] **Step 6: Run unit tests**

Run: `make run_tests`
Expected: all 75 tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/core/discoal_multipop.c
git commit -m "Update popShape on 'n' event firing

When -en or other size-change events fire, the new size is now
written to both currentSize[popID] and popShape[popID]. Future
sizeAt() reads at any time after the event get the new value.
Bit-equal regression preserved across the configuration grid."
```

---

## Task 8: Final regression sweep + tag

Run the full regression and unit-test suites end-to-end, document the milestone, and tag the branch.

**Files:** none modified.

- [ ] **Step 1: Full regression**

```bash
make discoal discoal_pre_phase3
./test/parity/phase3_bit_equality.sh
```
Expected: `PASS: all 8 configurations are byte-equal.`

- [ ] **Step 2: Full unit test suite**

```bash
make run_tests
```
Expected: `75 Tests 0 Failures 0 Ignored OK`.

- [ ] **Step 3: Smoke discoal main binary**

```bash
./build/discoal 6 5 1000 -t 5 -r 5 -p 2 3 3 -m 0 1 0.3 -ed 1.0 0 1 -d 12345 67890 | head -20
```
Expected: ms-format output, no crashes, multi-pop simulation completes.

- [ ] **Step 4: Confirm the back-derivation hack still works (regression for unrelated demes path)**

The hack at `discoalFunctions.c:222-261` is unchanged in Phase 3; demes-based runs should still work the same. Run:

```bash
ls config_examples/*.demes.yaml 2>/dev/null | head -1 | xargs -I {} ./build/discoal -Y {} 6 1 100 -t 1 -d 12345 67890 2>&1 | head -10
```

If a demes example exists, it should produce ms-format output; if not, this step is a no-op.

- [ ] **Step 5: Tag the milestone**

```bash
git tag -a phase3-wiring-complete -m "Phase 3 of issue-82 design complete

Shape state (popShape[], migShape[][]) is now seeded by
initializeShapesFromGlobals() at the start of each replicate from
the current parameter globals (currentSize[], migMatConst[][]).
The neutral phase reads sizes via sizeAt() and migrations via
migAt() instead of the sizeRatio[] parameter and migMat[][] global.
The 'n' event handler updates popShape alongside currentSize.

All shapes are SHAPE_CONSTANT in this phase, so the math is
mathematically identical to pre-Phase-3 behavior and the simulation
output is bit-equal across the regression grid (8 configurations
covering single-pop, two-pop with size changes / migration / split,
three-pop combos, admixture, and ancient samples).

The Exp(total_rate) sampler is unchanged. The per-component NHPP
draw is deferred to Phase 4 when SHAPE_EXPONENTIAL is added and
the rate-varies-within-interval problem actually requires it.
Sweep-phase code, importer, new event types, and the back-derivation
hack are unchanged in this phase."
```

- [ ] **Step 6: List tags**

Run: `git tag --list 'phase*'`
Expected: `phase1-shapes-complete`, `phase2-q1-complete`, `phase3-wiring-complete`.

(Push not done yet; user typically pushes when ready.)

---

## Self-Review Checklist (run before declaring this plan complete)

- [ ] Every task has explicit file paths.
- [ ] Every code step contains the actual code an engineer needs.
- [ ] Every test step includes the expected output of `make` / `./test_*`.
- [ ] No "TBD" markers anywhere.
- [ ] Function names match across tasks: `initializeShapesFromGlobals`, `sizeAt`, `migAt`, `popShape`, `migShape`.
- [ ] `Shape` struct fields used everywhere: `type`, `anchor_value`, `rate_param`, `anchor_time`.
- [ ] Commit messages contain no Claude/AI references and no emojis.
- [ ] All work lands on `feature/issue-82-time-varying-demography`.
- [ ] Bit-equality is the verification criterion at every regression checkpoint.

## What's Next After This Plan

Phases 4-7 of the design (per the spec's Phase 5 section):

- **Plan: Phase 4** — Add `SHAPE_EXPONENTIAL` support. This is when the per-component NHPP draw becomes necessary (because rates vary within an interval); the `Exp(total)` sampler can no longer compose across heterogeneous shapes. msprime parity for single-pop and 2-pop exp growth.

- **Plan: Phase 5** — Sweep accessor wiring using the closed-form `detSweepFreqGeneral` per the spec revision (no Euler). Removes the `--det-sweep-mode` flag and `detSweepFreqEuler` function as side effect.

- **Plan: Phase 6** — Add `SHAPE_LINEAR` support, msprime parity.

- **Plan: Phase 7** — Migration shape changes via `'em'` events, importer rewrite, removal of `discoalFunctions.c:222-261` back-derivation hack. **Closes issue #82.**

- **Plan: Phase 8** — Documentation, full-vocabulary parity sweep, CHANGELOG, PR prep.
