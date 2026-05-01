# Phase 7: Migration Shape Events + Importer Rewrite (closes issue #82)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close issue #82 by replacing the broken paired-`'M'`-event emission in the demes importer with a proper interval-based `'em'`-event scheme, removing the back-derivation hack at `discoalFunctions.c:222-261`, and adding the CLI surface (`-em`/`-eM`) and runtime handler for the new event type. Also lift the importer's exp/linear-epoch rejection so demes graphs with growth events can be simulated. Validate via msprime parity on the issue #82 fixture.

**Architecture:** Mirrors the Phase 4 `'g'` event pattern for the new `'em'` migration event. The CLI parser emits `'em'` events; the main dispatch handler updates `migShape[src][dst]` to `(SHAPE_CONSTANT, new_rate, 0, t_event)` (matching the `'n'` size-event pattern but for migration). The importer rewrite computes piecewise-constant migration matrices per disjoint time interval implied by the demes graph's migration windows, writes the t=0 active matrix directly to `migMatConst`, and emits `'em'` events at each interval boundary for pairs whose rate changes. Exp/linear demes epochs (`size_function: exponential` or `linear`) are now handled by emitting `'g'`/`'l'` events with the appropriately-converted growth rate. The back-derivation hack at `discoalFunctions.c:222-261` is deleted.

**Tech Stack:** C99, Unity, GNU make, bash, Python+msprime+demes for parity tests.

**Spec:** `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md` §4.4 (event vocabulary), §4.7 (importer changes), §4.9 (what gets removed), §6.2 (msprime parity test design).

**Deliverables:**
- `-em time srcPop dstPop rate` and `-eM time rate` CLI flags emitting `'em'` events.
- `'em'` event handler in the main event-dispatch switch updating `migShape[src][dst]`.
- `-em` smoke test (analogous to `-eg` smoke).
- Importer rewrite: `convertDemesToEvents` in `src/core/demesInterface.c` no longer emits `'M'` events; instead computes piecewise-constant migration matrices and emits `'em'` events per interval boundary. Direct write to `migMatConst[i][j]` for the t=0 interval. Lifts the rejection of exp/linear epochs and emits `'g'`/`'l'` events instead.
- Deletion: the back-derivation hack at `discoalFunctions.c:222-261` (the fprintf-laced scan).
- msprime parity test for the issue #82 fixture (`config_examples/demes_example.demes.yaml`) demonstrating that discoal under the new importer matches msprime running the same demes graph.
- All Phase 3-6 regressions still PASS.
- Tag `phase7-issue-82-closed`.

**Out of scope:**
- Adding shape support (EXP/LIN) for migration rates, beyond constant-rate `'em'` (would extend `'em'` semantics to allow rate_param != 0; not needed to close issue #82 since demes' migration windows are all constant-rate within a window).
- Phase 8 documentation / final docs / CHANGELOG / PR prep.

**Convention reminder:** Never mention Claude or AI in commits/code/docs. Never use emojis. Stay on `feature/issue-82-time-varying-demography`. Do not push (push is a separate manual step at end).

---

## Task 1: CLI flags `-em` and `-eM`

**Files:**
- Modify: `src/core/discoal_multipop.c` (add to `case 'e':` parser dispatch, alongside `case 'g':` / `case 'G':` / `case 'l':` / `case 'L':`)

### Step 1: Find the existing `-eg` / `-el` parser blocks

```bash
grep -n "case 'g':\|case 'l':" src/core/discoal_multipop.c
```

Locate the inner `switch(argv[args][2])` of the outer `case 'e':`. We're adding two more cases.

### Step 2: Add `-em` and `-eM` parsers

Inside the inner switch, after the existing `case 'L':` block, add:

```c
					case 'm':
						/* -em time srcPop dstPop rate — change migration rate at this time */
						{
							ensureEventsCapacity();
							events[eventNumber].time = atof(argv[++args]) * 2.0;
							events[eventNumber].popID2 = atoi(argv[++args]);  /* source pop */
							events[eventNumber].popID = atoi(argv[++args]);   /* destination pop */
							events[eventNumber].popnSize = atof(argv[++args]); /* new rate stored in popnSize field */
							events[eventNumber].type = 'm';  /* lowercase 'm' for time-varying migration event */
							eventNumber++;
						}
						break;
					case 'M':
						/* -eM time rate — set all off-diagonal pairs to the given rate at this time */
						{
							double t = atof(argv[++args]) * 2.0;
							double rate_val = atof(argv[++args]);
							for (int src = 0; src < npops; src++) {
								for (int dst = 0; dst < npops; dst++) {
									if (src == dst) continue;
									ensureEventsCapacity();
									events[eventNumber].time = t;
									events[eventNumber].popID2 = src;
									events[eventNumber].popID = dst;
									events[eventNumber].popnSize = rate_val;
									events[eventNumber].type = 'm';
									eventNumber++;
								}
							}
						}
						break;
```

(Event field convention: `popID` = destination, `popID2` = source — matches the existing 'M' event convention used by the demes importer at line 432-433. The runtime handler will use these to index `migShape[src][dst]`. We use lowercase `'m'` as the event type letter to avoid colliding with the legacy `'M'` events emitted by the unrewritten part of the importer; once Task 5 rewrites the importer, no `'M'` events are emitted any longer.)

### Step 3: Build and smoke-check parsing

```bash
make discoal
./build/discoal 8 1 1000 -t 5 -r 5 -p 2 4 4 -em 0.5 0 1 0.3 -d 12345 67890 | head -5
```

Expected: ms-format output. The simulation may not yet show distinguishable behavior because Task 2 hasn't added the runtime handler — `'m'` events fall through silently. That's fine.

### Step 4: Bit-equality regression

```bash
./test/parity/phase3_bit_equality.sh
./test/parity/phase5_sweep_bit_equality.sh
```

Expected: PASS for both.

### Step 5: Commit

```bash
git add src/core/discoal_multipop.c
git commit -m "$(cat <<'EOF'
Add -em and -eM CLI flags for time-varying migration events

-em time srcPop dstPop rate emits one 'm' event setting
migShape[src][dst] to SHAPE_CONSTANT with the given migration
rate at the given time (in 2N units, multiplied by 2.0 for
internal 4N representation matching -en / -eg / -el convention).

-eM time rate emits one 'm' event per off-diagonal pair, all
with the same rate. -p must come before -eM for npops to be set.

Event field convention: events[].popID is the destination
population, events[].popID2 is the source — matches the legacy
'M' event convention used by the demes importer. The runtime
handler in Task 2 indexes migShape[src][dst] accordingly.

Lowercase 'm' avoids colliding with the legacy 'M' events emitted
by the unrewritten parts of the importer; once Task 5 rewrites
the importer, no 'M' events are emitted any longer.
EOF
)"
```

---

## Task 2: `'m'` event handler in main event-dispatch switch

**Files:**
- Modify: `src/core/discoal_multipop.c`

### Step 1: Find the `'g'` / `'l'` handlers in the main switch

```bash
grep -n "case 'g':\|case 'l':" src/core/discoal_multipop.c | head -5
```

The blocks added in Phase 4 Task 3 and Phase 6 Task 2. The structure: set popShape, then call neutralPhase or sweepPhase.

### Step 2: Add `case 'm':` after `case 'l':`

```c
			case 'm':
				currentTime = events[j].time;
				migShape[events[j].popID2][events[j].popID].type = SHAPE_CONSTANT;
				migShape[events[j].popID2][events[j].popID].anchor_value = events[j].popnSize;
				migShape[events[j].popID2][events[j].popID].rate_param = 0.0;
				migShape[events[j].popID2][events[j].popID].anchor_time = currentTime;
				/* Run the inter-event interval as 'n' / 'g' / 'l' do. */
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

(Note: `events[j].popID2` is source, `events[j].popID` is destination. This matches the CLI parser's storage from Task 1 and the existing 'M' event convention. `migShape[src][dst]` is indexed as `migShape[popID2][popID]`.)

### Step 3: Build and verify

```bash
make discoal
./build/discoal 8 1 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -em 0.5 0 1 0.0 -d 12345 67890 | head -5
```

Expected: ms-format output. The `-em 0.5 0 1 0.0` should turn off the migration from pop 0 to pop 1 at time 0.5 (CLI units).

### Step 4: Bit-equality regression

```bash
./test/parity/phase3_bit_equality.sh
./test/parity/phase5_sweep_bit_equality.sh
make run_tests
```

Expected: all PASS. Regression configs don't use `-em`.

### Step 5: Commit

```bash
git add src/core/discoal_multipop.c
git commit -m "$(cat <<'EOF'
Add 'm' event handler for time-varying migration

When an 'm' event fires, migShape[src][dst] is set to
(SHAPE_CONSTANT, new_rate, 0, t_event). src is events[j].popID2,
dst is events[j].popID — matches the CLI parser convention from
Task 1 and the legacy 'M' event convention.

The case body otherwise mirrors 'g' / 'l' — runs the inter-event
interval through neutralPhase or sweepPhase as appropriate.

The neutral-phase sampler already reads migAt(src, dst, t) which
indexes migShape, so the new rate takes effect immediately for
subsequent inter-event waiting-time draws.
EOF
)"
```

---

## Task 3: `-em` smoke test

**Files:**
- Create: `test/parity/phase7_em_smoke.sh`
- Create: `test/parity/phase7_em_smoke_results/.gitignore`

### Step 1: Write the test

```bash
#!/usr/bin/env bash
# Phase 7 smoke: confirm -em produces distributionally different output
# from comparable constant-migration configs.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"

[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

NREPS=20
configs=(
  "no_migration|8 $NREPS 1000 -t 5 -r 5 -p 2 4 4 -ed 1.0 0 1 -d 12345 67890"
  "constant_migration|8 $NREPS 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -m 1 0 0.5 -ed 1.0 0 1 -d 12345 67890"
  "with_em_off|8 $NREPS 1000 -t 5 -r 5 -p 2 4 4 -m 0 1 0.5 -m 1 0 0.5 -em 0.3 0 1 0.0 -em 0.3 1 0 0.0 -ed 1.0 0 1 -d 12345 67890"
)

OUT="$HERE/phase7_em_smoke_results"
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

if diff <(sed '1d' "$OUT/no_migration.ms") <(sed '1d' "$OUT/constant_migration.ms") > /dev/null; then
  echo "FAIL: no-migration and constant-migration configs are byte-identical (migration appears not to be exercised)"
  exit 1
fi

if diff <(sed '1d' "$OUT/constant_migration.ms") <(sed '1d' "$OUT/with_em_off.ms") > /dev/null; then
  echo "FAIL: constant-migration and -em-off configs are byte-identical (-em appears to be a no-op)"
  exit 1
fi

echo
echo "PASS: -em produces distributionally different output from constant-migration baseline."
```

Make executable.

### Step 2: Add gitignore

```
*.ms
*.err
```

### Step 3: Run

```bash
./test/parity/phase7_em_smoke.sh
```

Expected: PASS.

### Step 4: Commit

```bash
git add test/parity/phase7_em_smoke.sh test/parity/phase7_em_smoke_results/.gitignore
git commit -m "$(cat <<'EOF'
Add -em smoke test to confirm migShape mutation is exercised

Three configs at the same RNG seed: no migration, constant
migration via -m, and -em turning off migration mid-simulation.
The smoke test confirms -em produces distributionally different
output from the constant-migration baseline — sanity check that
the new event type actually mutates migShape and the inner-loop
sampler sees the change.
EOF
)"
```

---

## Task 4: Importer rewrite — migration windows to `'em'` events

**Files:**
- Modify: `src/core/demesInterface.c`

This is the core of the issue #82 fix. The current code at `demesInterface.c:402-443` (or thereabouts) emits paired `'M'` events that the runtime back-derivation hack misuses. Replace with proper interval-based `'em'` event emission.

### Step 1: Inspect the current migration handling

```bash
grep -n "Process migrations\|graph->migrations" src/core/demesInterface.c
```

The block should be around `// Process migrations` followed by a loop over `graph->migrations` and emission of paired 'M' events for each migration.

### Step 2: Algorithm design

For each pair `(src, dst)` of populations, demes guarantees that migration windows are non-overlapping in time. So we can build, per pair, a list of `(start_time_internal, end_time_internal, rate_internal)` tuples, then walk them in time order to compute when the rate changes.

For each migration in `graph->migrations`:

```c
int sourceID = findPopulationIndex(graph, mig->source->name);
int destID = findPopulationIndex(graph, mig->dest->name);
double scaledRate = 4.0 * N * mig->rate;  /* discoal's 4Nm convention */
double t_start_internal = demesTimeToCoalTime(mig->start_time, graph->generation_time, N);
double t_end_internal = demesTimeToCoalTime(mig->end_time, graph->generation_time, N);
```

(`t_start_internal` is the older time; `t_end_internal` is the more recent time. demes "start_time" is when the migration begins going forward, which is the EARLIER (more ancient) discoal time, so `t_start_internal > t_end_internal` in discoal's backward-time convention.)

Algorithm:
1. Initialize `migMatConst[src][dst] = 0` for all (src, dst) pairs.
2. For each migration window `(start_time_internal, end_time_internal, rate)` covering t=0 (i.e., end_time_internal == 0 in discoal time): set `migMatConst[src][dst] += scaledRate`.
3. For each pair (src, dst), sort migration entries by `end_time_internal` ascending (most recent first).
4. For each entry whose `end_time_internal > 0`: emit a `'m'` event at `end_time_internal` setting the rate to whatever value is active going further into the past at that point.
5. For each entry whose `start_time_internal > 0` (window opens going backward): emit a `'m'` event at `start_time_internal` setting the rate to whatever value is active going further into the past at that point.

A simpler implementation: collect ALL boundary times for the pair, sort, walk, and at each boundary compute the active rate by checking which windows cover the time slightly older than the boundary. Emit a `'m'` event with that rate.

```c
/* Per-pair piecewise-constant migration construction */
for each pair (src, dst) where src != dst:
    collect_windows = [];
    for mig in graph->migrations:
        if mig->source matches src and mig->dest matches dst:
            collect_windows.append((t_start_internal, t_end_internal, scaledRate));

    /* Determine the t=0 active rate */
    rate_at_t0 = 0;
    for win in collect_windows:
        if win.t_end_internal == 0:
            rate_at_t0 += win.rate;
    migMatConst[src][dst] = rate_at_t0;

    /* Find boundaries: end_times and start_times, sorted ascending */
    boundaries = sorted({win.t_end_internal for win in windows} | {win.t_start_internal for win in windows});

    for boundary in boundaries:
        if boundary <= 0: continue;
        /* Compute rate just past this boundary (going further into the past) */
        rate_after = sum(win.rate for win in windows if win.t_end_internal <= boundary < win.t_start_internal);
        emit 'm' event at time=boundary with rate=rate_after;
```

(For each unique boundary time, emit one `'m'` event per pair. The runtime processes them in time order.)

### Step 3: Implement

This is the substantive code change. The existing block at `demesInterface.c:402-443` (paired-'M' emission) is fully replaced. Plus the t=0 active matrix is written directly to `migMatConst[src][dst]`.

Use the existing helpers `findPopulationIndex`, `demesTimeToCoalTime`. The factor 4*N for migration rate scaling is the discoal convention (4Nm), already used in the legacy code.

The `events[].popID2` is the source, `events[].popID` is the destination — matches the CLI from Task 1 and the existing 'M' convention (which the legacy importer uses).

For each migration window, the emission logic:
- Symmetric demes migrations (`demes: [A, B]`) produce TWO migration objects in `graph->migrations` (A→B and B→A). We handle each as a unidirectional migration. (The existing code does this.)
- Asymmetric (`source: A, dest: B`) produces ONE migration object.

### Step 4: Build and run unit tests + regressions

```bash
make discoal
make run_tests
./test/parity/phase3_bit_equality.sh
./test/parity/phase5_sweep_bit_equality.sh
```

Expected: all PASS. The regression configs don't use demes import.

### Step 5: Commit

```bash
git add src/core/demesInterface.c
git commit -m "$(cat <<'EOF'
Importer: emit interval-based 'm' events for demes migration windows

Replaces the broken paired-'M' event emission (issue #82) with a
proper interval-based scheme. For each unidirectional migration
window in graph->migrations:

- The active rate at t=0 (the present) is summed across windows
  covering t=0 and written directly to migMatConst[src][dst].
- Windows that end at t > 0 (going backward into the past) trigger
  'm' events at the boundary, setting the rate to whatever is
  active going further back.
- Windows that start at t > 0 trigger 'm' events similarly.

The runtime now correctly tracks time-varying migration rates,
including the issue #82 fixture (multi-window migration where
some demes don't exist across the full timespan).

The legacy 'M'-event back-derivation hack at
discoalFunctions.c:222-261 is no longer needed; removal in Task 6.
EOF
)"
```

---

## Task 5: Importer — emit `'g'`/`'l'` for exp/linear demes epochs

**Files:**
- Modify: `src/core/demesInterface.c`

The current importer rejects `size_function: exponential` and `size_function: linear` epochs (around line 323-345 per the spec). Lift the rejection and emit `'g'`/`'l'` events.

### Step 1: Find the rejection block

```bash
grep -n "exponential\|linear\|size_function" src/core/demesInterface.c
```

### Step 2: Replace rejection with event emission

For each demes epoch with `size_function` other than `CONSTANT`:

- Compute the per-generation forward-time growth rate from `start_size` and `end_size`:
  - For exponential: `alpha_per_gen = log(start_size / end_size) / (start_time - end_time)` (forward in real time the population grew from `end_size` at the more-ancient `start_time` to `start_size` at the more-recent `end_time`; backward, it shrinks).
  - For linear: `gamma_per_gen = (start_size - end_size) / (start_time - end_time)` (similar).
- Convert to the 4N-scaled internal alpha: `alpha_internal = alpha_per_gen * 2 * N` (where N is the haploid effective; matches Phase 4b Convention B's effective haploid = 2 * EFFECTIVE_POPN_SIZE in msprime parity).

  Actually: the conversion math from Phase 4b parity is `alpha_per_gen = alpha_internal / (2 * N_e)`. Inverted: `alpha_internal = alpha_per_gen * 2 * N_e`.

- Emit a 'g' event (for exponential) or 'l' event (for linear) at the more-recent epoch boundary `end_time` (in discoal internal time):

```c
double alpha_per_gen = log(start_size / end_size) / (start_time - end_time);  /* in demes time units */
/* Convert from per-generation to per-internal-time-unit (alpha_internal = alpha_per_gen * 2 * N) */
double alpha_internal = alpha_per_gen * 2.0 * N;
double t_internal = demesTimeToCoalTime(end_time, graph->generation_time, N);

events[eventNumber].time = t_internal;
events[eventNumber].popID = popID;
events[eventNumber].popnSize = alpha_internal;  /* alpha stored in popnSize field, matching -eg semantics */
events[eventNumber].type = 'g';  /* or 'l' for linear */
eventNumber++;
```

Verify the conversion direction empirically: load a demes graph with a known growth pattern, run discoal under the new importer, and compare against the same demes graph in msprime via `msprime.Demography.from_demes()`. If the simulations agree, the conversion is right.

### Step 3: Lift the rejection

Delete the rejection block (originally `if (size_function == EXPONENTIAL) { error; return -1; }` etc.). Replace with the emission described above.

### Step 4: Build and verify

```bash
make discoal
make run_tests
```

Expected: all PASS.

Run a smoke test with a demes file that has exponential growth:

```bash
# Create a minimal demes file with exponential growth
cat > /tmp/exp_growth.yaml <<'YAMLEOF'
description: Test exponential growth
time_units: generations
defaults:
  epoch:
    start_size: 10000
demes:
  - name: pop0
    epochs:
      - start_size: 10000
        end_size: 100000
        end_time: 0
        size_function: exponential
        start_time: 1000000
YAMLEOF
./build/discoal -D /tmp/exp_growth.yaml 6 1 1000 -t 5 -r 5 -d 12345 67890 2>&1 | head -10
```

(Adjust based on actual `-D` flag syntax. May need `-Y` for a wrapping config; verify by reading discoal docs or `getParameters`.)

Expected: discoal accepts the demes file with exponential growth and produces ms-format output.

### Step 5: Commit

```bash
git add src/core/demesInterface.c
git commit -m "$(cat <<'EOF'
Importer: emit 'g'/'l' events for exponential / linear demes epochs

Lifts the rejection of size_function: exponential / linear epochs
in convertDemesToEvents. Computes the per-generation forward-time
growth rate from (start_size, end_size, duration) and converts to
the discoal internal alpha (= alpha_per_gen * 2 * N, where N is
the discoal effective population size).

Emits 'g' events for exponential and 'l' events for linear
epochs at the more-recent boundary (end_time). The neutral and
sweep phases already handle these shapes (Phases 4 and 6).

Demes graphs with growth-mode epochs can now be simulated by
discoal without manual approximation.
EOF
)"
```

---

## Task 6: Remove the back-derivation hack

**Files:**
- Modify: `src/core/discoalFunctions.c` (delete lines 222-261, the fprintf-laced scan)

The hack reconstructs `migMat` at t=0 by scanning all 'M' events and finding the most recent one for each pair. Since the importer now writes `migMatConst` directly for the t=0 interval (Task 4), this scan is unnecessary AND wrong (the new importer no longer emits 'M' events, only 'm' events).

### Step 1: Find and delete the block

```bash
sed -n '218,265p' src/core/discoalFunctions.c
```

The block spans approximately lines 222-261 with comments like "Process migration events from demes" and a loop scanning `events[i].type == 'M'`. Delete the entire block, leaving the surrounding code (which initializes migMat from migMatConst, lines ~218-225) intact.

The result: `initialize()` continues to set `migMat[i][j] = migMatConst[i][j]` from the parser-time global, then proceeds to the rest of the per-rep init. No event scanning, no debug fprintf.

### Step 2: Build and run regressions

```bash
make discoal
./test/parity/phase3_bit_equality.sh
./test/parity/phase5_sweep_bit_equality.sh
./test/parity/phase4_eg_smoke.sh
./test/parity/phase5_sweep_exp_smoke.sh
./test/parity/phase6_el_smoke.sh
./test/parity/phase7_em_smoke.sh
make run_tests
```

Expected: all PASS. The regression configs don't use demes import (they use CLI flags), so the hack removal doesn't affect them.

### Step 3: Commit

```bash
git add src/core/discoalFunctions.c
git commit -m "$(cat <<'EOF'
Remove back-derivation hack from initialize()

The hack at discoalFunctions.c:222-261 scanned 'M' events to
reconstruct the t=0 migration matrix, working around the
importer's broken paired-'M' emission (issue #82). With Task 4's
importer rewrite, migMatConst[src][dst] is written directly at
parser time for the t=0 active matrix; the runtime no longer
needs to back-derive anything.

Deletes the entire fprintf-laced scanning block. initialize()
now does only the standard migMat = migMatConst copy and
proceeds.

Closes the structural part of issue #82. The msprime parity
test in Task 7 confirms behavioral closure on the issue's
fixture.
EOF
)"
```

---

## Task 7: msprime parity test for the issue #82 fixture

**Files:**
- Create: `test/parity/phase7_msprime/test_issue82_fixture.py`
- Modify: `test/parity/phase7_msprime/parity_utils.py` (or reuse Phase 4b's via import)

### Step 1: Build the test

The issue #82 fixture is `config_examples/demes_example.demes.yaml`. Run discoal under that demes file and msprime under the same demes graph; compare summary statistics via Bonferroni-corrected KS.

```python
"""test_issue82_fixture.py — multi-window migration parity (closes issue #82).

Runs the issue #82 demes fixture through both discoal (via -D / -Y demes
import) and msprime (via Demography.from_demes), then compares summary
statistics.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import demes
import msprime
import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent.parent / "phase4b_msprime"))
from parity_utils import (
    ACTIVE_CONVENTION,
    DEFAULT_NE,
    collect_stats,
    msprime_to_ms,
    parse_ms,
    run_discoal,
    validate_active_convention,
)


def main(reps: int = 1000):
    print(f"=== Phase 7: issue #82 fixture parity (REPS={reps}) ===")
    validate_active_convention(reps=200)

    fixture = Path(__file__).parent.parent.parent.parent / "config_examples" / "demes_example.demes.yaml"
    if not fixture.exists():
        raise RuntimeError(f"Fixture not found: {fixture}")
    print(f"Fixture: {fixture}")

    /* discoal run via -D flag (or -Y; check discoal docs for the right flag). */
    /* The fixture has 3 pops; sample size and seeds determined by what works.
     * Verify the right CLI form via: `./build/discoal --help` or by inspecting
     * config_examples/demes_example.yaml for an `# Equivalent to: discoal …` header. */

    /* Convert to discoal command line: */
    /* TODO: Determine n, theta, rho, nsites for the fixture. The demes fixture
     * specifies population sizes and times but not sample size or theta; these
     * must come from outside (likely via the YAML config or from a default in
     * the test). */

    n_total = 9   # 3 pops, 3 samples each
    L = 10000
    theta = 10.0
    rho = 10.0

    discoal_args = [
        str(n_total), str(reps), str(L),
        "-t", str(theta), "-r", str(rho),
        "-D", str(fixture),
        "-d", "12345", "67890",
    ]
    print(f"discoal: {' '.join(['./build/discoal'] + discoal_args)}")
    t0 = time.time()
    discoal_out = run_discoal(discoal_args)
    print(f"  discoal: {time.time()-t0:.1f}s")
    discoal_reps = parse_ms(discoal_out)

    /* msprime run via Demography.from_demes */
    print("msprime: matched demography from demes graph")
    graph = demes.load(str(fixture))
    demography = msprime.Demography.from_demes(graph)
    /* Determine sample sizes per pop matching discoal */
    samples = {}
    for pop in demography.populations:
        if pop.name in {"A", "B", "C"}:  /* fixture-specific pop names */
            samples[pop.name] = 3 / ACTIVE_CONVENTION.ploidy   /* 3 haploid samples / ploidy */
    /* Convert mu and r */
    /* Note: discoal's theta/rho conventions need conversion, similar to Phase 4b.
     * Use the same parity_utils helpers. */
    ...

    msp_reps = []
    t0 = time.time()
    for seed in range(1, reps + 1):
        ts = msprime.sim_ancestry(
            samples=samples,
            sequence_length=L,
            recombination_rate=rho / (4 * DEFAULT_NE * L),
            demography=demography,
            ploidy=ACTIVE_CONVENTION.ploidy,
            random_seed=seed,
        )
        ts = msprime.sim_mutations(
            ts, rate=theta / (4 * DEFAULT_NE * L),
            model=msprime.BinaryMutationModel(),
            discrete_genome=False,
            random_seed=seed,
        )
        msp_reps.append(msprime_to_ms(ts))
    print(f"  msprime: {time.time()-t0:.1f}s")

    s_discoal = collect_stats(discoal_reps, n_total)
    s_msp = collect_stats(msp_reps, n_total)

    /* Bonferroni-corrected KS on ss/pi/td/wtheta/hapdiv/nhap. */
    comparisons = []
    for stat_name in ("ss", "pi", "td", "wtheta", "hapdiv", "nhap"):
        D, p = stats.ks_2samp(s_discoal[stat_name], s_msp[stat_name])
        comparisons.append((f"issue82_{stat_name}", D, p))
    n_comp = len(comparisons)
    bonf_alpha = 0.01 / n_comp
    rejected = sum(1 for _, _, p in comparisons if p < bonf_alpha)

    print(f"\n{n_comp} comparisons, Bonferroni alpha = {bonf_alpha:.2e}")
    for name, D, p in comparisons:
        flag = "  **" if p < bonf_alpha else ""
        print(f"  {name:<24} {D:>10.4f} {p:>10.2e}{flag}")

    if rejected == 0:
        print(f"\nPASS: all {n_comp} comparisons within Bonferroni-corrected p > {bonf_alpha:.2e}")
        sys.exit(0)
    else:
        print(f"\nFAIL: {rejected}/{n_comp} comparisons reject equality")
        sys.exit(1)


if __name__ == "__main__":
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    main(reps=reps)
```

(The `...` blocks are intentional placeholders to be filled in based on inspection of `config_examples/demes_example.yaml` to determine the right discoal CLI form — `-D` vs `-Y`, sample size, theta/rho values. The implementer should resolve these from the actual fixture content.)

### Step 2: Smoke and full runs

```bash
python3 test/parity/phase7_msprime/test_issue82_fixture.py 100
python3 test/parity/phase7_msprime/test_issue82_fixture.py 1000
```

Expected: PASS. The smallest p across the comparison table should be > Bonferroni alpha (~1.7e-3 for 6 comparisons).

If FAIL: the importer-rewrite or the back-derivation removal may have introduced a bug. Investigate: is the discoal output reasonable (sane segsites, pi)? Does the msprime configuration match the discoal one? Are the migration rates being applied at the correct times?

### Step 3: Commit

```bash
git add test/parity/phase7_msprime/
git commit -m "$(cat <<'EOF'
Add msprime parity test for the issue #82 fixture

Loads config_examples/demes_example.demes.yaml in both discoal
(via the -D demes import) and msprime (via Demography.from_demes),
runs 1000 replicates each, and compares summary statistics via
Bonferroni-corrected KS at alpha = 0.01.

This is the structural validation that issue #82 is closed:
the fixture has multi-window migration that the legacy importer
mishandled, producing simulation output inconsistent with the
demes spec. With the Phase 7 importer rewrite (interval-based
'm' events, direct migMatConst writes for t=0, removed back-
derivation hack), discoal under the fixture should match
msprime running the same demes graph.
EOF
)"
```

---

## Task 8: Final regression sweep + tag

- [ ] **Step 1: Full regression**

```bash
make discoal discoal_pre_phase3 discoal_pre_phase5
./test/parity/phase3_bit_equality.sh
./test/parity/phase4_eg_smoke.sh
./test/parity/phase5_sweep_bit_equality.sh
./test/parity/phase5_sweep_exp_smoke.sh
./test/parity/phase6_el_smoke.sh
./test/parity/phase7_em_smoke.sh
python3 test/parity/phase7_msprime/test_issue82_fixture.py 1000
make run_tests
```

Expected: all PASS.

- [ ] **Step 2: Tag**

```bash
git tag -a phase7-issue-82-closed -m "$(cat <<'TAGEOF'
Phase 7 of issue-82 design complete: closes issue #82

The demes importer no longer emits broken paired-'M' migration
events. Instead, convertDemesToEvents:

- Computes per-pair piecewise-constant migration matrices over
  the disjoint time intervals implied by graph->migrations.
- Writes the t=0 active matrix directly to migMatConst[src][dst].
- Emits one 'm' event per pair per interval boundary (t > 0)
  where the rate changes.
- Lifts the rejection of exponential / linear demes epochs and
  emits 'g' / 'l' events instead.

The runtime back-derivation hack at discoalFunctions.c:222-261
is removed; the runtime simply copies migMat from migMatConst
at sim start as before, with no event scanning.

CLI surface added: -em time srcPop dstPop rate (single pair)
and -eM time rate (matrix-wide). 'm' event handler in the main
event-dispatch switch updates migShape[src][dst].

msprime parity test on the issue #82 fixture
(config_examples/demes_example.demes.yaml) PASSes Bonferroni-
corrected KS at alpha = 0.01 on summary statistics, validating
that discoal under the new importer matches msprime running
the same demes graph.

Out of scope for Phase 7: SHAPE_EXP / SHAPE_LIN migration shapes
(constant-rate 'm' events suffice for demes' migration windows
which are constant-rate within their bounds), Phase 8 docs.
TAGEOF
)"
```

- [ ] **Step 3: List tags**

Expected: 8 tags including `phase7-issue-82-closed`.

(Push not done yet; user pushes when ready.)

---

## Self-Review Checklist

- [ ] All tasks have explicit file paths.
- [ ] CLI flags `-em` / `-eM` mirror `-eg` / `-eG` exactly.
- [ ] `'m'` event handler mirrors `'g'` / `'l'` with `migShape` instead of `popShape`.
- [ ] Importer no longer emits paired `'M'` events.
- [ ] Importer writes `migMatConst[src][dst]` directly for the t=0 interval.
- [ ] Importer emits `'g'` / `'l'` events for exp / linear epochs.
- [ ] Back-derivation hack at `discoalFunctions.c:222-261` deleted.
- [ ] Bit-equality preserved for SHAPE_CONSTANT regression configs (Phase 3 + Phase 5).
- [ ] msprime parity passes for the issue #82 fixture.
- [ ] No emojis, no Claude/AI references.

## What's Next After This Plan

- **Plan: Phase 8** — Documentation (population_structure.rst, demes_integration.md, yaml_configuration.md), full-vocabulary parity sweep across all shape types, CHANGELOG, PR prep. May include msprime parity for LINEAR via discretization if time permits. **This is the final plan in the issue #82 work.**
