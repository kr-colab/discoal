# Phase 6: SHAPE_LINEAR Engine Wiring

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Wire `SHAPE_LINEAR` into the simulation engine via CLI flags `-el time popID gamma` and `-eL time gamma`, plus the `'l'` event type that updates `popShape[popID]` to a linear shape at the event time. Bit-equality preserved for `SHAPE_CONSTANT`-only configs. Smoke test confirms `-el` produces output distributionally distinct from `-en` (constant size change) and `-eg` (exponential growth).

**Architecture:** Mirrors the Phase 4 `'g'` event pattern. The CLI parser emits `'l'` events; the main dispatch handler sets `popShape[popID]` to `(SHAPE_LINEAR, sizeAt(popID, t), gamma, t)` (Convention 1 continuity); `proposeTrajectory`'s events walk also tracks `'l'` events when iterating forward. The neutral phase, sweep phase, and accessor functions already handle `SHAPE_LINEAR` (math primitives wired in Phase 1; engine accessors use shapes uniformly).

**Tech Stack:** C99, Unity, GNU make, bash for parity tests.

**Spec:** `docs/superpowers/specs/2026-04-30-time-varying-demographic-parameters-design.md` §3.2 (LINEAR math), §4.3 (`'l'` event), §4.8 (CLI surface — note: the spec acknowledges linear may stay YAML-only at first cut; this plan adds CLI flags for parity with -eg/-eG and to exercise the engine path).

**Deliverables:**
- `-el` and `-eL` CLI flags in `discoal_multipop.c`.
- `'l'` event handler in the main event-dispatch switch.
- `'l'` handling in `proposeTrajectory`'s events walk.
- `test/parity/phase6_el_smoke.sh` confirming `-el` is exercised.
- All Phase 3/4/5 regressions still PASS.
- Tag `phase6-linear-complete`.

**Out of scope:**
- msprime parity for LINEAR (msprime lacks native linear growth; parity would require discretization approximation, deferred to Phase 8 if needed).
- YAML/demes import for linear epochs (Phase 7 — importer rewrite handles all shapes including LINEAR).
- Removal of the back-derivation hack (Phase 7).

**Convention reminder:** Never mention Claude or AI in commits/code/docs. Never use emojis. Stay on `feature/issue-82-time-varying-demography`. Do not push (push is a separate manual step at end).

---

## Task 1: CLI flags `-el` and `-eL`

**Files:**
- Modify: `src/core/discoal_multipop.c` (add to `case 'e':` parser dispatch, alongside `case 'g':` / `case 'G':`)

### Step 1: Find the existing `-eg` parser

```bash
grep -n "case 'g':" src/core/discoal_multipop.c
```

Expected: shows the parser block added in Phase 4 Task 4. The structure is:

```c
					case 'g':
						/* -eg time popID alpha */
						ensureEventsCapacity();
						events[eventNumber].time = atof(argv[++args]) * 2.0;
						events[eventNumber].popID = atoi(argv[++args]);
						events[eventNumber].popnSize = atof(argv[++args]);
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

### Step 2: Add `-el` and `-eL` parsers immediately after

Add inside the inner `switch(argv[args][2])`, after the `case 'G':` block:

```c
					case 'l':
						/* -el time popID gamma — linear-growth event */
						ensureEventsCapacity();
						events[eventNumber].time = atof(argv[++args]) * 2.0;
						events[eventNumber].popID = atoi(argv[++args]);
						events[eventNumber].popnSize = atof(argv[++args]);  /* gamma stored in popnSize field */
						events[eventNumber].type = 'l';
						eventNumber++;
						break;
					case 'L':
						/* -eL time gamma — applies to all populations */
						{
							double t = atof(argv[++args]) * 2.0;
							double gamma_val = atof(argv[++args]);
							for (int p = 0; p < npops; p++) {
								ensureEventsCapacity();
								events[eventNumber].time = t;
								events[eventNumber].popID = p;
								events[eventNumber].popnSize = gamma_val;
								events[eventNumber].type = 'l';
								eventNumber++;
							}
						}
						break;
```

(Same time multiplier `* 2.0` as `-eg` — CLI is in 2N units, internal in 4N units. `gamma` is the per-generation forward-time linear growth rate, scaled identically to alpha in the EXP case.)

### Step 3: Build and smoke-check parsing

```bash
make discoal
./build/discoal 6 1 1000 -t 5 -r 5 -el 0.5 0 0.1 -d 12345 67890 | head -5
```

Expected: ms-format output. Don't worry about the simulation result yet — Task 3 wires the runtime handler. At this point the `'l'` event is parsed and added to `events[]`, but the main dispatch falls through to nothing (Task 3 adds the case). The event MIGHT be silently ignored; the simulation should still complete.

If discoal crashes or hangs: the `'l'` event with no main-dispatch case may cause undefined behavior. Skip this smoke step and proceed to Task 2 / Task 3 to add the handler before testing.

### Step 4: Bit-equality regression

```bash
./test/parity/phase3_bit_equality.sh
./test/parity/phase5_sweep_bit_equality.sh
```

Expected: PASS for both. None of the regression configs use `-el`, so the new parser doesn't affect them.

### Step 5: Commit

```bash
git add src/core/discoal_multipop.c
git commit -m "$(cat <<'EOF'
Add -el and -eL CLI flags for linear-growth events

-el time popID gamma emits one 'l' event setting popShape[popID]
to SHAPE_LINEAR with the given linear growth rate at the given
time (in 2N units, multiplied by 2.0 for internal 4N representation
matching -en / -eg convention).

-eL time gamma emits one 'l' event per active population. -p must
come before -eL for npops to be set.

gamma is the per-generation forward-time linear rate of size
change, msprime/forward-time convention; positive gamma means
the population grew forward (was smaller in the past).

The 'l' event handler in the main event-dispatch switch is
added in Task 2; until then 'l' events are parsed and added
to events[] but silently fall through.
EOF
)"
```

---

## Task 2: `'l'` event handler in main event-dispatch switch

**Files:**
- Modify: `src/core/discoal_multipop.c` (add `case 'l':` alongside `case 'g':` in the main event loop around line 230)

### Step 1: Find the `'g'` handler in the main switch

```bash
grep -n "case 'g':" src/core/discoal_multipop.c
```

The relevant block (added in Phase 4 Task 3):

```c
			case 'g':
				currentTime = events[j].time;
				popShape[events[j].popID].type = SHAPE_EXPONENTIAL;
				popShape[events[j].popID].anchor_value = sizeAt(events[j].popID, currentTime);
				popShape[events[j].popID].rate_param = events[j].popnSize;
				popShape[events[j].popID].anchor_time = currentTime;
				/* Run the inter-event interval as 'n' does. */
				if(activeSweepFlag == 0){
					...
				}
				else{
					...
				}
				break;
```

### Step 2: Add `case 'l':` immediately after

Mirror the `'g'` block but use `SHAPE_LINEAR`:

```c
			case 'l':
				currentTime = events[j].time;
				popShape[events[j].popID].type = SHAPE_LINEAR;
				popShape[events[j].popID].anchor_value = sizeAt(events[j].popID, currentTime);
				popShape[events[j].popID].rate_param = events[j].popnSize;
				popShape[events[j].popID].anchor_time = currentTime;
				/* Run the inter-event interval as 'n' / 'g' does. */
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

(Verbatim copy of the `'g'` body, only the `popShape[*].type` line differs.)

### Step 3: Build and verify

```bash
make discoal
./build/discoal 6 1 1000 -t 5 -r 5 -el 0.5 0 0.1 -d 12345 67890 | head -5
```

Expected: ms-format output, simulation completes.

### Step 4: Bit-equality regression

```bash
./test/parity/phase3_bit_equality.sh
./test/parity/phase5_sweep_bit_equality.sh
make run_tests
```

Expected: all PASS.

### Step 5: Commit

```bash
git add src/core/discoal_multipop.c
git commit -m "$(cat <<'EOF'
Add 'l' event handler for SHAPE_LINEAR transitions

When an 'l' event fires, popShape[popID] is set to (LINEAR,
sizeAt(popID, t), gamma, t). Continuity preserved across the
shape boundary: anchor_value uses sizeAt at the event time,
matching the 'g'-event pattern from Phase 4.

The case body otherwise mirrors 'g' verbatim — runs the
inter-event interval through neutralPhase or sweepPhase as
appropriate based on activeSweepFlag and recurSweepMode.
EOF
)"
```

---

## Task 3: `'l'` handling in `proposeTrajectory` events walk

**Files:**
- Modify: `src/core/discoalFunctions.c` (events walk in `proposeTrajectory`)

### Step 1: Find the existing `'g'` handling

```bash
awk '/^double proposeTrajectory/,/^double sweepPhaseEvents/' src/core/discoalFunctions.c | grep -n "type == '"
```

The walk contains (after Phase 5 Task 5):

```c
		if(events[i].type == 'n'){
			popShape[events[i].popID].type = SHAPE_CONSTANT;
			popShape[events[i].popID].anchor_value = events[i].popnSize;
			popShape[events[i].popID].rate_param = 0.0;
			popShape[events[i].popID].anchor_time = events[i].time;
			currentSizeRatio = events[i].popnSize;
			N = floor(N_0 * events[i].popnSize);
			if(currentSizeRatio > Nmax) Nmax = currentSizeRatio;
		}
		if(events[i].type == 'g'){
			popShape[events[i].popID].type = SHAPE_EXPONENTIAL;
			popShape[events[i].popID].anchor_value = sizeAt(events[i].popID, events[i].time);
			popShape[events[i].popID].rate_param = events[i].popnSize;
			popShape[events[i].popID].anchor_time = events[i].time;
			double sr_now = sizeAt(events[i].popID, events[i].time);
			if(sr_now > Nmax) Nmax = sr_now;
		}
```

### Step 2: Add `'l'` handling after `'g'`

```c
		if(events[i].type == 'l'){
			popShape[events[i].popID].type = SHAPE_LINEAR;
			popShape[events[i].popID].anchor_value = sizeAt(events[i].popID, events[i].time);
			popShape[events[i].popID].rate_param = events[i].popnSize;
			popShape[events[i].popID].anchor_time = events[i].time;
			double sr_now = sizeAt(events[i].popID, events[i].time);
			if(sr_now > Nmax) Nmax = sr_now;
		}
```

(Identical structure to `'g'`, only the `.type` differs.)

### Step 3: Build and verify

```bash
make discoal
./build/discoal 6 1 1000 -t 5 -r 5 -wd 0.05 -a 200 -x 0.5 -el 0.5 0 0.1 -d 12345 67890 | head -5
```

Expected: ms-format output, simulation completes (sweep + LINEAR shape).

### Step 4: Bit-equality regression

```bash
./test/parity/phase3_bit_equality.sh
./test/parity/phase5_sweep_bit_equality.sh
make run_tests
```

Expected: all PASS. None of the regression configs use `-el`, so the new walk handler doesn't affect them.

### Step 5: Commit

```bash
git add src/core/discoalFunctions.c
git commit -m "$(cat <<'EOF'
Track 'l' events in proposeTrajectory events walk

Mirrors the 'g'-event tracking added in Phase 5 Task 5. When
proposeTrajectory walks future events to predict the trajectory
shape, 'l' events set popShape to SHAPE_LINEAR. The same
save/restore pattern (memcpy at entry/exit) keeps the global
state unchanged for the caller.
EOF
)"
```

---

## Task 4: `-el` smoke test

**Files:**
- Create: `test/parity/phase6_el_smoke.sh`

### Step 1: Write the test

Confirm `-el` produces distributionally different output from baseline (no demography), `-en` (constant size change), and `-eg` (exponential growth). Same pattern as `test/parity/phase4_eg_smoke.sh`.

```bash
#!/usr/bin/env bash
# Phase 6 smoke: confirm -el produces distributionally different output
# from comparable -en and -eg configurations. Sanity check that the
# SHAPE_LINEAR path is exercised.

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
DISCOAL="$ROOT/build/discoal"

[[ -x "$DISCOAL" ]] || { echo "FAIL: build/discoal missing"; exit 1; }

NREPS=20
configs=(
  "no_demography|6 $NREPS 1000 -t 5 -r 5 -d 12345 67890"
  "with_el_gamma_pos|6 $NREPS 1000 -t 5 -r 5 -el 0.5 0 0.5 -d 12345 67890"
  "with_en_size_change|6 $NREPS 1000 -t 5 -r 5 -en 0.5 0 0.5 -d 12345 67890"
  "with_eg_alpha|6 $NREPS 1000 -t 5 -r 5 -eg 0.5 0 50 -d 12345 67890"
)

OUT="$HERE/phase6_el_smoke_results"
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

# Confirm with_el_gamma_pos differs from no_demography
if diff <(sed '1d' "$OUT/no_demography.ms") <(sed '1d' "$OUT/with_el_gamma_pos.ms") > /dev/null; then
  echo "FAIL: -el config produced byte-identical output to no-demography (el appears to be a no-op)"
  exit 1
fi

# Confirm with_el_gamma_pos differs from with_en_size_change
if diff <(sed '1d' "$OUT/with_el_gamma_pos.ms") <(sed '1d' "$OUT/with_en_size_change.ms") > /dev/null; then
  echo "FAIL: -el and -en configs produced byte-identical output (el may be silently treated as en)"
  exit 1
fi

# Confirm with_el_gamma_pos differs from with_eg_alpha
if diff <(sed '1d' "$OUT/with_el_gamma_pos.ms") <(sed '1d' "$OUT/with_eg_alpha.ms") > /dev/null; then
  echo "FAIL: -el and -eg configs produced byte-identical output (el may be silently treated as eg)"
  exit 1
fi

echo
echo "PASS: -el produces distributionally different output from no-demography, -en, and -eg."
```

Make executable.

### Step 2: Add `.gitignore`

Create `test/parity/phase6_el_smoke_results/.gitignore`:

```
*.ms
*.err
```

### Step 3: Run

```bash
./test/parity/phase6_el_smoke.sh
```

Expected: PASS. The three diffs should each show substantive differences (different segsites counts, different positions, different haplotypes). If any pair is byte-identical, investigate which step (CLI parsing, event handling, runtime dispatch) is failing.

### Step 4: Commit

```bash
git add test/parity/phase6_el_smoke.sh test/parity/phase6_el_smoke_results/.gitignore
git commit -m "$(cat <<'EOF'
Add -el smoke test to confirm SHAPE_LINEAR is exercised

Four configs at the same RNG seed: no demography, -el gamma=0.5
at t=0.5, -en size change at t=0.5, and -eg alpha=50 at t=0.5.
The smoke test confirms -el produces distributionally different
output from all three reference configs — sanity check that the
LINEAR shape path is wired into the simulation engine and not
silently treated as constant or exponential.
EOF
)"
```

---

## Task 5: Final regression sweep + tag

- [ ] **Step 1: Full regression**

```bash
make discoal discoal_pre_phase3 discoal_pre_phase5
./test/parity/phase3_bit_equality.sh
./test/parity/phase4_eg_smoke.sh
./test/parity/phase5_sweep_bit_equality.sh
./test/parity/phase5_sweep_exp_smoke.sh
./test/parity/phase6_el_smoke.sh
make run_tests
```

Expected: all PASS.

- [ ] **Step 2: Tag**

```bash
git tag -a phase6-linear-complete -m "$(cat <<'TAGEOF'
Phase 6 of issue-82 design complete

SHAPE_LINEAR is now exercised end-to-end in the simulation engine:

- CLI flags -el time popID gamma and -eL time gamma emit 'l'
  events.
- Main event-dispatch handler sets popShape[popID] to (LINEAR,
  sizeAt(popID, t), gamma, t). Continuity preserved across the
  shape boundary.
- proposeTrajectory tracks 'l' events alongside 'n' and 'g'.
- The neutral phase, sweep phase, and accessor functions already
  handled SHAPE_LINEAR via the math primitives wired in Phase 1
  and the engine accessor pattern from Phases 3-5.

Smoke test confirms -el output differs distributionally from
baseline, -en, and -eg — the LINEAR path is genuinely exercised.

All Phase 3/4/5 regressions still PASS (bit-equality preserved
for SHAPE_CONSTANT-only configs).

Out of scope for Phase 6: msprime parity for LINEAR (msprime
lacks native linear growth; would require discretization
approximation), YAML/demes import for LINEAR (Phase 7).
TAGEOF
)"
```

- [ ] **Step 3: List tags**

```bash
git tag --list 'phase*'
```

Expected: 7 tags including `phase6-linear-complete`.

(Push not done yet; user pushes when ready.)

---

## Self-Review Checklist

- [ ] All tasks have explicit file paths.
- [ ] CLI flags `-el` and `-eL` mirror `-eg` / `-eG` exactly except for the event type letter.
- [ ] `'l'` event handler mirrors `'g'` exactly except for the shape type.
- [ ] proposeTrajectory tracks `'l'` events with the same structure as `'g'`.
- [ ] Bit-equality preserved for SHAPE_CONSTANT regression configs.
- [ ] Smoke test verifies `-el` differs from `-en`, `-eg`, AND no-demography.
- [ ] No emojis, no Claude/AI references.

## What's Next After This Plan

- **Plan: Phase 7** — `'em'` migration shape events, importer rewrite, removal of back-derivation hack. **Closes issue #82.**
- **Plan: Phase 8** — Documentation, full-vocabulary parity sweep, CHANGELOG, PR prep. May include msprime parity for LINEAR via discretization if time permits.
