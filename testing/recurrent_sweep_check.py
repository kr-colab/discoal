#!/usr/bin/env python3
"""Check that recurrent sweeps (-R, -L) follow the demographic model
between sweeps.

With a recurrent-sweep rate so low that no sweep happens, the
recurrent-sweep sampler must give the same output as the neutral
sampler: the same events, drawn from the same random numbers. A model
with migration changes (-eM, -em) checks that the recurrent-sweep sampler
reads the migration rates that these events set.

Needs build/discoal. Exit code 0 if every check passes.
"""

import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DISCOAL = os.path.join(REPO, "build", "discoal")
BASE = ["20", "200", "1000", "-t", "10", "-r", "5", "-p", "2", "10", "10",
        "-ed", "3.0", "1", "0", "-d", "3", "4"]
MODELS = {
    "-eM (all migration rates change)": ["-M", "2.0", "-eM", "0.05", "0.01"],
    "-em (one migration rate changes)": ["-m", "0", "1", "1", "-m", "1", "0",
                                         "0.5", "-em", "0.3", "0", "1", "4"],
    "constant migration": ["-M", "1.0"],
}
NO_SWEEPS = ["-R", "0.000001", "-a", "500"]


def run(args):
    out = subprocess.run([DISCOAL] + args, capture_output=True, text=True,
                         check=True).stdout
    return out.split("\n", 2)[2]  # without the command line and the seeds


def main():
    ok = True
    for name, model in MODELS.items():
        same = run(BASE + model + NO_SWEEPS) == run(BASE + model)
        ok &= same
        print(f"{'PASS' if same else 'FAIL'}  recurrent sweeps at a negligible "
              f"rate give the neutral output: {name}", flush=True)
    print("\nAll checks pass." if ok else "\nOne or more checks FAILED.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
