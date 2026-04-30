#!/usr/bin/env python3
"""Compare two tskit tree-sequence files for structural equality.

Used by yaml_validation_suite.sh to parity-check the per-replicate
.trees output of `discoal -Y <fixture>` against the equivalent CLI
invocation declared in the fixture's header.

Provenance is ignored because discoal writes a timestamp into each
tree-sequence's provenance table; without that exclusion every
pair would differ on the timestamp alone. Everything else
(nodes, edges, sites, mutations, populations, sequence_length,
top-level metadata) must match byte-for-byte.

Exit codes:
  0  files are structurally equal
  1  files differ
  2  usage error or load failure (missing file, corrupt store, ...)
"""
import sys
import tskit


def main():
    if len(sys.argv) != 3:
        print("usage: cmp_trees.py <a.trees> <b.trees>", file=sys.stderr)
        return 2
    try:
        a = tskit.load(sys.argv[1]).tables
        b = tskit.load(sys.argv[2]).tables
    except Exception as e:
        print(f"cmp_trees: load failed: {e}", file=sys.stderr)
        return 2
    return 0 if a.equals(b, ignore_provenance=True) else 1


if __name__ == "__main__":
    sys.exit(main())
