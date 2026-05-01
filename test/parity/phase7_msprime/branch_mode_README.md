# Branch-mode parity diagnostic

`test_issue82_branch_mode.py` complements the site-mode test
(`test_issue82_fixture.py`) by computing tskit `mode="branch"` statistics
directly on tree sequences output by both simulators. This bypasses the
mutation process entirely.

## Why both tests

The site-mode test reports a residual divergence at REPS=500: the
demographic backbone (pi, haplotype diversity, n haplotypes) matches, but
mutation-distribution stats (segregating sites, Watterson theta, Tajima's D)
reject Bonferroni-corrected equality. The shape of the divergence (more
singletons in discoal, ~8% larger total ss) is consistent with either:

1. A mutation-placement convention difference (e.g., theta-to-mu conversion).
2. A coalescent simulator difference in tree shape (more terminal branches).

Branch-mode stats discriminate:

- If branch-mode stats MATCH -> the trees match, and the divergence is in mutations (case 1).
- If branch-mode stats DIFFER -> the trees themselves differ (case 2).

## Running

```bash
python3 test/parity/phase7_msprime/test_issue82_branch_mode.py 500
```

## Notes on units

discoal's tree-sequence output uses internal coalescent time units
(1 unit = 4 * Ne_diploid generations = 20000 generations for the issue #82
fixture's Ne_ref=10000 haploid). msprime's tree sequences are in generations.
Branch-mode `diversity` and `segregating_sites` therefore differ by a fixed
scale factor of ~20000.

To compare distribution SHAPE (not absolute values), each replicate's value is
divided by the across-replicate mean before KS comparison. Tajima's D is
scale-invariant (it's a ratio) and is compared directly.
