#ifndef FENWICK_H
#define FENWICK_H

#include <stdint.h>

/* A Fenwick (binary indexed) tree of non-negative integer weights, for
 * drawing an index with probability proportional to its weight. Setting a
 * weight and finding the index that holds a cumulative position both take
 * O(log n). The weights are integers (numbers of sites), so the sums are
 * exact and do not drift as weights are added and removed. */
typedef struct {
    int64_t *tree;     /* tree[1..capacity]: partial sums */
    int64_t *weight;   /* weight[0..capacity-1]: the weights themselves */
    int64_t total;     /* sum of all weights */
    int capacity;
} Fenwick;

void fenwickInit(Fenwick *f);
void fenwickFree(Fenwick *f);
/* Set every weight to 0; keep the memory. */
void fenwickClear(Fenwick *f);
/* Set the weight of index i (0-based); grows the tree as needed. */
void fenwickSet(Fenwick *f, int i, int64_t w);
int64_t fenwickGet(const Fenwick *f, int i);
static inline int64_t fenwickTotal(const Fenwick *f) { return f->total; }
/* The index i whose cumulative range holds `position`, for 0 <= position <
 * total: the sum of the weights before i is <= position, and with i's
 * weight it is > position. *offset receives position minus the sum of the
 * weights before i, which is uniform on [0, weight(i)) when position is
 * uniform on [0, total). */
int fenwickFind(const Fenwick *f, int64_t position, int64_t *offset);

#endif
