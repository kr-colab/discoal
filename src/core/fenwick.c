#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fenwick.h"

void fenwickInit(Fenwick *f) {
    f->tree = NULL;
    f->weight = NULL;
    f->total = 0;
    f->capacity = 0;
}

void fenwickFree(Fenwick *f) {
    free(f->tree);
    free(f->weight);
    fenwickInit(f);
}

void fenwickClear(Fenwick *f) {
    if (f->capacity > 0) {
        memset(f->tree, 0, (f->capacity + 1) * sizeof(int64_t));
        memset(f->weight, 0, f->capacity * sizeof(int64_t));
    }
    f->total = 0;
}

/* Grow to at least `need` entries. The new tree is rebuilt from the
 * weights, in O(capacity). */
static void fenwickGrow(Fenwick *f, int need) {
    int cap = f->capacity ? f->capacity : 16;
    while (cap < need) cap *= 2;
    int64_t *weight = realloc(f->weight, cap * sizeof(int64_t));
    int64_t *tree = weight ? realloc(f->tree, (cap + 1) * sizeof(int64_t)) : NULL;
    if (!weight || !tree) {
        fprintf(stderr, "Error: out of memory for recombination weights\n");
        exit(1);
    }
    memset(weight + f->capacity, 0, (cap - f->capacity) * sizeof(int64_t));
    f->weight = weight;
    f->tree = tree;
    f->capacity = cap;
    /* Rebuild: each node adds itself to its parent. */
    memset(f->tree, 0, (cap + 1) * sizeof(int64_t));
    for (int i = 1; i <= cap; i++) {
        f->tree[i] += f->weight[i - 1];
        int parent = i + (i & -i);
        if (parent <= cap) f->tree[parent] += f->tree[i];
    }
}

void fenwickSet(Fenwick *f, int i, int64_t w) {
    if (i >= f->capacity) fenwickGrow(f, i + 1);
    int64_t delta = w - f->weight[i];
    if (delta == 0) return;
    f->weight[i] = w;
    f->total += delta;
    for (int k = i + 1; k <= f->capacity; k += k & -k) f->tree[k] += delta;
}

int64_t fenwickGet(const Fenwick *f, int i) {
    return (i < f->capacity) ? f->weight[i] : 0;
}

int fenwickFind(const Fenwick *f, int64_t position, int64_t *offset) {
    /* Descend from the largest power of two: k is the number of leading
     * entries whose sum is <= position. */
    int k = 0;
    int step = 1;
    while (step * 2 <= f->capacity) step *= 2;
    for (; step > 0; step /= 2) {
        if (k + step <= f->capacity && f->tree[k + step] <= position) {
            k += step;
            position -= f->tree[k];
        }
    }
    if (offset) *offset = position;
    return k;  /* 0-based index of the entry that holds the position */
}
