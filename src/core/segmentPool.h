#ifndef SEGMENT_POOL_H
#define SEGMENT_POOL_H

#include <stddef.h>
#include <stdio.h>
#include "ancestrySegment.h"

// Pool statistics structure
typedef struct {
    size_t total_allocations;
    size_t total_frees;
    size_t current_in_use;
    size_t high_water_mark;
    size_t pool_capacity;
    size_t resets;
    size_t expansions;
} PoolStats;

// Initialize the global segment pool
// initial_size: number of segments to pre-allocate (0 for default)
void initSegmentPool(size_t initial_size);

// Destroy the pool and free all memory
void destroySegmentPool(void);

// Reset pool for next replicate (O(1) operation)
// This is much faster than freeing individual segments
void resetSegmentPool(void);

// Allocate a segment from the pool
// Returns zeroed segment ready for use
AncestrySegment* allocSegmentFromPool(void);

// Return a segment to the pool
void freeSegmentToPool(AncestrySegment *seg);

// Calculate recommended initial pool size based on simulation parameters
size_t calculateInitialPoolSize(int sampleSize, double rho, int nSites);

// Get current pool statistics
PoolStats getPoolStatistics(void);

// Print pool statistics to file
void printPoolStatistics(FILE *fp);

// Get current memory usage in bytes
size_t getPoolMemoryUsage(void);

// Optional: trim pool to reduce memory usage
void trimSegmentPool(size_t target_size);

// Debug support
#ifdef DEBUG_SEGMENT_POOL
#define SEGMENT_MAGIC 0xDEADBEEF
void validateSegmentPool(void);
#endif

#endif // SEGMENT_POOL_H