#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "segmentPool.h"
#include "ancestrySegmentAVL.h"

#define DEFAULT_POOL_SIZE 10000
#define MIN_POOL_SIZE 1000
#define MAX_POOL_SIZE 10000000
#define GROWTH_FACTOR 2

// Global pool instance (single-threaded application)
typedef struct SegmentPool {
    // Memory blocks
    AncestrySegment **blocks;     // Array of allocated blocks
    size_t *block_sizes;          // Size of each block
    size_t num_blocks;            // Number of blocks allocated
    size_t blocks_capacity;       // Capacity of blocks array
    
    // Free list (LIFO for cache locality)
    AncestrySegment **free_list;  // Stack of available segments
    size_t free_count;            // Number of segments in free list
    size_t free_capacity;         // Capacity of free list
    
    // All segments tracking (for reset operation)
    size_t total_segments;        // Total segments allocated across all blocks
    
    // Statistics
    PoolStats stats;
    
    // State
    int initialized;
} SegmentPool;

static SegmentPool g_pool = {0};

// Internal: expand the pool when we run out of segments
static void expandPool(size_t min_segments_needed) {
    if (!g_pool.initialized) return;
    
    // Calculate new block size
    size_t new_block_size = g_pool.total_segments > 0 ? 
                           g_pool.total_segments : DEFAULT_POOL_SIZE;
    
    // Ensure we allocate enough for the request
    if (new_block_size < min_segments_needed) {
        new_block_size = min_segments_needed;
    }
    
    // Apply growth factor
    new_block_size = new_block_size * GROWTH_FACTOR;
    
    // Enforce maximum
    if (new_block_size > MAX_POOL_SIZE) {
        new_block_size = MAX_POOL_SIZE;
    }
    
    // Expand blocks array if needed
    if (g_pool.num_blocks >= g_pool.blocks_capacity) {
        size_t new_capacity = g_pool.blocks_capacity * 2;
        if (new_capacity == 0) new_capacity = 4;
        
        AncestrySegment **new_blocks = realloc(g_pool.blocks, 
                                              new_capacity * sizeof(AncestrySegment*));
        size_t *new_sizes = realloc(g_pool.block_sizes,
                                   new_capacity * sizeof(size_t));
        
        if (!new_blocks || !new_sizes) {
            fprintf(stderr, "Failed to expand pool block arrays\n");
            exit(1);
        }
        
        g_pool.blocks = new_blocks;
        g_pool.block_sizes = new_sizes;
        g_pool.blocks_capacity = new_capacity;
    }
    
    // Allocate new block
    AncestrySegment *new_block = calloc(new_block_size, sizeof(AncestrySegment));
    if (!new_block) {
        fprintf(stderr, "Failed to allocate segment pool block of size %zu\n", new_block_size);
        exit(1);
    }
    
    // Store block info
    g_pool.blocks[g_pool.num_blocks] = new_block;
    g_pool.block_sizes[g_pool.num_blocks] = new_block_size;
    g_pool.num_blocks++;
    g_pool.total_segments += new_block_size;
    
    // Expand free list capacity if needed
    size_t required_free_capacity = g_pool.free_count + new_block_size;
    if (required_free_capacity > g_pool.free_capacity) {
        size_t new_free_capacity = required_free_capacity * 2;
        AncestrySegment **new_free_list = realloc(g_pool.free_list,
                                                  new_free_capacity * sizeof(AncestrySegment*));
        if (!new_free_list) {
            fprintf(stderr, "Failed to expand free list\n");
            exit(1);
        }
        g_pool.free_list = new_free_list;
        g_pool.free_capacity = new_free_capacity;
    }
    
    // Add all segments from new block to free list
    for (size_t i = 0; i < new_block_size; i++) {
        g_pool.free_list[g_pool.free_count++] = &new_block[i];
        // Debug magic field would go here if added to struct
    }
    
    g_pool.stats.expansions++;
}

void initSegmentPool(size_t initial_size) {
    if (g_pool.initialized) return;
    
    memset(&g_pool, 0, sizeof(g_pool));
    
    // Determine initial size
    if (initial_size == 0) initial_size = DEFAULT_POOL_SIZE;
    if (initial_size < MIN_POOL_SIZE) initial_size = MIN_POOL_SIZE;
    if (initial_size > MAX_POOL_SIZE) initial_size = MAX_POOL_SIZE;
    
    // Initialize block arrays
    g_pool.blocks_capacity = 4;
    g_pool.blocks = calloc(g_pool.blocks_capacity, sizeof(AncestrySegment*));
    g_pool.block_sizes = calloc(g_pool.blocks_capacity, sizeof(size_t));
    
    if (!g_pool.blocks || !g_pool.block_sizes) {
        fprintf(stderr, "Failed to initialize segment pool\n");
        exit(1);
    }
    
    // Initialize free list
    g_pool.free_capacity = initial_size * 2;
    g_pool.free_list = calloc(g_pool.free_capacity, sizeof(AncestrySegment*));
    if (!g_pool.free_list) {
        fprintf(stderr, "Failed to allocate free list\n");
        exit(1);
    }
    
    g_pool.initialized = 1;
    
    // Allocate initial block
    expandPool(initial_size);
}

void destroySegmentPool(void) {
    if (!g_pool.initialized) return;
    
    // Free all blocks
    for (size_t i = 0; i < g_pool.num_blocks; i++) {
        free(g_pool.blocks[i]);
    }
    
    // Free management structures
    free(g_pool.blocks);
    free(g_pool.block_sizes);
    free(g_pool.free_list);
    
    memset(&g_pool, 0, sizeof(g_pool));
}

void resetSegmentPool(void) {
    if (!g_pool.initialized) return;
    
    // O(1) reset - just rebuild the free list with all segments
    g_pool.free_count = 0;
    
    // Add all segments from all blocks back to free list
    for (size_t b = 0; b < g_pool.num_blocks; b++) {
        AncestrySegment *block = g_pool.blocks[b];
        size_t block_size = g_pool.block_sizes[b];
        
        for (size_t i = 0; i < block_size; i++) {
            // Clear the segment
            memset(&block[i], 0, sizeof(AncestrySegment));
            // Debug magic field would go here if added to struct
            
            // Add to free list
            g_pool.free_list[g_pool.free_count++] = &block[i];
        }
    }
    
    // Update statistics
    g_pool.stats.current_in_use = 0;
    g_pool.stats.resets++;
}

AncestrySegment* allocSegmentFromPool(void) {
    if (!g_pool.initialized) {
        initSegmentPool(DEFAULT_POOL_SIZE);
    }
    
    // Expand pool if necessary
    if (g_pool.free_count == 0) {
        expandPool(1000);  // Allocate at least 1000 more
    }
    
    // Pop from free list (LIFO)
    AncestrySegment *seg = g_pool.free_list[--g_pool.free_count];
    
    // Debug assertion would go here if magic field added to struct
    
    // The segment should already be zeroed
    // Just update statistics
    g_pool.stats.total_allocations++;
    g_pool.stats.current_in_use++;
    
    if (g_pool.stats.current_in_use > g_pool.stats.high_water_mark) {
        g_pool.stats.high_water_mark = g_pool.stats.current_in_use;
    }
    
    return seg;
}

void freeSegmentToPool(AncestrySegment *seg) {
    if (!seg || !g_pool.initialized) return;
    
    #ifdef DEBUG_SEGMENT_POOL
    // Check for double-free
    for (size_t i = 0; i < g_pool.free_count; i++) {
        assert(g_pool.free_list[i] != seg);
    }
    #endif
    
    // Clear segment data (important for safety)
    void *avl_tree = seg->avlTree;  // Preserve AVL tree pointer for cleanup
    memset(seg, 0, sizeof(AncestrySegment));
    
    // AVL tree needs special cleanup if present
    if (avl_tree) {
        freeAVLTree((AVLTree*)avl_tree);
    }
    
    // Debug magic field would go here if added to struct
    
    // Return to free list
    g_pool.free_list[g_pool.free_count++] = seg;
    
    // Update statistics
    g_pool.stats.total_frees++;
    g_pool.stats.current_in_use--;
}

size_t calculateInitialPoolSize(int sampleSize, double rho, int nSites) {
    // Base allocation for samples and their initial coalescences
    size_t base = sampleSize * 10;
    
    // rho is the total recombination rate for the region (from command line)
    // Each recombination can create ~2-4 new segments on average
    size_t recomb_segments = (size_t)(rho * 3);
    
    // Total initial size
    size_t initial = base + recomb_segments;
    
    // Apply bounds
    if (initial < DEFAULT_POOL_SIZE) initial = DEFAULT_POOL_SIZE;
    if (initial > 1000000) initial = 1000000;  // Cap at 1M for initial allocation
    
    return initial;
}

PoolStats getPoolStatistics(void) {
    g_pool.stats.pool_capacity = g_pool.total_segments;
    return g_pool.stats;
}

void printPoolStatistics(FILE *fp) {
    if (!g_pool.initialized) {
        fprintf(fp, "Segment pool not initialized\n");
        return;
    }
    
    PoolStats stats = getPoolStatistics();
    
    fprintf(fp, "\n=== Segment Pool Statistics ===\n");
    fprintf(fp, "Total allocations:    %zu\n", stats.total_allocations);
    fprintf(fp, "Total frees:          %zu\n", stats.total_frees);
    fprintf(fp, "Currently in use:     %zu\n", stats.current_in_use);
    fprintf(fp, "High water mark:      %zu\n", stats.high_water_mark);
    fprintf(fp, "Pool capacity:        %zu\n", stats.pool_capacity);
    fprintf(fp, "Pool expansions:      %zu\n", stats.expansions);
    fprintf(fp, "Resets:               %zu\n", stats.resets);
    fprintf(fp, "Memory usage:         %.2f MB\n", getPoolMemoryUsage() / (1024.0 * 1024.0));
    fprintf(fp, "Efficiency:           %.1f%%\n", 
            stats.pool_capacity > 0 ? 
            (100.0 * stats.high_water_mark / stats.pool_capacity) : 0.0);
}

size_t getPoolMemoryUsage(void) {
    if (!g_pool.initialized) return 0;
    
    size_t total = 0;
    
    // Segment storage
    total += g_pool.total_segments * sizeof(AncestrySegment);
    
    // Management structures
    total += g_pool.blocks_capacity * sizeof(AncestrySegment*);
    total += g_pool.blocks_capacity * sizeof(size_t);
    total += g_pool.free_capacity * sizeof(AncestrySegment*);
    total += sizeof(SegmentPool);
    
    return total;
}

void trimSegmentPool(size_t target_size) {
    if (!g_pool.initialized) return;
    
    // This is a placeholder for future implementation
    // Would free unused blocks to reduce memory usage
    // For now, we don't shrink the pool
}

#ifdef DEBUG_SEGMENT_POOL
void validateSegmentPool(void) {
    if (!g_pool.initialized) return;
    
    size_t total_validated = 0;
    
    // Check all blocks
    for (size_t b = 0; b < g_pool.num_blocks; b++) {
        AncestrySegment *block = g_pool.blocks[b];
        size_t block_size = g_pool.block_sizes[b];
        
        for (size_t i = 0; i < block_size; i++) {
            // Debug assertion would go here
            total_validated++;
        }
    }
    
    assert(total_validated == g_pool.total_segments);
    assert(g_pool.free_count + g_pool.stats.current_in_use == g_pool.total_segments);
}
#endif