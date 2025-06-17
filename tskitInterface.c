// tskitInterface.c
// Implementation of tskit tree sequence recording for discoal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "tskitInterface.h"
#include "ranlib.h"

// Global tskit table collection
tsk_table_collection_t *tsk_tables = NULL;

// Mapping between discoal node IDs and tskit node IDs
// We'll use the allNodes array index as the key
tsk_id_t *node_id_map = NULL;
int node_id_map_capacity = 0;

// Initialize tskit recording
int tskit_initialize(double sequence_length) {
    int ret;
    
    // Allocate table collection
    tsk_tables = malloc(sizeof(tsk_table_collection_t));
    if (tsk_tables == NULL) {
        return -1;
    }
    
    // Initialize table collection
    ret = tsk_table_collection_init(tsk_tables, 0);
    if (ret != 0) {
        free(tsk_tables);
        tsk_tables = NULL;
        return ret;
    }
    
    // Set sequence length
    tsk_tables->sequence_length = sequence_length;
    
    // Add population records for all populations in discoal
    // We need to add enough populations to cover the maximum population ID used
    extern int npops;  // From discoal globals
    for (int i = 0; i < npops; i++) {
        ret = tsk_population_table_add_row(&tsk_tables->populations, NULL, 0);
        if (ret < 0) {
            tsk_table_collection_free(tsk_tables);
            free(tsk_tables);
            tsk_tables = NULL;
            return ret;
        }
    }
    
    // Initialize node ID mapping
    node_id_map_capacity = 1000;  // Start with reasonable size
    node_id_map = malloc(sizeof(tsk_id_t) * node_id_map_capacity);
    if (node_id_map == NULL) {
        tsk_table_collection_free(tsk_tables);
        free(tsk_tables);
        tsk_tables = NULL;
        return -1;
    }
    
    // Initialize all mappings to TSK_NULL
    for (int i = 0; i < node_id_map_capacity; i++) {
        node_id_map[i] = TSK_NULL;
    }
    
    return 0;
}

// Ensure node_id_map has enough capacity
static void ensure_node_map_capacity(int required_size) {
    if (required_size >= node_id_map_capacity) {
        int new_capacity = node_id_map_capacity;
        while (new_capacity <= required_size) {
            new_capacity *= 2;
        }
        
        tsk_id_t *new_map = realloc(node_id_map, sizeof(tsk_id_t) * new_capacity);
        if (new_map == NULL) {
            fprintf(stderr, "Error: Failed to reallocate node_id_map\n");
            exit(1);
        }
        
        // Initialize new entries to TSK_NULL
        for (int i = node_id_map_capacity; i < new_capacity; i++) {
            new_map[i] = TSK_NULL;
        }
        
        node_id_map = new_map;
        node_id_map_capacity = new_capacity;
    }
}

// Record a new node
tsk_id_t tskit_add_node(double time, int population, int is_sample) {
    tsk_flags_t flags = is_sample ? TSK_NODE_IS_SAMPLE : 0;
    tsk_id_t node_id;
    
    if (tsk_tables == NULL) {
        return TSK_NULL;
    }
    
    // tskit doesn't allow negative population IDs, so map -1 to TSK_NULL
    tsk_id_t tsk_population = (population < 0) ? TSK_NULL : population;
    
    node_id = tsk_node_table_add_row(&tsk_tables->nodes, 
                                     flags, 
                                     time / 2.0,  // Scale time by 0.5 to match msprime convention
                                     tsk_population, 
                                     TSK_NULL,  // individual
                                     NULL, 0);  // metadata
    
    // Debug output (commented out for production)
    // if (is_sample) {
    //     fprintf(stderr, "Added sample node %lld: time=%f (scaled=%f), pop=%d\n", 
    //             (long long)node_id, time, time / 2.0, tsk_population);
    // }
    
    return node_id;
}

// Record edges during coalescence
int tskit_add_edges(tsk_id_t parent, tsk_id_t child, double left, double right) {
    int ret;
    
    if (tsk_tables == NULL) {
        return -1;
    }
    
    ret = tsk_edge_table_add_row(&tsk_tables->edges,
                                 left, right,
                                 parent, child,
                                 NULL, 0);  // metadata
    
    // Debug output (commented out for production)
    // if (right - left < 0.0001) {
    //     fprintf(stderr, "WARNING: Very small edge interval: parent=%lld, child=%lld, left=%.10f, right=%.10f, length=%.10f\n", 
    //             (long long)parent, (long long)child, left, right, right - left);
    // }
    
    return ret;
}

// Record a site
tsk_id_t tskit_add_site(double position, const char *ancestral_state) {
    tsk_id_t site_id;
    
    if (tsk_tables == NULL) {
        return TSK_NULL;
    }
    
    site_id = tsk_site_table_add_row(&tsk_tables->sites,
                                     position,
                                     ancestral_state,
                                     strlen(ancestral_state),
                                     NULL, 0);  // metadata
    
    return site_id;
}

// Record a mutation
int tskit_add_mutation(tsk_id_t site, tsk_id_t node, const char *derived_state) {
    int ret;
    
    if (tsk_tables == NULL) {
        return -1;
    }
    
    ret = tsk_mutation_table_add_row(&tsk_tables->mutations,
                                     site,
                                     node,
                                     TSK_NULL,  // parent mutation
                                     TSK_UNKNOWN_TIME,  // time
                                     derived_state,
                                     strlen(derived_state),
                                     NULL, 0);  // metadata
    
    return ret;
}

// Finalize and write tree sequence
int tskit_finalize(const char *filename) {
    int ret;
    
    if (tsk_tables == NULL) {
        return -1;
    }
    
    
    // Sort tables (required for valid tree sequences)
    ret = tsk_table_collection_sort(tsk_tables, NULL, 0);
    if (ret != 0) {
        fprintf(stderr, "Error sorting tables: %s\n", tsk_strerror(ret));
        return ret;
    }
    
    // Note: We're not simplifying because with recombination, it's natural
    // to have multiple roots - different genomic regions may not fully coalesce
    
    // Build index after simplification
    ret = tsk_table_collection_build_index(tsk_tables, 0);
    if (ret != 0) {
        fprintf(stderr, "Error building index: %s\n", tsk_strerror(ret));
        return ret;
    }
    
    // Write to file
    ret = tsk_table_collection_dump(tsk_tables, filename, 0);
    if (ret != 0) {
        fprintf(stderr, "Error writing tree sequence: %s\n", tsk_strerror(ret));
        return ret;
    }
    
    return 0;
}

// Clean up tskit resources
void tskit_cleanup(void) {
    if (tsk_tables != NULL) {
        tsk_table_collection_free(tsk_tables);
        free(tsk_tables);
        tsk_tables = NULL;
    }
    
    if (node_id_map != NULL) {
        free(node_id_map);
        node_id_map = NULL;
        node_id_map_capacity = 0;
    }
}

// Map discoal node to tskit node ID
tsk_id_t get_tskit_node_id(rootedNode *node) {
    // In tskit-only mode, tskit ID is stored directly in the node
    return node->tskit_node_id;

}

// Store mapping from discoal node to tskit node ID using index
void set_tskit_node_id_at_index(int index, tsk_id_t tsk_id) {
    ensure_node_map_capacity(index);
    node_id_map[index] = tsk_id;
    // fprintf(stderr, "Stored mapping at index %d -> tskit ID %lld\n", 
    //         index, (long long)tsk_id);
}

// Store mapping from discoal node to tskit node ID
void set_tskit_node_id(rootedNode *node, tsk_id_t tsk_id) {
    // In tskit-only mode, store tskit ID directly in the node
    node->tskit_node_id = tsk_id;
    return;

}

// Record all mutations after they've been placed on the tree
int tskit_record_mutations(void) {
    if (tsk_tables == NULL) {
        return -1;
    }
    
    // TODO: Implement tskit-only mutation recording
    // In tskit-only mode, mutations should be recorded directly during simulation
    return 0;

}

// Structure for temporary mutation storage
typedef struct {
    double position;
    tsk_id_t node_id;
} temp_mutation_t;

// Comparison function for sorting mutations by position
static int compare_mutations_by_position(const void *a, const void *b) {
    const temp_mutation_t *mut_a = (const temp_mutation_t *)a;
    const temp_mutation_t *mut_b = (const temp_mutation_t *)b;
    
    if (mut_a->position < mut_b->position) return -1;
    if (mut_a->position > mut_b->position) return 1;
    return 0;
}

// Place mutations using edge-based algorithm (similar to msprime's design)
int tskit_place_mutations_edge_based(double theta) {
    if (tsk_tables == NULL) {
        return -1;
    }
    
    extern int nSites;
    
    // First pass: collect all mutation information
    temp_mutation_t mutations[MAXMUTS];
    int mutation_count = 0;
    
    // Iterate through all edges and place mutations independently
    for (tsk_id_t edge_id = 0; edge_id < tsk_tables->edges.num_rows; edge_id++) {
        // Get edge information
        double left = tsk_tables->edges.left[edge_id];
        double right = tsk_tables->edges.right[edge_id];
        tsk_id_t parent_id = tsk_tables->edges.parent[edge_id];
        tsk_id_t child_id = tsk_tables->edges.child[edge_id];
        
        // Calculate branch length (in tskit time units)
        double parent_time = tsk_tables->nodes.time[parent_id];
        double child_time = tsk_tables->nodes.time[child_id];
        double branch_length = parent_time - child_time;
        
        // Calculate genomic span
        double genomic_span = right - left;
        
        // Calculate expected number of mutations on this edge
        // theta is the population-scaled mutation rate for the entire sequence
        // For a fraction of the sequence: mutations ~ Poisson(branch_length * (genomic_span / sequence_length) * theta)
        double expected_mutations = branch_length * (genomic_span / tsk_tables->sequence_length) * theta;
        
        if (expected_mutations > 0.0) {
            // Draw number of mutations from Poisson distribution
            int num_mutations = ignpoi(expected_mutations);
            
            // Collect mutations for later sorting
            for (int i = 0; i < num_mutations; i++) {
                if (mutation_count >= MAXMUTS) {
                    fprintf(stderr, "Error: Too many mutations (limit: %d)\n", MAXMUTS);
                    return -1;
                }
                
                // Generate position uniformly within this edge's interval
                double mut_pos = genunf(left, right);
                
                // Store mutation info for later processing
                mutations[mutation_count].position = mut_pos;
                mutations[mutation_count].node_id = child_id;
                mutation_count++;
            }
        }
    }
    
    // Second pass: sort mutations by position and add to tskit
    if (mutation_count > 0) {
        // Sort mutations by position using efficient qsort - O(M log M) instead of O(M²)
        qsort(mutations, mutation_count, sizeof(mutations[0]), compare_mutations_by_position);
        
        // Add sites and mutations in sorted order
        for (int i = 0; i < mutation_count; i++) {
            // Create a new site for each mutation
            tsk_id_t site_id = tskit_add_site(mutations[i].position, "0");  // ancestral state is "0"
            if (site_id < 0) {
                fprintf(stderr, "Error: Failed to add site at position %f\n", mutations[i].position);
                return -1;
            }
            
            // Add the mutation at this site on the appropriate node
            int ret = tskit_add_mutation(site_id, mutations[i].node_id, "1");  // derived state is "1"
            if (ret < 0) {
                fprintf(stderr, "Error: Failed to add mutation at site %lld on node %lld\n", 
                        (long long)site_id, (long long)mutations[i].node_id);
                return -1;
            }
        }
    }
    
    return 0;
}

// Original function name kept for backward compatibility
int tskit_place_mutations_directly(double theta) {
    // Use edge-based algorithm
    return tskit_place_mutations_edge_based(theta);
}

// Populate discoal mutation arrays from tskit data (for ms output compatibility)
int tskit_populate_discoal_mutations(void) {
    if (tsk_tables == NULL) {
        return -1;
    }
    
    // TODO: Implement tskit-only mutation population
    // In tskit-only mode, use genotype matrix directly for output
    return 0;

}