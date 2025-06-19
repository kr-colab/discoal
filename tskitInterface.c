// tskitInterface.c
// Implementation of tskit tree sequence recording for discoal

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/utsname.h>
#include "tskitInterface.h"
#include "ranlib.h"
#include "version.h"

// Global tskit table collection
tsk_table_collection_t *tsk_tables = NULL;

// Mapping between discoal node IDs and tskit node IDs
// We'll use the totNodeNumber as the key
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
    
    // Set time units to indicate coalescent scaling
    const char *time_units_str = "coalescent units (2N generations)";
    ret = tsk_table_collection_set_time_units(tsk_tables, time_units_str, strlen(time_units_str));
    if (ret != 0) {
        tsk_table_collection_free(tsk_tables);
        free(tsk_tables);
        tsk_tables = NULL;
        return ret;
    }
    
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
                                     time,  
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
// Global variable to store command line
static char *global_command_line = NULL;

// Function to store command line (called from main)
void tskit_store_command_line(int argc, const char **argv) {
    if (global_command_line != NULL) {
        free(global_command_line);
    }
    
    // Calculate total length needed
    size_t total_len = 0;
    for (int i = 0; i < argc; i++) {
        total_len += strlen(argv[i]) + 1; // +1 for space
    }
    
    global_command_line = malloc(total_len);
    if (global_command_line == NULL) {
        return;
    }
    
    // Build command string
    global_command_line[0] = '\0';
    for (int i = 0; i < argc; i++) {
        if (i > 0) strcat(global_command_line, " ");
        strcat(global_command_line, argv[i]);
    }
}

// Add provenance record
static int add_provenance(void) {
    char timestamp[64];
    char provenance_json[8192];
    struct utsname sys_info;
    int ret;
    
    // Get current time in ISO 8601 format
    time_t now = time(NULL);
    struct tm *utc_tm = gmtime(&now);
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%dT%H:%M:%SZ", utc_tm);
    
    // Get system information
    uname(&sys_info);
    
    // Get simulation parameters from global variables
    extern int sampleSize;
    extern int nSites;
    extern double theta;
    extern double rho;
    extern int npops;
    extern int seed1, seed2;
    extern double *currentSize;  // Current population sizes
    extern double alpha;         // Selection coefficient
    extern double f0;            // Starting frequency
    extern double sweepSite;     // Sweep position
    extern char sweepMode;       // Type of sweep
    extern int activeSweepFlag;  // Whether sweep is active
    extern double migMat[MAXPOPS][MAXPOPS];  // Migration matrix
    
    // Get tskit version
    char tskit_version[32] = "C-API";  // tskit C API doesn't expose version easily
    
    // Use stored command line or reconstruct from globals
    const char *command_str = global_command_line ? global_command_line : "discoal [command line not available]";
    
    // Build parameters object
    char params_str[4096];
    snprintf(params_str, sizeof(params_str),
        "\"sample_size\": %d,"
        "\"num_sites\": %d,"
        "\"theta\": %.6f,"
        "\"rho\": %.6f,"
        "\"num_populations\": %d,"
        "\"random_seeds\": [%d, %d]",
        sampleSize, nSites, theta, rho, npops, seed1, seed2);
    
    // Add population sizes if multiple populations
    if (npops > 1 && currentSize != NULL) {
        char pop_sizes[512];
        strcat(params_str, ",\"population_sizes\": [");
        for (int i = 0; i < npops; i++) {
            char size_str[64];
            snprintf(size_str, sizeof(size_str), "%.6f%s", currentSize[i], 
                     (i < npops - 1) ? ", " : "]");
            strcat(params_str, size_str);
        }
    }
    
    // Add selection parameters if applicable
    if (alpha != 0.0) {
        char selection_str[512];
        snprintf(selection_str, sizeof(selection_str),
            ",\"selection\": {"
            "\"coefficient\": %.6f,"
            "\"sweep_site\": %.6f,"
            "\"starting_frequency\": %.6f,"
            "\"sweep_mode\": \"%c\""
            "}",
            alpha, sweepSite, f0, sweepMode);
        strcat(params_str, selection_str);
    }
    
    // Build complete provenance record
    snprintf(provenance_json, sizeof(provenance_json),
        "{"
        "\"schema_version\": \"1.0.0\","
        "\"software\": {"
            "\"name\": \"discoal\","
            "\"version\": \"" DISCOAL_VERSION "\""
        "},"
        "\"parameters\": {"
            "%s,"
            "\"command\": \"%s\""
        "},"
        "\"environment\": {"
            "\"os\": {"
                "\"system\": \"%s\","
                "\"node\": \"%s\","
                "\"release\": \"%s\","
                "\"version\": \"%s\","
                "\"machine\": \"%s\""
            "},"
            "\"libraries\": {"
                "\"tskit\": {\"version\": \"%s\"}"
            "}"
        "}"
        "}",
        params_str,
        command_str,
        sys_info.sysname,
        sys_info.nodename,
        sys_info.release,
        sys_info.version,
        sys_info.machine,
        tskit_version);
    
    // Add to provenance table
    ret = tsk_provenance_table_add_row(&tsk_tables->provenances, 
        timestamp, strlen(timestamp),
        provenance_json, strlen(provenance_json));
    
    if (ret < 0) {
        fprintf(stderr, "Error adding provenance: %s\n", tsk_strerror(ret));
        return ret;
    }
    
    return 0;
}

int tskit_finalize(const char *filename) {
    int ret;
    
    if (tsk_tables == NULL) {
        return -1;
    }
    
    // Add provenance record
    ret = add_provenance();
    if (ret != 0) {
        fprintf(stderr, "Warning: Failed to add provenance record\n");
        // Continue anyway - provenance is not critical
    }
    
    // Sort tables (required for valid tree sequences)
    ret = tsk_table_collection_sort(tsk_tables, NULL, 0);
    if (ret != 0) {
        fprintf(stderr, "Error sorting tables: %s\n", tsk_strerror(ret));
        return ret;
    }
    
    // Note: We simplify before mutation placement to ensure mutations are only
    // placed on edges ancestral to samples
    
    // Build index
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
    
    // First pass: estimate mutation count and allocate dynamically
    // Estimate based on total branch length * theta * sequence length
    double total_branch_length = 0.0;
    for (tsk_id_t edge_id = 0; edge_id < tsk_tables->edges.num_rows; edge_id++) {
        tsk_id_t parent_id = tsk_tables->edges.parent[edge_id];
        tsk_id_t child_id = tsk_tables->edges.child[edge_id];
        double parent_time = tsk_tables->nodes.time[parent_id];
        double child_time = tsk_tables->nodes.time[child_id];
        double branch_length = parent_time - child_time;
        double genomic_span = tsk_tables->edges.right[edge_id] - tsk_tables->edges.left[edge_id];
        total_branch_length += branch_length * (genomic_span / tsk_tables->sequence_length);
    }
    
    // Estimate expected mutations with buffer
    int estimated_mutations = (int)(total_branch_length * theta) + 50;
    if (estimated_mutations > MAXMUTS) estimated_mutations = MAXMUTS;
    if (estimated_mutations < 10) estimated_mutations = 10;  // Minimum allocation
    
    // Allocate mutations array dynamically
    temp_mutation_t *mutations = malloc(sizeof(temp_mutation_t) * estimated_mutations);
    if (mutations == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for mutations array\\n");
        return -1;
    }
    
    int mutation_count = 0;
    int mutations_capacity = estimated_mutations;
    
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
        // Note: branch_length is in units of 2N generations, but theta=4Nμ assumes time in generations
        // So we need to multiply by 0.5 to account for the time scaling
        double expected_mutations = branch_length * (genomic_span / tsk_tables->sequence_length) * theta * 0.5;
        
        if (expected_mutations > 0.0) {
            // Draw number of mutations from Poisson distribution
            int num_mutations = ignpoi(expected_mutations);
            
            // Collect mutations for later sorting
            for (int i = 0; i < num_mutations; i++) {
                // Check if we need to expand capacity
                if (mutation_count >= mutations_capacity) {
                    int new_capacity = mutations_capacity * 2;
                    if (new_capacity > MAXMUTS) new_capacity = MAXMUTS;
                    
                    if (mutation_count >= MAXMUTS) {
                        fprintf(stderr, "Error: Too many mutations (limit: %d)\\n", MAXMUTS);
                        free(mutations);
                        return -1;
                    }
                    
                    temp_mutation_t *new_mutations = realloc(mutations, sizeof(temp_mutation_t) * new_capacity);
                    if (new_mutations == NULL) {
                        fprintf(stderr, "Error: Failed to expand mutations array\\n");
                        free(mutations);
                        return -1;
                    }
                    mutations = new_mutations;
                    mutations_capacity = new_capacity;
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
        // Sort mutations by position 
        qsort(mutations, mutation_count, sizeof(mutations[0]), compare_mutations_by_position);
        
        // Add sites and mutations in sorted order
        for (int i = 0; i < mutation_count; i++) {
            // Create a new site for each mutation
            tsk_id_t site_id = tskit_add_site(mutations[i].position, "0");  // ancestral state is "0"
            if (site_id < 0) {
                fprintf(stderr, "Error: Failed to add site at position %f\\n", mutations[i].position);
                free(mutations);
                return -1;
            }
            
            // Add the mutation at this site on the appropriate node
            int ret = tskit_add_mutation(site_id, mutations[i].node_id, "1");  // derived state is "1"
            if (ret < 0) {
                fprintf(stderr, "Error: Failed to add mutation at site %lld on node %lld\\n", 
                        (long long)site_id, (long long)mutations[i].node_id);
                free(mutations);
                return -1;
            }
        }
    }
    
    // Clean up allocated memory
    free(mutations);
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

// Record sweep mutations for nodes that carry them
int tskit_record_sweep_mutations(double sweepSite) {
    tsk_id_t site_id = TSK_NULL;
    int mutation_count = 0;
    int i;
    
    if (tsk_tables == NULL || sweepSite < 0.0) {
        return -1;
    }
    
    // First, check if any nodes carry the sweep mutation
    int has_sweep_carriers = 0;
    for (i = 0; i < alleleNumber; i++) {
        if (nodes[i]->carriesSweepMutation) {
            has_sweep_carriers = 1;
            break;
        }
    }
    
    if (!has_sweep_carriers) {
        return 0;  // No sweep mutations to add
    }
    
    // Add the site for the sweep mutation
    site_id = tskit_add_site(sweepSite, "0");
    if (site_id < 0) {
        fprintf(stderr, "Error: Failed to add sweep site at position %f\n", sweepSite);
        return -1;
    }
    
    // Add mutations for all nodes that carry the sweep mutation
    for (i = 0; i < alleleNumber; i++) {
        if (nodes[i]->carriesSweepMutation) {
            tsk_id_t node_id = get_tskit_node_id(nodes[i]);
            if (node_id >= 0) {
                if (tskit_add_mutation(site_id, node_id, "1") >= 0) {
                    mutation_count++;
                }
            }
        }
    }
    
    return mutation_count;
}