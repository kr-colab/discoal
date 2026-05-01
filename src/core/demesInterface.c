#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <yaml.h>
#include "demesInterface.h"
#include "discoal.h"
#include "discoalFunctions.h"
#include "demes.h"

// Helper function to ensure events array has enough capacity
static void ensureDemesEventsCapacity(event **events, int *eventsCapacity, int needed) {
    if (needed >= *eventsCapacity) {
        int newCapacity = needed * 2;
        struct event *newEvents = (struct event*) realloc(*events, newCapacity * sizeof(struct event));
        if (newEvents == NULL) {
            fprintf(stderr, "Error: Failed to reallocate events array\n");
            exit(1);
        }
        *events = newEvents;
        *eventsCapacity = newCapacity;
        // Initialize new elements to zero
        memset(&(*events)[*eventsCapacity/2], 0, (newCapacity - *eventsCapacity/2) * sizeof(struct event));
    }
}

// Helper function to find population index by name
static int findPopulationIndex(struct demes_graph *graph, const char *name) {
    for (int i = 0; i < graph->n_demes; i++) {
        if (strcmp((char*)graph->demes[i].name, name) == 0) {
            return i;
        }
    }
    return -1;
}

// Convert demes time to coalescent time. demes_time is in the units the
// graph specifies; for time_units="generations" generation_time is 1 (demes
// validates this), and for time_units="years" generation_time is the years
// per generation, so dividing yields generations either way.
static double demesTimeToCoalTime(double demes_time, double generation_time, double N) {
    return demes_time / generation_time / (2.0 * N);
}

// Convert demes size to coalescent size (relative to ancestral N)
static double demesSizeToCoalSize(double demes_size, double N) {
    return demes_size / N;
}

// Compare events by time for sorting (descending order - most recent first)
static int compareEventsByTime(const void *a, const void *b) {
    const struct event *ea = (const struct event *)a;
    const struct event *eb = (const struct event *)b;
    if (ea->time < eb->time) return 1;
    if (ea->time > eb->time) return -1;
    return 0;
}

// Compare doubles ascending (used for migration boundary deduplication)
static int compareDoublesAscending(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return  1;
    return 0;
}

// Per-pair migration window collected from graph->migrations.
// Times are in discoal internal coalescent units. Demes spec says
// start_time > end_time (start_time is older, end_time is more recent),
// so after conversion t_start > t_end as well.
struct mig_window {
    int src, dst;
    double t_start;  // older boundary, larger internal time
    double t_end;    // more-recent boundary, smaller internal time
    double rate;     // 4*N*demes_rate
};

int loadDemesFile(const char *filename, event **events, int *eventNumber, int *eventsCapacity, 
                  double *currentSize, int *npops, int *sampleSizes, double N) {
    // External declarations for discoal globals
    extern int migFlag;
    extern double tDiv;
    struct demes_graph *graph;
    int ret;
    
    // Load the demes file
    ret = demes_graph_load(filename, &graph);
    if (ret != 0) {
        fprintf(stderr, "Error loading demes file '%s': ", filename);
        switch(ret) {
            case DEMES_ERR_MEMORY:
                fprintf(stderr, "Memory allocation failed\n");
                break;
            case DEMES_ERR_IO:
                fprintf(stderr, "File I/O error\n");
                break;
            case DEMES_ERR_YAML:
                fprintf(stderr, "YAML parsing error\n");
                break;
            case DEMES_ERR_TYPE:
                fprintf(stderr, "Type mismatch in YAML\n");
                break;
            case DEMES_ERR_KEY:
                fprintf(stderr, "Invalid field in YAML\n");
                break;
            case DEMES_ERR_VALUE:
                fprintf(stderr, "Invalid value in YAML\n");
                break;
            case DEMES_ERR_MISSING_REQUIRED:
                fprintf(stderr, "Required field missing\n");
                break;
            case DEMES_ERR_PROPORTIONS:
                fprintf(stderr, "Invalid proportions\n");
                break;
            default:
                fprintf(stderr, "Unknown error code %d\n", ret);
        }
        return ret;
    }
    
    // First, find the present-day size of the first PRESENT-DAY population
    // This will be used as the reference N for all scaling
    double referenceN = N;  // default if we can't find a present-day population
    
    // Find the first population that exists at time 0 (present day)
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        if (deme->n_epochs > 0) {
            struct demes_epoch *last_epoch = &deme->epochs[deme->n_epochs - 1];
            // Check if this population exists at present day (end_time == 0)
            if (last_epoch->end_time == 0.0) {
                // This is a present-day population, use its size as reference
                if (last_epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL) {
                    referenceN = last_epoch->end_size;
                } else {
                    referenceN = last_epoch->start_size;
                }
                fprintf(stderr, "Using population '%s' (index %d) as reference with N=%g\n", 
                        deme->name, i, referenceN);
                break;  // Use the first present-day population
            }
        }
    }
    
    // Convert the graph to discoal events using the reference N
    ret = convertDemesToEvents(graph, events, eventNumber, eventsCapacity, 
                              currentSize, npops, sampleSizes, referenceN);
    
    if (ret != 0) {
        demes_graph_free(graph);
        return ret;
    }
    
    // Validate that populations can coalesce
    // For multi-population models, we need either migration or population mergers
    // Count only present-day populations (those that exist at time 0)
    int presentDayPops = 0;
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        if (deme->n_epochs > 0) {
            // Check if this deme exists at time 0 (present day)
            struct demes_epoch *last_epoch = &deme->epochs[deme->n_epochs - 1];
            if (last_epoch->end_time == 0.0) {
                presentDayPops++;
            }
        }
    }
    
    if (presentDayPops > 1) {
        fprintf(stderr, "\nValidating population connectivity (%d present-day populations)...\n", presentDayPops);
        // Check if there's any migration
        int hasMigration = 0;
        for (int i = 0; i < graph->n_migrations; i++) {
            if (graph->migrations[i].rate > 0) {
                hasMigration = 1;
                break;
            }
        }
        fprintf(stderr, "Has migration: %s\n", hasMigration ? "yes" : "no");
        
        // Check if all populations eventually merge
        // We need to check if all present-day populations trace back to a common ancestor
        int hasCommonAncestor = 0;
        
        // Check specifically if present-day populations have ancestors
        for (int i = 0; i < graph->n_demes; i++) {
            struct demes_deme *deme = &graph->demes[i];
            // Only check demes that exist at present
            if (deme->n_epochs > 0 && deme->epochs[deme->n_epochs - 1].end_time == 0.0) {
                if (deme->n_ancestors > 0) {
                    hasCommonAncestor = 1;
                    break;
                }
            }
        }
        fprintf(stderr, "Present-day populations have ancestors: %s\n", hasCommonAncestor ? "yes" : "no");
        
        // If there's no migration and no apparent mergers, issue an error
        if (!hasMigration && !hasCommonAncestor) {
            fprintf(stderr, "\nError: The demographic model has %d present-day populations but no migration or merger events.\n", presentDayPops);
            fprintf(stderr, "This would result in infinite coalescent time as populations cannot coalesce.\n");
            fprintf(stderr, "Please add migration between populations or ensure populations merge.\n");
            demes_graph_free(graph);
            return -1;
        }
        
        // Even with ancestors, we should verify all populations can eventually coalesce
        // This is a more thorough check
        if (!hasMigration && hasCommonAncestor) {
            // Just provide information, not a warning, since merged populations can coalesce
            fprintf(stderr, "\nNote: Populations will coalesce through common ancestors.\n");
        }
    }
    
    // Set flags for discoal
    if (graph->n_migrations > 0) {
        migFlag = 1;  // Enable migration flag if there are migrations
    }
    
    // If there are population splits/merges, set tDiv to a non-666 value
    for (int i = 0; i < graph->n_demes; i++) {
        if (graph->demes[i].n_ancestors > 0) {
            tDiv = 1.0;  // Set to non-666 value to indicate splits exist
            break;
        }
    }
    
    // Free the graph
    demes_graph_free(graph);
    
    return ret;
}

int convertDemesToEvents(struct demes_graph *graph, event **events, int *eventNumber, 
                        int *eventsCapacity, double *currentSize, int *npops, 
                        int *sampleSizes, double N) {
    
    // Set number of populations
    *npops = graph->n_demes;
    if (*npops > MAXPOPS) {
        fprintf(stderr, "Error: Demes file has %d populations, but discoal supports at most %d\n", 
                *npops, MAXPOPS);
        return -1;
    }
    
    // Debug output for conversion parameters
    fprintf(stderr, "\n=== DEMES TO DISCOAL CONVERSION DEBUG ===\n");
    fprintf(stderr, "Reference N (present-day pop0 size): %g\n", N);
    fprintf(stderr, "Number of populations: %d\n", *npops);
    fprintf(stderr, "Generation time: %g\n", graph->generation_time);
    
    // Initialize population sizes and collect all events
    int totalEvents = *eventNumber;  // Start from existing events
    
    // First pass: count events needed
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        
        // Each epoch transition is an event (size change)
        totalEvents += deme->n_epochs;
        
        // Population split/merge events
        if (deme->n_ancestors > 0) {
            totalEvents += 1;  // One event for the split/merge
        }
    }
    
    // Add migration events (start and end of each migration)
    totalEvents += graph->n_migrations * 2;
    
    // Add pulse (admixture) events
    totalEvents += graph->n_pulses;
    
    // Ensure capacity
    ensureDemesEventsCapacity(events, eventsCapacity, totalEvents + 10);
    
    // Second pass: create events
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        int popID = i;
        
        fprintf(stderr, "\nPopulation %d: '%s'\n", popID, deme->name);
        fprintf(stderr, "  Start time: %g generations\n", deme->start_time);
        fprintf(stderr, "  Number of epochs: %zu\n", deme->n_epochs);
        if (deme->n_ancestors > 0) {
            fprintf(stderr, "  Ancestors: ");
            for (size_t j = 0; j < deme->n_ancestors; j++) {
                fprintf(stderr, "'%s' ", deme->ancestors[j]->name);
            }
            fprintf(stderr, "\n");
        }
        
        // Process epochs (population size changes)
        for (int j = 0; j < deme->n_epochs; j++) {
            struct demes_epoch *epoch = &deme->epochs[j];
            
            // Check for unsupported features
            if (epoch->selfing_rate > 0.0) {
                fprintf(stderr, "\nError: Selfing is not supported in discoal.\n");
                fprintf(stderr, "Population '%s' epoch %d has selfing_rate = %g.\n", 
                        deme->name, j, epoch->selfing_rate);
                fprintf(stderr, "Please remove selfing from your demes model.\n");
                return -1;
            }
            
            if (epoch->cloning_rate > 0.0) {
                fprintf(stderr, "\nError: Cloning is not supported in discoal.\n");
                fprintf(stderr, "Population '%s' epoch %d has cloning_rate = %g.\n", 
                        deme->name, j, epoch->cloning_rate);
                fprintf(stderr, "Please remove cloning from your demes model.\n");
                return -1;
            }
            
            // Get epoch boundaries
            double startTime = (j == 0) ? deme->start_time : deme->epochs[j-1].end_time;
            double endTime = epoch->end_time;
            
            // For constant size epochs, we need an event at the start of each epoch
            // (except the first epoch which extends to infinity in the past)
            if (j > 0) {
                double coalTime = demesTimeToCoalTime(startTime, graph->generation_time, N);
                // Use the size from the PREVIOUS epoch since we're transitioning FROM that size
                struct demes_epoch *prev_epoch = &deme->epochs[j-1];
                double prevSize = (prev_epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL) 
                                 ? prev_epoch->end_size : prev_epoch->start_size;
                double coalSize = demesSizeToCoalSize(prevSize, N);
                
                fprintf(stderr, "  Epoch %d->%d transition at time %g:\n", j-1, j, startTime);
                fprintf(stderr, "    Previous epoch size: %g -> Discoal relative size: %g\n", prevSize, coalSize);
                fprintf(stderr, "    Demes time: %g -> Discoal coalescent time: %g\n", startTime, coalTime);
                fprintf(stderr, "    Creating size change event: time=%g, pop=%d, size=%g\n", coalTime, popID, coalSize);
                
                (*events)[*eventNumber].time = coalTime;
                (*events)[*eventNumber].popID = popID;
                (*events)[*eventNumber].popnSize = coalSize;
                (*events)[*eventNumber].type = 'n';  // Size change
                (*eventNumber)++;
            }
            
            // Handle exponential growth within epoch.
            // Emit a 'g' event at the more-recent (smaller-internal-time) boundary
            // of the epoch with the forward-time per-generation growth rate
            // converted to discoal's internal alpha (4N-scaled).  Phase 4 wires
            // 'g' to SHAPE_EXPONENTIAL.  See -eg CLI parsing for the convention.
            if (epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL &&
                epoch->start_size != epoch->end_size) {
                /* Forward-time per-generation rate.  start_time is the older
                 * boundary (forward-time start), end_time is the more recent
                 * (forward-time end).  end_size > start_size => positive alpha
                 * (forward growth, past was smaller).  Divide the time delta
                 * by generation_time to land in generations regardless of
                 * the graph's time_units. */
                double dt_gens = (startTime - endTime) / graph->generation_time;
                double alpha_per_gen = log(epoch->end_size / epoch->start_size)
                                       / dt_gens;
                /* Convert to discoal internal alpha (4N-scaled).  Mirrors the
                 * Phase 4b parity-validated alpha_per_gen = alpha_internal/(2N)
                 * inverse. */
                double alpha_internal = alpha_per_gen * 2.0 * N;
                double t_internal = demesTimeToCoalTime(epoch->end_time,
                                                        graph->generation_time, N);

                ensureDemesEventsCapacity(events, eventsCapacity, *eventNumber + 1);
                (*events)[*eventNumber].time = t_internal;
                (*events)[*eventNumber].popID = popID;
                (*events)[*eventNumber].popnSize = alpha_internal;  /* alpha stored in popnSize, matches -eg */
                (*events)[*eventNumber].type = 'g';
                (*eventNumber)++;

                fprintf(stderr, "  Epoch %d exponential: pop %d, %g -> %g, alpha=%g per gen (alpha_internal=%g) at t_internal=%g\n",
                        j, popID, epoch->start_size, epoch->end_size,
                        alpha_per_gen, alpha_internal, t_internal);
            }

            // Linear growth: any non-EXPONENTIAL size_function with differing
            // start_size and end_size.  demes-c may not expose a LINEAR enum
            // but the demes spec recognises this case.  Emit an 'l' event at
            // the more-recent boundary; Phase 6 wires 'l' to SHAPE_LINEAR.
            if (epoch->size_function != DEMES_SIZE_FUNCTION_EXPONENTIAL &&
                epoch->start_size != epoch->end_size) {
                /* Forward-time linear growth rate per generation: positive when
                 * end_size > start_size (forward growth, past was smaller).
                 * (start_time - end_time) is positive (older - more recent).
                 * Divide by generation_time to handle non-generation time_units. */
                double dt_gens = (startTime - endTime) / graph->generation_time;
                double gamma_per_gen = (epoch->end_size - epoch->start_size)
                                       / dt_gens;
                double gamma_internal = gamma_per_gen * 2.0 * N;
                double t_internal = demesTimeToCoalTime(epoch->end_time,
                                                        graph->generation_time, N);

                ensureDemesEventsCapacity(events, eventsCapacity, *eventNumber + 1);
                (*events)[*eventNumber].time = t_internal;
                (*events)[*eventNumber].popID = popID;
                (*events)[*eventNumber].popnSize = gamma_internal;  /* gamma stored in popnSize, matches -el */
                (*events)[*eventNumber].type = 'l';
                (*eventNumber)++;

                fprintf(stderr, "  Epoch %d linear: pop %d, %g -> %g, gamma=%g per gen (gamma_internal=%g) at t_internal=%g\n",
                        j, popID, epoch->start_size, epoch->end_size,
                        gamma_per_gen, gamma_internal, t_internal);
            }
        }
        
        // Handle population splits/merges
        if (deme->n_ancestors > 0) {
            double coalTime = demesTimeToCoalTime(deme->start_time, graph->generation_time, N);
            
            if (deme->n_ancestors == 1) {
                // Simple split from one ancestor
                int ancestorID = findPopulationIndex(graph, (char*)deme->ancestors[0]->name);
                if (ancestorID < 0) {
                    fprintf(stderr, "Error: Could not find ancestor population %s\n", 
                            deme->ancestors[0]->name);
                    return -1;
                }
                
                fprintf(stderr, "Creating split event: time=%g, new pop=%d joins ancestor pop=%d\n",
                        coalTime, popID, ancestorID);
                (*events)[*eventNumber].time = coalTime;
                (*events)[*eventNumber].popID = popID;  // New population
                (*events)[*eventNumber].popID2 = ancestorID;  // Source population
                (*events)[*eventNumber].type = 'p';  // Population join (backwards in time)
                (*eventNumber)++;
            } else {
                // Admixture from multiple ancestors
                // Validate proportions sum to 1.0 (demes-c should enforce this, but be safe)
                double propSum = 0.0;
                for (int k = 0; k < deme->n_ancestors; k++) {
                    propSum += deme->proportions[k];
                }
                if (fabs(propSum - 1.0) > 1e-9) {
                    fprintf(stderr, "Error: Ancestry proportions for population '%s' sum to %g, not 1.0\n", 
                            deme->name, propSum);
                    return -1;
                }
                
                for (int k = 0; k < deme->n_ancestors; k++) {
                    int ancestorID = findPopulationIndex(graph, (char*)deme->ancestors[k]->name);
                    if (ancestorID < 0) {
                        fprintf(stderr, "Error: Could not find ancestor population %s\n", 
                                deme->ancestors[k]->name);
                        return -1;
                    }
                    
                    (*events)[*eventNumber].time = coalTime;
                    (*events)[*eventNumber].popID = popID;  // New population
                    (*events)[*eventNumber].popID2 = ancestorID;  // Source population
                    (*events)[*eventNumber].admixProp = deme->proportions[k];
                    (*events)[*eventNumber].type = 'A';  // Admixture
                    (*eventNumber)++;
                }
            }
        }
    }
    
    // Process migrations: emit interval-based 'm' events (issue #82).
    //
    // Each unidirectional migration window in graph->migrations is a tuple
    // (src, dst, t_start, t_end, rate) with t_start > t_end (both in the
    // past, demes convention: end_time is the more-recent boundary).
    //
    // For each ordered pair (src, dst) we treat the rate as the sum of all
    // windows that cover a given time t (this matches demes semantics; in
    // practice a well-formed demes graph has at most one window per
    // ordered pair active at a time, but we sum to be safe).
    //
    // Strategy:
    //   A) collect all windows in their internal time units.
    //   B) write the active-at-t=0 rate directly into migMatConst[src][dst]
    //      so the runtime picks it up via its standard migMatConst -> migMat
    //      copy at the start of dropChunks. This replaces the old back-
    //      derivation hack that scanned 'M' events.
    //   C) for every distinct boundary time > 0 (start or end of any window),
    //      determine the active rate for each ordered pair just-after the
    //      boundary (i.e., for time slightly older than the boundary, going
    //      further into the past). If that differs from the rate just-before
    //      the boundary, emit one 'm' event setting the post-boundary rate.
    fprintf(stderr, "\nMigrations: %zu\n", graph->n_migrations);
    if (graph->n_migrations > 0) {
        struct mig_window *windows = (struct mig_window *)malloc(
            graph->n_migrations * sizeof(struct mig_window));
        if (windows == NULL) {
            fprintf(stderr, "Error: Failed to allocate migration windows\n");
            return -1;
        }
        int n_windows = 0;
        for (int i = 0; i < (int)graph->n_migrations; i++) {
            struct demes_migration *mig = &graph->migrations[i];
            int sourceID = findPopulationIndex(graph, (char*)mig->source->name);
            int destID   = findPopulationIndex(graph, (char*)mig->dest->name);
            if (sourceID < 0 || destID < 0) {
                fprintf(stderr, "Error: Could not find population for migration\n");
                free(windows);
                return -1;
            }
            double t_start = demesTimeToCoalTime(mig->start_time, graph->generation_time, N);
            double t_end   = demesTimeToCoalTime(mig->end_time,   graph->generation_time, N);
            double rate    = 4.0 * N * mig->rate;
            fprintf(stderr, "Migration %d: %s (pop%d) -> %s (pop%d)\n",
                    i, mig->source->name, sourceID, mig->dest->name, destID);
            fprintf(stderr, "  Time interval (gens): %g to %g\n", mig->start_time, mig->end_time);
            fprintf(stderr, "  Internal interval:    %g to %g\n", t_start, t_end);
            fprintf(stderr, "  Demes rate: %g per gen, scaled (4Nm): %g\n", mig->rate, rate);
            windows[n_windows].src     = sourceID;
            windows[n_windows].dst     = destID;
            windows[n_windows].t_start = t_start;
            windows[n_windows].t_end   = t_end;
            windows[n_windows].rate    = rate;
            n_windows++;
        }

        // Step B: write the t=0 active matrix directly into migMatConst.
        // The active rate at t=0 for (src, dst) is the sum of rates from
        // windows whose more-recent boundary is at t=0 (i.e., they cover
        // the present).
        for (int w = 0; w < n_windows; w++) {
            if (windows[w].t_end == 0.0) {
                migMatConst[windows[w].src][windows[w].dst] += windows[w].rate;
                fprintf(stderr, "  Active at t=0: pop%d -> pop%d = %g (added to migMatConst)\n",
                        windows[w].src, windows[w].dst, windows[w].rate);
            }
        }

        // Step C: collect boundary times > 0 (start_time and end_time of
        // every window) and dedupe.
        double *boundaries = (double *)malloc(2 * n_windows * sizeof(double));
        if (boundaries == NULL) {
            fprintf(stderr, "Error: Failed to allocate migration boundaries\n");
            free(windows);
            return -1;
        }
        int n_bounds = 0;
        for (int w = 0; w < n_windows; w++) {
            if (windows[w].t_start > 0.0) boundaries[n_bounds++] = windows[w].t_start;
            if (windows[w].t_end   > 0.0) boundaries[n_bounds++] = windows[w].t_end;
        }
        qsort(boundaries, n_bounds, sizeof(double), compareDoublesAscending);
        int n_unique = 0;
        for (int i = 0; i < n_bounds; i++) {
            if (i == 0 || boundaries[i] != boundaries[i-1]) {
                boundaries[n_unique++] = boundaries[i];
            }
        }

        // Worst case: at each boundary we emit one 'm' event per ordered
        // (src, dst) pair. Reserve capacity accordingly.
        int worst_extra = n_unique * (*npops) * (*npops);
        ensureDemesEventsCapacity(events, eventsCapacity, *eventNumber + worst_extra + 10);

        // Step D: for each boundary, for each ordered pair, emit an 'm'
        // event when the post-boundary rate differs from the pre-boundary
        // rate.
        //
        // Convention: discoal time grows going backward into the past.
        // A window with internal interval [t_end, t_start) is active for
        // times t in [t_end, t_start). Just-after a boundary (going further
        // into the past) means t slightly larger than the boundary.
        for (int b = 0; b < n_unique; b++) {
            double t = boundaries[b];
            for (int src = 0; src < *npops; src++) {
                for (int dst = 0; dst < *npops; dst++) {
                    if (src == dst) continue;
                    double rate_after  = 0.0;  // for t' slightly > t (older)
                    double rate_before = 0.0;  // for t' slightly < t (more recent)
                    for (int w = 0; w < n_windows; w++) {
                        if (windows[w].src != src || windows[w].dst != dst) continue;
                        // Window covers [t_end, t_start). Just-after t means
                        // t in [t_end, t_start), i.e. t_end <= t < t_start.
                        if (windows[w].t_end <= t && t < windows[w].t_start) {
                            rate_after += windows[w].rate;
                        }
                        // Just-before t means t' < t with t' in [t_end, t_start),
                        // i.e. t_end < t and t <= t_start (window's recent
                        // boundary is older than t', and t' < t_start).
                        if (windows[w].t_end < t && t <= windows[w].t_start) {
                            rate_before += windows[w].rate;
                        }
                    }
                    if (rate_after != rate_before) {
                        (*events)[*eventNumber].time     = t;
                        (*events)[*eventNumber].popID2   = src;
                        (*events)[*eventNumber].popID    = dst;
                        (*events)[*eventNumber].popnSize = rate_after;
                        (*events)[*eventNumber].type     = 'm';
                        (*eventNumber)++;
                        fprintf(stderr, "  Boundary t=%g: pop%d -> pop%d rate %g -> %g ('m' event emitted)\n",
                                t, src, dst, rate_before, rate_after);
                    }
                }
            }
        }

        free(boundaries);
        free(windows);
    }
    
    // Process pulses (admixture events)
    fprintf(stderr, "\nPulses: %zu\n", graph->n_pulses);
    for (int i = 0; i < graph->n_pulses; i++) {
        struct demes_pulse *pulse = &graph->pulses[i];
        int destID = findPopulationIndex(graph, (char*)pulse->dest->name);
        
        if (destID < 0) {
            fprintf(stderr, "Error: Could not find destination population '%s' for pulse\n", 
                    pulse->dest->name);
            return -1;
        }
        
        // Validate pulse proportions sum to <= 1.0
        double propSum = 0.0;
        for (int j = 0; j < pulse->n_sources; j++) {
            propSum += pulse->proportions[j];
        }
        if (propSum > 1.0 + 1e-9) {
            fprintf(stderr, "Error: Pulse proportions sum to %g, which exceeds 1.0\n", propSum);
            return -1;
        }
        
        fprintf(stderr, "Pulse %d at time %g into %s (pop%d)\n", 
                i, pulse->time, pulse->dest->name, destID);
        
        for (int j = 0; j < pulse->n_sources; j++) {
            int sourceID = findPopulationIndex(graph, (char*)pulse->sources[j]->name);
            if (sourceID < 0) {
                fprintf(stderr, "Error: Could not find source population '%s' for pulse\n",
                        pulse->sources[j]->name);
                return -1;
            }
            
            fprintf(stderr, "  Source: %s (pop%d), proportion: %g\n", 
                    pulse->sources[j]->name, sourceID, pulse->proportions[j]);
            
            (*events)[*eventNumber].time = demesTimeToCoalTime(pulse->time, graph->generation_time, N);
            (*events)[*eventNumber].popID = destID;
            (*events)[*eventNumber].popID2 = sourceID;
            (*events)[*eventNumber].admixProp = pulse->proportions[j];
            (*events)[*eventNumber].type = 'a';  // Pulse (lowercase for single generation)
            (*eventNumber)++;
        }
    }
    
    // Sort events by time (most recent first)
    qsort(&(*events)[1], *eventNumber - 1, sizeof(struct event), compareEventsByTime);
    
    // Set initial population sizes
    fprintf(stderr, "\nPresent-day population sizes:\n");
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        if (deme->n_epochs > 0) {
            // Check if this population exists at present day
            struct demes_epoch *last_epoch = &deme->epochs[deme->n_epochs-1];
            if (last_epoch->end_time != 0.0) {
                // This population doesn't exist at present day
                fprintf(stderr, "  Pop %d (%s): not present at time 0 (ends at %g)\n", 
                        i, deme->name, last_epoch->end_time);
                currentSize[i] = 0.0;  // Will be handled by events
                continue;
            }
            
            // Use the end size of the last epoch (most recent time)
            double presentSize = (last_epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL) 
                                ? last_epoch->end_size : last_epoch->start_size;
            currentSize[i] = demesSizeToCoalSize(presentSize, N);
            fprintf(stderr, "  Pop %d (%s): %g (relative to N=%g: %g)\n", 
                    i, deme->name, presentSize, N, currentSize[i]);
        } else {
            currentSize[i] = 1.0;  // Default
            fprintf(stderr, "  Pop %d (%s): default size 1.0\n", i, deme->name);
        }
    }
    
    fprintf(stderr, "\nTotal events created: %d\n", *eventNumber - 1);
    
    // Debug: Show what currentSize array contains
    fprintf(stderr, "\nDISCOAL currentSize array after conversion:\n");
    for (int i = 0; i < *npops; i++) {
        fprintf(stderr, "  currentSize[%d] = %f\n", i, currentSize[i]);
    }
    
    fprintf(stderr, "=== END DEMES CONVERSION DEBUG ===\n\n");
    
    return 0;
}