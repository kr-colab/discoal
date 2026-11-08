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

// Index of a deme in graph->demes. demes-c stores ancestors, pulse sources
// and migration endpoints as pointers into that array.
static int demeIndex(const struct demes_graph *graph, const struct demes_deme *deme) {
    return (int)(deme - graph->demes);
}

// A deme exists at time 0 if its most recent epoch ends there. demes-c
// gives every deme at least one epoch.
static int demeIsPresent(const struct demes_deme *deme) {
    return deme->epochs[deme->n_epochs - 1].end_time == 0.0;
}

// Append a zeroed event and return it, growing the array as needed.
static event *appendEvent(event **events, int *eventNumber, int *eventsCapacity,
                          double time, char type, int popID) {
    ensureDemesEventsCapacity(events, eventsCapacity, *eventNumber + 1);
    event *e = &(*events)[(*eventNumber)++];
    memset(e, 0, sizeof *e);
    e->time = time;
    e->type = type;
    e->popID = popID;
    return e;
}

// Move the lineages of popID to n target demes with the given proportions,
// as a chain of events at one time. An 'a' event moves a lineage to one of
// two populations, so target j gets an 'a' event that moves a lineage there
// with probability props[j] / (share not yet assigned) and else keeps it in
// popID; msprime uses the same sequential proportions. If lastIsMerge (the
// proportions sum to 1), the last target takes the remaining lineages with a
// 'p' event. The event sort is stable, so the chain runs in this order.
static void appendSplitChain(event **events, int *eventNumber, int *eventsCapacity,
                             const struct demes_graph *graph, double time, int popID,
                             struct demes_deme **targets, const double *props,
                             size_t n, int lastIsMerge) {
    double remaining = 1.0;
    for (size_t j = 0; j < n; j++) {
        int target = demeIndex(graph, targets[j]);
        if (lastIsMerge && j == n - 1) {
            appendEvent(events, eventNumber, eventsCapacity, time, 'p', popID)->popID2 = target;
            break;
        }
        event *e = appendEvent(events, eventNumber, eventsCapacity, time, 'a', popID);
        e->popID2 = target;
        e->popID3 = popID;
        e->admixProp = props[j] / remaining;
        remaining -= props[j];
    }
}

// Convert demes time to discoal internal time (units of 2N generations).
// demes_time is in the units the graph specifies; for
// time_units="generations" generation_time is 1 (demes validates this), and
// for time_units="years" generation_time is the years per generation, so
// dividing yields generations either way.
static double demesTimeToCoalTime(double demes_time, double generation_time, double N) {
    return demes_time / generation_time / (2.0 * N);
}

// Convert a per-generation rate to a rate per unit of internal time
// (2N generations), the unit of Shape.rate_param.
static double demesRateToInternal(double rate_per_gen, double N) {
    return rate_per_gen * 2.0 * N;
}

// Convert demes size to coalescent size (relative to ancestral N)
static double demesSizeToCoalSize(double demes_size, double N) {
    return demes_size / N;
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
// src and dst are backward-time directions, as in discoal's migMat[src][dst]:
// lineages in src move to dst. A demes migration from source to dest
// (forward in time) is therefore src = dest, dst = source.
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
    
    // deme_sample_size is indexed by the order of demes in the file, extinct
    // demes included. Samples in a deme that does not exist at time 0 would
    // all coalesce at once (its size is 0), so reject them.
    for (int i = 0; i < MAXPOPS; i++) {
        if (sampleSizes[i] <= 0) continue;
        if (i >= (int)graph->n_demes) {
            fprintf(stderr, "Error: deme_sample_size[%d] = %d, but the demes file "
                    "has only %zu demes\n", i, sampleSizes[i], graph->n_demes);
            demes_graph_free(graph);
            return -1;
        }
        struct demes_deme *deme = &graph->demes[i];
        if (!demeIsPresent(deme)) {
            fprintf(stderr, "Error: deme_sample_size[%d] = %d puts samples in deme "
                    "'%s', which ends %g time units ago and does not exist at "
                    "time 0. deme_sample_size follows the order of demes in the "
                    "demes file, extinct demes included; give such demes 0.\n",
                    i, sampleSizes[i], deme->name,
                    deme->epochs[deme->n_epochs - 1].end_time);
            demes_graph_free(graph);
            return -1;
        }
    }

    // The reference N for all scaling is the present-day size of the first
    // deme that exists at time 0.
    double referenceN = N;  // default if no deme exists at time 0
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        if (demeIsPresent(deme)) {
            referenceN = deme->epochs[deme->n_epochs - 1].end_size;
            VERBOSE_LOG("Using population '%s' (index %d) as reference with N=%g\n",
                    deme->name, i, referenceN);
            break;
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
        if (demeIsPresent(&graph->demes[i])) presentDayPops++;
    }
    
    if (presentDayPops > 1) {
        VERBOSE_LOG("\nValidating population connectivity (%d present-day populations)...\n", presentDayPops);
        // Check if there's any migration
        int hasMigration = 0;
        for (int i = 0; i < graph->n_migrations; i++) {
            if (graph->migrations[i].rate > 0) {
                hasMigration = 1;
                break;
            }
        }
        VERBOSE_LOG("Has migration: %s\n", hasMigration ? "yes" : "no");
        
        // Check if all populations eventually merge
        // We need to check if all present-day populations trace back to a common ancestor
        int hasCommonAncestor = 0;
        
        // Check specifically if present-day populations have ancestors
        for (int i = 0; i < graph->n_demes; i++) {
            struct demes_deme *deme = &graph->demes[i];
            if (demeIsPresent(deme) && deme->n_ancestors > 0) {
                hasCommonAncestor = 1;
                break;
            }
        }
        VERBOSE_LOG("Present-day populations have ancestors: %s\n", hasCommonAncestor ? "yes" : "no");
        
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
            VERBOSE_LOG("\nNote: Populations will coalesce through common ancestors.\n");
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

// Convert a demes graph into discoal events. Times are internal (units of
// 2 N_ref generations), sizes are relative to N_ref, and migration rates are
// 4 N_ref m.
//
// Events at the same time run in the order emitted here (sortEventArray is
// stable). The order follows msprime's Demography.from_demes:
//   1. For each deme, for each epoch that ends before the present: the end
//      size ('n'), then the epoch's growth ('g' or 'l'), which anchors at
//      that size. Then the deme's split from a single ancestor ('p'). demes
//      lists ancestors before descendants.
//   2. Migration-rate changes ('m'), after the splits.
//   3. Pulses, in reverse file order: demes applies same-time pulses in file
//      order forward in time, so backward in time they run reversed.
//   4. Demes with several ancestors (a chain of 'a' events and a final 'p'),
//      after the pulses, so a pulse at the same time runs first.
int convertDemesToEvents(struct demes_graph *graph, event **events, int *eventNumber,
                        int *eventsCapacity, double *currentSize, int *npops,
                        int *sampleSizes, double N) {
    (void)sampleSizes;
    double gt = graph->generation_time;

    *npops = graph->n_demes;
    if (*npops > MAXPOPS) {
        fprintf(stderr, "Error: Demes file has %d populations, but discoal supports at most %d\n",
                *npops, MAXPOPS);
        return -1;
    }

    VERBOSE_LOG("\n=== DEMES TO DISCOAL CONVERSION DEBUG ===\n");
    VERBOSE_LOG("Reference N (first present-day deme): %g\n", N);
    VERBOSE_LOG("Number of populations: %d\n", *npops);
    VERBOSE_LOG("Generation time: %g\n", gt);

    // 1. Epoch sizes and growth, then single-ancestor splits.
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];

        VERBOSE_LOG("\nPopulation %d: '%s'\n", i, deme->name);
        VERBOSE_LOG("  Start time: %g generations\n", deme->start_time);
        VERBOSE_LOG("  Number of epochs: %zu\n", deme->n_epochs);
        if (deme->n_ancestors > 0) {
            VERBOSE_LOG("  Ancestors: ");
            for (size_t j = 0; j < deme->n_ancestors; j++) {
                VERBOSE_LOG("'%s' ", deme->ancestors[j]->name);
            }
            VERBOSE_LOG("\n");
        }

        for (int j = 0; j < deme->n_epochs; j++) {
            struct demes_epoch *epoch = &deme->epochs[j];

            if (epoch->selfing_rate > 0.0 || epoch->cloning_rate > 0.0) {
                fprintf(stderr, "\nError: Selfing and cloning are not supported in discoal.\n");
                fprintf(stderr, "Population '%s' epoch %d has selfing_rate = %g, cloning_rate = %g.\n",
                        deme->name, j, epoch->selfing_rate, epoch->cloning_rate);
                return -1;
            }

            double startTime = (j == 0) ? deme->start_time : deme->epochs[j-1].end_time;
            double endTime = epoch->end_time;
            double tEnd = demesTimeToCoalTime(endTime, gt, N);

            // Going backward in time, the size jumps to the epoch's end size
            // where the epoch ends. The present-day epoch (end_time 0) takes
            // its size from currentSize instead.
            if (endTime > 0.0) {
                double coalSize = demesSizeToCoalSize(epoch->end_size, N);
                VERBOSE_LOG("  Epoch %d ends at time %g: size %g (relative %g) at discoal time %g\n",
                            j, endTime, epoch->end_size, coalSize, tEnd);
                appendEvent(events, eventNumber, eventsCapacity, tEnd, 'n', i)->popnSize = coalSize;
            }

            // Growth over the epoch, as a forward-time rate per generation:
            // positive means the past was smaller. Dividing by generation_time
            // gives generations for any time_units. Exponential epochs give a
            // 'g' event (rate of log size); any other non-constant epoch gives
            // an 'l' event (rate of size relative to N). demes-c accepts only
            // constant and exponential epochs, so 'l' does not occur today.
            if (epoch->start_size != epoch->end_size) {
                double dt_gens = (startTime - endTime) / gt;
                int exponential = epoch->size_function == DEMES_SIZE_FUNCTION_EXPONENTIAL;
                double rate_per_gen = exponential
                    ? log(epoch->end_size / epoch->start_size) / dt_gens
                    : demesSizeToCoalSize(epoch->end_size - epoch->start_size, N) / dt_gens;
                double rate_internal = demesRateToInternal(rate_per_gen, N);
                appendEvent(events, eventNumber, eventsCapacity, tEnd,
                            exponential ? 'g' : 'l', i)->popnSize = rate_internal;
                VERBOSE_LOG("  Epoch %d %s: %g -> %g, rate %g per gen (%g internal) at discoal time %g\n",
                            j, exponential ? "exponential" : "linear", epoch->start_size,
                            epoch->end_size, rate_per_gen, rate_internal, tEnd);
            }
        }

        if (deme->n_ancestors == 1) {
            double coalTime = demesTimeToCoalTime(deme->start_time, gt, N);
            int ancestorID = demeIndex(graph, deme->ancestors[0]);
            VERBOSE_LOG("Creating split event: time=%g, new pop=%d joins ancestor pop=%d\n",
                        coalTime, i, ancestorID);
            appendEvent(events, eventNumber, eventsCapacity, coalTime, 'p', i)->popID2 = ancestorID;
        }
    }

    // 2. Migration windows.
    //
    // Each unidirectional window is (src, dst, t_start, t_end, rate) with
    // t_start > t_end. For each ordered pair the rate at time t is the sum of
    // the windows that cover t.
    //   A) Collect all windows in internal time units.
    //   B) Write the rate active at t = 0 into migMatConst[src][dst].
    //   C) At every distinct boundary time > 0, emit an 'm' event for each
    //      ordered pair whose rate just older than the boundary differs from
    //      the rate just more recent than it.
    VERBOSE_LOG("\nMigrations: %zu\n", graph->n_migrations);
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
            int sourceID = demeIndex(graph, mig->source);
            int destID   = demeIndex(graph, mig->dest);
            double t_start = demesTimeToCoalTime(mig->start_time, gt, N);
            double t_end   = demesTimeToCoalTime(mig->end_time,   gt, N);
            double rate    = 4.0 * N * mig->rate;
            VERBOSE_LOG("Migration %d: %s (pop%d) -> %s (pop%d)\n",
                    i, mig->source->name, sourceID, mig->dest->name, destID);
            VERBOSE_LOG("  Time interval (gens): %g to %g\n", mig->start_time, mig->end_time);
            VERBOSE_LOG("  Internal interval:    %g to %g\n", t_start, t_end);
            VERBOSE_LOG("  Demes rate: %g per gen, scaled (4Nm): %g\n", mig->rate, rate);
            windows[n_windows].src     = destID;
            windows[n_windows].dst     = sourceID;
            windows[n_windows].t_start = t_start;
            windows[n_windows].t_end   = t_end;
            windows[n_windows].rate    = rate;
            n_windows++;
        }

        for (int w = 0; w < n_windows; w++) {
            if (windows[w].t_end == 0.0) {
                migMatConst[windows[w].src][windows[w].dst] += windows[w].rate;
                VERBOSE_LOG("  Active at t=0: pop%d -> pop%d = %g (added to migMatConst)\n",
                        windows[w].src, windows[w].dst, windows[w].rate);
            }
        }

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

        // A window with internal interval [t_end, t_start) covers t when
        // t_end <= t < t_start. Just older than boundary t, it is active if
        // t_end <= t < t_start; just more recent, if t_end < t <= t_start.
        for (int b = 0; b < n_unique; b++) {
            double t = boundaries[b];
            for (int src = 0; src < *npops; src++) {
                for (int dst = 0; dst < *npops; dst++) {
                    if (src == dst) continue;
                    double rate_after  = 0.0;
                    double rate_before = 0.0;
                    for (int w = 0; w < n_windows; w++) {
                        if (windows[w].src != src || windows[w].dst != dst) continue;
                        if (windows[w].t_end <= t && t < windows[w].t_start) {
                            rate_after += windows[w].rate;
                        }
                        if (windows[w].t_end < t && t <= windows[w].t_start) {
                            rate_before += windows[w].rate;
                        }
                    }
                    if (rate_after != rate_before) {
                        event *e = appendEvent(events, eventNumber, eventsCapacity, t, 'm', dst);
                        e->popID2   = src;
                        e->popnSize = rate_after;
                        VERBOSE_LOG("  Boundary t=%g: pop%d -> pop%d rate %g -> %g ('m' event emitted)\n",
                                t, src, dst, rate_before, rate_after);
                    }
                }
            }
        }

        free(boundaries);
        free(windows);
    }

    // 3. Pulses, in reverse order. Lineages that stay in the destination keep
    // their place (no final merge).
    VERBOSE_LOG("\nPulses: %zu\n", graph->n_pulses);
    for (int i = graph->n_pulses - 1; i >= 0; i--) {
        struct demes_pulse *pulse = &graph->pulses[i];
        int destID = demeIndex(graph, pulse->dest);
        VERBOSE_LOG("Pulse %d at time %g into %s (pop%d)\n",
                i, pulse->time, pulse->dest->name, destID);
        for (size_t j = 0; j < pulse->n_sources; j++) {
            VERBOSE_LOG("  Source: %s (pop%d), proportion: %g\n", pulse->sources[j]->name,
                        demeIndex(graph, pulse->sources[j]), pulse->proportions[j]);
        }
        appendSplitChain(events, eventNumber, eventsCapacity, graph,
                         demesTimeToCoalTime(pulse->time, gt, N), destID,
                         pulse->sources, pulse->proportions, pulse->n_sources, 0);
    }

    // 4. Demes with several ancestors. demes-c checks that the proportions
    // sum to 1, so the chain ends with a merge into the last ancestor.
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        if (deme->n_ancestors < 2) continue;
        appendSplitChain(events, eventNumber, eventsCapacity, graph,
                         demesTimeToCoalTime(deme->start_time, gt, N), i,
                         deme->ancestors, deme->proportions, deme->n_ancestors, 1);
    }

    // Present-day sizes. A deme that does not exist at time 0 gets size 0;
    // its size is set by the 'n' event where its last epoch ends.
    VERBOSE_LOG("\nPresent-day population sizes:\n");
    for (int i = 0; i < graph->n_demes; i++) {
        struct demes_deme *deme = &graph->demes[i];
        struct demes_epoch *last = &deme->epochs[deme->n_epochs - 1];
        if (!demeIsPresent(deme)) {
            VERBOSE_LOG("  Pop %d (%s): not present at time 0 (ends at %g)\n",
                    i, deme->name, last->end_time);
            currentSize[i] = 0.0;
            continue;
        }
        currentSize[i] = demesSizeToCoalSize(last->end_size, N);
        VERBOSE_LOG("  Pop %d (%s): %g (relative to N=%g: %g)\n",
                i, deme->name, last->end_size, N, currentSize[i]);
    }

    VERBOSE_LOG("\nTotal events created: %d\n", *eventNumber - 1);
    VERBOSE_LOG("=== END DEMES CONVERSION DEBUG ===\n\n");
    return 0;
}
