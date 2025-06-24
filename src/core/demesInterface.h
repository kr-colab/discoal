#ifndef DEMES_INTERFACE_H
#define DEMES_INTERFACE_H

#include "discoal.h"

// Forward declaration
struct demes_graph;

// Function to load demographic events from a demes file
// Returns 0 on success, non-zero on error
int loadDemesFile(const char *filename, event **events, int *eventNumber, int *eventsCapacity, 
                  double *currentSize, int *npops, int *sampleSizes, double N);

// Function to convert a demes graph to discoal events
// Populates the events array and updates eventNumber
int convertDemesToEvents(struct demes_graph *graph, event **events, int *eventNumber, 
                        int *eventsCapacity, double *currentSize, int *npops, 
                        int *sampleSizes, double N);

#endif // DEMES_INTERFACE_H