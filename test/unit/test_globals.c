/*
 * test_globals.c
 * 
 * Provides minimal global variables required by discoal functions
 * when they are tested in isolation. These globals are normally
 * defined in discoal_multipop.c but need to be provided for unit tests.
 */

#include <stdlib.h>
#include "discoal.h"

// Random number generator seeds
long seed1 = 12345;
long seed2 = 67890;

// Current population size (used by some functions)
double *currentSize = NULL;

// Event management globals
event *events = NULL;
int eventNumber = 0;
int eventsCapacity = 0;

// Stub implementation of ensureEventsCapacity for testing
void ensureEventsCapacity() {
    if (events == NULL) {
        eventsCapacity = 100;
        events = (event *)malloc(eventsCapacity * sizeof(event));
    } else if (eventNumber >= eventsCapacity - 1) {
        eventsCapacity *= 2;
        events = (event *)realloc(events, eventsCapacity * sizeof(event));
    }
}