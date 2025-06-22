/*
 * test_globals.c
 * 
 * Provides minimal global variables required by discoal functions
 * when they are tested in isolation. These globals are normally
 * defined in discoal_multipop.c but need to be provided for unit tests.
 */

// Random number generator seeds
long seed1 = 12345;
long seed2 = 67890;

// Current population size (used by some functions)
int currentSize = 0;