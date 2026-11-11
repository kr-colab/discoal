/* version.h - Version information for discoal
 * 
 * This file defines the version of discoal using semantic versioning.
 * Update these values when making a new release.
 */

#ifndef DISCOAL_VERSION_H
#define DISCOAL_VERSION_H

#define DISCOAL_VERSION_MAJOR 2
#define DISCOAL_VERSION_MINOR 0
#define DISCOAL_VERSION_PATCH 0

/* Version string for display */
#define DISCOAL_VERSION_STRING "2.0.0"

/* A pre-release defines DISCOAL_VERSION_STATUS (for example "beta"), which
 * the version string then carries ("2.0.0-beta"); a release defines none. */

/* Build the full version string */
#ifdef DISCOAL_VERSION_STATUS
    #define DISCOAL_VERSION DISCOAL_VERSION_STRING "-" DISCOAL_VERSION_STATUS
#else
    #define DISCOAL_VERSION DISCOAL_VERSION_STRING
#endif

#endif /* DISCOAL_VERSION_H */