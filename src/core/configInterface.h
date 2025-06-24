#ifndef CONFIGINTERFACE_H
#define CONFIGINTERFACE_H

#include "discoal.h"
#include <yaml.h>

// Structure to hold complete simulation configuration
typedef struct {
    // Simulation settings
    int sample_size;
    int num_replicates; 
    int num_sites;
    int seed1, seed2;
    int has_seed;           // flag to track if seed was set
    
    // Population configuration
    int npops;
    int sample_sizes[MAXPOPS];
    double effective_popn_size;
    int has_populations;    // flag to track if populations were set
    
    // Genetics parameters
    double theta;           // mutation rate
    double rho;             // recombination rate
    double gamma;           // gene conversion rate
    int gc_mean;            // gene conversion tract length
    double gamma_co_ratio;  // crossover to gene conversion ratio
    int gamma_co_ratio_mode;
    int has_genetics;       // flag to track if genetics were set
    
    // Selection parameters
    double alpha;           // selection coefficient
    double sweep_site;      // position of selected site
    char sweep_mode;        // 'd', 's', or 'N'
    double tau;             // time of selection start
    double f0;              // soft sweep initial frequency
    double ua;              // beneficial mutation rate
    double partial_sweep_final_freq;
    int recur_sweep_mode;
    double recur_sweep_rate;
    int partial_sweep_mode;
    int soft_sweep_mode;
    int has_selection;      // flag to track if selection was set
    
    // Output settings
    char output_style;      // 'h' for haplotype, etc.
    int finite_output_flag;
    int hide_partial_snp;
    int tskit_output_mode;
    char tskit_output_filename[256];
    int minimal_tree_seq;
    int has_output;         // flag to track if output settings were set
    
    // Migration matrix (only used when not using demes)
    int mig_flag;
    double mig_mat_const[MAXPOPS][MAXPOPS];
    
    // Events (can be from demes file or explicit specification)
    char demes_file[256];   // path to demes file if specified
    int use_demes;          // flag indicating whether to use demes file
    
    // Explicit events list for when not using demes
    struct event explicit_events[1000];  // maximum number of explicit events
    int num_explicit_events;
    
    // Ancient sampling
    int anc_sample_flag;
    int anc_sample_size;
    
} SimulationConfig;

// Function to load configuration from YAML file
int loadConfigFile(const char *filename, SimulationConfig *config);

// Function to apply configuration to global variables (called after command line parsing)
int applyConfiguration(const SimulationConfig *config);

// Function to initialize default configuration
void initializeDefaultConfig(SimulationConfig *config);

#endif // CONFIGINTERFACE_H