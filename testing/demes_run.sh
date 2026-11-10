#!/bin/bash
# demes_run DEMES_FILE nSamples nReplicates nSites [options] [-p k n1 ... nk]
#
# Run discoal on a demes file through a YAML config (-Y). The config holds
# the simulation size, the demes file and the samples per deme (from -p;
# all samples in one deme without -p); the other options go on the
# command line after -Y, where they replace the YAML values (-t, -r, -d)
# or add events. The config is written next to the demes file.
# DISCOAL must name the discoal binary.
demes_run() {
    local demes=$1 nsamp=$2 nreps=$3 nsites=$4
    shift 4
    local sizes="" rest=() i k
    while [ $# -gt 0 ]; do
        if [ "$1" = "-p" ]; then
            k=$2
            shift 2
            for ((i = 0; i < k; i++)); do
                sizes="$sizes${sizes:+, }$1"
                shift
            done
        else
            rest+=("$1")
            shift
        fi
    done
    [ -n "$sizes" ] || sizes=$nsamp
    local config="$demes.config.yaml"
    # A single-quoted YAML string, so any path works ('' is a quote).
    local quoted=${demes//\'/\'\'}
    cat > "$config" <<CONFIG
simulation: {sample_size: $nsamp, num_replicates: $nreps, num_sites: $nsites}
genetics: {mutation_rate: 0, recombination_rate: 0}
demography: {deme_sample_size: [$sizes], demes_filename: '$quoted'}
CONFIG
    "$DISCOAL" -Y "$config" "${rest[@]}"
}
