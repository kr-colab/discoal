#!/usr/bin/env python3
"""
Convert demes YAML files to equivalent discoal command line arguments.

This tool reads a demes demographic model specification and outputs the 
equivalent discoal command line arguments that would produce the same
demographic history.

Usage:
    python demes_to_discoal.py model.yaml [--N0 1000000] [--samples 20 20]
"""

import sys
import argparse
import math
from pathlib import Path

# Add the demes-python directory to the path for importing
sys.path.insert(0, str(Path(__file__).parent.parent / "reference_code" / "demes-python"))

try:
    import demes
except ImportError:
    print("Error: demes-python library not found. Please check that reference_code/demes-python exists.")
    sys.exit(1)


def float_str(a: float) -> str:
    """
    Convert float to string, for use in command line arguments.
    Handles negative numbers properly for argparse.
    """
    if a < 0:
        return format(a, ".10f")
    else:
        return str(a)


def demes_to_discoal(demes_file: str, N0: float = 1000000, samples=None) -> str:
    """
    Convert a demes YAML file to discoal command line arguments.
    
    Args:
        demes_file: Path to the demes YAML file
        N0: Reference population size (default: 1000000, matching discoal's EFFECTIVE_POPN_SIZE)
        samples: List of sample sizes for each present-day population
        
    Returns:
        String containing discoal command line arguments
    """
    
    # Load the demes graph
    try:
        graph = demes.load(demes_file)
    except Exception as e:
        raise ValueError(f"Error loading demes file '{demes_file}': {e}")
    
    # Convert to generations if not already
    graph = graph.in_generations()
    
    # Build the command line arguments
    cmd_parts = []
    
    # Count present-day populations (those that exist at time 0)
    present_day_demes = [deme for deme in graph.demes if deme.end_time == 0]
    num_present_pops = len(present_day_demes)
    
    # Add population structure (-p flag)
    if num_present_pops > 1:
        if samples is None:
            # Use placeholder sample sizes - user must specify
            samples = [10] * num_present_pops
        elif len(samples) != num_present_pops:
            raise ValueError(f"Number of sample sizes ({len(samples)}) doesn't match "
                           f"number of present-day populations ({num_present_pops})")
        
        pop_spec = f"-p {num_present_pops} " + " ".join(str(s) for s in samples)
        cmd_parts.append(pop_spec)
    
    # Collect all demographic events
    events = []
    
    # Create a mapping from deme names to population indices (0-based)
    # Only for present-day populations
    deme_to_index = {}
    for i, deme in enumerate(present_day_demes):
        deme_to_index[deme.name] = i
    
    # Also map all demes for ancestry tracking
    all_deme_to_index = {}
    for i, deme in enumerate(graph.demes):
        all_deme_to_index[deme.name] = i
    
    # Process population size changes
    for deme in graph.demes:
        if deme.name not in deme_to_index:
            continue  # Skip non-present-day populations for now
            
        pop_index = deme_to_index[deme.name]
        
        # Process epochs from most recent to oldest
        for epoch in reversed(deme.epochs):
            if epoch.end_time > 0:  # Don't create event for present day
                time_discoal = epoch.end_time / (4 * N0)  # Convert to discoal time units
                size_relative = epoch.end_size / N0  # Relative to N0
                
                # Add size change event
                events.append((time_discoal, f"-en {float_str(time_discoal)} {pop_index} {float_str(size_relative)}"))
    
    # Process migrations
    for migration in graph.migrations:
        # Check if both source and dest are present-day populations
        if migration.source in deme_to_index and migration.dest in deme_to_index:
            source_idx = deme_to_index[migration.source]
            dest_idx = deme_to_index[migration.dest]
            
            # Convert migration rate from per-generation to 4Nm units
            rate_discoal = migration.rate * 4 * N0
            
            # Migration start (going backwards in time, this turns migration on)
            if not math.isinf(migration.start_time):
                time_discoal = migration.start_time / (4 * N0)
                events.append((time_discoal, f"-em {float_str(time_discoal)} {dest_idx} {source_idx} 0"))
            
            # Migration end (going backwards in time, this turns migration off) 
            time_discoal = migration.end_time / (4 * N0)
            events.append((time_discoal, f"-em {float_str(time_discoal)} {dest_idx} {source_idx} {float_str(rate_discoal)}"))
    
    # Process population splits/joins and admixture
    for deme in graph.demes:
        if len(deme.ancestors) == 1 and deme.name in deme_to_index:
            # Simple population split
            ancestor_name = deme.ancestors[0]
            if ancestor_name in deme_to_index:
                # This is a join event going backwards in time
                child_idx = deme_to_index[deme.name] 
                parent_idx = deme_to_index[ancestor_name]
                time_discoal = deme.start_time / (4 * N0)
                
                events.append((time_discoal, f"-ej {float_str(time_discoal)} {child_idx} {parent_idx}"))
                
                # Also add size change for the ancestor at this time if needed
                ancestor_deme = graph[ancestor_name]
                # Find the epoch that contains this time
                for epoch in ancestor_deme.epochs:
                    if epoch.end_time <= deme.start_time < epoch.start_time or math.isinf(epoch.start_time):
                        ancestor_size = epoch.start_size if hasattr(epoch, 'start_size') else epoch.end_size
                        size_relative = ancestor_size / N0
                        events.append((time_discoal, f"-en {float_str(time_discoal)} {parent_idx} {float_str(size_relative)}"))
                        break
        
        elif len(deme.ancestors) > 1:
            # Admixture event - this is more complex and may need -es/-ej combinations
            print(f"Warning: Admixture events not fully implemented yet. Deme '{deme.name}' has multiple ancestors.", file=sys.stderr)
    
    # Process pulses (admixture events)
    for pulse in graph.pulses:
        print(f"Warning: Pulse migrations not implemented yet. Pulse from {pulse.sources} to {pulse.dest}.", file=sys.stderr)
    
    # Sort events by time (most recent first, which is how discoal processes them)
    events.sort(key=lambda x: x[0])
    
    # Add events to command
    for time, event_str in events:
        cmd_parts.append(event_str)
    
    return " ".join(cmd_parts)


def main():
    parser = argparse.ArgumentParser(
        description="Convert demes YAML files to discoal command line arguments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Convert a simple two-population model
    python demes_to_discoal.py two_pop_split.yaml --samples 20 20
    
    # Convert with custom reference population size
    python demes_to_discoal.py model.yaml --N0 50000 --samples 10 10 10
    
    # Convert single population model
    python demes_to_discoal.py single_pop.yaml --samples 30
        """
    )
    
    parser.add_argument("demes_file", help="Path to the demes YAML file")
    parser.add_argument("--N0", type=float, default=1000000,
                        help="Reference population size (default: 1000000)")
    parser.add_argument("--samples", type=int, nargs='+',
                        help="Sample sizes for each present-day population")
    parser.add_argument("--output", "-o", help="Output file (default: stdout)")
    
    args = parser.parse_args()
    
    try:
        discoal_cmd = demes_to_discoal(args.demes_file, args.N0, args.samples)
        
        if args.output:
            with open(args.output, 'w') as f:
                f.write(discoal_cmd + '\n')
            print(f"Discoal command written to {args.output}")
        else:
            print(discoal_cmd)
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()