# Demes Features Not Currently Supported by discoal

This document lists demographic features that are part of the Demes specification but are not currently supported by discoal's demographic model.

## 1. Reproduction Mode Features

### Selfing Rate
- **Demes**: Supports per-epoch selfing rates (probability of self-fertilization)
- **discoal**: No support for selfing
- **Example**: `selfing_rate: 0.1` in an epoch

### Cloning Rate  
- **Demes**: Supports per-epoch cloning rates (asexual reproduction)
- **discoal**: No support for cloning
- **Example**: `cloning_rate: 0.5` in an epoch

## 2. Multiple Ancestors (Complex Admixture)

### Multi-way Admixture at Foundation
- **Demes**: A deme can have multiple ancestors with specified proportions at its founding time
- **discoal**: Only supports 2-way admixture events (-ea flag)
- **Example**: 
  ```yaml
  ancestors: [AFR, EUR, EAS]
  proportions: [0.167, 0.333, 0.5]
  ```

## 3. Pulse Events

### Instantaneous Gene Flow
- **Demes**: Supports pulse events - instantaneous migration of a proportion of one population into another
- **discoal**: No direct support for pulse events (would need to approximate with very high migration for short time)
- **Example**: 
  ```yaml
  pulses:
  - sources: [Neanderthal]
    dest: Modern_Human
    time: 1853.0
    proportions: [0.024]
  ```

### Multi-source Pulses
- **Demes**: Can have multiple source populations in a single pulse
- **discoal**: Not supported
- **Example**: Multiple populations contributing to admixture at the same time

## 4. Time-Bounded Migration

### Migration with Start and End Times
- **Demes**: Migrations can have both start_time and end_time
- **discoal**: Migration changes are instantaneous events (rate changes at specific times)
- **Example**:
  ```yaml
  migrations:
  - demes: [Pop1, Pop2]
    rate: 0.001
    start_time: 1000
    end_time: 500
  ```

## 5. Defaults System

### Global and Local Defaults
- **Demes**: Sophisticated defaults system for epochs and demes
- **discoal**: No defaults system
- **Impact**: More verbose configuration needed in discoal

## 6. Metadata Features

### DOI References
- **Demes**: Can include DOI references for citations
- **discoal**: No metadata support

### Descriptions
- **Demes**: Supports descriptions for demes and the overall model
- **discoal**: No description fields

## 7. Exponential Population Growth

### Continuous Exponential Growth
- **Demes**: Supports exponential growth within epochs via `size_function: exponential`
- **discoal**: NO SUPPORT for exponential growth - only discrete size changes via `-en` events
- **Impact**: Exponential growth must be approximated as a series of discrete size changes
- **Example**: 
  ```yaml
  epochs:
  - start_size: 1000
    end_size: 10000
    size_function: exponential  # NOT SUPPORTED - must approximate with steps
  ```

## 8. Time Units Flexibility

### Flexible Time Units
- **Demes**: Supports "generations" or "years" with generation_time conversion
- **discoal**: Works in coalescent time units (generations)
- **Status**: Partially supported in our converter

## 9. Continuous Migration Matrices

### Symmetric Migration Shortcuts
- **Demes**: Can specify symmetric migration between groups of demes easily
- **discoal**: Must specify each pairwise rate (though we added a helper)

## Testing Implications

When testing Demes models in discoal, we should:

1. **Skip or warn** about models using:
   - Selfing or cloning rates
   - Multi-way admixture (>2 ancestors)
   - Pulse events
   - Time-bounded migrations

2. **Approximate** where possible:
   - Convert pulses to brief high migration periods
   - Simplify complex admixture to sequential 2-way events

3. **Focus testing** on supported features:
   - Simple population splits
   - Two-way admixture
   - Population size changes
   - Constant migration rates
   - Exponential growth

4. **Create test suite** with:
   - Positive tests: Demes models that should work
   - Negative tests: Demes models that should be rejected
   - Warning tests: Models that work with limitations