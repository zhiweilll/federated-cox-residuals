# Federated Learning Residual Test Agreement Analysis

This project simulates and analyzes the agreement between site-specific Schoenfeld residual test results and full data results in a federated learning context, specifically focusing on Cox proportional hazards models.

## Project Overview

The project aims to: 1. Simulate survival data with various time-varying covariate effects 2. Split data into multiple sites to simulate a federated learning environment 3. Perform Schoenfeld residual tests on individual site data (site-specific analysis) 4. Compare these results with residual tests performed on the full combined dataset 5. Analyze the agreement between site-specific and full data results using meta-analysis

## Methodology

### Data Generation and Scenarios

-   Supports multiple simulation scenarios (1-35+) with different covariate configurations:
    -   Univariate scenarios (x1 only)
    -   Bivariate scenarios (x1, x2)
    -   Trivariate scenarios (x1, x2, x3)
-   Each scenario can have different time-varying patterns:
    -   Proportional hazards (PH)
    -   Step function
    -   Spline function
    -   Logarithmic function
    -   Linear function

### Simulation Process

1.  Data Generation:
    -   Generate covariates with specified time-varying patterns
    -   Simulate survival times using Cox proportional hazards model
    -   Split data into multiple sites
2.  Analysis:
    -   Site-specific Analysis:
        -   Fit Cox models on each site's data
        -   Perform Schoenfeld residual tests
        -   Collect p-values for each covariate
    -   Full Data Analysis:
        -   Combine all site data
        -   Fit Cox model on complete dataset
        -   Perform Schoenfeld residual test
3.  Meta-analysis:
    -   Aggregate site-specific p-values using Fisher's method
    -   Compare meta-analysis results with full data results
    -   Calculate agreement metrics

## Usage

### Basic Usage

``` r
# Run a single simulation
result <- run_single_simulation_scenario(
  scenario_id = 1,
  seed = 123,
  n_total = 1000,
  n_sites = 5
)

# Run multiple replications
replications <- replicate_simulation_scenario(
  scenario_id = 1,
  replicate_seed = 123,
  n_total = 1000,
  n_sites = 5,
  n_rep = 100
)
```

### Key Parameters

-   `scenario_id`: Defines the covariate structure and time-varying patterns
-   `n_total`: Total number of individuals to simulate
-   `n_sites`: Number of sites to split the data into
-   `n_rep`: Number of simulation replications (for stability analysis)

## Requirements

-   R packages:
    -   survival
    -   simsurv
    -   dplyr
    -   metap
    -   tibble

## Results

The analysis provides: 1. Comparison of p-values between full data and meta-analysis 2. Agreement metrics between site-specific and full data results 3. Statistical significance of differences 4. Detailed debugging information for each site

## Scenario Structure

The project supports various scenarios (1-35+) with different combinations of: - Number of covariates (1-3) - Time-varying patterns (PH, step, spline, log, linear) - Interaction patterns between covariates

See `generate_scenario_structure.R` for detailed scenario mapping.
