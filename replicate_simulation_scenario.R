#' Replicate a Simulation Scenario Multiple Times
#'
#' Repeats a simulation scenario multiple times with different random seeds to compare
#' full-data and meta-analysis Schoenfeld residual tests under varying data realizations.
#'
#' @param scenario_id Integer. ID defining covariate structure and time-varying patterns.
#' @param replicate_seed Integer. Base seed to control reproducibility across repetitions.
#' @param n_total Integer. Total number of individuals to simulate per replicate.
#' @param n_sites Integer. Number of sites to split the dataset into.
#' @param n_rep Integer. Number of replicate simulations to run. Default is 100.
#' @param intercept_range Numeric vector of length 2. Range to sample intercept from.
#' @param beta_fixed_range Numeric vector of length 2. Range to sample fixed-effect coefficients.
#' @param beta_tde_range Numeric vector of length 2. Range to sample time-varying coefficients.
#' @param basehaz_range Numeric vector of length 2. Range for the baseline hazard multiplier.
#' @param tdefun_library Named list of time-varying functions. Use \code{get_tdefun_library()} to generate.
#'
#' @return A list with
#' \describe{
#'   \item{summary}{Combined data-frame of p-values, agreement flags, and metadata.}
#'   \item{replicate_debug}{List (length = \code{n_rep}) — each element contains
#'         \code{site_fit}, \code{site_zph}, \code{site_data}, \code{maxt}, \code{seed}.}
#' }
#'
#' @examples
#' replicate_simulation_scenario(
#'   scenario_id    = 2,
#'   replicate_seed = 12,
#'   n_total        = 100,
#'   n_sites        = 3,
#'   n_rep          = 10)
#'
#' @export



replicate_simulation_scenario <- function(scenario_id,
                                          replicate_seed,
                                          n_total,
                                          n_sites,
                                          n_rep = 100,
                                          intercept_range = c(0.4, 0.6),
                                          beta_fixed_range = c(0.4, 0.6),
                                          beta_tde_range = c(1.5, 2.5),
                                          basehaz_range = c(0.01, 0.05),
                                          tdefun_library = get_tdefun_library()) {
  # ===================
  # Loop over Repetitions
  # ===================
  summary_list       <- vector("list", n_rep)   # p-value summary
  replicate_debug    <- vector("list", n_rep)   # raw site-level caches
  
  for (i in 1:n_rep) {
    sim_result <- run_single_simulation_scenario(
      scenario_id = scenario_id,
      seed = i + replicate_seed,
      # control by a new seed
      n_total = n_total,
      n_sites = n_sites,
      intercept_range = intercept_range,
      beta_fixed_range = beta_fixed_range,
      beta_tde_range = beta_tde_range,
      basehaz_range = basehaz_range,
      tdefun_library = tdefun_library
    )
    
    # ---- store p-value comparison (one row per term) ----
    summary_list[[i]] <- sim_result$zph_compare %>%
      mutate(
        rep = i,
        scenario_id    = scenario_id,
        scenario_label = sim_result$scenario_label,
        replicate_seed = replicate_seed,
        scenario_seed  = sim_result$seed,
        maxt           = sim_result$maxt,
        n_total        = n_total,
        n_sites        = n_sites
      )
    
    # ---- keep raw site-level objects for debugging ----
    replicate_debug[[i]] <- list(
      site_fit  = sim_result$site_fit,
      site_zph  = sim_result$site_zph,
      site_data = sim_result$site_data
    )
    
  }
  
  # ===================
  # Combine Results
  # ===================
  return(list(
    summary = dplyr::bind_rows(summary_list),
    replicate_debug  = replicate_debug
  ))
  
}
