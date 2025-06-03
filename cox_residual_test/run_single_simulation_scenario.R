#' Run a Single Simulation Scenario
#'
#' This function runs a single simulation under a defined scenario using the Cox
#' proportional hazards model. It compares Schoenfeld residual tests from full data
#' and meta-analyzed site-specific results using Fisher’s method.
#'
#' @param scenario_id Integer. ID of the scenario structure, determining covariates and their time-dependency.
#' @param seed Integer. Random seed for generating covariates, betas, hazard, and time range.
#' @param n_total Integer. Total number of individuals to simulate.
#' @param n_sites Integer. Number of sites to split the data into for meta-analysis.
#' @param intercept_range Numeric vector. Range of intercept values.
#' @param beta_fixed_range Numeric vector. Range of time-fixed beta coefficients.
#' @param beta_tde_range Numeric vector. Range of time-varying beta slopes.
#' @param basehaz_range Numeric vector. Range of baseline hazard multipliers.
#' @param tdefun_library Named list. Mapping of time-varying forms to expression strings.
#'
#' @return A list containing:
#' \describe{
#'   \item{zph_compare}{Data frame comparing full-data and meta p-values and agreement.}
#'   \item{scenario_label}{String label of the scenario.}
#'   \item{haz_str}{Hazard function as a character string.}
#'   \item{haz_expr}{Evaluated hazard function (for use in `simsurv`).}
#'   \item{maxt}{Maximum follow-up time used in the simulation.}
#'   \item{seed}{Seed used.}
#'   \item{n_total}{Total number of samples.}
#'   \item{n_sites}{Number of sites simulated.}
#' }
#'
#' @export
#'
run_single_simulation_scenario <- function(scenario_id, seed, n_total, n_sites,
                                           intercept_range, beta_fixed_range, beta_tde_range, basehaz_range,
                                           tdefun_library) {
  set.seed(seed)
  
  # ===================
  # Generate Full Data
  # ===================
  # Randomly sample maxt within a reasonable range, e.g. between 8 and 15
  maxt <- as.integer(runif(1, min = 8, max = 15))

  # Generate covariate data
  data_full <- data.frame(id = 1:n_total)
  vars <- generate_scenario_structure(scenario_id)$vars
  for (v in vars) data_full[[paste0("x", v)]] <- rnorm(n_total)

  # ===================
  # Generate Scenario
  # ===================
  scenario_info <- generate_scenario_structure(scenario_id) # no randomness

  # ===================
  # Generate Betas and Hazard Function
  # ===================
  beta_haz <- generate_beta_and_haz(
    seed = seed,
    n_total = n_total,
    vars = scenario_info$vars,
    forms = scenario_info$forms,
    intercept_range = intercept_range,
    beta_fixed_range = beta_fixed_range,
    beta_tde_range = beta_tde_range,
    basehaz_range = basehaz_range,
    tdefun_library = tdefun_library
  )
  betas    <- beta_haz$betas
  haz_str  <- beta_haz$haz_str
  haz_expr <- beta_haz$haz_expr

  # ===================
  # Simulate Survival Data
  # ===================
  simdat <- simsurv(
    hazard = haz_expr,
    x = data_full,
    betas = betas,
    maxt = maxt
  )

  data_full <- merge(simdat, data_full, by = "id")

  # ===================
  # Cox PH Model and Schoenfeld Test (Full Data)
  # ===================
  formula <- as.formula(Surv(eventtime, status) ~ . - id)
  fit_full <- coxph(formula, data = data_full)
  zph_full <- cox.zph(fit_full, transform = "rank")
  zph_full <- tibble::rownames_to_column(as.data.frame(zph_full$table), var = "term")

  # ===================
  # Site-specific Cox + Meta Schoenfeld Tests
  # ===================
  split_ids <- split(1:n_total, rep(1:n_sites, length.out = n_total))
  site_data <- lapply(split_ids, function(ids) data_full[ids, ])
  
  site_zph <- list()
  site_fit <- list()
  
  # for (i in 1:n_sites) {
  #   df <- site_data[[i]]  
  #   
  #   fit <- tryCatch(
  #     coxph(formula, data = df),
  #     error = function(e) return(NULL)
  #   )
  #   site_fit[[i]] <- fit
  #   
  #   if (!is.null(fit)) {
  #     zph <- tryCatch(
  #       as.data.frame(cox.zph(fit, transform = "rank")$table),
  #       error = function(e) return(NULL)
  #     )
  #     site_zph[[i]] <- zph
  #   } else {
  #     site_zph[[i]] <- NULL
  #   }
  # }
  
  for (i in 1:n_sites) {
    df <- site_data[[i]]
    
    fit <- tryCatch(
      coxph(formula, data = df),
      error = function(e) {
        message(sprintf("[Site %d] coxph() failed: %s", i, e$message))
        return(NULL)
      }
    )
    site_fit[[i]] <- fit
    
    zph <- tryCatch(
      if (!is.null(fit)) as.data.frame(cox.zph(fit, transform = "rank")$table) else NULL,
      error = function(e) {
        message(sprintf("[Site %d] cox.zph() failed: %s", i, e$message))
        return(NULL)
      }
    )
    site_zph[[i]] <- zph
  }
  
  site_zph_summary <- bind_rows(
    lapply(site_zph, function(df) tibble::rownames_to_column(df, var = "term")),
    .id = "site"
  )
  
  # ===================
  # Meta-analysis Aggregation of Schoenfeld P-values
  # ===================
  zph_meta <- site_zph_summary %>%
    group_split(term) %>%
    map_dfr(~ {
      if (all(is.na(.x$p))) return(NULL)
      res <- metap::sumlog(.x$p)
      tibble(
        term = unique(.x$term),
        chisq_meta = res$chisq,
        df_meta = res$df,
        p_meta = res$p
      )
    }) %>%
    mutate(term = factor(term, levels = unique(site_zph_summary$term))) %>%
    arrange(term)


  # ===================
  # Combine Full and Meta-analysis Results
  # ===================
  zph_compare <- left_join(zph_full, zph_meta, by = "term") %>%
    mutate(
      agreement = (p_meta < 0.05) == (p < 0.05),
      diff = p_meta - p,
      log_diff = log(p_meta + 1e-10) - log(p + 1e-10)
    )

  return(list(
    zph_compare = zph_compare,
    scenario_label = scenario_info$scenario_label,
    haz_str = haz_str,
    haz_expr = haz_expr,
    maxt = maxt,
    seed = seed,
    n_total = n_total,
    n_sites = n_sites,
    site_fit = site_fit,
    site_zph = site_zph,
    site_data = site_data
  ))
}
