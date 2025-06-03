#' Run one simulation – **sites generated first**
#'
#' This helper creates an *independent* dataset for each site, fits a Cox model
#' per-site, performs the Schoenfeld residual test, aggregates the p-values
#' with Fisher’s method, **then** binds all sites into a single “gold-standard”
#' dataset and runs the full–data test.  The function returns the two sets of
#' p-values plus complete site-level caches for debugging.
#'
#' @param scenario_id   Integer. ID produced by `generate_scenario_structure()`.
#' @param seed          Integer. Master seed controlling scenario-wide β/λ and
#'                      the base for per-site seeds.
#' @param n_sites       Integer. Number of federated sites.
#' @param n_per_site    Integer. Sample size *per* site (total = `n_sites*n_per_site`).
#' @param intercept_range,beta_fixed_range,beta_tde_range,basehaz_range
#'                      Length-2 numeric vectors giving the sampling ranges for
#'                      intercept, fixed effects, time-varying slopes, and
#'                      baseline-hazard multiplier.
#' @param tdefun_library Named list created by `get_tdefun_library()` that maps
#'                      `"linear"`, `"log"`, `"step"`, `"spline"` → character
#'                      expressions of *t*.
#'
#' @return A list with  
#'   \describe{
#'     \item{\code{compare}}{data-frame: full-data vs meta p-values and agreement.}
#'     \item{\code{site_data}}{list of site‐level data .frames.}
#'     \item{\code{site_fit}}{list of \code{coxph} objects (or errors).}
#'     \item{\code{site_zph}}{list of Schoenfeld residual tables (or \code{NULL}).}
#'     \item{\code{site_seed_vec}}{integer vector, seed actually used for each site.}
#'     \item{\code{data_full}}{row-bound full dataset used for gold-standard fit.}
#'     \item{\code{scenario_label}}{human-readable description of the scenario.}
#'   }
#' @export
run_single_simulation_scenario_2 <- function(
    scenario_id, seed, n_sites, n_per_site,
    intercept_range, beta_fixed_range, beta_tde_range,
    basehaz_range, tdefun_library
) {
  set.seed(seed)
  
  # ===================
  # Scenario-wide β / hazard
  # ===================
  scen <- generate_scenario_structure(scenario_id)
  beta_haz <- generate_beta_and_haz(
    seed             = seed,
    n_total          = n_per_site,
    vars             = scen$vars,
    forms            = scen$forms,
    intercept_range  = intercept_range,
    beta_fixed_range = beta_fixed_range,
    beta_tde_range   = beta_tde_range,
    basehaz_range    = basehaz_range,
    tdefun_library   = tdefun_library
  )
  
  # ===================
  # Containers
  # ===================
  site_data      <- vector("list", n_sites)
  site_fit       <- vector("list", n_sites)
  site_zph       <- vector("list", n_sites)
  site_seed_vec  <- integer(n_sites)
  
  site_base_seed <- seed * 10  # deterministic offset for per-site seeds
  
  # ============================
  # Generate / fit each site
  # ============================
  for (s in seq_len(n_sites)) {
    
    ## ---- site-specific RNG seed ----
    site_seed <- site_base_seed + s
    site_seed_vec[s] <- site_seed
    set.seed(site_seed)
    
    ## ---- (a) covariates ----
    df_cov <- data.frame(id = 1:n_per_site)
    for (v in scen$vars)
      df_cov[[paste0("x", v)]] <- rnorm(n_per_site)
    
    ## ---- (b) survival times via simsurv ----
    simdat <- simsurv(
      hazard = beta_haz$haz_expr,
      x      = df_cov,
      betas  = beta_haz$betas,
      maxt   = 10
    )
    site_df <- dplyr::left_join(simdat, df_cov, by = "id") |>
      dplyr::mutate(site = s, site_seed = site_seed)
    
    ## ---- (c) Cox fit & Schoenfeld test ----
    fml <- Surv(eventtime, status) ~ . - id - site - site_seed
    fit <- try(coxph(fml, data = site_df), silent = TRUE)
    
    zph <- if (inherits(fit, "coxph"))
      try(as.data.frame(cox.zph(fit, "rank")$table),
          silent = TRUE) else NULL
    
    site_data[[s]] <- site_df
    site_fit[[s]]  <- fit
    site_zph[[s]]  <- if (inherits(zph, "data.frame")) zph else NULL
  }
  
  # ===================
  # Fisher meta-analysis
  # ===================
  site_zph_summary <- dplyr::bind_rows(
    lapply(site_zph, function(z)
      if (!is.null(z))
        tibble::rownames_to_column(z, "term") else NULL),
    .id = "site"
  )
  
  zph_meta <- site_zph_summary |>
    dplyr::group_split(term) |>
    purrr::map_dfr(~{
      res <- metap::sumlog(.x$p)
      tibble(term = unique(.x$term),
             chisq_meta = res$chisq,
             df_meta    = res$df,
             p_meta     = res$p)
    })
  
  # ===================
  # Full-data gold-standard
  # ===================
  data_full <- dplyr::bind_rows(site_data)
  fit_full  <- coxph(Surv(eventtime, status) ~ . - id - site - site_seed,
                     data = data_full)
  zph_full  <- tibble::rownames_to_column(
    as.data.frame(cox.zph(fit_full, "rank")$table), "term")
  
  # ===================
  # Combine Results
  # ===================
  compare <- dplyr::left_join(zph_full, zph_meta, by = "term") |>
    dplyr::mutate(agreement = (p_meta < 0.05) == (p < 0.05))
  
  # ---- keep raw site-level objects for debugging ----
  list(
    compare        = compare,
    site_data      = site_data,
    site_fit       = site_fit,
    site_zph       = site_zph,
    site_seed_vec  = site_seed_vec,
    data_full      = data_full,
    scenario_label = scen$scenario_label
  )
}