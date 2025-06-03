debug_simulation_output <- function(sim_result) {
  message("=== Simulation Debug Summary ===")
  
  site_fit <- sim_result$site_fit
  site_zph <- sim_result$site_zph
  safe_fit_check <- function(fit) {
    tryCatch({
      if (!inherits(fit, "coxph")) return(TRUE)
      if (is.null(fit$converged) || is.na(fit$converged)) return(TRUE)
      return(!fit$converged)
    }, error = function(e) TRUE)
  }
  
  # Check Cox model convergence
  failed_converge_sites <- which(sapply(site_fit, safe_fit_check))
  if (length(failed_converge_sites) > 0) {
    message("❌ Sites where coxph() did NOT converge or errored: ", paste(failed_converge_sites, collapse = ", "))
  } else {
    message("✅ All site-level coxph() models converged.")
  }
  
  # Check failed zph
  failed_zph_sites <- which(sapply(site_zph, function(x) is.null(x) || inherits(x, "try-error")))
  if (length(failed_zph_sites) > 0) {
    message("❌ Sites where cox.zph() FAILED: ", paste(failed_zph_sites, collapse = ", "))
  } else {
    message("✅ All site-level cox.zph() succeeded.")
  }
  
  # Check zph_compare NA
  if ("zph_compare" %in% names(sim_result)) {
    zph_compare <- sim_result$zph_compare
    na_terms <- zph_compare %>%
      filter(is.na(p) | is.na(p_meta)) %>%
      pull(term) %>%
      unique()
    
    if (length(na_terms) > 0) {
      message("⚠️ Terms with NA in zph_compare: ", paste(na_terms, collapse = ", "))
    } else {
      message("✅ No NA in zph_compare p-values.")
    }
  }
  
  # Check NA in site-level p
  site_zph_summary <- tryCatch({
    bind_rows(
      lapply(site_zph, function(df) {
        if (is.null(df)) return(NULL)
        df <- as.data.frame(df)
        tibble::rownames_to_column(df, var = "term")
      }),
      .id = "site"
    )
  }, error = function(e) NULL)
  
  if (!is.null(site_zph_summary)) {
    na_term_site <- site_zph_summary %>%
      filter(is.na(p)) %>%
      select(site, term)
    
    if (nrow(na_term_site) > 0) {
      message("⚠️ Specific site+term with NA p-values:")
      print(na_term_site)
    } else {
      message("✅ No NA in site-level p-values.")
    }
  } else {
    message("⚠️ Could not summarize site-level zph results.")
  }
}