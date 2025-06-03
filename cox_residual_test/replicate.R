# ===================
# Load Library
# ===================
library(survival)
library(simsurv)
library(ggplot2)
library(metap)
library(purrr)
library(reshape2)
library(dplyr)
library(tibble)


run_single_simulation_scenario2 <- function(n_total, n_sites, seed) {
# ===================
# Configuration
# ===================
set.seed(seed)

data_full <- data.frame(
  id = 1:n_total,
  x1 = rnorm(n_total),
  x2 = rnorm(n_total),
  x3 = rnorm(n_total)
)

# ===================
# Set Time-Varying Beta(t)
# ===================
# x1: proportional hazards (PH), constant beta
# x2: log-time violation
# x3: sinusoidal violation

betas <- data.frame(
  intercept        = rep(0.01, n_total),  # baseline hazard rate
  beta1            = rep(0.5,  n_total),  # x1 constant effect
  beta2_intercept  = rep(0.5,  n_total),  # x2 baseline effect
  beta2_tde        = rep(2,    n_total),  # x2 time-dependent effect
  beta3_intercept  = rep(0.5,  n_total),  # x3 baseline effect
  beta3_tde        = rep(2,    n_total)   # x3 time-dependent effect
)

# ------------ Define Hazard Function ------------ #
haz <- function(t, x, betas, ...) {
  beta0 <- betas[["intercept"]]
  beta1 <- betas[["beta1"]]
  beta2 <- betas[["beta2_intercept"]] + betas[["beta2_tde"]] * log(t + 1)
  beta3 <- betas[["beta3_intercept"]] + betas[["beta3_tde"]] * sin(pi * t / 5)
  
  lp <- beta0 + beta1 * x[["x1"]] + beta2 * x[["x2"]] + beta3 * x[["x3"]]
  basehaz <- 0.01
  
  return(basehaz * exp(lp))
}

# ===================
# Simulate Survival Data
# ===================
simdat <- simsurv(
  hazard = haz,
  x = data_full,
  betas = betas,
  maxt = 10
)

data_full <- merge(simdat, data_full, by = "id")

# ===================
# Cox PH Model and Schoenfeld Test (Full Data)
# ===================
fit_full <- coxph(Surv(eventtime, status) ~ x1 + x2 + x3, data = data_full)
zph_full <- cox.zph(fit_full, transform = "rank")
zph_full <- tibble::rownames_to_column(as.data.frame(zph_full$table), var = "term")

# Visualize Schoenfeld residuals
# plot(cox.zph(fit_full, transform = "rank"))

# ===================
# Split Data by Site and Run Site-specific PH Tests
# ===================
split_ids <- split(1:n_total, rep(1:n_sites, length.out = n_total))
site_data <- lapply(split_ids, function(ids) data_full[ids, ])

site_zph <- list()

for (i in 1:n_sites) {
  df <- site_data[[i]]
  fit <- coxph(Surv(eventtime, status) ~ x1 + x2 + x3, data = df)
  zph <- as.data.frame(cox.zph(fit, transform = "rank")$table)
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
    agreement = (p_meta < 0.05) == (p < 0.05) ,
    diff = p_meta - p,
    log_diff = log(p_meta) - log(p)
  )

return(zph_compare)
}



# ============
# Batch Simulation Grid
# ============
replicate_simulations <- function(n_total, n_sites, n_rep) {
  results <- vector("list", n_rep)
  
  for (i in 1:n_rep) {
    results[[i]] <- run_single_simulation_scenario2(n_total, n_sites, seed = i)
  }
  
  result_df <- bind_rows(results, .id = "rep") %>%
    mutate(rep = as.integer(rep))
  
  # Summary statistics
  summary_df <- result_df %>%
    group_by(term) %>%
    summarise(
      agreement_rate = mean(agreement),
      sd_agreement = sd(agreement),
      mean_diff = mean(diff),
      sd_diff = sd(diff),
      mean_log_diff = mean(log_diff),
      sd_log_diff = sd(log_diff),
      mean_p = mean(p),
      mean_p_meta = mean(p_meta),
      .groups = "drop"'
    ) %>%
    mutate(n_total = n_total, n_sites = n_sites)
  
  return(list(result_df, summary_df))
}

out <- replicate_simulations(n_total = 1000, n_sites = 2, n_rep = 50)
