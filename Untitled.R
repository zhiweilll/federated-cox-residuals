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

# ===================
# Configuration
# ===================
set.seed(12)
n_total <- 1000
n_sites <- 2

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
# Visualize Beta(t)
# ===================

# ------------ Extract parameters from subject 1 ------------ #
beta1           <- betas[1, "beta1"]
beta2_intercept <- betas[1, "beta2_intercept"]
beta2_tde       <- betas[1, "beta2_tde"]
beta3_intercept <- betas[1, "beta3_intercept"]
beta3_tde       <- betas[1, "beta3_tde"]

# ------------ Generate time-varying beta(t) ------------ #
t_seq <- seq(1, 10, length.out = 300)

beta_t_df <- data.frame(
  t = t_seq,
  beta1 = rep(beta1, length(t_seq)),
  beta2 = beta2_intercept + beta2_tde * log(t_seq + 1),
  beta3 = beta3_intercept + beta3_tde * sin(pi * t_seq / 5)
)

beta_long <- melt(beta_t_df, id.vars = "t", variable.name = "Covariate", value.name = "beta")

# ------------ Plot beta(t) ------------ #
ggplot(beta_long, aes(x = t, y = beta, color = Covariate)) +
  geom_line(size = 1.2) +
  labs(title = "Time-varying Coefficients beta(t) for x1, x2, x3",
       x = "Time", y = expression(beta(t))) +
  theme_minimal() +
  scale_color_manual(values = c("darkgreen", "blue", "red"))

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
plot(cox.zph(fit_full, transform = "rank"))

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

print(zph_meta)

# ===================
# Combine Full and Meta-analysis Results
# ===================
zph_compare <- left_join(zph_full, zph_meta, by = "term") %>%
  mutate(
    p_agree = (p < 0.05) == (p_meta < 0.05),
    p_log_diff = log(p + 1e-10) - log(p_meta + 1e-10)
  )

print(zph_compare)

