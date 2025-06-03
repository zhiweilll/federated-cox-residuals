# =========================================
# Spline Function Cox PH Test (beta(t) ~ sin-like curve)
# =========================================

# ======== Load Libraries =========
library(survival)
library(simsurv)
library(ggplot2)
library(metap)
library(purrr)
library(tibble)
library(dplyr)
library(tibble)
# ======== Set Parameters =========
set.seed(12)
n_total <- 10000
n_sites <- 3

data_full <- data.frame(id = 1:n_total, x1 = rnorm(n_total))

# ======== Set Beta(t) =========
# Spline-like time-varying beta(t): beta(t) = 0.5 + 0.5 * sin(pi * t / 5)
betas_spline <- c(x1 = 0.0)
tde_spline   <- c(x1 = 1.0)
tdefun_spline <- function(t) {
  sin(pi * t / 5)
}

# ======== Visualize beta(t) =========
t_seq <- seq(0.1, 10, length.out = 300)
beta_t <- betas_spline + tde_spline * tdefun_spline(t_seq)

plot(data.frame(t = t_seq, beta = beta_t), type = "l", col = "blue", lwd = 2,
     main = "Spline Function: beta(t) = sin(pi * t / 5)", ylab = "beta(t)", xlab = "Time")

# ======== Simulate Time-to-Event Data =========
simdat <- simsurv(
  dist = "exponential",
  lambdas = 0.01,
  x = data_full,
  betas = betas_spline,
  tde = tde_spline,
  tdefunction = tdefun_spline,
  maxt = 10
)

data_full <- merge(simdat, data_full, by = "id")

# ======== Full Data Fit and PH Test =========
fit_full <- coxph(Surv(eventtime, status) ~ x1, data = data_full)
zph_full <- cox.zph(fit_full, transform = "rank")
plot(zph_full) 

zph_full <- rownames_to_column(as.data.frame(zph_full$table), var = "term") %>%
  mutate(method = "Full Data")


# ======== Split into Sites and Test PH Per Site =========
split_ids <- split(1:n_total, rep(1:n_sites, length.out = n_total))
site_data <- lapply(split_ids, function(ids) data_full[ids, ])

site_zph <- list()

for (i in 1:n_sites) {
  df <- site_data[[i]]
  fit <- coxph(Surv(eventtime, status) ~ x1, data = df)
  zph <- as.data.frame(cox.zph(fit, transform = "rank")$table)
  site_zph[[i]] <- zph
}

site_zph_summary <- bind_rows(
  lapply(site_zph, function(df) rownames_to_column(df, var = "term")),
  .id = "site"
)

# ======== Meta-Analysis Aggregation =========
zph_meta <- site_zph_summary %>%
  group_split(term) %>%
  map_dfr(~ {
    res <- metap::sumlog(.x$p)
    tibble(
      term = unique(.x$term),
      chisq = res$chisq,
      df = res$df,
      p = res$p
    )
  }) %>%
  mutate(
    term = factor(term, levels = unique(site_zph_summary$term)),
    method = "Meta-analysis"
  ) %>%
  arrange(term)

# ======== Store All P-values =========
zph_compare <- bind_rows(zph_full, zph_meta) %>%
  relocate(method, term, chisq, df, p)

print(zph_compare)

