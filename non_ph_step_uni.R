# =========================================
# Step Function Cox PH Test (0.8 → -0.8)
# =========================================

# ======== Load Libraries =========
library(survival)
library(simsurv)
library(ggplot2)
library(metap)
library(purrr)

# ======== Set Parameters =========
set.seed(12)
n_total <- 10000
n_sites <- 3

data_full <- data.frame(id = 1:n_total, x1 = rnorm(n_total))

# ======== Set Beta(t) =========
# Step function: beta(t) = 0.8 → -0.8 after t = 5
betas_step <- c(x1 = 0.8)
tde_step   <- c(x1 = -1.6)
tdefun_step <- function(t) {
  step <- numeric(length(t))
  step[t > 5] <- 1.0
  return(step)
}

# ======== Visualize beta(t) =========
t_seq <- seq(0.1, 10, length.out = 300)
beta_t <- betas_step + tde_step * tdefun_step(t_seq)

plot(data.frame(t = t_seq, beta = beta_t), type = "l", col = "red", lwd = 2,
     main = "Step Function: beta(t) = 0.8 → -0.8", ylab = "beta(t)", xlab = "Time")

# ======== Simulate Time-to-Event Data =========
simdat <- simsurv(
  dist = "exponential",
  lambdas = 0.01,
  x = data_full,
  betas = betas_step,
  tde = tde_step,
  tdefunction = tdefun_step,
  maxt = 10
)

data_full <- merge(simdat, data_full, by = "id")

# ======== Full Data Fit and PH Test =========
fit_full <- coxph(Surv(eventtime, status) ~ x1, data = data_full)
zph_full <- cox.zph(fit_full, transform = "rank")
zph_full = tibble::rownames_to_column(as.data.frame(zph_full$table), var = "term")

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

library(dplyr)
site_zph_summary <- bind_rows(
  lapply(site_zph, function(df) tibble::rownames_to_column(df, var = "term")),
  .id = "site"
)
# dplyr::filter(site_zph_summary, term == "x1")

# ======== Meta-Analysis Aggregation =========
# # Fisher's method Example
# ## Fisher's method using base R
# p_vals <- site_zph_summary |> filter(term == "x1") |> pull(p)
# stat_fisher <- -2 * sum(log(p_vals))
# df_fisher <- 2 * length(p_vals)
# p_fisher <- pchisq(stat_fisher, df_fisher, lower.tail = FALSE)
# print(c(stat_fisher,df_fisher,p_fisher))
# 
# ## Fisher's method using meta package
# meta_fisher <- metap::sumlog(p_vals)
# stat_fisher <- meta_fisher$chisq
# df_fisher <- meta_fisher$df
# p_fisher <- meta_fisher$p
# print(c(stat_fisher,df_fisher,p_fisher))

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
  })%>%
  mutate(term = factor(term, levels = unique(site_zph_summary$term))) %>%
  arrange(term)

print(zph_meta)


# ======== Store All P-values =========
zph_full <- zph_full %>%
  mutate(method = "Full Data")

zph_meta <- zph_meta %>%
  mutate(method = "Meta-analysis")

zph_compare <- bind_rows(zph_full, zph_meta) %>%
  relocate(method, term, chisq, df, p)

print(zph_compare)




