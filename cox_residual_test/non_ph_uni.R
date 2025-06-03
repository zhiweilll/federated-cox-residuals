# federated_ph_test_singlevar_nonph_simsurv.R (使用 simsurv + exponential baseline 模拟违反 PH 的情景)

library(survival)
library(simsurv)
library(ggplot2)

set.seed(42)

# ========== 理论注释 ==========
# 模拟模型为 h(t|x) = h0(t) * exp(beta(t) * x)
# beta(t) = beta0 + alpha * log(t + 1)，明确违反 PH 假设
# 使用 simsurv 提供的 tde + tdefunction 参数来设定 beta(t)

# ========== 模拟设置 ==========
n_total <- 500
n_sites <- 2

data_full <- data.frame(id = 1:n_total, x1 = rnorm(n_total))

# 设置 β0 和 time-varying α 系数
betas <- c(x1 = 0.5)
tde <- c(x1 = 0.5)
tdefun <- function(t) log(t + 1)

# 使用 simsurv 模拟具有时间变系数的生存数据
simdat <- simsurv(
  dist = "exponential",
  lambdas = 0.01,
  x = data_full,
  betas = betas,
  tde = tde,
  tdefunction = tdefun,
  maxt = 10
)

# 合并协变量
data_full <- merge(simdat, data_full, by = "id")
head(data_full)
# # 添加删失
# censor_flag <- rbinom(n_total, 1, 0.9)
# data_full$time <- ifelse(censor_flag == 1,
#                          data_full$eventtime,
#                          runif(n_total, 0, data_full$eventtime))
# data_full$status <- censor_flag

# ========== 拟合 Cox 模型 + PH 检验 ==========
fit_full <- coxph(Surv(eventtime, status) ~ x1, data = data_full)
zph_full <- cox.zph(fit_full, transform = "rank")  # 更稳健地检验时间趋势

cat("\n[Gold Standard - Full Data]\n")
print(zph_full)

# ========== 拆分站点 ==========
split_ids <- split(1:n_total, rep(1:n_sites, length.out = n_total))
site_data <- lapply(split_ids, function(ids) data_full[ids, ])

site_zph <- lapply(site_data, function(df) {
  fit <- coxph(Surv(eventtime, status) ~ x1, data = df)
  cox.zph(fit, transform = "rank")
})

# ========== Meta-analysis ==========
chisq_total <- sum(sapply(site_zph, function(z) z$table["x1", "chisq"]))
df_total <- sum(sapply(site_zph, function(z) z$table["x1", "df"]))
p_meta_chisq <- 1 - pchisq(chisq_total, df_total)

p_vals <- sapply(site_zph, function(z) z$table["x1", "p"])
fisher_stat <- -2 * sum(log(p_vals))
df_fisher <- 2 * length(p_vals)
p_meta_fisher <- 1 - pchisq(fisher_stat, df_fisher)

# ========== 输出结果 ==========
p_gold <- zph_full$table["x1", "p"]

cat("\n[Comparison]\n")
cat("Full data p-value:", p_gold, "\n")
cat("Meta-analysis p-value (Chi-square sum):", p_meta_chisq, "\n")
cat("Meta-analysis p-value (Fisher method):", p_meta_fisher, "\n")

# ========== 可视化 beta(t), S(t|x), F(t|x) ==========
t_seq <- seq(0.1, 10, length.out = 300)
beta0 <- 0.5
alpha <- 0.5
x_vals <- c(-2, 0, 2)

plot_data <- do.call(rbind, lapply(x_vals, function(x) {
  beta_t <- beta0 + alpha * log(t_seq + 1)
  lambda_t <- 0.01 * exp(beta_t * x)
  S <- exp(-lambda_t * t_seq)
  F <- 1 - S
  data.frame(t = t_seq, beta = beta_t, S = S, F = F, x = paste0("x = ", x))
}))

p_beta <- ggplot(plot_data, aes(x = t, y = beta, color = x)) +
  geom_line(size = 1.1) +
  labs(title = "Time-varying Coefficient beta(t)", y = "beta(t)", x = "Time") +
  theme_minimal()

p_surv <- ggplot(plot_data, aes(x = t, y = S, color = x)) +
  geom_line(size = 1.1) +
  labs(title = "Survival Function S(t|x)", y = "S(t|x)", x = "Time") +
  theme_minimal()

p_cdf <- ggplot(plot_data, aes(x = t, y = F, color = x)) +
  geom_line(size = 1.1) +
  labs(title = "CDF F(t|x)", y = "F(t|x)", x = "Time") +
  theme_minimal()

print(p_beta)
print(p_surv)
print(p_cdf)

# ========== 注释：为何 transform = "rank"？ ==========
# Schoenfeld 残差回归需要对时间建模，默认使用 log(time)
# 若时间分布极度偏态或有重复值，log(time) 可能不稳
# transform = "rank" 使用事件时间的秩作为时间变量 → 更鲁棒
# 特别适用于模拟或非 parametric 情况下进行 PH 检验
