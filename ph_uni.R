# federated_ph_test_singlevar_nonph.R (模拟一个变量违反 PH 假设的情形)

library(survival)
library(ggplot2)

set.seed(42)

# ========== 理论注释：如何构造一个变量违反 PH 假设的情形 ==========
# 我们仍使用 Cox 模型 h(t|x) = h0(t) * exp(beta(t) * x1)
# 与 PH 假设不同的是：beta 不再是常数，而是随时间变化
# 比如：beta(t) = beta0 + alpha * log(t + 1)，使得 log hazard 比随时间增长变化
# 这就违反了比例风险假设，因为风险比不再是恒定倍数关系

# ========== 模拟数据 ==========
n_total <- 500
n_sites <- 2

# 协变量 x1
data_full <- data.frame(
  id = 1:n_total,
  x1 = rnorm(n_total)
)

# 设置非PH的 beta(t): beta(t) = beta0 + alpha * log(t + 1)
beta0 <- 0.5
alpha <- 0.5  # 控制时间依赖程度
lambda0 <- 0.01

# 自定义生存时间模拟函数（违反 PH）
sim_nonph <- function(x1, lambda0, beta0, alpha) {
  # 反解 F(t) = u → t = -log(1 - u) / (lambda0 * exp(beta(t) * x1))
  # 但 beta(t) 非常复杂，不能解析 → 我们用 rejection sampling 近似模拟
  n <- length(x1)
  t <- numeric(n)
  for (i in 1:n) {
    repeat {
      t_try <- rexp(1, rate = lambda0)  # 提议分布
      beta_t <- beta0 + alpha * log(t_try + 1)
      hazard_ratio <- exp(beta_t * x1[i])
      accept_prob <- hazard_ratio / (1 + hazard_ratio)  # 缩放到 (0,1) 区间
      if (runif(1) < accept_prob) {
        t[i] <- t_try
        break
      }
    }
  }
  return(t)
}

data_full$eventtime <- sim_nonph(data_full$x1, lambda0, beta0, alpha)

# 添加删失
censor_flag <- rbinom(n_total, 1, 0.9)
data_full$time <- ifelse(censor_flag == 1,
                         data_full$eventtime,
                         runif(n_total, 0, data_full$eventtime))
data_full$status <- censor_flag

# ========== 全数据 Cox + PH 检验 ==========
fit_full <- coxph(Surv(time, status) ~ x1, data = data_full)
zph_full <- cox.zph(fit_full, transform = "rank")

cat("\n[Gold Standard - Full Data]\n")
print(zph_full)

# ========== 拆分站点数据 ==========
split_ids <- split(1:n_total, rep(1:n_sites, length.out = n_total))
site_data <- lapply(split_ids, function(ids) data_full[ids, ])

# ========== 每站点拟合 Cox + PH 检验 ==========
site_zph <- lapply(site_data, function(df) {
  fit <- coxph(Surv(time, status) ~ x1, data = df)
  cox.zph(fit, transform = "rank")
})

# ========== 方法 1：Chi-square 加和 ==========
chisq_total <- sum(sapply(site_zph, function(z) z$table["x1", "chisq"]))
df_total <- sum(sapply(site_zph, function(z) z$table["x1", "df"]))
p_meta_chisq <- 1 - pchisq(chisq_total, df_total)

# ========== 方法 2：Fisher ==========
p_vals <- sapply(site_zph, function(z) z$table["x1", "p"])
fisher_stat <- -2 * sum(log(p_vals))
df_fisher <- 2 * length(p_vals)
p_meta_fisher <- 1 - pchisq(fisher_stat, df_fisher)

# ========== 输出对比结果 ==========
p_gold <- zph_full$table["x1", "p"]

cat("\n[Comparison]\n")
cat("Full data p-value:", p_gold, "\n")
cat("Meta-analysis p-value (Chi-square sum):", p_meta_chisq, "\n")
cat("Meta-analysis p-value (Fisher method):", p_meta_fisher, "\n")

# ========== 可视化 S(t|x) 与非PH系数曲线（近似） ==========
t_seq <- seq(0.1, 100, length.out = 300)
x_vals <- c(-2, 0, 2)

plot_beta <- do.call(rbind, lapply(x_vals, function(x) {
  beta_t <- beta0 + alpha * log(t_seq + 1)
  lambda_x <- lambda0 * exp(beta_t * x)
  S <- exp(-lambda_x * t_seq)
  data.frame(t = t_seq, beta = beta_t, S = S, x = paste0("x = ", x))
}))

p_beta <- ggplot(plot_beta, aes(x = t, y = beta, color = x)) +
  geom_line(size = 1.1) +
  labs(title = "Time-varying beta(t)", y = "beta(t)", x = "Time") +
  theme_minimal()

p_surv <- ggplot(plot_beta, aes(x = t, y = S, color = x)) +
  geom_line(size = 1.1) +
  labs(title = "Survival Function S(t|x) under non-PH", y = "S(t|x)", x = "Time") +
  theme_minimal()

print(p_beta)
print(p_surv)
