# ============
# Cox PH Simulation Study Framework with Scenario Generator
# ============

library(survival)
library(simsurv)
library(dplyr)
library(tibble)
library(purrr)
library(metap)
library(pROC)

# =============================
# Time-varying Effect Library
# =============================
tdefun_library <- list(
  linear = "t",
  log = "log(t + 1)",
  step = "as.numeric(t > 5)",
  spline = "sin(pi * t / 5)"
)

# get_tdefun_library <- function(seed = 1) {
#   set.seed(seed)
# 
#   # ------- Linear: a * t + b -------
#   a <- sample(seq(0.5, 2, 0.1), 1)   # slope
#   b <- sample(seq(-1, 1, 0.2), 1)    # intercept
#   linear_expr <- sprintf("%.2f * t + %.2f", a, b)
# 
#   # ------- Log: log(t + shift) -------
#   shift <- sample(seq(0.5, 2.5, 0.1), 1)
#   log_expr <- sprintf("log(t + %.2f)", shift)
# 
#   # ------- Step: h * 1{t > t1} -------
#   h <- sample(seq(0.5, 2.0, 0.1), 1)
#   t <- sample(2:4, 1)
#   step_expr <- sprintf("%.2f * as.numeric(t > %d)", h, t)
# 
#   # ------- Spline: sin(pi * t / freq) -------
#   freq <- sample(seq(3, 7, 0.5), 1)
#   spline_expr <- sprintf("sin(pi * t / %.1f)", freq)
# 
#   return(list(
#     linear = linear_expr,
#     log = log_expr,
#     step = step_expr,
#     spline = spline_expr
#   ))
# }




# ============
# Generate Covariate DataFrame
# ============
get_covariate_x <- function(vars, n_total) {
  df <- data.frame(id = 1:n_total)
  for (v in vars) {
    df[[paste0("x", v)]] <- rnorm(n_total)
  }
  return(df)
}
x <- get_covariate_x(vars, n_total)






# ============
# Scenario Generator Function
# ============
generate_scenario <- function(scenario_id, 
                              n_total,
                              intercept_range = c(0, 1),
                              beta_fixed_range = c(0, 1),
                              beta_tde_range = c(0, 3),
                              basehaz_range = c(0, 0.5)) {

  tdefun_library <- list(
    linear = "t",
    log = "log(t + 1)",
    step = "as.numeric(t > 5)",
    spline = "sin(pi * t / 5)"
  )
  
  # Generate per-variable fixed and time-varying beta values
  get_beta_names <- function(vars, forms, n_total, scenario_id,
                             intercept_range,
                             beta_fixed_range,
                             beta_tde_range) {
    set.seed(scenario_id)
    
    intercept <- runif(1, intercept_range[1], intercept_range[2])
    beta_fixed <- runif(length(vars), beta_fixed_range[1], beta_fixed_range[2])
    beta_tde <- runif(length(vars), beta_tde_range[1], beta_tde_range[2])
    
    betas <- data.frame(intercept = rep(intercept, n_total))
    for (i in seq_along(vars)) {
      v <- vars[i]
      betas[[paste0("beta", v)]] <- rep(beta_fixed[i], n_total)
      if (!is.na(forms[i]) && forms[i] != "ph") {
        betas[[paste0("beta", v, "_tde")]] <- rep(beta_tde[i], n_total)
      }
    }
    return(betas)
  }
  
  # Generate basehazard from seed and range
  get_basehaz <- function(scenario_id, basehaz_range) {
    set.seed(scenario_id)
    runif(1, basehaz_range[1], basehaz_range[2])
  }
  
  # Build hazard expression string
  build_haz_expr <- function(vars, forms, basehaz, tdefun_library) {
    header <- c("beta0 <- betas[['intercept']]")
    terms <- c("beta0")
    for (i in seq_along(vars)) {
      v <- vars[i]
      if (forms[i] == "ph") {
        header <- c(header, sprintf("beta%s <- betas[['beta%s']]", v, v))
      } else {
        header <- c(header, sprintf(
          "beta%s <- betas[['beta%s']] + betas[['beta%s_tde']] * %s",
          v, v, v, tdefun_library[[forms[i]]]
        ))
      }
      terms <- c(terms, sprintf("beta%s * x[['x%s']]", v, v))
    }
    
    expr <- paste0("lp <- ", paste(terms, collapse = " + "))
    body <- c(header, expr, sprintf("return(%f * exp(lp))", basehaz))
    full <- paste(c("function(t, x, betas, ...) {", paste0("  ", body), "}"), collapse = "\n")
    return(full)
  }
  
  # =====================
  # Scenario Mapping
  # =====================
  # 1       : x1 PH (univariate)
  # 2-5     : x1 time-varying (step/spline/log/linear)
  # 6-9     : x1 PH, x2 time-varying
  # 10-19   : x1 and x2 both time-varying (pairwise combinations)
  # 20      : x1, x2, x3 all PH
  # 21-24   : x1,x2 PH, x3 time-varying
  # 25-34   : x1 PH, x2+x3 time-varying
  # 35+     : x1+x2+x3 all time-varying (triple combinations)
  if (scenario_id == 1) {
    vars <- c("1")
    forms <- c("ph")
    
    # 2-5     : x1 time-varying (step/spline/log/linear)
  } else if (scenario_id %in% 2:5) {
    vars <- c("1")
    forms <- c("step", "spline", "log", "linear")[scenario_id - 1]
    
    # 6-9     : x1 PH, x2 time-varying
  } else if (scenario_id %in% 6:9) {
    vars <- c("1", "2")
    forms <- c("ph", c("step", "spline", "log", "linear")[scenario_id - 5])
    
    # 10-19   : x1 and x2 both time-varying (pairwise combinations)
  } else if (scenario_id %in% 10:19) {
    vars <- c("1", "2")
    forms_list <- combn(c("step", "spline", "log", "linear"), 2, simplify = FALSE)
    forms <- forms_list[[scenario_id - 9]]
    
    # 20      : x1, x2, x3 all PH
  } else if (scenario_id == 20) {
    vars <- c("1", "2", "3")
    forms <- rep("ph", 3)
    
    # 21-24   : x1,x2 PH, x3 time-varying
  } else if (scenario_id %in% 21:24) {
    vars <- c("1", "2", "3")
    forms <- c("ph", "ph", c("step", "spline", "log", "linear")[scenario_id - 20])
    
    # 25-34   : x1 PH, x2+x3 time-varying
  } else if (scenario_id %in% 25:34) {
    vars <- c("1", "2", "3")
    forms_list <- combn(c("step", "spline", "log", "linear"), 2, simplify = FALSE)
    forms <- c("ph", forms_list[[scenario_id - 24]])
    
    # 35+     : x1+x2+x3 all time-varying (triple combinations)
  } else {
    vars <- c("1", "2", "3")
    forms_list <- combn(c("step", "spline", "log", "linear"), 3, simplify = FALSE)
    forms <- forms_list[[scenario_id - 34]]
  }
  
  # =============================
  # Output
  # =============================
  betas <- get_beta_names(vars, forms, n_total, scenario_id,
                          intercept_range,
                          beta_fixed_range,
                          beta_tde_range)
  basehaz <- get_basehaz(scenario_id, basehaz_range)
  haz_str <- build_haz_expr(vars, forms, basehaz, tdefun_library)
  haz_expr <- eval(parse(text = haz_str)) 
  
  return(list(
    betas = betas,
    haz_expr = haz_expr,
    basehaz = basehaz,
    vars = vars,
    forms = forms
  ))
}


s <- generate_scenario(scenario_id = 23, n_total = 10,
                       intercept_range = c(0, 1),
                       beta_fixed_range = c(0, 1),
                       beta_tde_range = c(0, 3),
                       basehaz_range = c(0, 0.5))

