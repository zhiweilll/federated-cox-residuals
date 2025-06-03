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

# ============
# Time-varying Effect Functions
# ============
tdefun_library <- list(
  linear = "t",
  log = "log(t + 1)",
  step = "as.numeric(t > 5)",
  spline = "sin(pi * t / 5)"
)

# ============
# Scenario Generator Function
# ============
generate_scenario <- function(scenario_id, n_total) {
  x <- data.frame(
    id = 1:n_total,
    x1 = rnorm(n_total),
    x2 = rnorm(n_total),
    x3 = rnorm(n_total)
  )
  
  get_beta_names <- function(vars, forms) {
    betas <- data.frame(intercept = rep(0.01, n_total))
    for (i in seq_along(vars)) {
      v <- vars[i]
      betas[[paste0("beta", v)]] <- rep(0.5, n_total)
      if (forms[i] != "ph") {
        betas[[paste0("beta", v, "_tde")]] <- rep(2, n_total)
      }
    }
    return(betas)
  }
  
  # build hazard expression string with beta extraction
  build_haz_expr <- function(vars, forms) {
    header <- c("beta0 <- betas[['intercept']]")
    terms <- c("beta0")
    for (i in seq_along(vars)) {
      v <- vars[i]
      if (forms[i] == "ph") {
        header <- c(header, sprintf("beta%s <- betas[['beta%s']]", v, v))
        terms <- c(terms, sprintf("beta%s * x[['x%s']]", v, v))
      } else {
        header <- c(header,
                    sprintf("beta%s <- betas[['beta%s']] + betas[['beta%s_tde']] * %s",
                            v, v, v, tdefun_library[[forms[i]]]))
        terms <- c(terms, sprintf("beta%s * x[['x%s']]", v, v))
      }
    }
    expr <- paste0("lp <- ", paste(terms, collapse = " + "))
    body <- c(header, expr, "return(0.01 * exp(lp))")
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
  
  betas <- get_beta_names(vars, forms)
  haz_str <- build_haz_expr(vars, forms)
  
  return(list(
    x = x,
    betas = betas,
    haz_expr = eval(parse(text = haz_str)) 
  ))
}


s <- generate_scenario(scenario_id = 23, n_total = 1000)

