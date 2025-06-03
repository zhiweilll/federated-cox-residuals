#' Generate Random Beta Coefficients and Hazard Function
#'
#' This function generates the time-fixed and time-varying coefficients for each covariate,
#' as well as a character string and parsed R expression for the hazard function.
#' It supports proportional hazards and multiple types of time-dependent effects.
#'
#' @param seed Integer. Seed used to generate all random values deterministically.
#' @param n_total Integer. Number of observations (used to repeat betas into a dataframe).
#' @param vars Character vector. Covariate variable indices (e.g., `"1"`, `"2"`, `"3"`).
#' @param forms Character vector. Time-dependency form for each variable (e.g., `"ph"`, `"log"`).
#' @param intercept_range Numeric vector of length 2. Range to draw the intercept coefficient.
#' @param beta_fixed_range Numeric vector of length 2. Range to draw the time-fixed coefficients.
#' @param beta_tde_range Numeric vector of length 2. Range to draw time-dependent effect coefficients.
#' @param basehaz_range Numeric vector of length 2. Range to sample the baseline hazard rate.
#' @param tdefun_library Named list. A list of time-transformation expressions (e.g., from `get_tdefun_library()`).
#'
#' @return A list with:
#' \describe{
#'   \item{betas}{A `data.frame` with intercept and per-variable coefficients.}
#'   \item{haz_str}{The character string for the constructed hazard function.}
#'   \item{haz_expr}{The parsed and evaluated hazard function expression.}
#'   \item{basehaz}{The sampled baseline hazard value.}
#' }
#'
#' @examples
#' scenario <- generate_scenario_structure(10)
#' beta_haz <- generate_beta_and_haz(
#'   seed = 12,
#'   n_total = 100,
#'   vars = scenario$vars,
#'   forms = scenario$forms,
#'   intercept_range = c(0.4, 0.6),
#'   beta_fixed_range = c(0.4, 0.6),
#'   beta_tde_range = c(1.5, 2.5),
#'   basehaz_range = c(0.01, 0.05),
#'   tdefun_library = get_tdefun_library()
#' )
#'
#' @export

generate_beta_and_haz <- function(seed, n_total, vars, forms,
                                  intercept_range, beta_fixed_range,
                                  beta_tde_range, basehaz_range,
                                  tdefun_library) {
  # ===================
  # Set Random Seed and Draw Parameters
  # ===================
  set.seed(seed)  
  
  intercept <- runif(1, intercept_range[1], intercept_range[2])
  beta_fixed <- runif(length(vars), beta_fixed_range[1], beta_fixed_range[2])
  beta_tde <- runif(length(vars), beta_tde_range[1], beta_tde_range[2])
  basehaz <- runif(1, basehaz_range[1], basehaz_range[2])
  
  
  # ===================
  # Construct Beta Data Frame
  # ===================
  betas <- data.frame(intercept = rep(intercept, n_total))
  for (i in seq_along(vars)) {
    v <- vars[i]
    betas[[paste0("beta", v)]] <- rep(beta_fixed[i], n_total)
    if (!is.na(forms[i]) && forms[i] != "ph") {
      betas[[paste0("beta", v, "_tde")]] <- rep(beta_tde[i], n_total)
    }
  }
  
  # ===================
  # Construct Hazard Function (as string and evaluated expression)
  # ===================
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
  haz_str <- paste(c("function(t, x, betas, ...) {", paste0("  ", body), "}"), collapse = "\n")
  haz_expr <- eval(parse(text = haz_str))
  
  return(list(
    betas = betas,
    haz_str = haz_str,
    haz_expr = haz_expr,
    basehaz = basehaz
  ))
}

beta_haz <- generate_beta_and_haz(
  seed = 12,
  n_total = 100,
  vars = scenario$vars,
  forms = scenario$forms,
  intercept_range = c(0.4, 0.6),
  beta_fixed_range = c(0.4, 0.6),
  beta_tde_range = c(1.5, 2.5),
  basehaz_range = c(0.01, 0.05),
  tdefun_library = get_tdefun_library())
