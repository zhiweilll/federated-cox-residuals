# ============ Scenario Label Generator ============
generate_scenario_labels <- function() {
  labels <- list()
  vars <- list()
  id <- 1
  
  # Univariate
  labels[[id]] <- "x1: PH"; vars[[id]] <- 1; id <- id + 1
  for (type in c("step", "spline", "log", "linear")) {
    labels[[id]] <- paste0("x1: ", type); vars[[id]] <- 1; id <- id + 1
  }
  
  # Bivariate: x1 PH, x2 non-PH
  for (type in c("step", "spline", "log", "linear")) {
    labels[[id]] <- paste0("x1: PH, x2: ", type); vars[[id]] <- 2; id <- id + 1
  }
  
  # Bivariate: x1 + x2 both non-PH
  combos2 <- combn(c("step", "spline", "log", "linear"), 2, simplify = FALSE)
  for (cmb in combos2) {
    labels[[id]] <- paste0("x1: ", cmb[1], ", x2: ", cmb[2]); vars[[id]] <- 2; id <- id + 1
  }
  
  # Trivariate: all PH
  labels[[id]] <- "x1: PH, x2: PH, x3: PH"; vars[[id]] <- 3; id <- id + 1
  
  # Trivariate: x1,x2 PH, x3 non-PH
  for (type in c("step", "spline", "log", "linear")) {
    labels[[id]] <- paste0("x1: PH, x2: PH, x3: ", type); vars[[id]] <- 3; id <- id + 1
  }
  
  # Trivariate: x1 PH, x2+x3 non-PH
  combos2 <- combn(c("step", "spline", "log", "linear"), 2, simplify = FALSE)
  for (cmb in combos2) {
    labels[[id]] <- paste0("x1: PH, x2: ", cmb[1], ", x3: ", cmb[2]); vars[[id]] <- 3; id <- id + 1
  }
  
  # Trivariate: x1+x2+x3 all non-PH
  combos3 <- combn(c("step", "spline", "log", "linear"), 3, simplify = FALSE)
  for (cmb in combos3) {
    labels[[id]] <- paste0("x1: ", cmb[1], ", x2: ", cmb[2], ", x3: ", cmb[3])
    vars[[id]] <- 3
    id <- id + 1
  }
  
  tibble(scenario_id = seq_along(labels),
         label = unlist(labels),
         n_vars = unlist(vars))
}


# generate_scenario_labels() %>% print(n = Inf)
