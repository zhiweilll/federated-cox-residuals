#' Generate Variable Structure for a Given Scenario
#'
#' This function maps a scenario ID to its corresponding variable configuration 
#' and time-dependency form (PH or time-varying). Each scenario defines how many 
#' covariates are included (x1, x2, x3) and what time-varying transformation (if any) 
#' is applied to each.
#'
#' ### Scenario Index Map:
#' - **1**       : x1 PH (univariate)
#' - **2–5**     : x1 time-varying (step, spline, log, linear)
#' - **6–9**     : x1 PH, x2 time-varying
#' - **10–19**   : x1 & x2 both time-varying (all 2-way combinations)
#' - **20**      : x1, x2, x3 all PH
#' - **21–24**   : x1 & x2 PH, x3 time-varying
#' - **25–34**   : x1 PH, x2 & x3 time-varying (all 2-way combinations)
#' - **35+**     : x1, x2, x3 all time-varying (all 3-way combinations)
#'
#' @param scenario_id Integer. A number representing the scenario configuration.
#'
#' @return A list with:
#' \describe{
#'   \item{vars}{Character vector of variable indices involved (e.g., `"1"`, `"2"`, `"3"`).}
#'   \item{forms}{Character vector indicating time-dependency for each variable (e.g., `"ph"`, `"log"`).}
#'   \item{scenario_label}{A descriptive string summarizing the scenario setting.}
#' }
#'
#' @examples
#' generate_scenario_structure(1)
#' generate_scenario_structure(23)
#'
#' @export
#' 
generate_scenario_structure <- function(scenario_id) {
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
  
  scenario_label <- paste0("Scenario ", scenario_id, " (", 
                           paste0("x", vars, "=", forms, collapse = ", "), ")")
  
  return(list(
    vars = vars,
    forms = forms,
    scenario_label = scenario_label
  ))
}


scenario = generate_scenario_structure(23)
