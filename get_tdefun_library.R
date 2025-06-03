#' Generate Time-Varying Function Expressions
#'
#' Creates a list of time-varying effect expressions to be used in hazard functions.
#' The returned expressions are character strings representing how each variable's
#' time-dependent effect (beta(t)) evolves over time.
#'
#' @param linear_slope Numeric. Slope of the linear time function. Default is 1.
#' @param linear_intercept Numeric. Intercept of the linear time function. Default is 0.
#' @param log_shift Numeric. A positive shift value added inside the logarithmic function (e.g., log(t + shift)). Default is 1.
#' @param step_height Numeric. Height (multiplier) of the step function. Default is 1.
#' @param step_cutoff Numeric. Cutoff time for the step function (i.e., time at which the function jumps). Default is 5.
#' @param spline_freq Numeric. Frequency denominator for the sine-based spline function (i.e., sin(pi * t / freq)). Default is 5.
#'
#' @return A named list of character strings representing each time-dependent function:
#' \itemize{
#'   \item \code{linear}: Linear function, e.g., "1.000 * t + 0.000"
#'   \item \code{log}: Logarithmic function, e.g., "log(t + 1.000)"
#'   \item \code{step}: Step function, e.g., "1.000 * as.numeric(t > 5.000)"
#'   \item \code{spline}: Sine-based function, e.g., "sin(pi * t / 5.000)"
#' }
#'
#' @examples
#' generate_tdefun_library()
#' generate_tdefun_library(linear_slope = 1.2,linear_intercept = -0.5, log_shift = 2,step_height = 1.5,step_cutoff = 4,spline_freq = 6)
#' @export

get_tdefun_library <- function(
    linear_slope = 1,
    linear_intercept = 0,
    log_shift = 1,
    step_height = 1,
    step_cutoff = 5,
    spline_freq = 5
) {
  list(
    linear = sprintf("%.3f * t + %.3f", linear_slope, linear_intercept),
    log = sprintf("log(t + %.3f)", log_shift),
    step = sprintf("%.3f * as.numeric(t > %.3f)", step_height, step_cutoff),
    spline = sprintf("sin(pi * t / %.3f)", spline_freq)
  )
}
