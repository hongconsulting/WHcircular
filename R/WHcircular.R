#' @useDynLib WHcircular, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Periodic B-spline basis for "HH:MM" times
#'
#' Generates the periodic B-spline basis matrix for times-of-day in "HH:MM" string
#' format with the specified interior knots and degree, evaluated at the values of x.
#' @param x String vector of times in "HH:MM" format.
#' @param knots Numeric vector of interior knot locations in hours 
#' (boundary knots are fixed at 0 and 24). Default = 6, 12, 18.
#' @param degree Degree of the piecewise polynomial. Default = 3.
#' @details Wrapper function for \code{splines2::bSpline()}.
#' @return A numeric matrix of \code{length(x)} rows and \code{length(knots)} columns.
#' @export
HHMM.Bspline <- function(x, knots = NULL, degree = 3) {
  if (!requireNamespace("splines2", quietly = TRUE)) {
    stop("[HHMM.Bspline] requires package 'splines2'")
  }
  HM <- do.call(rbind, strsplit(x, ":", fixed = TRUE))
  hrs <- (as.numeric(HM[,1])*60 + as.numeric(HM[,2])) / 60
  if (is.null(knots)) knots <- c(6, 12, 18)
  b <- splines2::bSpline(hrs, knots = knots, degree = degree,
                         Boundary.knots = c(0, 24), periodic = TRUE, intercept = FALSE)
  return(b)
}

#' Circular mean of times-of-day in "HH:MM" format
#'
#' Computes the circular mean of times-of-day in "HH:MM" format,
#' rounded to the nearest minute.
#' @param x String vector of times-of-day in "HH:MM" format.
#' @details Wrapper function for \code{WH_HHMM_to_rad()},
#' \code{WH_rad_mean()}, and \code{WH_rad_to_HHMM()}.
#' @return The circular mean time-of-day in "HH:MM" string format.
#' @examples
#' y <- c("01:00","01:15","01:30","01:45","02:00",
#'        "23:00","23:15","23:30","23:45","00:00")
#' print(HHMM.mean(y))
#' @export
HHMM.mean <- function(x) {
  return(WH_rad_to_HHMM(WH_rad_mean(WH_HHMM_to_rad(x))))
}

#' Circular-linear regression with "HH:MM" outcomes
#'
#' Fits a circular-linear regression¹ model for time-of-day outcomes in
#' "HH:MM" format.
#' @param y String vector of times-of-day in "HH:MM" format.
#' @param X.formula A model formula specifying the predictors.
#' @param name Optional string vector of variable names.
#' @details Wrapper function for \code{WH_HHMM_to_rad()},
#' \code{WH_rad_to_dHHMM()}, \code{WH_rad_to_HHMM()}, and 
#' \code{WH_reg_circlinear()}.
#' @return A string matrix with columns:
#' \itemize{
#'   \item "": Parameter names
#'   \item "Time (95%CI)": Parameter estimates and 95% confidence intervals
#'   \item "p": \emph{P}-values
#' }
#' The μ intercept is presented as a time-of-day whereas predictor coefficients
#' are presented as signed time-of-day differences.
#' @references
#' 1. Fisher, N.I. and Lee, A.J., 1992. Regression models for an angular response.
#' \emph{Biometrics}, pp. 665–677.
#' @examples
#' y <- c("01:00","01:15","01:30","01:45","02:00",
#'        "23:00","23:15","23:30","23:45","00:00")
#' X <- c(rep(0, 5), rep(1, 5))
#' print(HHMM.mean(y))
#' print(HHMM.lm(y, ~X))
#' @export
HHMM.lm <- function(y, X.formula, name = NULL) {
  X <- as.matrix(stats::model.matrix(X.formula)[, -1])
  if (is.null(name)) name <- c("(Intercept)", paste0("\u03b2", 1:(ncol(X))))
  fit <- WH_reg_circlinear(WH_HHMM_to_rad(y), X)
  coef <- fit$coef[2:length(fit$coef)]
  SE <- fit$SE[2:length(fit$SE)]
  p <- c(NA, fit$p)
  LL <- coef - stats::qnorm(0.975) * SE
  UL <- coef + stats::qnorm(0.975) * SE
  output <- matrix("", ncol = 3, nrow = length(coef) + 1)
  output[1, 2] <- "Time (95%CI)"
  output[1, 3] <- "p"
  f <- function(x) 2 * atan(x)
  f.signif <- function(x, digits = 2, threshold = 0.001) {
    if (x < threshold) return(paste0("< ", toString(threshold)))
    output <- formatC(signif(x, digits), format = "fg", digits = digits, flag = "#")
    output[x == 0] <- paste0("0.", strrep("0", digits - 1))
    return(sub("\\.$", "", output))
  }
  output[2, 1] <- name[1]
  output[2, 2] <- paste0(WH_rad_to_HHMM(coef[1]), " (", 
                         WH_rad_to_HHMM(LL[1]), " to ", 
                         WH_rad_to_HHMM(UL[1]), ")")
  for (i in 2:length(coef)) {
    output[1 + i, 1] <- name[i]
    output[1 + i, 2] <- paste0(WH_rad_to_dHHMM(f(coef[i])), " (", 
                               WH_rad_to_dHHMM(f(LL[i])), " to ", 
                               WH_rad_to_dHHMM(f(UL[i])), ")")
    output[1 + i, 3] <- f.signif(p[i])
  }
  return(output)
}
