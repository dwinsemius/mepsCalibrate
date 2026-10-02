#' Extract Cross-Partial Derivatives from Bivariate Tensor Hazards
#'
#' @description Computes the second-order cross-partial derivative 
#'   (d2f / dx dy) of an estimated multivariate smooth surface from a GAM or PAMM.
#'   Useful for evaluating departures from the homogeneity assumption over 
#'   bivariate spatial topologies (like height and weight).
#'
#' @param model A fitted \code{mgcv::gam} or PAMM survival object.
#' @param focal_x Character string indicating the first coordinate axis variable (e.g., "calibrated_height").
#' @param focal_y Character string indicating the second coordinate axis variable (e.g., "calibrated_weight").
#' @param data A data frame containing evaluation nodes at which to compute derivatives.
#' @param eps Numeric value setting the finite difference step size. Defaults to \code{1e-05}.
#'
#' @return A data frame identical to \code{data} with appended columns for the cross-partial 
#'   derivative point estimates (\code{derivative}) and standard errors (\code{se}).
#'
#' @importFrom stats predict vcov
#' @export
extract_cross_partials <- function(model, focal_x, focal_y, data, eps = 1e-05) {
  if (!inherits(model, "gam")) {
    stop("Model must be an object of class 'gam' from the mgcv package.")
  }
  
  if (!all(c(focal_x, focal_y) %in% names(data))) {
    stop("Focal coordinate variables must be present within evaluation data parameters.")
  }
  
  data_pp <- data
  data_pm <- data
  data_mp <- data
  data_mm <- data
  
  data_pp[[focal_x]] <- data_pp[[focal_x]] + eps
  data_pp[[focal_y]] <- data_pp[[focal_y]] + eps
  
  data_pm[[focal_x]] <- data_pm[[focal_x]] + eps
  data_pm[[focal_y]] <- data_pm[[focal_y]] - eps
  
  data_mp[[focal_x]] <- data_mp[[focal_x]] - eps
  data_mp[[focal_y]] <- data_mp[[focal_y]] + eps
  
  data_mm[[focal_x]] <- data_mm[[focal_x]] - eps
  data_mm[[focal_y]] <- data_mm[[focal_y]] - eps
  
  X_pp <- stats::predict(model, newdata = data_pp, type = "lpmatrix")
  X_pm <- stats::predict(model, newdata = data_pm, type = "lpmatrix")
  X_mp <- stats::predict(model, newdata = data_mp, type = "lpmatrix")
  X_mm <- stats::predict(model, newdata = data_mm, type = "lpmatrix")
  
  X_deriv <- (X_pp - X_pm - X_mp + X_mm) / (4 * eps^2)
  
  beta <- coef(model)
  derivative_estimates <- as.numeric(X_deriv %*% beta)
  
  V_beta <- stats::vcov(model)
  derivative_variances <- rowSums((X_deriv %*% V_beta) * X_deriv)
  derivative_se <- sqrt(pmax(0, derivative_variances))
  
  output <- data
  output$derivative <- derivative_estimates
  output$se         <- derivative_se
  
  return(output)
}
