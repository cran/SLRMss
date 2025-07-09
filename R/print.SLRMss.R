#' Print values of a SLRMss object
#' 
#' This function displays a succinct summary of the fitted model. It includes the mean and dispersion parameter estimates, corrected Akaike information criterion and Bayesian information criterion.
#'
#' @param x An object of class \code{SLRMss}.
#' @param ... Currently ignored.
#'
#' @return Coefficients, AICc and BIC extracted from a SLRMss object. 
#' @export
#'
#' @examples
#' data(orange)
#' fit <- SLRMss(emulsion ~ arabicgum + xanthangum + orangeoil, family="Student", xi=3,
#' testingbeta="xanthangum", statistic="LR", data=orange)
#' print(fit)
print.SLRMss <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nBeta Coefficients:\n")
  print(round(x$beta.coefficients[,1],4))
  cat("\nPhi:\n") 
  print(round(x$phi[,1],4))
  cat("\nAICc:\n")
  print(x$AIC)
  cat("\nBIC:\n")
  print(x$BIC)
  invisible(x)
}
