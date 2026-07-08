#' Extract Model Coefficients for SLRMss Objects
#' 
#' This function provides the coefficients of a \code{SLRMss} model.
#'
#' @param object An object of class \code{SLRMss}.
#' @param H0 If TRUE, choose the coefficients under the null hypothesis, if FALSE, choose the coefficients under alternative hypothesis (default).
#' @param ... Currently ignored.
#'
#' @return Coefficients extracted from the SLRMss object.
#' @export
#'
#' @examples
#' data(orange)
#' fit <- SLRMss(emulsion ~ arabicgum + xanthangum + orangeoil, family="Student",
#' xi=3, testingbeta="xanthangum", statistic="LR", data=orange)
#' coef(fit)
coef.SLRMss <- function(object,H0=FALSE, ...){
  out=NULL
  if(H0==FALSE){
  out$beta=object$beta.coefficients[,1]
  out$phi=object$phi[,1]
  }else{
  out$beta=object$beta.coefficients.h0[,1]
  out$phi=object$phi.h0[,1]
  }
  return(out)
} 
