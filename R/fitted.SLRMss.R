#' Extract model fitted values for a SLRMss object
#' 
#' This function provides the fitted values of a SLRMss model.
#'
#' @param object An object of class \code{SLRMss}.
#' @param H0 If TRUE, choose the fitted values under null hypothesis, if FALSE, choose the fitted values under alternative hypothesis (default).
#' @param ... Currently ignored.
#'
#' @return Fitted values extracted from the SLRMss object.
#' @export
#'
#' @examples
#' data(orange)
#' fit <- SLRMss(emulsion ~ arabicgum + xanthangum + orangeoil, family="Student", xi=3,
#' testingbeta="xanthangum", statistic="LR", data=orange)
#' fitted(fit)
fitted.SLRMss <- function(object,H0=FALSE, ...){
  if(H0==FALSE){
  fitted=object$y.fitted
  }else{
  fitted=object$X[,colnames(object$X)%in%rownames(object$beta.coefficients.h0)]%*%object$beta.coefficients.h0[,1]
  }
  return(fitted)
 }
