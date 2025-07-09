#' Extract residuals model for SLRMss objects
#' 
#' This function provides the residuals of a SLRMss model.
#'
#' @param object An object of class \code{SLRMss}.
#' @param H0 If TRUE, choose the residuals under null hypothesis, if FALSE, choose the residuals under alternative hypothesis (default).
#' @param std If TRUE, choose the standardized residuals, if FALSE, choose the non-standardized residuals (default).
#' @param ... Currently ignored.
#'
#' @return Residuals extracted from a SLRMss object.
#' @export
#'
#' @examples
#' data(orange)
#' fit <- SLRMss(emulsion ~ arabicgum + xanthangum + orangeoil, family="Student",
#' xi=3, testingbeta="xanthangum", statistic="LR", data=orange)
#' residuals(fit)
residuals.SLRMss <- function(object,H0=FALSE,std=FALSE, ...){
   res=(object$y-fitted.SLRMss(object,H0=H0))
   if(std==TRUE){
   return(res/coef.SLRMss(object,H0=H0)$phi)  
   }else{
   return(res)  
   }
}
