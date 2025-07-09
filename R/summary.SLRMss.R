#' SLRMss object summaries
#' 
#' This function displays the summary of the fitted model. It includes parameter estimates under both, null and alternative hypothesis, corrected Akaike information criterion, bayesian information criterion and choosed statistics.
#'
#' @param object An object of class \code{SLRMss}.
#' @param ... Currently ignored.
#'
#' @return A selected components extracted from a SLRMss object.
#' @export
#' 
#' @examples
#' data(orange)
#' fit <- SLRMss(emulsion ~ arabicgum + xanthangum + orangeoil, family="Student",
#' xi=3, testingbeta="xanthangum", statistic="LR", data=orange)
#' summary(fit)
summary.SLRMss <- function(object, ...){
    out=NULL
    out$call=match.call()
    out$beta=round(object$beta.coefficients,4)
    out$phi=round(object$phi,4)
    out$beta.h0=round(object$beta.coefficients.h0,4)
    out$phi.h0=round(object$phi.h0,4)
    out$null.hypotesis=object$nul 
    out$statistics=round(object$statistics,4)
    out$statistic.distribution=object$statistic.distribution
    out$df=object$df
    out$AICc=object$AIC
    out$BIC=object$BIC
    return(out)
  }
