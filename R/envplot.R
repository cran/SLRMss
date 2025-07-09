#' Quantile-quantile plots with simulated envelope of residuals for SLRMss objects
#' 
#' This function provides an envelope plot of a fitted SLRMss model.
#'
#' @param object An object of class \code{SLRMss}.
#' @param J The number of Monte Carlo replications. 100 by default.
#' @param conf The confidence level. 0.95 by default.
#' @param seed An optional seed for the simulation.
#' @param H0 If TRUE, choose the fitted values under null hypothesis, if FALSE, choose the fitted values under alternative hypothesis (default).
#' @param colors A vector with one or two characters. If it has one character, that represents the color of the plotted points. If it has two characters, the first one represents the color of the points out of the limits and the second one represents the color of the points under the limits. Red and green by default.
#' @param pch A vector with one or two numbers. If it has one numeric, that represents the plot \code{pch}. If it has two numbers, the first one represents the \code{pch} of the points out of the limits and the second one represents the \code{pch} of the points under the limits. 16 by default.
#' @param lty A vector with one or two numbers. If it has one number, that represents the \code{lty} of all lines. If it has two numbers, the first one represents the \code{lty} of the middle line and the second one represents the lty of the limits line. 2 by default.
#' @param xlab A title for the x axis.
#' @param ylab A title for the y axis.
#' @param main A title for the plot.
#'
#' @return Quantile-quantile plot with simulated envelope for a SLRMss object.
#' @export
#' 
#' @importFrom stats rnorm rlogis qqnorm rt as.formula quantile
#' @importFrom normalp rnormp
#' @importFrom graphics lines
#'
#' @examples
#' data(orange)
#' fit <- SLRMss(emulsion ~ arabicgum + xanthangum + orangeoil, family="Student",
#' xi=3, testingbeta="xanthangum", statistic="LR", data=orange)
#' envplot(fit)
envplot <- function (object, J=100, conf = 0.95, seed = NULL, H0=FALSE,colors=c("red","green"),pch=16,lty=2,xlab,ylab,main) 
{   
    if(!inherits(object,"SLRMss")){
    stop("Object class must be 'SLRMss'.")
    }
    fit=object
    family = fit$family
    if (family == "Normal") {
        rfam = function(x, mu, sigma) {
            stats::rnorm(x, mu, sigma)
        }
    }
    else {
        if (family == "Student") {
            rfam = function(x, mu, sigma) {
                sigma * stats::rt(x, df = fit$xi) + mu
            }
        }
        else {
            if (family == "Logistic") {
                rfam = function(x, mu, sigma) {
                  stats::rlogis(x, mu, sigma)
                }
            }
            else {
                rfam = function(x, mu, sigma) {
                  normalp::rnormp(x, mu, sigma, p = 2/(fit$xi + 1))
                }
            }
        }
    }
    set.seed(seed)
    y = fit$y
    n = length(y)
    fv = fitted.SLRMss(fit,H0=H0)
    sige = coef.SLRMss(fit,H0=H0)$phi
    rqobs <- residuals.SLRMss(fit,H0=H0,std=TRUE)
    mrq <- matrix(NA, J, n)
    for (j in 1:J) {
        Yj <- rfam(n, fv, sige)
        form = stats::as.formula(paste("Yj ~ ", paste(colnames(fit$X)[-1], 
            collapse = "+")))
        mj <- SLRMss(form, data = data.frame(Yj, fit$X), statistic = "Wald", 
            testingbeta = fit$testingbeta, family = family, 
            xi = fit$xi)
        mrq[j, ] <- residuals.SLRMss(mj,H0=H0,std=TRUE)
        mrq[j, ] <- sort(mrq[j, ])
    }
    infsup <- apply(mrq, 2, stats::quantile, probs = c((1 - conf)/2, 
        (1 + conf)/2), type = 6)
    media <- colMeans(mrq)
    faixay <- range(infsup)
    if(missingArg(xlab)) xlab="Standard Normal Quantiles"
    if(missingArg(ylab)) ylab="Standardized Residuals"
    if(missingArg(main)) main=paste0(100*conf,"% Confidence Envelope Plot")
    ylim=c(min(faixay,min(rqobs)),max(faixay,max(rqobs)))
    if(length(colors)>=2){
        if(length(colors)>2) warning("Only the first two colors entries were used.")
        colors=ifelse(sort(rqobs)<infsup[1,] | sort(rqobs) > infsup[2,],colors[1],colors[2])
    }
    if(length(pch)>=2){
        if(length(pch)>2) warning("Only the first two 'pch' entries were used.")
        pch=ifelse(sort(rqobs)<infsup[1,] | sort(rqobs) > infsup[2,],pch[1],pch[2])
    }
    qq0 <- stats::qqnorm(sort(rqobs), main = main, xlab = xlab,ylab=ylab, 
        col=colors,pch=pch, ylim = faixay)
    eixox <- sort(qq0$x)
    if(length(lty)>2)warning("Only the first two 'lty' entries were used")
    if(length(lty)==1) lty[2]=lty
    graphics::lines(eixox, media,lty=lty[1])
    graphics::lines(eixox, infsup[1, ],lty=lty[2])
    graphics::lines(eixox, infsup[2, ],lty=lty[2])
    out=NULL
    out$x=qq0$x
    out$y=qq0$y
    invisible(out)
}
