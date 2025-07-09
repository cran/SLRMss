#' Diagnostic plots for SLRMss objects
#' 
#' This function provides four plot for residual analysis. The first plot shows the standardized residuals against the fitted values. The second one shows the standardized residuals against by their index. The third one presents QQ-normal plot of them and the last one shows their density estimate.
#'
#' @param x An object of class \code{SLRMss}.
#' @param H0 If TRUE, plot the graphics under null hypothesis, if FALSE, plot the graphics under alternative hypothesis (default).
#' @param xlab A vector containing the four x-axis titles.
#' @param ylab A vector containing the four y-axis titles.
#' @param main A vector containing the four main plot titles.
#' @param ... Currently ignored.
#'
#' @return Four diagnostic plots extracted from a SLRMss object.
#' @export
#' 
#' @importFrom stats qqnorm density
#' @importFrom graphics par abline points
#' @importFrom grDevices devAskNewPage
#'
#' @examples
#' data(orange)
#' fit <- SLRMss(emulsion ~ arabicgum + xanthangum + orangeoil, family="Student",
#' xi=3, testingbeta="xanthangum", statistic="LR", data=orange)
#' plot(fit)
plot.SLRMss <- function(x,H0=FALSE,
             xlab=c("Fitted Values","Index","Theoretical Quantiles","Standardized Residuals"),
             ylab=c("Standardized Residuals", "Standardized Residuals","Standardized Residuals","Density"),
             main=c("Residuals Against Fitted Values"," Residuals Against Index","Normal Q-Q Plot","Density Estimate"), ...){
    resid=residuals.SLRMss(x,std=TRUE,H0=H0)
    op <- grDevices::devAskNewPage(TRUE)
    on.exit(devAskNewPage(op))
    plot(fitted.SLRMss(x),resid,
         xlab=xlab[1],ylab=ylab[1],
         main=main[1],pch=16)
    graphics::abline(h=0,lty=2,col="red")
    plot(resid,xlab=xlab[2],
         ylab=ylab[2],main=main[2],pch=16)
    graphics::abline(h=0,lty=2,col="red")
    stats::qqnorm(resid,xlab=xlab[3],ylab=ylab[3],main=main[3],pch=16)
    graphics::abline(0,1,lty=2,col="red")
    plot(stats::density(resid),xlab=xlab[4],
         ylab=ylab[4],main=main[4])
    graphics::points(resid,rep(0,length(resid)),pch=16)
    invisible()
  }
