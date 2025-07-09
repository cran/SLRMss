#' Symmetric Linear Regression Models for small samples
#' 
#' Computes Wald, Likelihood Ratio, Score, or Gradient statistics for symmetric linear regression models. Also computes modified versions of the Likelihood Ratio, Score, and Gradient tests for small sample sizes.
#'
#' @param formula An object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param family A description of the error distribution to be used in the model. There are four supported families, Normal, t-Student, Power Exponential and Logistic ("Normal", "Student", "Powerexp" and "Logistic", respectively)
#' @param xi An extra parameter of some specified error distribution. For t-Student is a positive value and for Power Exponential is a real number between -1 and 1/3.
#' @param statistic The statistic which will be used. It includes "Wald", "LR", "Score" or "Gradient".
#' @param testingbeta A vector containing the names of the variables to be testing.
#' @param data An optional data frame containing the variables in the model.
#'
#' @return A list with the following components
#' @return \item{beta.coefficients}{A matrix with the estimated position parameters under alternative hypothesis.}
#' @return \item{phi}{A numeric value with the estimated precision paramater under alternative hypothesis.}
#' @return \item{beta.coefficients.h0}{A matrix with the estimated position parameters under null hypothesis.}
#' @return \item{phi.h0}{A numeric value with the estimated precision paramater under null hypothesis.}
#' @return \item{y.fitted}{A vector with the fitted values of the model.}
#' @return \item{null.hypothesis}{The description of the null hypothesis.}
#' @return \item{statistics}{A matrix with the selected statistics and theirs p-values. The corrected statistic is marked with an asterisk.}
#' @return \item{statistic.distribution}{The name of the statistics' distribution used to test null hypothesis. It always return "Chi-Squared".}
#' @return \item{df}{The degrees of freedom of the statistics' distribution. It's the length of the testingbeta vector.}
#' @return \item{residuals}{The difference among the real y values and the fitted y.}
#' @return \item{std.residuals}{The residuals divided by the precision parameter}
#' @return \item{AICc}{The corrected Akaike Information Criterion for small samples.}
#' @return \item{BIC}{The Bayesian Information Criterion.}
#' 
#' @export
#' 
#' @importFrom ssym ssym.l
#' @importFrom stats nlminb optim nlm sd model.matrix pchisq
#' @importFrom methods missingArg
#' 
#' @references Medeiros, F. M. C and Ferrari, S. L. P. (2017). Small-sample testing inference in symmetric and log-symmetric linear regression models, Statistica Neerlandica. \doi{10.1111/stan.12107}.
#' 
#' @examples
#' data(orange)
#' fit1 <- SLRMss(emulsion ~ arabicgum + xanthangum + orangeoil, family="Student", 
#' xi=3, testingbeta="xanthangum", statistic="LR", data=orange)
#' print(fit1)
#' 
#' data(cheese)
#' fit2 <- SLRMss(cohe ~ fat + xangum + sodcase, family="Normal", 
#' testingbeta=c("xangum","sodcase"), statistic="Gradient", data=cheese)
#' print(fit2)
SLRMss <- function(formula, family, xi, statistic, testingbeta,data){
    if(!is.character(testingbeta)){
      stop("'testingbeta' must be in character.")
    }
    
    if(family=="Normal"|family=="Student"|family=="Powerexp"){
      
      f1 = toString(formula[2])
      f2 = toString(formula[3])
      if(methods::missingArg(data)) a <- ssym::ssym.l(formula, family=family, xi=xi)
      else a <- ssym::ssym.l(formula, family=family, xi=xi, data=data)
      y= a$y
      n= length(y)
      X= a$model.matrix.mu
      phichapeu = exp(a$theta.phi/2)
      betachapeu = a$theta.mu
      echapeu= y-X%*%betachapeu
      betash0 = which(colnames(X)%in%testingbeta)
      if(sum(colnames(X)%in%testingbeta)==0){
        stop("No variable in 'testingbeta' is in 'formula'.")
      }
      q = sum(colnames(X)%in%testingbeta)
      if(length(testingbeta)>q){
        warning(paste0("Variable(s) ",paste(testingbeta[!(testingbeta %in% colnames(X))],collapse=",")," ignored."))
      }
      testingbeta=testingbeta[(testingbeta %in% colnames(X))]
      p = length(betachapeu)
      m = p+1
      
      if(q < p & all(betash0<=p)){
        X1 = as.matrix(X[,(colnames(X)%in%testingbeta)])
        colnames(X1)=colnames(X)[colnames(X)%in%testingbeta]
        X2 = as.matrix(X[,!colnames(X)%in%testingbeta])
        colnames(X2)=colnames(X)[!colnames(X)%in%testingbeta]
        R= X1 - X2%*%solve(crossprod(X2))%*%crossprod(X2,X1)
        if(colnames(X2)[1]=="(Intercept)"){
          if(ncol(X2)==1){
            formula2=formula(paste(f1,"~","1"))
            }else{ 
             formula2=formula(paste(f1,"~",paste0(colnames(X2)[-1],collapse = "+")))
            }
        }else{
          formula2=formula(paste(f1,"~","-1 + ",paste0(colnames(X2),collapse = "+")))
        }
        b=ssym.l(formula2, data=data,family=family,xi=xi)
        betatil=b$theta.mu
        phitil=exp(b$theta.phi/2)
        etil=y-X2%*%betatil
        beta1 = betachapeu[betash0]
        beta2 = betachapeu[-betash0]
        Z2 = X2%*%solve(crossprod(X2))%*%t(X2)
        
      }else{
        if (q == p & all(betash0<=p)){
          X1 = X
          X2 = 0
          R = X
          beta1= betachapeu
          beta2=0
          betatil = 0
          phitil = stats::sd(y)
          etil= y
          Z2 = diag(0,n)
        }else{
          if(q>p){
            stop("More parameters to be testing than variables in 'formula'.")
          }
          else{
            stop("There is no corresponding variable.")
          }
        }
      }
      
      switch(family,
             Normal={
               xi=NULL
               h= function(u){
                 exp(-u/2)/sqrt(2*pi)
               }
               Wchapeu=diag(1,n)
               Wtil=diag(1,n)
               delta00010=0;delta20000 = 1;delta01002=-1;delta00103=0;
               delta00101=0;delta01000=-1;delta00012=0;delta11001=1;
               delta21000=-1;delta20002=3;delta30001=-3;delta40002=15;delta21002=-3
             },
             
             Student={
               h= function(u){
                 xi^(xi/2)/beta(1/2,xi/2)*(xi+u)^(-(xi+1)/2)
               }
               Wchapeu=diag(c((xi+1)/(xi+(echapeu/phichapeu)^2)),nrow=n)
               Wtil=diag(c((xi+1)/(xi+(etil/phitil)^2)),nrow=n)
               delta00010=6*(xi+1)*(xi+2)/(xi*(xi+5)*(xi+7));delta20000 = (xi+1)/(xi+3)
               delta01002=(3-xi)/(xi+3);delta00103=6*(3*xi-5)/((xi+3)*(xi+5))
               delta00101=6*(xi+1)/((xi+3)*(xi+5));delta01000=-(xi+1)/(xi+3)
               delta00012=6*(xi^2-12*xi-13)/((xi+3)*(xi+5)*(xi+7));delta11001=(xi+1)*(xi-1)/((xi+3)*(xi+5))
               delta21000=-(xi+1)^3*(xi+2)/(xi*(xi+3)*(xi+5)*(xi+7));delta20002 = 3*(xi+1)/(xi+3)
               delta30001=-3*(xi+1)^2/((xi+3)*(xi+5));delta40002=15*(xi+1)^3/((xi+3)*(xi+5)*(xi+7))
               delta21002=3*(xi+1)^2*(3-xi)/((xi+3)*(xi+5)*(xi+7))
             },
             
             Powerexp={
               h=function(u){
                 C=gamma(1+(1+xi)/2)*2^(1+(1+xi)/2)
                 return(exp(-u^(1/(1+xi))/2)/C)
               }
               Wchapeu=diag(c(1/((1+xi)*(((echapeu/phichapeu)^2)^(xi/(1+xi))))),nrow=n)
               Wtil=diag(c(1/((1+xi)*(((etil/phitil)^2)^(xi/(1+xi))))),nrow=n)
               if(-1<xi & xi<1/3){
                 aux1 = (1+xi);aux2= gamma(aux1/2)
                 delta00010=2^(1-2*xi)*xi*(1-xi)*gamma((1-3*xi)/2)/(aux1^4*aux2)
                 delta20000=2^(1-xi)*gamma((3-xi)/2)/(aux1^2*aux2)
                 delta01002=-(1-xi)/aux1; delta00103=2*xi*(1-xi)/aux1^2
                 delta00101=2^(2-xi)*xi*gamma((3-xi)/2)/(aux1^3*aux2)
                 delta01000=-2^(1-xi)*gamma((3-xi)/2)/(aux1^2*aux2)
                 delta00012=-2^(2-xi)*xi*(1+3*xi)*gamma((3-xi)/2)/(aux1^4*aux2)
                 delta11001=2^(1-xi)*(1-xi)*gamma((3-xi)/2)/(aux1^3*aux2)
                 delta21000=-2^(1-2*xi)*(1-xi)*gamma((3-3*xi)/2)/(aux1^4*aux2)
                 delta20002=(3+xi)/aux1
                 delta30001=-2^(2-xi)*gamma((5-xi)/2)/(aux1^3*aux2)
                 delta40002=2^(3-xi)*gamma((7-xi)/2)/(aux1^4*aux2)
                 delta21002=2^(2-xi)*(xi-1)*gamma((5-xi)/2)/(aux1^4*aux2)
               }else{
                  stop("'xi' must be in (-1, 1/3).")  
               }
             }
             
      )
      
      llog = function(phi,beta,X){
        if(all(beta==0)){
          z=y/phi
        }else{
          z=(y-X%*%beta)/phi
        }
        g= log(h(z^2))
        - n * log (phi) + sum(g)
      }
      
    }else{
      if(family=="Logistic"){
        
        f1=toString(formula[2])
        fo=formula(paste("~",f1,"-1"))
        if(methods::missingArg(data)){
          X=stats::model.matrix(formula)
          y=stats::model.matrix(fo)
        }else{
          X= stats::model.matrix(formula,data=data)
          y= stats::model.matrix(fo,data=data)
        }
        n=length(y)
        betash0 = which(colnames(X)%in%testingbeta)
        if(sum(colnames(X)%in%testingbeta)==0){
          stop("No variable in 'testingbeta' is in 'formula'.")
        }
        q = sum(colnames(X)%in%testingbeta)
        if(length(testingbeta)>q){
          warning(paste0("Variable(s) ",paste(testingbeta[!(testingbeta %in% colnames(X))],collapse=",")," ignored."))
        }
        testingbeta=testingbeta[(testingbeta %in% colnames(X))]
        beta0=solve(crossprod(X))%*%crossprod(X,y)
        phi0=as.numeric(sqrt(crossprod(y-X%*%beta0)/n))
        m=length(beta0)+1
        xi=NULL
        h=function(u){
          exp(sqrt(u))/(1+exp(sqrt(u)))^2
        }
        W=function(beta,phi,X){
          z=abs((y-X%*%beta)/phi)
          diag(c((exp(z)-1)/(z*(1+exp(z)))))
        }
        delta00010=1/15;delta20000=1/3;delta01002=-0.43;delta00103=0.645
        delta00101=1/6;delta01000=-1/3;delta00012=-0.114;delta11001=1/6
        delta21000=-1/15;delta20002=2.43;delta30001=-2/3;delta40002=1.9913
        delta21002=-0.2193
        
        llog = function(phi,beta,X){
          if(all(beta==0)){
            z=y/phi
          }else{
            z=(y-X%*%beta)/phi
          }
          g= log(h(z^2))
          - n * log (phi) + sum(g)
        }
        
        llog1=function(x){
          beta<-x[-m]
          phi<-x[m]
          -llog(phi,beta,X)
        }
        
        bp=try(stats::nlminb(c(beta0,phi0),llog1),silent=T)
        if(!is.list(bp)){
          bp=try(stats::nlm(llog1,c(beta0,phi0)),silent=T)
          if(!is.list(bp)){
            bp=stats::optim(c(beta0,phi0),llog1,method="BFGS")
            betachapeu= bp$par[-m]
            phichapeu= bp$par[m]
          }else{
            betachapeu= bp$estimate[-m]
            phichapeu= bp$estimate[m]
          }
        }else{
          betachapeu= bp$par[-m]
          phichapeu= bp$par[m]
        }
        
        
        echapeu= y-X%*%betachapeu
        Wchapeu=W(betachapeu,phichapeu,X)
        p = length(betachapeu)
        
        
        if(q < p & all(betash0<=p)){
          X1 = as.matrix(X[,(colnames(X)%in%testingbeta)])
          colnames(X1)=colnames(X)[colnames(X)%in%testingbeta]
          X2 = as.matrix(X[,!colnames(X)%in%testingbeta])
          colnames(X2)=colnames(X)[!colnames(X)%in%testingbeta]
          R= X1 - X2%*%solve(crossprod(X2))%*%crossprod(X2,X1)
          beta1 = betachapeu[betash0]
          beta2 = betachapeu[-betash0]
          beta02=solve(crossprod(X2))%*%crossprod(X2,y)
          phi02=as.numeric(sqrt(crossprod(y-X2%*%beta02)/n))
          m2=length(beta02)+1
          
          llog2 = function(x){
            beta<-x[-m2]
            phi<-x[m2]
            -llog(phi,beta,X2)
          }
          
          bp2=try(stats::nlminb(c(beta02,phi02),llog2),silent=T)
          if(!is.list(bp2)){
            bp2=try(stats::nlm(llog2,c(beta02,phi02)),silent=T)
            if(!is.list(bp2)){
              bp2=stats::optim(c(beta02,phi02),llog2,method="BFGS")
              betatil=bp2$par[-m2]
              phitil=bp2$par[m2]
            }else{
              betatil=bp2$estimate[-m2]
              phitil=bp2$estimate[m2]
            }
          }else{
            betatil=bp2$par[-m2]
            phitil=bp2$par[m2]
          }
          
          
          etil=y-X2%*%betatil
          Z2 = X2%*%solve(crossprod(X2))%*%t(X2)
          Wtil=W(betatil,phitil,X2)
          
        }else{
          if (q == p & all(betash0<=p)){
            X1 = X
            X2 = 0
            R = X
            beta1= betachapeu
            beta2=0
            betatil = 0
            phitil = sd(y)
            etil= y
            Z2 = diag(0,n)
            Wtil=W(betatil,phitil,X2)
          }else{
            if(q>p){
              stop("More parameters to be testing than variables in 'formula'.")
            }else{
              stop("There is no corresponding variable.")
            }
          }
        }
        
        
      }else{
        stop("Family is not supported.")
      }
    }
    
    Z= X%*%solve(crossprod(X))%*%t(X)
    Zd = diag(c(diag(Z)),ncol=n)
    Z2d = diag(c(diag(Z2)),ncol=n)
    rozz=n*sum(diag(Zd%*%Zd))
    roz2z2=n*sum(diag(Z2d%*%Z2d))
    rozz2=n*sum(diag(Zd%*%Z2d))
    
    
    if(statistic=="Wald"){
      st = delta20000*t(beta1)%*%crossprod(R)%*%beta1/phichapeu^2
      stcor = st
    }else{
      
      if(statistic=="LR"){
        st = 2*(llog(phichapeu,betachapeu,X) - llog(phitil,betatil,X2))
        
        m1=delta01002 -1; m2= 4 - delta00103 - 6*delta01002;
        m3=(delta00101 + 2*delta01000)/delta20000;m4=(delta00012-6*delta11001)/delta20000
        
        d0=delta00010/(4*delta20000^2);d2=-m3^2/(2*m1)
        d1=-(m2*m3)/(2*m1^2) - (2*m3+m3^2+m4)/(2*m1)
        
        ALR= d0/(n*q)*(rozz - roz2z2)
        ALRbp= d1/n + (d2*(2*p-q))/(2*n)
        aLR=ALR + ALRbp
        stcor = st*(1-aLR)
        
      }else{
        if(statistic=="Score"){
          st = crossprod(etil,Wtil)%*%X1%*%solve(crossprod(R))%*%crossprod(X1,Wtil)%*%etil/(delta20000*phitil^2)
          
          b0=delta21000/(delta20000^2) + 1
          b1= delta11001*(delta11001 - delta01000)/(delta20000^2*(delta20002- 1))
          b2=(2*delta11001*(2*delta01002 + delta00103)+(delta20002 -1)*(4*delta30001 +delta40002 + delta21002 -        2 *delta01000))/(delta20000  *(delta20002-1)^2)
          b3=delta11001^2/(delta20000^2*(delta20002-1))
          AR1 = 12*b0/n*(rozz2 - roz2z2)
          AR2= -9*b0/n*(rozz - 2*rozz2 + roz2z2)
          AR1bp = 12*b1/n * q * (p-q) -6*b2/n*q
          AR2bp = -12*b3/n *q *(q+2)
          AR11 = AR1 + AR1bp
          AR22 = AR2 + AR2bp
          bR = AR22/(12*q*(q+2))
          cR= (AR11-AR22)/(12*q)
          stcor = st * (1 - (cR + bR*st))
        }else{
          if(statistic=="Gradient"){
            st = crossprod(etil,Wtil)%*%X1%*%beta1/phitil^2
            
            
            m1=delta01002 -1; m2= 4 - delta00103 - 6*delta01002;
            m3=(delta00101 + 2*delta01000)/delta20000;m4=(delta00012-6*delta11001)/delta20000
            
            c0= delta00010/delta20000^2
            c1= - m3^2/m1
            c2= - (m2*m3+2*m1*m3)/m1^2 - m4/m1
            
            
            AT1= 6*c0/n*(rozz2 - roz2z2)
            AT2= -3*c0/n*(rozz - 2*rozz2 + roz2z2)
            AT1bp = 6*c1/n*q*(p-q)+6*c2/n*q
            AT2bp = -3*c1/n*q*(q+2)
            AT11 = AT1 + AT1bp
            AT22 = AT2 + AT2bp
            bT = (AT22)/(12*q*(q+2))
            cT = (AT11 - AT22)/(12*q)
            stcor = st*(1-(cT+bT*st))
          }else{
            stop("The available statistics are: 'Wald', 'LR','Score' or 'Gradient'.")
          }
        }
      }
    }
    
    pst = 1-stats::pchisq(st,q)
    pstcor = 1-stats::pchisq(stcor,q)
    if(st < 0 | stcor < 0){
      stop("Did not converge.")
    }
    
    aic = -2*llog(phichapeu,betachapeu,X) + 2*m+2*m*(m+1)/(n-m-1)
    bic = -2*llog(phichapeu,betachapeu,X) + m *log(n)
    
    out=NULL
    
    if(statistic!="Wald"){
      ests = c(st,stcor)
      vs = c(pst,pstcor)
      estatisticas = matrix(c(ests,vs),ncol=2)
      rownames(estatisticas)=c(statistic,paste(statistic,"*",sep=""))
      colnames(estatisticas)=c("Observed value","p-value")
      estatisticas = as.table(t(estatisticas))
    }else{
      estatisticas = matrix(c(st,pst),ncol=2)
      rownames(estatisticas)=statistic
      colnames(estatisticas)=c("Observed value","p-value")
      estatisticas = as.table(t(estatisticas))
    }
    
    kbinv = sqrt(diag(phichapeu^2 * solve(crossprod(X))/delta20000))
    mu = matrix(c(betachapeu,kbinv),ncol=2)
    rownames(mu)= colnames(X)
    colnames(mu)=c("Estimate","Std. error")
    
    kbinv2 = sqrt(diag(phitil^2 * solve(crossprod(X2))/delta20000))
    mu0 = matrix(c(betatil,kbinv2),ncol=2)
    rownames(mu0)= colnames(X2)
    colnames(mu0)=c("Estimate","Std. error")
    
    kbinvp = sqrt(phichapeu^2/(n*(delta20002 - 1)))
    phi = matrix(c(phichapeu,kbinvp),nrow = 1)
    colnames(phi) = c("Estimate","Std. error")
    rownames(phi) = "phi"
    
    kbinvp2 = sqrt(phitil^2/(n*(delta20002 - 1)))
    phi0 = matrix(c(phitil,kbinvp2),nrow = 1)
    colnames(phi0) = c("Estimate","Std. error")
    rownames(phi0) = "phi"
    
    betatest = testingbeta
    text1= matrix(paste(betatest,"= 0"),nrow=q)
    colnames(text1)=""
    rownames(text1)=rep("",q)
    
    
    distribution = "Chi-Squared"
    
    df= q
    
    residuals=y-X%*%betachapeu
    std.residuals=residuals/(phichapeu)
    
    
    call=match.call()
    out=list(formula=formula,residuals=residuals,std.residuals=std.residuals,beta.coefficients=mu,
             phi=phi,beta.coefficients.h0=mu0,phi.h0=phi0,X=X,y=y,
             y.fitted=X%*%betachapeu,null.hypothesis=text1,statistics=estatisticas,
             statistic.distribution=distribution,df=df,family=family,
             xi=xi,testingbeta=testingbeta,AICc=aic,BIC=bic,call=call)
    class(out)="SLRMss"
    return(out)
  }
