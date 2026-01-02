##############################################################
#
# Created by Thiago Tavares Lopes (thiago.tavares@acad.ufsm.br)
# Reviewed by Renata Rojas Guerra (renata.r.guerra@ufsm.br)
# December/2025
#
##############################################################
simu.ugoarma <- function(n, phi = 0.2, theta = 0.4, alpha = 1, sigma = 6, 
                         tau = 0.5, freq=12, beta = 0, X = NA,
                         link="logit"){
  source("ugo-functions.R")
  
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }else{
    ar<-0
    phi<-0
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }else{
    ma<-0
    theta<-0
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\"
               and \"cloglog\""))
  }
  
  link <- structure(list(link = linktemp,
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv ))
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  
  
  beta <- as.matrix(beta, 1, 2)
  
  ##### X definitions
  if (is.na(X)) {
    X<-matrix(0, c(n,1))
  }else{
    if(X=="cos"){
      X=as.matrix(cos(2*pi*(1:n)/12))
      if(beta==0) stop("Inform the value of beta")
    }
  }
  
  
  
  
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    error<-rep(0,n+m)
    eta<- y <- rep(NA,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha  + as.numeric(phi%*%ynew[i-ar]) + 
        as.numeric(theta%*%error[i-ma]) + X[i-m,] %*% beta
      # print(as.numeric(theta%*%error[i-ma]))
      # print(eta[i])
      mu[i]   <- linkinv(eta[i])
      y[i]    <- rUGO(1,mu[i],sigma,tau) #mudar aqui
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]
      
    }
    
    
    return( ts(y[(m+1):(n+m)],frequency=freq) )
  }
}




n = 1000
y <- simu.ugoarma(n, 
                  phi = 0.96, 
                  theta = NA,
                  alpha= -0.09,
                  sigma = 25, 
                  tau=0.5,
                  freq=12, 
                  beta = 0.03, link="logit", X ="cos")

X_cos <- matrix(cos(2*pi*(1:length(y))/12))

source("ugo_fit.R")

fit <- uGoarma.fit(y, ma = NA, ar = 1, X = X_cos)
fit$model


