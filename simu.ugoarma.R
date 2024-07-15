simu.ugoarma <- function(n,phi=0.2,theta=0.4, alpha=1,sigma=6, tau=0.5,freq=12,
                         link="logit")
{
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
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
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

# 
# plot(simu.ugoarma(100))
# 
#y<-simu.ugoarma(10000)
#plot(simu.ugoarma(1000))
