
#rm(list = ls())


uGoarma.fit<-function(y, ar = NA, ma = NA, tau = .5, link = "logit", h = 1, 
                      diag = 0, X = NA, X_hat = NA)
{
  source("ugo-functions.R")
  # adicionar o teste lógico para y
  # adicionar o teste lógico para a série temporal
  
  z<-c()
  maxit1<-50
  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p, q, na.rm = TRUE) 
  y1 <- y[(m+1):n]
  p1 <- length(ar)
  q1 <- length(ma)
  error <- rep(0,n)
  eta <- rep(NA,n)
  
  y_prev <- c(rep(NA, (n+h)))
  
  linktemp <- substitute(link)
  
  # verifica se o que foi escrito para a funcao de ligação é um valor válido
  
  if(!is.character(linktemp)){
    linktemp <- deparse(linktemp)
    if (linktemp == 'link'){
      linktemp <- eval(link)
    }
  }
  
  valid_links<-c("logit", "probit", "cloglog")
  if (linktemp %in% valid_links) {
    stats <- make.link(linktemp)
  } else {
    stop(paste(linktemp, "link not available, available links are", 
               paste(valid_links, collapse = ", ")))
  }
  
  link = linktemp 
  linkfun = stats$linkfun
  linkinv = stats$linkinv 
  mu.eta = stats$mu.eta 
  diflink = function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  ynew = linkfun(y) #g (y) avaliado em y - funcao de ligacao
  ynew_ar <- suppressWarnings(matrix(ynew,(n-1),max(p,1,na.rm=T)))
  
  
  ## --------
  
  # Ajuste das variaveis e dados para um modelo autoregressivo
  
  # Verifica se há elemntos NA em `ar`
  
  if(any(is.na(ar)) == F) {
    names_phi <- c(paste("phi", ar, sep = "")) 
    Z <- suppressWarnings(matrix(ynew, (n-1), p1)[m:(n-1),])} else {
      ar = p1<-0; Z <- NA  
    } 
  
  # Verifica se há elemntos NA em `am`
  if(any(is.na(ma)) == F) {
    names_theta <- c(paste("theta", ma, sep = ""))
  } else ma = q1 <- 0 
  
  # Verifica se há elemntos NA em `X` - Covariaveis
  if(any(is.na(X)) == F){
    names_beta<-c(paste("beta", 1 : ncol(as.matrix(X)), sep = ""))
    Xm <- X[(m+1):n, ]     
    k = ncol(X)
  } else {
    k = 0 
    X <- matrix(rep(0,n), nrow = n)
    Xm <- NA
  }
  
  # Recorrences
  q_1 <- max(q1, 1)
  R <- matrix(rep(NA, (n-m)*q_1), ncol = q_1)  
  k_i <- q1/q_1                                 
  
  deta.dalpha <- rep(0, n)
  deta.dbeta <- matrix(0, ncol=max(k,1), nrow=n)
  deta.dphi <- matrix(0, ncol=p1, nrow=n)
  deta.dtheta <- matrix(0, ncol=q_1, nrow=n)
  
  Xstart <- (cbind(rep(1, (n-m)), Xm, Z))         
  Xstart <- matrix(apply(Xstart, 1, na.omit),nrow = (n-m),byrow = T)
  ols <- lm.fit(Xstart, ynew[(m+1) : n])$coef
  initial <- c(rep(0, k+p1+q1+1),1)
  initial[1 : (k+p1+1)] <- ols
  
  
  loglik <- function(z) 
  {
    alpha <- z[1]
    if(k==0)  beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else  theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    sigma <- z[length(z)]
    
    #if (is.na(sigma) || sigma <= 0) stop("problemas com sigma") 
    
    Xbeta <- X %*% beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    
    for(i in (m+1):n)
    {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar]) %*% phi + 
        t(theta) %*% error[i-ma]
      error[i] <- ynew[i] - eta[i] 
    }
    mu <- linkinv(eta[(m+1):n])
    x <- y[(m+1):n]
    
    ll <- log(log(tau)/(1 - mu^-sigma)) + log(sigma) +
      (-1 - sigma) * log(x) + (log(tau)/(1 - mu^-sigma)) * (1 - x^-sigma)
    sum(ll)
  } 
  #############################################################################
  
  # VETOR SCORE 
  escore.UGoarma <- function(z)
  {
    alpha <- z[1]
    if(k==0)  beta = as.matrix(0) else beta = as.matrix(z[2:(k+1)])
    if(p1==0) {phi = as.matrix(0);ar=1} else phi = as.matrix(z[(k+2):(k+p1+1)]) 
    if(q1==0) theta = as.matrix(0) else  theta = as.matrix(z[(k+p1+2):(k+p1+q1+1)])
    #c_par <- z[length(z)]
    
    Xbeta <- X%*%beta
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    for(i in (m+1):n)
    {
      eta[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1),ar] - Xbeta_ar[(i-1),ar])%*%phi + t(theta)%*%error[i-ma]
      error[i] <- ynew[i] - eta[i] 
    }
    
    # TENHO QUE VERIFICAR 
    sigma <- z[length(z)]
    # q_t <- linkinv(eta[(m+1):n])
    mu <- linkinv(eta[(m+1):n])
    # x  <- y[(m+1):n] 
    x<-y1
    Xbeta <- X%*%beta
    for(i in 1:(n-m)){
      R[i,] <- error[i+m-ma]*k_i}
    
    for(i in (m+1):n)
    {
      deta.dalpha[i] <- 1 - deta.dalpha[i-ma]%*%theta
      deta.dbeta[i,] <- X[i,] - t(phi)%*%X[i-ar,] - t(theta)%*%deta.dbeta[i-ma,]
      deta.dphi[i,] <- ynew_ar[i-ar]- Xbeta[i-ar] - t(theta)%*%deta.dphi[i-ma,]
      deta.dtheta[i,] <- R[(i-m),] - t(theta)%*%deta.dtheta[i-ma,]
    }
    
    v <- deta.dalpha[(m+1):n]
    rM <- deta.dbeta[(m+1):n,]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    
    mT <- diag(mu.eta(eta[(m+1):n]))
   

    #ell_q
    
    a_t <- -sigma*mu^(-sigma-1)/(1-mu^-sigma)-log(tau)*(1-x^-sigma)*sigma*mu^(sigma-1)/(1-mu^sigma)^2
      
      
    #ell_c_par
    y_sust <- as.vector(1/sigma + log(mu) / (1 - mu^sigma) -log(x) +(mu^sigma * log(tau) * x^(-sigma) * 
                                                                       ((mu^sigma - 1) * log(x) - log(mu) * (x^sigma - 1))) / (mu^sigma - 1)^2)
    
    Ualpha <- t(v) %*% mT %*% a_t
    Ubeta <- t(rM) %*% mT %*% a_t
    Uphi <-   t(rP) %*% mT %*% a_t
    Utheta <- t(rR) %*% mT %*% a_t
    Uc <- sum(y_sust)
    
    rval <- c(Ualpha,Ubeta,Uphi,Utheta,Uc)
    return(rval[rval!=0])
  }
  
  
  ##############################################################################
  #ATENCAO AQUI - USO DO VETOR SCORE
  opt<-optim(initial, loglik,
             escore.UGoarma , #mudar aqui
             method = "BFGS", hessian = TRUE,
             control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
  
  
  
  if (opt$conv != 0)
  {
    warning("FUNCTION DID NOT CONVERGE WITH ANALITICAL GRADIENT!")
    opt<-optim(initial, loglik, 
               method = "BFGS", hessian = TRUE,
               control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12))
    if (opt$conv != 0)
    {
      warning("FUNCTION DID NOT CONVERGE NEITHER WITH NUMERICAL GRADIENT!")
    }else{
      warning("IT WORKS WITH NUMERICAL GRADIENT!")
    }
  }
  
  z$conv <- opt$conv
  coef <- (opt$par)[1:(p1+q1+k+2)]
  alpha <- coef[1]
  if(k==0) beta=names_beta=NULL else z$beta <- coef[2:(k+1)]
  if(p1==0) phi=names_phi=NULL else z$phi <- coef[(k+2):(k+p1+1)]
  if(q1==0) theta=names_theta=NULL else z$theta <- coef[(k+p1+2):(k+p1+q1+1)]
  z$sigma <- coef[length(coef)]
  
  names_par <- c("alpha",names_beta,names_phi,names_theta,"sigma")
  names(coef)<-names_par
  z$coeff <- coef
  J_inv <- solve(-(opt$hessian))
  z$stderror<-sqrt(diag(J_inv))
  z$zstat <- z$coef/z$stderror
  z$pvalues <- 2*(1 - pnorm(abs(z$zstat)) )
  z$loglik <- opt$value
  z$counts <- as.numeric(opt$counts[1])
  

  if(any(is.na(X)==F))
  {
    z$k<- (p1+q1+k+2)
    z$aic <- -2*(z$loglik)+2*(z$k)
    z$bic <- -2*(z$loglik)+log(n)*(z$k)
    z$hq <- -2*(z$loglik)+log(log(n))*(z$k)
  }else{
    z$k<- (p1+q1+2)
    z$aic <- -2*(z$loglik)+2*(z$k)
    z$bic <- -2*(z$loglik)+log(n)*(z$k)
    z$hq <- -2*(z$loglik)+log(log(n))*(z$k)
  }
  
  model_presentation <- cbind(round(z$coef,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")
  z$model <- model_presentation
  
  # Predicted values   (NO REGRESSORS)
  if(k==0){                                     
    alpha <- as.numeric(coef[1])
    phi <- as.numeric(coef[2:(p1+1)])
    theta <- as.numeric(coef[(p1+2):(p1+q1+1)])
    sigma <- as.numeric(coef[p1+q1+2])  
    
    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    z$sigma <- sigma
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    if(p1==0) {phi = as.matrix(0);ar=1}
    if(q1==0) {theta = as.matrix(0);ma=1}
    for(i in (m+1):n)
    {
      etahat[i]<-alpha + ynew_ar[(i-1),ar]%*%as.matrix(phi) + (theta%*%errorhat[i-ma])
      errorhat[i]<- ynew[i]-etahat[i] # predictor scale
    }
    
    q_hat <- linkinv(etahat[(m+1):n])   # fitted values 
    
    z$fitted <- ts(c(rep(NA,m),q_hat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    
    # Forecasting
    ynew_prev <- c(ynew,rep(NA,h))
    y_prev[1:n] <- z$fitted
    
    for(i in 1:h)
    {
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0 
    }
    
    z$forecast <- y_prev[(n+1):(n+h)]
  } else{                              # with REGRESSORS
    X_hat <- as.matrix(X_hat)
    
    alpha <- as.numeric(coef[1])
    beta <- as.numeric(coef[2:(k+1)])
    phi <- as.numeric(coef[(k+2):(k+p1+1)])
    theta <- as.numeric(coef[(k+p1+2):(k+p1+q1+1)])
    sigma <- as.numeric(coef[length(coef)])  
    
    z$alpha <- alpha
    z$beta <- beta
    z$phi <- phi
    z$theta <- theta
    z$sigma <- sigma
    
    errorhat<-rep(0,n) # E(error)=0 
    etahat<-rep(NA,n)
    
    if(p1==0) {phi = as.matrix(0);ar=1}
    if(q1==0) {theta = as.matrix(0);ma=1}
    
    Xbeta <- X%*%as.matrix(beta)
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p, 1, na.rm = T)))
    
    for(i in (m+1):n)
    {
      etahat[i] <- alpha + Xbeta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar])%*%as.matrix(phi) +
        t(as.matrix(theta))%*%errorhat[i-ma]
      errorhat[i] <- ynew[i] - etahat[i]
    }
    q_hat <- linkinv(etahat[(m+1):n])   # fitted values 
    
    z$fitted <- ts(c(rep(NA,m),q_hat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat
    
    
    # Forecasting
    ynew_prev <- c(ynew,rep(NA,h))
    y_prev[1:n] <- z$fitted
    
    X_prev <- rbind(X,X_hat)
    
    for(i in 1:h)
    {
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) +
        (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta))) +
        (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }
    
    z$forecast <- y_prev[(n+1):(n+h)]
    
  }
  
  
  # Quantile residuals 
  z$residuals <- as.vector(qnorm(pUGo(y[(m+1):n],z$fitted[(m+1):n],z$sigma)))  #mudar aqui
  residc <- z$residuals 
  
  # # GRAPHICS  ---- Comentar 
  # 
  # if(diag>0)
  # {
  #   
  #   print(model_presentation)
  #   print(" ",quote=F)
  #   print(c("Log-likelihood:",round(z$loglik,4)),quote=F)
  #   print(c("Number of iterations in BFGS optim:",z$counts),quote=F)
  #   print(c("AIC:",round(z$aic,4)," SIC:",round(z$bic,4)," HQ:",round(z$hq,4)),quote=F)
  #   
  #   print("Residuals:",quote=F)
  #   print(summary(residc))
  #   
  #   par(mfrow=c(1,1))
  #   plot(y,type="l",ylab="Serie",xlab="Time",ylim=c(min(y),max(y)))
  #   lines(z$fitted,col="blue",lty=2)
  #   legend("topright",c("Observed data","Predicted median"),
  #          pt.bg="white", lty=c(1,2), bty="n",col=c(1,"blue"))
  #   
  #   
  #   w1<-5
  #   h1<-4
  #   
  #   if(diag>1)
  #   {
  #     postscript(file = "resid_v_ind.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1))
  #       par(mgp=c(1.7, 0.45, 0))
  #       plot(residc,main=" ",xlab="Index",ylab="Residuals", pch = "+",
  #            ylim=c(-4,4))
  #       lines(t,rep(-3,n+h+6),lty=2,col=1)
  #       lines(t,rep(3,n+h+6),lty=2,col=1)
  #       lines(t,rep(-2,n+h+6),lty=3,col=1)
  #       lines(t,rep(2,n+h+6),lty=3,col=1)
  #     }
  #     dev.off()
  #     
  #     postscript(file = "adjusted.eps",horizontal=F,paper="special",width = w1, height = h1,family = "Times")
  #     {
  #       par(mfrow=c(1,1))
  #       par(mar=c(2.8, 2.7, 1, 1))  
  #       par(mgp=c(1.7, 0.45, 0))
  #       plot(y,type="l",ylab="Serie",xlab="Time")
  #       lines(z$fitted,col=2,lty=2)
  #       legend("topright",c("Observed data","Predicted median"),
  #              pt.bg="white", lty=c(1,2), bty="n",col=c(1,2), cex=.8)
  #       
  #     }
  #     dev.off()
  #     
  #   }    
  # }  
  # 
  return(z)
  
}



# PARA TESTAR 

# set.seed(2)
# 
# source("simu.ugoarma.R")
# 
# y<-simu.ugoarma(100,phi=0.2,theta=0.4, alpha=1,sigma=6, tau=0.5,freq=12,link="logit")
# 
# fit<-uGoarma.fit(y, ma=1, ar=1)
# 
# fit$model





