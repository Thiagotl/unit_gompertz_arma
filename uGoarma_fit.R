
uGoarma.fit<-function(y, ar = NA, ma = NA, tau = .5, link = "logit", h = 0, 
                      diag = 0, X = NA, X_hat = NA)
{
  source("ugo-functions.r")
  # adicionar o teste lógico para y
  # adicionar o teste lógico para a série temporal
  
  z<-c()
  max_it1<-50
  p <- max(ar)
  q <- max(ma)
  n <- length(x)
  m <- max(p, q, na.rm = TRUE) 
  y1 <- y[(m+1):n]
  p1 <- length(ar)
  q1 <- length(ma)
  error <- rep(0,n)
  eta <- rep(NA,n)
  
  y_prev <- c(rep(NA, (n+h)))
  
  linktemp <- substitute(link)
  
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
  ynew = linkfun(y) 
  ynew_ar <- suppressWarnings(matrix(ynew,(n-1),max(p,1,na.rm=T)))
  
  
  ## --------
  
  if(any(is.na(ar)) == F) {
    names_phi <- c(paste("phi", ar, sep = ""))
    Z <- suppressWarnings(matrix(ynew, (n-1), p1)[m:(n-1),])} else {
      ar = p1<-0; Z <- NA  
    } 
  
  if(any(is.na(ma)) == F) {
    names_theta <- c(paste("theta", ma, sep = ""))
  } else ma = q1 <- 0 
  
  if(any(is.na(X)) == F){
    names_beta<-c(paste("beta", 1 : ncol(as.matrix(X)), sep = ""))
    Xm <- X[(m+1):n, ]     
    k = ncol(X)
  } else {
    k = 0 
    X <- matrix(rep(0,n), nrow = n)
    Xm <- NA
  }
  
  
  
}











