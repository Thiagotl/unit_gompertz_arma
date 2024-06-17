
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
  
  if(!is.character(linktemp))
  {
    linktemp<-deparse
  }
  
  link = linktemp


}