# PROBABILITY DENSITY FUNCTION

dUGo <- function(x, mu=.5, sigma=1.2, tau=.5)
{
  
  fx<- log(tau)/(1 - mu^(-sigma)) * sigma * x^(-(1 + sigma)) * exp((log(tau) / (1 - mu^(-sigma))) *
                                                                     (1 - x^(-sigma)))
  return(fx)
  
}

integrate(dUGo,0.01, 0.99) 

# CUMULATIVE DISTRIBUTION FUNCTION 

pUGo <-  function(q, mu=.5, sigma=1.2, tau=.5) 
{
  cdf <- exp(((log(tau)) / (1 - mu^(-sigma))) * (1 - q^(-sigma)))
  return(cdf)
}

pUGo(.25)
integrate(dUGo, 0, .25)

# QUANTILE FUNCTION
qUGo<-function(u, mu=.5, sigma=1.2, tau=.5)
{
  q<-(1 - (log(u) * (1 - mu^(-sigma)) / log(tau)))^(-1 / sigma)
  q
}

u=pUGo(.82)
qUGo(u,mu=.5, sigma=1.2)


# INVERSION METHOD FOR RANDOM GENERATION

rUGO <- function(n, mu=.5, sigma=1.2, tau=.5) {
  u <- runif(n)
  y <- (1 - (log(u) * (1 - mu^(-sigma)) / log(tau)))^(-1 / sigma)
  return(y)
}

#rUGO(100)