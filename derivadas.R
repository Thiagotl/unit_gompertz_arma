tau=.5;mu=.85;sigma=2;x=.63
ll<-expression(
  log(log(tau)/(1-mu^-sigma))+log(sigma)+(-1-sigma)*log(x)+(log(tau)/(1-mu^-sigma))*(1-x^-sigma)
  # 1/(1-mu^-sigma)
  )
lmu<-D(ll,"mu")
eval(lmu)

-sigma*mu^(-sigma-1)/(1-mu^-sigma)-
  log(tau)*(1-x^-sigma)*sigma*mu^(sigma-1)/(1-mu^sigma)^2




ll<-expression(
  log(log(tau)/(1-mu^-sigma))+log(sigma)+(-1-sigma)*log(x)+(log(tau)/(1-mu^-sigma))*(1-x^-sigma)
  
)


lsigma<-D(ll,"sigma")
eval(lsigma)

1/sigma + log(mu) / (1 - mu^sigma) -log(x) +(mu^sigma * log(tau) * x^(-sigma) * 
                                                ((mu^sigma - 1) * log(x) - log(mu) * (x^sigma - 1))) / (mu^sigma - 1)^2









