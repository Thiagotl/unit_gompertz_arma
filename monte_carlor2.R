rm(list = ls())


source("simu.ugoarma.R")
source("ugo_fit.R")


phi = 0.2 # AR
theta = 0.4 # MA
alpha = 1
sigma = 6
tau = 0.5
true_values = c(1, 0.2, 0.4, 6)

# Tamanhos da amostra ajustados
vn = c(30, 35, 40, 50) 
z = 1.96
R = 100


set.seed(2024)

for (n in vn){
  #matriz de resultados
  estim <-ICi<-ICs<- err <- matrix(NA, nrow = R, ncol = length(true_values))
  #contadores
  calpha<-cphi<-ctheta<-csigma<-0 # talvez preciso colocar o tau
  i<-0
  
  for(i in 1:R){
    
    y <- simu.ugoarma(n, phi = phi, theta = theta, alpha = alpha, sigma = sigma, tau = tau, freq = 12, link = "logit")
    fit1 <- uGoarma.fit(y)
    
    #i<-i+1
    print(c("i=",i))
    
    estim[i,]<-fit1$model[,1]
    err[i,]<-fit1$model[,2]
    
    #IC
    ICi[i,]<-estim[i,]-(z*err[i,])
    ICs[i,]<-estim[i,]+ (z*err[i,])
    
    if (ICi[i,1]<=alpha && ICi[i,1]>=alpha) 
    {
      calpha<-calpha+1
    }
    
    if (ICi[i,2]<= phi && ICs[i,2]>=phi) 
    {
      cphi<-cphi+1
    }
    
    if (ICi[i,3]<=theta && ICs[i,3]>=theta) 
    {
      ctheta<-ctheta+1
    }
    
    if (ICi[i,4]<=sigma && ICs[i,4]>=sigma) 
    {
      csigma<-csigma+1
    }
    
  }
  
  ### mean
  m <- apply(estim, 2, mean)
  ### bias
  bias <- (true_values-m)
  ### relative percentage bias
  biasP <- bias/true_values *100
  ### SD
  erro <- apply(estim, 2, sd)
  ### MSE
  MSE <- apply(estim, 2, var)+bias^2
  ### 
  TC<-c(calpha,
        cphi,
        ctheta,
        csigma)/R
  ## final results
  results <- rbind(m, bias, biasP, erro, MSE,TC)
  rownames(results) <- c("Mean", "Bias","RB%", "SE", "MSE","TC")
  print(c("Tamanho da Amostra:",n))
  print(round(results,4))
  
  
}

