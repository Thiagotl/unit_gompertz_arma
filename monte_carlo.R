# # # Checking the results
# # Case 2: with regressors
# 
# #initial parametres
# 
# phi=0.2 # AR
# theta=0.4 #MA
# alpha=1
# sigma=6 
# tau=0.5
# true_values=c(1,0.2,0.4,6)
# 
# #sample size
# n=c(70, 150, 300, 500)
# z=1.96
# R=100
# 
# mu_result<-matrix(NA,R,4)
# 
# for (i in 1:R) {
# 
#   y<-simu.ugoarma(n,phi=phi,theta=theta, alpha=alpha,sigma=sigma, tau=tau,freq=12,
#                     link="logit")
#   fit1<-uGoarma.fit(y)
#   mu_result[i,]<-fit1$coeff
# }
# 
# 
# mean_values<-c(apply(mu_result,2,mean))
# b_values<-(true_values-mean_values)/true_values*100
# eqm_values<-c(apply(mu_result,2,var))+(true_values-mean_values)^2
# result1<- cbind(true_values,
#                 mean_values,
#                 b_values,
#                 eqm_values
# )
# colnames(result1)<-c("true value","mean","bias","eqm")
# rownames(result1)<-c("alpha","phi","theta","sigma") 
# print(round(result1,5))
#########

rm(list = ls())


source("simu.ugoarma.R")
source("ugo_fit.R")

phi = 0.2
theta = 0.4
alpha = 1
sigma = 6
tau = 0.5
true_values = c(1, 0.2, 0.4, 6)
vn = c(30,35,40,50 )#150, 300, 500
R = 10
z = 1.96

results = list()
coverage = list()

for (n in vn) {
  mu_result = matrix(NA, R, 4)
  ci_lower = matrix(NA, R, 4)
  ci_upper = matrix(NA, R, 4)
  
  for (i in 1:R) {
    #print(c("i=",i))
    y <- simu.ugoarma(n, phi = phi, theta = theta, alpha = alpha, sigma = sigma, tau = tau, freq = 12, link = "logit")
    fit1 <- uGoarma.fit(y)
    mu_result[i, ] <- fit1$coeff
    
    # CI 
   
    ci_lower[i, ] <- fit1$coeff - z * fit1$stderror
    ci_upper[i, ] <- fit1$coeff + z * fit1$stderror
  }
  
  mean_values <- apply(mu_result, 2, mean)
  b_values <- (true_values - mean_values) / true_values * 100
  eqm_values <- apply(mu_result, 2, var) + (true_values - mean_values)^2
  
  coverage_rate = colMeans((true_values >= ci_lower) & (true_values <= ci_upper)) * 100
  
  result = cbind(true_values, mean_values, b_values, eqm_values, coverage_rate)
  colnames(result) <- c("true value", "mean", "bias", "eqm", "coverage")
  rownames(result) <- c("b1", "b2", "g1", "g2")
  
  results[[as.character(n)]] = result
}

# Impressão dos resultados
for (n in vn) {
  cat("\nSample size:", n, "\n")
  print(round(results[[as.character(n)]], 5))
}












fit1$zstat





