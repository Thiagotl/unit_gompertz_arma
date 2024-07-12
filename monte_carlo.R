# # Checking the results
# Case 2: with regressors

phi=0.2;theta=0.4; alpha=1;sigma=6; tau=0.5
true_values=c(1,0.2,0.4,6)
n=50
R<-100
mu_result<-matrix(NA,R,4)
for (i in 1:R) {

  y<-simu.ugoarma(n,phi=phi,theta=theta, alpha=alpha,sigma=sigma, tau=tau,freq=12,
                    link="logit")
  fit1<-uGoarma.fit(y)
  mu_result[i,]<-fit1$coeff
}


mean_values<-c(apply(mu_result,2,mean))
b_values<-(true_values-mean_values)/true_values*100
eqm_values<-c(apply(mu_result,2,var))+(true_values-mean_values)^2
result1<- cbind(true_values,
                mean_values,
                b_values,
                eqm_values
)
colnames(result1)<-c("true value","mean","bias","eqm")
rownames(result1)<-c("alpha","phi","theta","sigma") 
print(round(result1,5))
#########


# Parâmetros iniciais
phi = 0.2
theta = 0.4
alpha = 1
sigma = 6
tau = 0.5
true_values = c(1, 0.2, 0.4, 6)

# Tamanhos das amostras
sample_sizes = c(70, 150, 300, 500)
R = 100  # Número de simulações

# Função para calcular a cobertura do intervalo de confiança de 95%
calculate_coverage <- function(estimates, true_value, alpha = 0.05) {
  lower_bound <- apply(estimates, 2, function(x) quantile(x, probs = alpha / 2))
  upper_bound <- apply(estimates, 2, function(x) quantile(x, probs = 1 - alpha / 2))
  coverage <- mean((true_value >= lower_bound) & (true_value <= upper_bound))
  return(coverage)
}

# Matriz para armazenar resultados
results <- list()

# Loop para diferentes tamanhos de amostra
for (n in sample_sizes) {
  mu_result <- matrix(NA, R, 4)
  
  for (i in 1:R) {
    y <- simu.ugoarma(n, phi = phi, theta = theta, alpha = alpha, sigma = sigma, tau = tau, freq = 12, link = "logit")
    fit1 <- uGoarma.fit(y)
    mu_result[i, ] <- fit1$coeff
  }
  
  mean_values <- apply(mu_result, 2, mean)
  b_values <- (true_values - mean_values) / true_values * 100
  eqm_values <- apply(mu_result, 2, var) + (true_values - mean_values)^2
  cr95_values <- sapply(1:4, function(i) calculate_coverage(mu_result[, i], true_values[i]))
  
  result <- cbind(true_values, mean_values, b_values, eqm_values, cr95_values)
  colnames(result) <- c("true value", "mean", "bias", "eqm", "cr95")
  rownames(result) <- c("alpha", "phi", "theta", "sigma")
  
  results[[paste0("n=", n)]] <- result
}

# Exibindo os resultados
for (n in sample_sizes) {
  print(paste("Sample size:", n))
  print(round(results[[paste0("n=", n)]], 5))
}

