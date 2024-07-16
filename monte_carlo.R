rm(list = ls())

source("simu.ugoarma.R")
source("ugo_fit.R")

alpha = 1
phi = 0.2
theta = 0.4
sigma = 6
tau = 0.5
true_values = c(1, 0.2, 0.4, 6) # alpha, phi, theta, sigma
vn = c(70, 150, 300, 500) # 150, 300, 500
R = 10000
z = 1.96

for (n in vn) {
  # matriz de resultados
  estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
  # contadores
  calpha <- cphi <- ctheta <- csigma <- 0 # para guardar os resultados
  i <- 0
  bug <- 0 # inicializa o contador de bugs
  
  for (i in 1:R) {
    y <- simu.ugoarma(n, phi = phi, theta = theta, alpha = alpha, sigma = sigma, tau = tau, freq = 12, link = "logit")
    fit1 <- try(uGoarma.fit(y), silent = TRUE)
    
    if (class(fit1) == "try-error" || fit1$conv != 0) {
      bug <- bug + 1
    } else {
      estim[i, ] <- fit1$model[, 1]
      err[i, ] <- fit1$model[, 2]
      
      if (!any(is.na(estim[i, ])) && !any(is.na(err[i, ]))) {
        # IC
        ICi[i, ] <- estim[i, ] - (z * err[i, ])
        ICs[i, ] <- estim[i, ] + (z * err[i, ])
        
        if (ICi[i, 1] <= alpha && ICs[i, 1] >= alpha) {
          calpha <- calpha + 1
        }
        
        if (ICi[i, 2] <= phi && ICs[i, 2] >= phi) {
          cphi <- cphi + 1
        }
        
        if (ICi[i, 3] <= theta && ICs[i, 3] >= theta) {
          ctheta <- ctheta + 1
        }
        
        if (ICi[i, 4] <= sigma && ICs[i, 4] >= sigma) {
          csigma <- csigma + 1
        }
      }
    }
  }
  
  ### mean
  m <- apply(estim, 2, mean, na.rm = TRUE)
  ### bias
  bias <- (true_values - m)
  ### relative percentage bias
  biasP <- bias / true_values * 100
  ### SD
  erro <- apply(estim, 2, sd, na.rm = TRUE)
  ### MSE
  MSE <- apply(estim, 2, var, na.rm = TRUE) + bias^2
  ###
  TC <- c(calpha, cphi, ctheta, csigma) / R
  ## final results
  results <- rbind(m, bias, biasP, erro, MSE, TC)
  rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
  colnames(results) <- c("alpha", "phi", "theta", "sigma")
  print(c("Tamanho da Amostra:", n))
  print(round(results, 4))
  
  # Exibir avisos, se houver
  print(warnings())
}
