# Limpar o ambiente
rm(list = ls())

# Carregar funções necessárias
source("simu.ugoarma.R")
source("ugo_fit.R")

# Definições de parâmetros
phi = 0.2 # AR
theta = 0.4 # MA
alpha = 1
sigma = 6
tau = 0.5
true_values = c(1, 0.2, 0.4, 6)

# Tamanhos da amostra ajustados
vn = c(70,100,150,200,350,400) 
z = 1.96
R = 5000

set.seed(2024)

# Função para simulação
simulacao <- function(n) {
  tryCatch({
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
    calpha <- cphi <- ctheta <- csigma <- 0
    
    # Inicializar a barra de progresso
    pb <- txtProgressBar(min = 0, max = R, style = 3)
    
    for (i in 1:R) {
      y <- simu.ugoarma(n, phi = phi, theta = theta, alpha = alpha, sigma = sigma, tau = tau, freq = 12, link = "logit")
      fit1 <- uGoarma.fit(y)
      
      estim[i,] <- fit1$model[,1]
      err[i,] <- fit1$model[,2]
      
      ICi[i,] <- estim[i,] - (z * err[i,])
      ICs[i,] <- estim[i,] + (z * err[i,])
      
      if (ICi[i,1] <= alpha && ICs[i,1] >= alpha) calpha <- calpha + 1
      if (ICi[i,2] <= phi && ICs[i,2] >= phi) cphi <- cphi + 1
      if (ICi[i,3] <= theta && ICs[i,3] >= theta) ctheta <- ctheta + 1
      if (ICi[i,4] <= sigma && ICs[i,4] >= sigma) csigma <- csigma + 1
      
      # Atualizar a barra de progresso
      setTxtProgressBar(pb, i)
    }
    
    close(pb) # Fechar a barra de progresso
    
    m <- apply(estim, 2, mean)
    bias <- (true_values - m)
    biasP <- bias / true_values * 100
    erro <- apply(estim, 2, sd)
    MSE <- apply(estim, 2, var) + bias^2
    TC <- c(calpha, cphi, ctheta, csigma) / R
    
    results <- rbind(m, bias, biasP, erro, MSE, TC)
    rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
    
    list(tamanho = n, resultados = round(results, 4))
  }, error = function(e) {
    message(paste("Erro ao processar tamanho de amostra", n, ":", e$message))
    NULL
  })
}

# Executar simulações em paralelo
library(parallel)
num_cores <- 31
resultados <- mclapply(vn, simulacao, mc.cores = num_cores)

# Remover resultados NULL (caso alguma simulação tenha falhado)
resultados <- Filter(Negate(is.null), resultados)

# Exibir resultados
for (resultado in resultados) {
  print(paste("Tamanho da Amostra:", resultado$tamanho))
  print(resultado$resultados)
}
