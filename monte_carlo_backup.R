#MONTE CARLO SIMULATION

# Created by Thiago Tavares Lopes (thiago.tavares@acad.ufsm.br), 
# Reviewed by Renata Rojas Guerra (renata.r.guerra@ufsm.br), december/2025

# simu - ARMA(1,1) - 

source("simu.ugoarma.R")
source("ugo_fit.R")

alpha = -0.14
phi =  0.94     #0.2 #AR
theta = 0.10  # 0.4 #MA
sigma = 25   #6
tau =  0.5   #0.5
true_values = c(-0.14, 0.94, 0.10, 24.59 ) # alpha, phi, theta, sigma / phi= 0.2
vn = c(70,150, 300, 500, 1000) # 70,150, 300, 500
R = 10
z = 1.96

ar1=1
ma1=1

MC_out <- list()

start_time <- Sys.time()

system.time({
  for (n in vn) {
  # result matrix 
  estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
  

  calpha <- cphi <- ctheta <- csigma <- 0 # para guardar os resultados cphi 
  i <- 0
  bug <- 0 # error + convergence fail
  
  grad_analytical <- 0
  grad_numerical  <- 0
  conv_fail       <- 0  # falhas de convergência (optim$convergence != 0)
  error_count     <- 0  # try-error
  
  # Inicializa a barra de progresso
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  for (i in 1:R) {
    y <- simu.ugoarma(n, phi = phi, theta = theta, alpha = alpha, sigma = sigma, tau = tau, freq = 12, link = "logit")
    fit1 <- try(uGoarma.fit(y, ma=ma1, ar=ar1), silent = TRUE) #, ma=ma1, ar=ar1
    
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
    # Atualiza a barra de progresso
    setTxtProgressBar(pb, i)
  }
  
  # Fecha a barra de progresso
  close(pb)
  
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
  TC <- c(calpha, cphi, ctheta, csigma) / R #, cphi
  ## final results
  results <- rbind(m, bias, biasP, erro, MSE, TC)
  rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
  colnames(results) <- c("alpha","phi", "theta", "sigma")
  print(c("Tamanho da Amostra:", n))
  print(round(results, 4))
  
  # Exibir avisos, se houver
  print(warnings())
}
  
})  


end_time <- Sys.time()

# Calcula e imprime o tempo total de execução
execution_time <- end_time - start_time
print(paste("Total time execution:", round(as.numeric(execution_time, units = "secs")), "segundos"))

###########44
###########44
##############################################################
# MONTE CARLO SIMULATION - UGo-ARMA(1,1)
#
# Created by Thiago Tavares Lopes (thiago.tavares@acad.ufsm.br)
# Reviewed by Renata Rojas Guerra (renata.r.guerra@ufsm.br)
# December/2025
#
##############################################################

set.seed(2025)

source("simu.ugoarma.R")
source("ugo_fit.R")

alpha <- 1
phi   <-  0.2   # AR
theta <-  0.4   # MA
sigma <- 0.6    # sigma
tau   <-  0.5

true_values <- c(alpha, phi, theta, sigma)
vn <- c(70, 150, 300, 500, 1000)
R  <- 100
z  <- 1.96
ar1 <- 1
ma1 <- 1

# Object to all results for "n" 
MC_out <- list()

start_time <- Sys.time()

system.time({
  for (n in vn) {
    
    cat("\n============================\n")
    cat("Start simulation to n =", n, "\n")
    cat("============================\n")
    
    # matrizes com R linhas, que serão preenchidas apenas com réplicas bem-sucedidas
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
    
    # contadores de cobertura
    calpha <- cphi <- ctheta <- csigma <- 0
    
    # contadores de problemas
    bug         <- 0   # error + convergence fail 
    conv_fail   <- 0   # convergence fail (optim$conv != 0)
    error_count <- 0   # try-error 
    
    # contadores de uso de gradiente
    grad_analytical <- 0   # convergência com gradiente analítico (escore)
    grad_numerical  <- 0   # convergência usando gradiente numérico
    
    # contadores de réplicas
    success <- 0          # quantas réplicas CONVERGIRAM (vai até R)
    attempt <- 0          # quantas tentativas no total (convergidas + bugs)
    
    pb <- txtProgressBar(min = 0, max = R, style = 3)
    
    # ---- LOOP PRINCIPAL: CONTINUA ATÉ TER R SUCESSOS ----
    while (success < R) {
      
      attempt <- attempt + 1
      
      # simula série
      y <- simu.ugoarma(
        n     = n,
        phi   = phi,
        theta = theta,
        alpha = alpha,
        sigma = sigma,
        tau   = tau,
        freq  = 12,
        link  = "logit"
      )
      
      # tenta ajustar o modelo
      fit1 <- try(
        uGoarma.fit(y, ma = ma1, ar = ar1),
        silent = TRUE
      )
      
      # error + convergence 
      if (inherits(fit1, "try-error")) {
        bug         <- bug + 1
        error_count <- error_count + 1
        next  # não conta como sucesso, volta pro while
      } 
      
      if (is.null(fit1$conv) || fit1$conv != 0) {
        # convergence fail 
        bug       <- bug + 1
        conv_fail <- conv_fail + 1
        next  # também não conta como sucesso
      }
      
      # Se chegou aqui: CONVERGIU (conv == 0)
      success <- success + 1
      i <- success           # índice da réplica bem-sucedida (1..R)
      
      # qual gradiente foi usado?
      if (!is.null(fit1$grad_used)) {
        if (fit1$grad_used == "analytical") grad_analytical <- grad_analytical + 1
        if (fit1$grad_used == "numerical")  grad_numerical  <- grad_numerical  + 1
      }
      
      # Estimativas e erros-padrão dos parâmetros (alpha, phi, theta, sigma)
      estim[i, ] <- fit1$model[, 1]
      err[i, ]   <- fit1$model[, 2]
      
      if (!any(is.na(estim[i, ])) && !any(is.na(err[i, ]))) {
        # Intervalos de confiança de 95%
        ICi[i, ] <- estim[i, ] - (z * err[i, ])
        ICs[i, ] <- estim[i, ] + (z * err[i, ])
        
        # Cobertura para cada parâmetro
        if (ICi[i, 1] <= alpha && ICs[i, 1] >= alpha) calpha <- calpha + 1
        if (ICi[i, 2] <= phi   && ICs[i, 2] >= phi)   cphi   <- cphi   + 1
        if (ICi[i, 3] <= theta && ICs[i, 3] >= theta) ctheta <- ctheta + 1
        if (ICi[i, 4] <= sigma && ICs[i, 4] >= sigma) csigma <- csigma + 1
      }
      
      # atualiza barra de progresso com número de SUCESSOS
      setTxtProgressBar(pb, success)
    }  # fim do while(success < R)
    
    close(pb)
    
    cat("Total de tentativas (sucessos + bugs) para n =", n, ":", attempt, "\n")
    cat("Sucessos (convergências)                   :", success, "\n")
    cat("Bugs (falha + erro)                        :", bug, "\n\n")
    
    # === A PARTIR DAQUI, SUCCESS = R (temos R réplicas convergidas) ===
    
    # Média das estimativas
    m <- apply(estim, 2, mean, na.rm = TRUE)
    
    # Viés absoluto
    bias <- (true_values - m)
    
    # Viés relativo percentual (RB%)
    biasP <- bias / true_values * 100
    
    # Desvio-padrão das estimativas
    erro <- apply(estim, 2, sd, na.rm = TRUE)
    
    # MSE = Var(est) + bias^2
    MSE <- apply(estim, 2, var, na.rm = TRUE) + bias^2
    
    # Taxa de cobertura dos ICs (TC) — AGORA SOBRE R SUCESSOS
    TC <- c(calpha, cphi, ctheta, csigma) / R
    
    results <- rbind(m, bias, biasP, erro, MSE, TC)
    rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
    colnames(results) <- c("alpha", "phi", "theta", "sigma")
    
    cat("Tamanho da amostra:", n, "\n")
    print(round(results, 4))
    
    # Resumo das contagens de gradiente e falhas
    cat("\nResumo de convergência para n =", n, "\n")
    cat("Total de tentativas (sucesso + bug)       :", attempt,        "\n")
    cat("Sucessos (convergências)                  :", success,        "\n")
    cat("Convergências com gradiente analítico     :", grad_analytical, "\n")
    cat("Convergências com gradiente numérico      :", grad_numerical,  "\n")
    cat("Falhas de convergência (conv != 0)        :", conv_fail,      "\n")
    cat("Erros (try-error)                         :", error_count,    "\n")
    cat("Total de bugs (erro + falha)              :", bug,            "\n\n")
    
    print(warnings())
    
    # Result list  n ,  MC_out
    MC_out[[as.character(n)]] <- list(
      n               = n,
      true_values     = true_values,
      estim           = estim,
      err             = err,
      ICi             = ICi,
      ICs             = ICs,
      results         = results,
      R               = R,          # número de réplicas BEM SUCEDIDAS
      attempts        = attempt,    # total de tentativas
      bug_total       = bug,
      conv_fail       = conv_fail,
      error_count     = error_count,
      grad_analytical = grad_analytical,
      grad_numerical  = grad_numerical
    )
    
  } 
})

end_time <- Sys.time()

# Tempo total de execução
execution_time <- end_time - start_time
cat("Tempo total de execução:",
    round(as.numeric(execution_time, units = "secs")),
    "segundos\n")

# save(MC_out, file = "MC_UGoARMA_ARMA11_MC_10k_success.RData")

######### nao apagar

source("simu.ugoarma.R")
source("ugo_fit.R")
alpha = 1
phi = 0.2 #AR
theta = NA #MA
sigma = 6
tau = 0.5
true_values = c(1, 0.2, 6) # alpha, phi, theta, sigma 
vn = c(70,150, 300, 500, 1000) # 

R = 10000
z = 1.96

ar1=1
ma1=NA
for (n in vn) {
  # matriz de resultados
  estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
  # contadores
  calpha <-  cphi <- csigma <- 0 # para guardar os resultados cphi <-
  i <- 0
  bug <- 0 # inicializa o contador de bugs
  
  # Inicializa a barra de progresso
  pb <- txtProgressBar(min = 0, max = R, style = 3)
  
  for (i in 1:R) {
    y <- simu.ugoarma(n, phi = phi, theta = theta, alpha = alpha, sigma = sigma, tau = tau, freq = 12, link = "logit")
    fit1 <- try(uGoarma.fit(y, ma=ma1, ar=ar1), silent = TRUE) #, ma=ma1, ar=ar1
    
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
        
        # if (ICi[i, 2] <= theta && ICs[i, 2] >= theta) {
        #   ctheta <- ctheta + 1
        # }
        # 
        if (ICi[i, 3] <= sigma && ICs[i, 3] >= sigma) {
          csigma <- csigma + 1
        }
      }
    }
    # Atualiza a barra de progresso
    setTxtProgressBar(pb, i)
  }
  
  # Fecha a barra de progresso
  close(pb)
  
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
  TC <- c(calpha, cphi, csigma) / R #, cphi
  ## final results
  results <- rbind(m, bias, biasP, erro, MSE, TC)
  rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
  colnames(results) <- c("alpha", "phi", "sigma") #, "phi"
  print(c("Tamanho da Amostra:", n))
  print(round(results, 4))
  
  # Exibir avisos, se houver
  #print(warnings())
}





