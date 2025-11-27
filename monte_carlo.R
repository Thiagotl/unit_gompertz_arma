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


alpha <- -0.14
phi   <-  0.94   # AR
theta <-  0.10   # MA
sigma <- 25      # sigma
tau   <-  0.5
true_values <- c(alpha, phi, theta, 24.59)
vn <- c(30, 50, 70, 100)   #c(70, 150, 300, 500, 1000)
R <- 10
z <- 1.96
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
    
   
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
    
    
    calpha <- cphi <- ctheta <- csigma <- 0
    
    bug         <- 0   # error + convergence fail 
    conv_fail   <- 0   # convergence fail (optim$conv != 0)
    error_count <- 0   # try-error 
    
    # Contadores de uso de gradiente
    grad_analytical <- 0   # convergência com gradiente analítico (escore)
    grad_numerical  <- 0   # convergência usando gradiente numérico
    
 
    pb <- txtProgressBar(min = 0, max = R, style = 3)
    
    for (i in 1:R) {
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
      fit1 <- try(
        uGoarma.fit(y, ma = ma1, ar = ar1),
        silent = TRUE
      )
      # error + convergence 
      if (inherits(fit1, "try-error")) {
        bug         <- bug + 1
        error_count <- error_count + 1
        
      } else if (is.null(fit1$conv) || fit1$conv != 0) {
        # convergence fail 
        bug       <- bug + 1
        conv_fail <- conv_fail + 1
        
      } else {
        # 
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
      }
      
      
      setTxtProgressBar(pb, i)
    } 
    
    close(pb)
    
    
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
    
    # Taxa de cobertura dos ICs (TC)
    TC <- c(calpha, cphi, ctheta, csigma) / R
    
    results <- rbind(m, bias, biasP, erro, MSE, TC)
    rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
    colnames(results) <- c("alpha", "phi", "theta", "sigma")
    
    cat("Tamanho da amostra:", n, "\n")
    print(round(results, 4))
    
    # Resumo das contagens de gradiente e falhas
    cat("\nResumo de convergência para n =", n, "\n")
    cat("Convergências com gradiente analítico :", grad_analytical, "\n")
    cat("Convergências com gradiente numérico  :", grad_numerical,  "\n")
    cat("Falhas de convergência (conv != 0)    :", conv_fail,       "\n")
    cat("Erros (try-error)                     :", error_count,     "\n")
    cat("Total de bugs (erro + falha)         :", bug,             "\n\n")
    
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
      R               = R,
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


#save(MC_out, file = "MC_UGoARMA_ARMA11_parallel.RData")