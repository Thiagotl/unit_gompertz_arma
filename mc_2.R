set.seed(2025)

source("simu.ugoarma.R")
source("ugo_fit.R")



alpha <- -0.09
phi   <- 0.96 #AR
theta <- NA # MA
sigma <- 25

beta  <- 0.03 #0.03
tau   <- 0.5


true_values <- c(alpha, phi, sigma, beta)
vn <- c(70, 150, 300, 500, 1000)
R  <- 100
z  <- 1.96
ar1 <- 1
ma1 <- NA

# Object to all results for "n" 
MC_out <- list()

start_time <- Sys.time()

system.time({
  for (n in vn) {
    
    cat("\n============================\n")
    cat("Start simulation to n =", n, "\n")
    cat("============================\n")
    
    
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
    
    calpha <- cphi <- csigma <- cbeta <- 0
    
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
    
    
    while (success < R) {
      
      attempt <- attempt + 1
      
      #X_cos <- matrix(cos(2 * pi * (1:n) / 12))
      # simula série
      y <- simu.ugoarma(
        n     = n,
        phi   = phi,
        theta = theta,
        alpha = alpha,
        sigma = sigma,
        tau   = tau,
        freq  = 12,
        link  = "logit",
        X     = 'cos',
        beta  = beta
      )
      
      X_cas = as.matrix(cos(2*pi*(1:n)/12))
      
      
      # tenta ajustar o modelo
      fit1 <- try(
        uGoarma.fit(y, ma = ma1, ar = ar1,  X = X_cos),
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
      
     
      estim[i, ] <- fit1$model[, 1]
      err[i, ]   <- fit1$model[, 2]
      
      if (!any(is.na(estim[i, ])) && !any(is.na(err[i, ]))) {
       
        ICi[i, ] <- estim[i, ] - (z * err[i, ])
        ICs[i, ] <- estim[i, ] + (z * err[i, ])
        
       
        if (ICi[i, 1] <= alpha && ICs[i, 1] >= alpha) calpha <- calpha + 1
        if (ICi[i, 2] <= phi   && ICs[i, 2] >= phi)   cphi   <- cphi   + 1
        #if (ICi[i, 3] <= theta && ICs[i, 3] >= theta) ctheta <- ctheta + 1
        if (ICi[i, 3] <= sigma && ICs[i, 3] >= sigma) csigma <- csigma + 1
        if (ICi[i, 4] <= beta  && ICs[i, 4] >= beta) cbeta <- cbeta + 1
      }
      
     
      setTxtProgressBar(pb, success)
    } 
    
    close(pb)
    
    cat("Total de tentativas (sucessos + bugs) para n =", n, ":", attempt, "\n")
    cat("Sucessos (convergências)                   :", success, "\n")
    cat("Bugs (falha + erro)                        :", bug, "\n\n")
    
   
    
    m <- apply(estim, 2, mean, na.rm = TRUE)
    
    bias <- (true_values - m)
    
    biasP <- bias / true_values * 100
    
    erro <- apply(estim, 2, sd, na.rm = TRUE)
    
    MSE <- apply(estim, 2, var, na.rm = TRUE) + bias^2
    
    TC <- c(calpha, cphi, csigma, cbeta) / R
    
    results <- rbind(m, bias, biasP, erro, MSE, TC)
    rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
    colnames(results) <- c("alpha", "phi", "sigma","beta")
    
    cat("Tamanho da amostra:", n, "\n")
    print(round(results, 4))
    
    cat("\nResumo de convergência para n =", n, "\n")
    cat("Total de tentativas (sucesso + bug)       :", attempt,        "\n")
    cat("Sucessos (convergências)                  :", success,        "\n")
    cat("Convergências com gradiente analítico     :", grad_analytical, "\n")
    cat("Convergências com gradiente numérico      :", grad_numerical,  "\n")
    cat("Falhas de convergência (conv != 0)        :", conv_fail,      "\n")
    cat("Erros (try-error)                         :", error_count,    "\n")
    cat("Total de bugs (erro + falha)              :", bug,            "\n\n")
    
    #print(warnings())
    
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