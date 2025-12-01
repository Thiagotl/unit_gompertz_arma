##############################################################
# MONTE CARLO SIMULATION - UGo-ARMA(1,1) - paralelo
#
# Created by Thiago Tavares Lopes (thiago.tavares@acad.ufsm.br)
# Reviewed by Renata Rojas Guerra (renata.r.guerra@ufsm.br)
# December/2025
#
##############################################################

# rm(list = ls())
set.seed(123)

library(parallel)

source("simu.ugoarma.R")
source("ugo_fit.R")

alpha <- 0.3 #1
phi   <- 0.9 #0.2   # AR
theta <- 0.14 #0.4   # MA
sigma <- 13 #0.6   # sigma
tau   <- 0.5

true_values <- c(alpha, phi,  theta, sigma)
vn <- c(70, 150, 300, 500, 1000)
R  <- 10000          
z  <- 1.96
ar1 <- 1 # p phi
ma1 <- 1 # q theta


n_cores <- max(1, detectCores() - 3)

# Object to all results for "n" 
MC_out <- list()

start_time <- Sys.time()

system.time({
  for (n in vn) {
    
    cat("\n============================\n")
    cat("Start simulation to n =", n, "\n")
    cat("Using", n_cores, "cores\n")
    cat("============================\n")
    
    n_par <- length(true_values)
    
    # matrizes com R linhas, que serão preenchidas apenas com réplicas bem-sucedidas
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = n_par)
    
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
    
    # -----------------------------------------------------------------
    # Função de UMA tentativa de Monte Carlo (roda em paralelo)
    # -----------------------------------------------------------------
    one_rep <- function(dummy) {
      res <- list(
        bug             = 0,
        conv_fail       = 0,
        error           = 0,
        grad_analytical = 0,
        grad_numerical  = 0,
        estim           = rep(NA_real_, n_par),
        err             = rep(NA_real_, n_par),
        ICi             = rep(NA_real_, n_par),
        ICs             = rep(NA_real_, n_par),
        calpha          = 0,
        cphi            = 0,
        ctheta          = 0,
        csigma          = 0
      )
      
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
      
      # erro de execução
      if (inherits(fit1, "try-error")) {
        res$bug   <- 1
        res$error <- 1
        return(res)
      }
      
      # falha de convergência
      if (is.null(fit1$conv) || fit1$conv != 0) {
        res$bug       <- 1
        res$conv_fail <- 1
        return(res)
      }
      
      # Se chegou aqui: CONVERGIU (conv == 0)
      if (!is.null(fit1$grad_used)) {
        if (fit1$grad_used == "analytical") res$grad_analytical <- 1
        if (fit1$grad_used == "numerical")  res$grad_numerical  <- 1
      }
      
      # Estimativas e erros-padrão dos parâmetros (alpha, phi, theta, sigma)
      est_i <- fit1$model[, 1]
      err_i <- fit1$model[, 2]
      
      # segurança: se o comprimento não bater, marca como bug
      if (length(est_i) != n_par || length(err_i) != n_par) {
        res$bug   <- 1
        res$error <- 1
        return(res)
      }
      
      res$estim <- est_i
      res$err   <- err_i
      
      if (!any(is.na(est_i)) && !any(is.na(err_i))) {
        # Intervalos de confiança de 95%
        ICi_i <- est_i - (z * err_i)
        ICs_i <- est_i + (z * err_i)
        
        res$ICi <- ICi_i
        res$ICs <- ICs_i
        
        # Cobertura para cada parâmetro
        if (ICi_i[1] <= alpha && ICs_i[1] >= alpha) res$calpha <- 1
        if (ICi_i[2] <= phi   && ICs_i[2] >= phi)   res$cphi   <- 1
        if (ICi_i[3] <= theta && ICs_i[3] >= theta) res$ctheta <- 1
        if (ICi_i[4] <= sigma && ICs_i[4] >= sigma) res$csigma <- 1
      }
      
      return(res)
    }
    
    # -----------------------------------------------------------------
    # LOOP PRINCIPAL EM BLOCOS PARA GARANTIR R SUCESSOS (em paralelo)
    # -----------------------------------------------------------------
    while (success < R) {
      remaining <- R - success
      
      # roda um bloco paralelo de "remaining" tentativas
      res_list <- mclapply(1:remaining, one_rep, mc.cores = n_cores)
      
      for (j in 1:remaining) {
        ri <- res_list[[j]]
        attempt <- attempt + 1
        
        # contadores globais de problemas e gradiente
        bug         <- bug         + ri$bug
        conv_fail   <- conv_fail   + ri$conv_fail
        error_count <- error_count + ri$error
        grad_analytical <- grad_analytical + ri$grad_analytical
        grad_numerical  <- grad_numerical  + ri$grad_numerical
        
        # se houve bug/falha/erro, não conta como sucesso
        if (ri$bug == 1 || ri$conv_fail == 1 || ri$error == 1) {
          next
        }
        
        # SUCESSO
        success <- success + 1
        i <- success
        
        estim[i, ] <- ri$estim
        err[i, ]   <- ri$err
        ICi[i, ]   <- ri$ICi
        ICs[i, ]   <- ri$ICs
        
        calpha <- calpha + ri$calpha
        cphi   <- cphi   + ri$cphi
        ctheta <- ctheta + ri$ctheta
        csigma <- csigma + ri$csigma
        
        setTxtProgressBar(pb, success)
        
        if (success >= R) break
      }
    }
    
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
    
    # Taxa de cobertura dos ICs (TC) — SOBRE R SUCESSOS
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

#save(MC_out, file = "MC_UGoARMA_ARMA11_MC_10k_success_parallel_2.RData")
