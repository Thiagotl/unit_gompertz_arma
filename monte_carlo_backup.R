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

# rm(list = ls())
set.seed(2025)

library(parallel)  # for mclapply()

source("simu.ugoarma.R")
source("ugo_fit.R")

alpha <- -0.14
phi   <-  0.94   # AR
theta <-  0.10   # MA
sigma <- 25      # sigma
tau   <-  0.5

true_values <- c(alpha, phi, theta, 24.59)

vn <- c(30, 50, 70, 100)   # c(70, 150, 300, 500, 1000)
R  <- 5000
z  <- 1.96
ar1 <- 1
ma1 <- 1

# Number of cores for parallel execution
n_cores <- max(1, detectCores() - 2)

# Object to store all results for each n
MC_out <- list()

start_time <- Sys.time()

system.time({
  for (n in vn) {
    
    cat("\n============================\n")
    cat("Start simulation to n =", n, "\n")
    cat("Using", n_cores, "cores\n")
    cat("============================\n")
    
    # Run R replications in parallel
    # Each element of res_list is a list with:
    #   status     = "ok", "error", or "conv_fail"
    #   grad_used  = "analytical", "numerical", or NA
    #   estim      = numeric vector (alpha, phi, theta, sigma) or NA
    #   err        = numeric vector of standard errors or NA
    res_list <- mclapply(
      X = 1:R,
      FUN = function(i) {
        # 1) Simulate series
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
        
        # 2) Fit model
        fit1 <- try(
          uGoarma.fit(y, ma = ma1, ar = ar1),
          silent = TRUE
        )
        
        # 3) Handle errors and convergence status
        if (inherits(fit1, "try-error")) {
          # Execution error inside uGoarma.fit
          return(list(
            status    = "error",
            grad_used = NA_character_,
            estim     = rep(NA_real_, length(true_values)),
            err       = rep(NA_real_, length(true_values))
          ))
        }
        
        if (is.null(fit1$conv) || fit1$conv != 0) {
          # Convergence failure
          gu <- if (!is.null(fit1$grad_used)) fit1$grad_used else NA_character_
          return(list(
            status    = "conv_fail",
            grad_used = gu,
            estim     = rep(NA_real_, length(true_values)),
            err       = rep(NA_real_, length(true_values))
          ))
        }
        
        # Success: converged (conv == 0)
        gu <- if (!is.null(fit1$grad_used)) fit1$grad_used else NA_character_
        
        return(list(
          status    = "ok",
          grad_used = gu,
          estim     = fit1$model[, 1],
          err       = fit1$model[, 2]
        ))
      },
      mc.cores = n_cores
    )
    
    # ========================
    # Combine results
    # ========================
    status_vec <- vapply(res_list, `[[`, character(1), "status")
    grad_vec   <- vapply(res_list, `[[`, character(1), "grad_used")
    
    estim_mat  <- t(vapply(
      res_list,
      FUN = `[[`,
      FUN.VALUE = numeric(length(true_values)),
      "estim"
    ))
    
    err_mat    <- t(vapply(
      res_list,
      FUN = `[[`,
      FUN.VALUE = numeric(length(true_values)),
      "err"
    ))
    
    estim <- estim_mat
    err   <- err_mat
    
    # Counters for errors and convergence
    error_count <- sum(status_vec == "error")
    conv_fail   <- sum(status_vec == "conv_fail")
    bug         <- sum(status_vec != "ok")  # error + convergence fail
    
    # Gradient usage
    grad_analytical <- sum(status_vec == "ok" & grad_vec == "analytical")
    grad_numerical  <- sum(status_vec == "ok" & grad_vec == "numerical")
    
    # ========================
    # Performance measures
    # ========================
    
    # Mean of estimates
    m <- apply(estim, 2, mean, na.rm = TRUE)
    
    # Absolute bias
    bias <- (true_values - m)
    
    # Relative percentage bias (RB%)
    biasP <- bias / true_values * 100
    
    # Standard deviation of estimates
    erro <- apply(estim, 2, sd, na.rm = TRUE)
    
    # MSE = Var(est) + bias^2
    MSE <- apply(estim, 2, var, na.rm = TRUE) + bias^2
    
    # Confidence intervals and coverage
    ICi <- estim - (z * err)
    ICs <- estim + (z * err)
    
    # Valid rows (no NA in estim or err)
    ok_row <- !apply(is.na(estim) | is.na(err), 1, any)
    
    calpha <- sum(ok_row & ICi[, 1] <= alpha & ICs[, 1] >= alpha)
    cphi   <- sum(ok_row & ICi[, 2] <= phi   & ICs[, 2] >= phi)
    ctheta <- sum(ok_row & ICi[, 3] <= theta & ICs[, 3] >= theta)
    csigma <- sum(ok_row & ICi[, 4] <= sigma & ICs[, 4] >= sigma)
    
    # Coverage rate (TC)
    TC <- c(calpha, cphi, ctheta, csigma) / R
    
    results <- rbind(m, bias, biasP, erro, MSE, TC)
    rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
    colnames(results) <- c("alpha", "phi", "theta", "sigma")
    
    cat("Tamanho da amostra:", n, "\n")
    print(round(results, 4))
    
    # Convergence summary
    cat("\nResumo de convergência para n =", n, "\n")
    cat("Convergências com gradiente analítico :", grad_analytical, "\n")
    cat("Convergências com gradiente numérico  :", grad_numerical,  "\n")
    cat("Falhas de convergência (conv != 0)    :", conv_fail,       "\n")
    cat("Erros (try-error)                     :", error_count,     "\n")
    cat("Total de bugs (erro + falha)         :", bug,             "\n\n")
    
    #print(warnings())
    
    # Store results for this n in MC_out
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
    
  }  # end for n
})   # end system.time

end_time <- Sys.time()

# Total execution time
execution_time <- end_time - start_time
cat("Tempo total de execução:",
    round(as.numeric(execution_time, units = "secs")),
    "segundos\n")

# (Opcional) salvar tudo em um RData:
# save(MC_out, file = "MC_UGoARMA_parallel_results.RData")



###########################



##############################################################
# MONTE CARLO SIMULATION - UGo-ARMA(1,1)
#
# Created by Thiago Tavares Lopes (thiago.tavares@acad.ufsm.br)
# Reviewed by Renata Rojas Guerra (renata.r.guerra@ufsm.br)
# December/2025
#
##############################################################

# rm(list = ls())
set.seed(2025)

library(parallel)

source("simu.ugoarma.R")
source("ugo_fit.R")

alpha <- -0.14
phi   <-  0.94   # AR
theta <-  NA #0.10   # MA
sigma <- 25      # sigma
tau   <-  0.5

true_values <- c(alpha, phi, sigma)

vn <- c(30, 50, 70, 100)   # c(70, 150, 300, 500, 1000)
R  <- 100
z  <- 1.96
ar1 <- NA
ma1 <- 1

# número de cores
n_cores <- max(1, detectCores() - 1)

# Objeto para guardar todos os resultados
MC_out <- list()

start_time <- Sys.time()

system.time({
  for (n in vn) {
    
    cat("\n============================\n")
    cat("Start simulation to n =", n, "\n")
    cat("Using", n_cores, "cores\n")
    cat("============================\n")
    
    # matrizes principais
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
    
    # contadores de cobertura
    calpha <- cphi <-  csigma <- 0
    
    # contadores de problemas
    bug         <- 0   # error + convergence fail 
    conv_fail   <- 0   # convergence fail (optim$conv != 0)
    error_count <- 0   # try-error 
    
    # contadores de uso de gradiente
    grad_analytical <- 0   # convergência com gradiente analítico (escore)
    grad_numerical  <- 0   # convergência usando gradiente numérico
    
    # -----------------------------------------------------------------
    # Função para UMA réplica (vai rodar em paralelo via mclapply)
    # -----------------------------------------------------------------
    one_rep <- function(i) {
      # estrutura padrão com tudo inicializado
      res <- list(
        bug            = 0,
        conv_fail      = 0,
        error          = 0,
        grad_analytical = 0,
        grad_numerical  = 0,
        estim          = rep(NA_real_, length(true_values)),
        err            = rep(NA_real_, length(true_values)),
        ICi            = rep(NA_real_, length(true_values)),
        ICs            = rep(NA_real_, length(true_values)),
        calpha         = 0,
        cphi           = 0,
        #ctheta         = 0,
        csigma         = 0
      )
      
      # gerar série
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
      
      # ajustar modelo
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
      
      # sucesso: convergiu
      if (!is.null(fit1$grad_used)) {
        if (fit1$grad_used == "analytical") res$grad_analytical <- 1
        if (fit1$grad_used == "numerical")  res$grad_numerical  <- 1
      }
      
      # estimativas e EP
      est_i <- fit1$model[, 1]
      err_i <- fit1$model[, 2]
      
      res$estim <- est_i
      res$err   <- err_i
      
      # se não tem NA, calcula IC e cobertura
      if (!any(is.na(est_i)) && !any(is.na(err_i))) {
        ICi_i <- est_i - (z * err_i)
        ICs_i <- est_i + (z * err_i)
        
        res$ICi <- ICi_i
        res$ICs <- ICs_i
        
        if (ICi_i[1] <= alpha && ICs_i[1] >= alpha) res$calpha <- 1
        if (ICi_i[2] <= phi   && ICs_i[2] >= phi)   res$cphi   <- 1
        if (ICi_i[3] <= theta && ICs_i[3] >= theta) res$ctheta <- 1
        if (ICi_i[4] <= sigma && ICs_i[4] >= sigma) res$csigma <- 1
      }
      
      return(res)
    }
    
    # -----------------------------------------------------------------
    # Aqui entra o paralelismo: roda as R réplicas em paralelo
    # -----------------------------------------------------------------
    res_list <- mclapply(1:R, one_rep, mc.cores = n_cores)
    
    # -----------------------------------------------------------------
    # Agora agregamos os resultados na MESMA estrutura que você já tinha
    # -----------------------------------------------------------------
    for (i in 1:R) {
      ri <- res_list[[i]]
      
      # contadores globais
      bug         <- bug         + ri$bug
      conv_fail   <- conv_fail   + ri$conv_fail
      error_count <- error_count + ri$error
      grad_analytical <- grad_analytical + ri$grad_analytical
      grad_numerical  <- grad_numerical  + ri$grad_numerical
      
      calpha <- calpha + ri$calpha
      cphi   <- cphi   + ri$cphi
      ctheta <- ctheta + ri$ctheta
      csigma <- csigma + ri$csigma
      
      # matrizes
      estim[i, ] <- ri$estim
      err[i, ]   <- ri$err
      ICi[i, ]   <- ri$ICi
      ICs[i, ]   <- ri$ICs
    }
    
    # -----------------------------------------------------------------
    # Daqui pra baixo é exatamente a mesma lógica de antes
    # -----------------------------------------------------------------
    
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
    
    #print(warnings())
    
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
    
  }  # fim do for (n in vn)
})   # fim do system.time

end_time <- Sys.time()

# Tempo total de execução
execution_time <- end_time - start_time
cat("Tempo total de execução:",
    round(as.numeric(execution_time, units = "secs")),
    "segundos\n")

# Se quiser salvar:
# save(MC_out, file = "MC_UGoARMA_ARMA11_parallel.RData")

#########


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