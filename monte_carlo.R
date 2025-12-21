##############################################################
# MONTE CARLO SIMULATION - UGo-ARMA(1,1) 
#
# Created by Thiago Tavares Lopes (thiago.tavares@acad.ufsm.br)
# Reviewed by Renata Rojas Guerra (renata.r.guerra@ufsm.br)
# December/2025
#
##############################################################


set.seed(123)

library(parallel)

source("simu.ugoarma.R")
source("ugo_fit.R")


alpha <- 0.3   
phi   <- 0.9   # AR (put NA for ARMA(0,1))
theta <- 0.14  # MA (put NA for ARMA(1,0))
sigma <- 13    
tau   <- 0.5

vn <- c(70, 150, 300, 500, 1000)
R  <- 10000          
z  <- 1.96

# Ordem do modelo a ser ajustado
# ARMA(1,0): ar1 <- 1; ma1 <- NA
# ARMA(0,1): ar1 <- NA; ma1 <- 1
# ARMA(1,1): ar1 <- 1; ma1 <- 1
ar1 <- 1 
ma1 <- 1 

# ------------------------------------------------------------------
# Escolha automática de quais parâmetros entram no vetor verdadeiro
# conforme phi/theta e ar1/ma1 (podem ser NA)
# ------------------------------------------------------------------
include_phi   <- !is.na(phi)   && !is.na(ar1)
include_theta <- !is.na(theta) && !is.na(ma1)

if (include_phi && include_theta) {
  par_names   <- c("alpha", "phi", "theta", "sigma")
  true_values <- c(alpha, phi, theta, sigma)
} else if (include_phi && !include_theta) {
  par_names   <- c("alpha", "phi", "sigma")
  true_values <- c(alpha, phi, sigma)
} else if (!include_phi && include_theta) {
  par_names   <- c("alpha", "theta", "sigma")
  true_values <- c(alpha, theta, sigma)
} else {
  par_names   <- c("alpha", "sigma")
  true_values <- c(alpha, sigma)
}

n_cores <- max(1, detectCores() - 3)


MC_out <- list()

start_time <- Sys.time()

system.time({
  for (n in vn) {
    
    cat("\n============================\n")
    cat("Start simulation to n =", n, "\n")
    cat("Using", n_cores, "cores\n")
    cat("============================\n")
    
    n_par <- length(true_values)
    
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = n_par)
    
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
      
      if (inherits(fit1, "try-error")) {
        res$bug   <- 1
        res$error <- 1
        return(res)
      }
      
      if (is.null(fit1$conv) || fit1$conv != 0) {
        res$bug       <- 1
        res$conv_fail <- 1
        return(res)
      }
      
      if (!is.null(fit1$grad_used)) {
        if (fit1$grad_used == "analytical") res$grad_analytical <- 1
        if (fit1$grad_used == "numerical")  res$grad_numerical  <- 1
      }
      
      est_i <- fit1$model[, 1]
      err_i <- fit1$model[, 2]
      
      if (length(est_i) != n_par || length(err_i) != n_par) {
        res$bug   <- 1
        res$error <- 1
        return(res)
      }
      
      res$estim <- est_i
      res$err   <- err_i
      
      if (!any(is.na(est_i)) && !any(is.na(err_i))) {
        ICi_i <- est_i - (z * err_i)
        ICs_i <- est_i + (z * err_i)
        
        res$ICi <- ICi_i
        res$ICs <- ICs_i
        
        idx_alpha <- which(par_names == "alpha")
        if (length(idx_alpha) == 1 &&
            ICi_i[idx_alpha] <= alpha && ICs_i[idx_alpha] >= alpha) {
          res$calpha <- 1
        }
        
        idx_phi <- which(par_names == "phi")
        if (length(idx_phi) == 1 &&
            ICi_i[idx_phi] <= phi && ICs_i[idx_phi] >= phi) {
          res$cphi <- 1
        }
        
        idx_theta <- which(par_names == "theta")
        if (length(idx_theta) == 1 &&
            ICi_i[idx_theta] <= theta && ICs_i[idx_theta] >= theta) {
          res$ctheta <- 1
        }
        
        idx_sigma <- which(par_names == "sigma")
        if (length(idx_sigma) == 1 &&
            ICi_i[idx_sigma] <= sigma && ICs_i[idx_sigma] >= sigma) {
          res$csigma <- 1
        }
      }
      
      return(res)
    }
    
    while (success < R) {
      remaining <- R - success
      
      res_list <- mclapply(1:remaining, one_rep, mc.cores = n_cores)
      
      for (j in 1:remaining) {
        ri <- res_list[[j]]
        attempt <- attempt + 1
        
        bug            <- bug            + ri$bug
        conv_fail      <- conv_fail      + ri$conv_fail
        error_count    <- error_count    + ri$error
        grad_analytical <- grad_analytical + ri$grad_analytical
        grad_numerical  <- grad_numerical  + ri$grad_numerical
        
        if (ri$bug == 1 || ri$conv_fail == 1 || ri$error == 1) {
          next
        }
        
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
    
    
    m <- apply(estim, 2, mean, na.rm = TRUE)
    
    bias <- (true_values - m)
    
    biasP <- bias / true_values * 100
    
    erro <- apply(estim, 2, sd, na.rm = TRUE)
    
    MSE <- apply(estim, 2, var, na.rm = TRUE) + bias^2
    
    TC_vec <- numeric(n_par)
    names(TC_vec) <- par_names
    
    if ("alpha" %in% par_names) TC_vec[par_names == "alpha"] <- calpha / R
    if ("phi"   %in% par_names) TC_vec[par_names == "phi"]   <- cphi   / R
    if ("theta" %in% par_names) TC_vec[par_names == "theta"] <- ctheta / R
    if ("sigma" %in% par_names) TC_vec[par_names == "sigma"] <- csigma / R
    
    results <- rbind(m, bias, biasP, erro, MSE, TC_vec)
    rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
    colnames(results) <- par_names
    
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

execution_time <- end_time - start_time
cat("Tempo total de execução:",
    round(as.numeric(execution_time, units = "secs")),
    "segundos\n")

# save(MC_out, file = "MC_UGoARMA_ARMA11_MC_10k_success_parallel_2.RData")
