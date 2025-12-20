best_ugo_2 <- function(y, h=6, pmax=6, qmax=6, nbest=8,
                     tau=0.5, link="logit", X=NA, X_hat=NA) {
  
  if (!exists("uGoarma.fit")) {
    source("ugo_fit.R")
  }
  
  # Inicializa o critério AIC
  fit <- uGoarma.fit(y, ma=1, diag=0, link=link)
  #fit <- uGoarma.fit(y, ma=1, diag=0, link=link, h=h, tau=tau, X=X, X_hat=X_hat)
  aicmin <- fit$aic
  
  cat("Initial AIC:", aicmin, "\n")
  
  #
  melhores <- data.frame(
    p = integer(nbest),
    q = integer(nbest),
    AIC = rep(Inf, nbest),
    stringsAsFactors = FALSE
  )
  
  tot <- 0  
  bug <- 0  
  best_model_aic <- NULL
  
  for (p in 0:pmax) { 
    for (q in 0:qmax) { 
      ar1 <- if (p == 0) NA else 1:p
      ma1 <- if (q == 0) NA else 1:q
      
      if (!all(is.na(c(ar1, ma1)))) {   
        cat("Testing p =", p, "q =", q, "\n")
        
        fituGo <- uGoarma.fit(y, ar=ar1, ma=ma1, tau=tau, link=link, h=h, X=X, X_hat=X_hat, diag=1)
        tot <- tot + 1
        
        if (fituGo$conv != 0) {  
          cat("No convergence for p =", p, "q =", q, "\n")
          bug <- bug + 1
          next 
        }          
        
        if (fituGo$aic < aicmin) {  
          aicmin <- fituGo$aic
          best_model_aic <- fituGo$model
          cat("New best AIC:", aicmin, "\n")
        }
        
        if (fituGo$aic < max(melhores$AIC)) {
          max_idx <- which.max(melhores$AIC)
          melhores[max_idx, ] <- c(p, q, fituGo$aic)
        }
      }
    }
  }
  
  
  melhores <- melhores[order(melhores$AIC),]
  
  cat("\nSelected model from AIC:\n")
  print(best_model_aic)
  cat("\nTotal tested models:", tot, "\n")
  cat("Total errors in estimation:", bug, "\n")
  cat("\nBest models:\n")
  print(melhores)
  
  
  return(melhores)
}
