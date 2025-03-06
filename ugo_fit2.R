#rm(list = ls())

uGoarma.fit <- function(y, ar = 1, ma = 1, tau = .5, link = "logit", h = 1, 
                        diag = 0, X = NA, X_hat = NA)
{
  # Se necessário, ajuste o caminho:
  # source("ugo-functions.R")  
  
  #----------------------------------------------------------------------------
  # 0) Preparação básica
  #----------------------------------------------------------------------------
  
  z <- list()
  maxit1 <- 50
  
  p  <- max(ar)
  q  <- max(ma)
  n  <- length(y)
  m  <- max(p, q, na.rm = TRUE) 
  p1 <- length(ar)
  q1 <- length(ma)
  
  # vetores auxiliares
  error <- rep(0, n)
  eta   <- rep(NA, n)
  
  # Para previsões
  y_prev <- rep(NA, (n + h))
  
  #----------------------------------------------------------------------------
  # 1) Função de ligação
  #----------------------------------------------------------------------------
  
  linktemp <- substitute(link)
  
  if(!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
    if (linktemp == "link"){
      linktemp <- eval(link)
    }
  }
  
  valid_links <- c("logit", "probit", "cloglog")
  if (linktemp %in% valid_links) {
    stats <- make.link(linktemp)
  } else {
    stop(paste(linktemp, "link not available. Possible links:", 
               paste(valid_links, collapse = ", ")))
  }
  
  linkfun  <- stats$linkfun
  linkinv  <- stats$linkinv 
  mu.eta   <- stats$mu.eta
  diflink  <- function(t) 1/(stats$mu.eta(stats$linkfun(t)))
  
  ynew <- linkfun(y) # g(y)
  
  # matriz para AR (transformada via link)
  ynew_ar <- suppressWarnings(matrix(ynew, (n-1), max(p,1,na.rm=TRUE)))
  
  #----------------------------------------------------------------------------
  # 2) Definição de ar, ma e covariáveis
  #----------------------------------------------------------------------------
  
  # AR
  if(!any(is.na(ar))) {
    names_phi <- paste0("phi", ar)
    Z <- suppressWarnings(matrix(ynew, (n-1), p1)[m:(n-1), ])
  } else {
    ar <- p1 <- 0
    Z <- NA
  }
  
  # MA
  if(!any(is.na(ma))) {
    names_theta <- paste0("theta", ma)
  } else {
    ma <- q1 <- 0
  }
  
  # Covariáveis
  if(!any(is.na(X))) {
    X <- as.matrix(X)
    names_beta <- paste0("beta", 1:ncol(X))
    Xm <- X[(m+1):n, ]
    k  <- ncol(X)
  } else {
    k <- 0 
    X <- matrix(0, nrow=n, ncol=0)  # sem colunas
    Xm <- NA
  }
  
  #----------------------------------------------------------------------------
  # 3) Ajuste das recursões
  #----------------------------------------------------------------------------
  
  q_1 <- max(q1, 1)
  
  # Em vez de max(k,1), usamos k mesmo (e p1, q1 etc.)
  deta.dalpha  <- rep(0, n)
  deta.dbeta   <- matrix(0, nrow=n, ncol=k)
  deta.dphi    <- matrix(0, nrow=n, ncol=p1)
  deta.dtheta  <- matrix(0, nrow=n, ncol=q1)
  
  # Vetor para o termo MA
  R <- matrix(NA, nrow=(n - m), ncol=q_1)
  k_i <- q1 / q_1  # pode ser 1 se q1=1, ou 0 se q1=0
  
  #----------------------------------------------------------------------------
  # 4) Estimação via log-verossimilhança
  #----------------------------------------------------------------------------
  
  # a) Valores iniciais (OLS)
  
  # Montamos Xstart apenas se temos X e AR
  Xstart <- cbind(rep(1, (n - m)), Xm, Z)  # intercepto + covariáveis + AR
  # remover NAs por linha
  Xstart <- t(apply(Xstart, 1, function(xx) xx[!is.na(xx)]))
  
  ols <- lm.fit(Xstart, ynew[(m+1):n])$coef
  
  # Número total de parâmetros:
  # alpha + (k betas) + (p1 phis) + (q1 thetas) + sigma
  n_par <- k + p1 + q1 + 2
  
  initial <- rep(0, n_par)
  # Preenche: alpha + betas + phis
  len_ols <- min(length(ols), (k + p1 + 1))
  initial[1:len_ols] <- ols[1:len_ols]
  initial[n_par] <- 1  # valor inicial para sigma
  
  # b) Função de log-verossimilhança
  loglik <- function(z) 
  {
    alpha <- z[1]
    
    if(k>0)  beta  <- z[2:(1+k)] else beta <- numeric(0)
    i_phi <- (1 + k + 1) : (1 + k + p1)
    if(p1>0) phi <- z[i_phi] else phi <- numeric(0)
    i_theta <- (1 + k + p1 + 1) : (1 + k + p1 + q1)
    if(q1>0) theta <- z[i_theta] else theta <- numeric(0)
    
    sigma <- z[length(z)]
    
    Xbeta <- as.vector(X %*% beta)
    if(any(p1>0)) {
      Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p,1,na.rm=TRUE)))
    } else {
      Xbeta_ar <- rep(0, (n-1))  # se sem AR
    }
    
    # Recursão de eta e erro
    for(i in (m+1):n){
      eta[i] <- alpha
      if(k>0)  eta[i] <- eta[i] + Xbeta[i]
      if(p1>0) {
        eta[i] <- eta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar]) %*% phi
      }
      if(q1>0) {
        eta[i] <- eta[i] + t(theta) %*% error[i - ma]
      }
      error[i] <- ynew[i] - eta[i]
    }
    
    mu <- linkinv(eta[(m+1):n])
    x  <- y[(m+1):n]
    
    # Log-verossimilhança UGo (exemplo adaptado)
    ll <-  log( log(tau)/(1 - mu^-sigma ) ) + 
      log(sigma) + 
      (-1 - sigma)*log(x) + 
      (log(tau)/(1 - mu^-sigma))*(1 - x^-sigma)
    
    sum(ll)
  }
  
  # c) Score (gradiente)
  
  escore.UGoarma <- function(z)
  {
    alpha <- z[1]
    
    if(k>0)  beta  <- z[2:(1+k)] else beta <- numeric(0)
    i_phi <- (1 + k + 1) : (1 + k + p1)
    if(p1>0) phi <- z[i_phi] else phi <- numeric(0)
    i_theta <- (1 + k + p1 + 1) : (1 + k + p1 + q1)
    if(q1>0) theta <- z[i_theta] else theta <- numeric(0)
    
    # sigma deve ser o último
    sigma <- z[length(z)]
    
    # Recalcular eta/error a cada iteração
    Xbeta <- as.vector(X %*% beta)
    if(p1>0) {
      Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p,1,na.rm=TRUE)))
    } else {
      Xbeta_ar <- rep(0, (n-1))
    }
    
    for(i in (m+1):n){
      eta[i] <- alpha
      if(k>0)  eta[i] <- eta[i] + Xbeta[i]
      if(p1>0) {
        eta[i] <- eta[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar]) %*% phi
      }
      if(q1>0) {
        eta[i] <- eta[i] + t(theta) %*% error[i - ma]
      }
      error[i] <- ynew[i] - eta[i]
    }
    
    # Precisamos de mu e x para as derivadas
    mu <- linkinv(eta[(m+1):n])
    x  <- y[(m+1):n]
    
    # Montar R p/ termo MA
    if(q1>0) {
      for(i in 1:(n-m)) {
        R[i,] <- error[i+m - ma] * k_i
      }
    }
    
    # Zeramos derivadas
    deta.dalpha[]  <- 0
    if(k>0)  deta.dbeta[,]  <- 0
    if(p1>0) deta.dphi[,]   <- 0
    if(q1>0) deta.dtheta[,] <- 0
    
    for(i in (m+1):n) {
      # d(eta_i)/d(alpha)
      deta.dalpha[i] <- 1
      if(q1>0 && (i-ma)>0) {
        deta.dalpha[i] <- deta.dalpha[i] - (deta.dalpha[i-ma] %*% theta)
      }
      
      # d(eta_i)/d(beta)
      if(k>0){
        temp_dbeta <- X[i, ]
        # subtrair a parte AR?
        if(p1>0) {
          temp_dbeta <- temp_dbeta - t(phi) %*% X[i-ar, ]
        }
        # parte MA:
        if(q1>0 && (i-ma)>0) {
          temp_dbeta <- temp_dbeta - t(theta) %*% deta.dbeta[i-ma, ]
        }
        deta.dbeta[i, ] <- temp_dbeta
      }
      
      # d(eta_i)/d(phi)
      if(p1>0){
        temp_dphi <- (ynew_ar[i - ar] - Xbeta[i - ar])
        if(q1>0 && (i-ma)>0) {
          temp_dphi <- temp_dphi - t(theta) %*% deta.dphi[i-ma, ]
        }
        deta.dphi[i, ] <- temp_dphi
      }
      
      # d(eta_i)/d(theta)
      if(q1>0){
        temp_dtheta <- R[i-m, ]
        if((i-ma)>0){
          temp_dtheta <- temp_dtheta - t(theta) %*% deta.dtheta[i-ma, ]
        }
        deta.dtheta[i, ] <- temp_dtheta
      }
    }
    
    # Pegamos as partes (m+1):n
    v  <- deta.dalpha[(m+1):n]
    rM <- if(k>0)  deta.dbeta[(m+1):n, , drop=FALSE]  else matrix(0,0,0)
    rP <- if(p1>0) deta.dphi[(m+1):n, , drop=FALSE]   else matrix(0,0,0)
    rR <- if(q1>0) deta.dtheta[(m+1):n, , drop=FALSE] else matrix(0,0,0)
    
    # diag(mu.eta(eta_i))
    mT <- diag(mu.eta(eta[(m+1):n]))
    
    # Derivada no preditor: a_t = d(logLik)/d(eta_i)
    # (segue fórmula da UGo; ajuste se necessário)
    a_t <- -sigma * mu^(-sigma - 1) / (1 - mu^-sigma) -
      log(tau) * (1 - x^-sigma) * sigma * mu^(sigma-1) / (1 - mu^sigma)^2
    
    # Derivada no sigma:
    y_sust <- as.vector(
      1/sigma + log(mu)/(1 - mu^sigma) - log(x) +
        ( mu^sigma * log(tau) * x^(-sigma) *
            ((mu^sigma - 1)*log(x) - log(mu)*(x^sigma - 1))
        ) / (mu^sigma - 1)^2
    )
    
    # Montar cada componente do score:
    Ualpha <- sum( v * (mT %*% a_t) )
    
    Ubeta  <- if(k>0) {
      colSums( rM * as.vector(mT %*% a_t) )
    } else numeric(0)
    
    Uphi   <- if(p1>0) {
      colSums( rP * as.vector(mT %*% a_t) )
    } else numeric(0)
    
    Utheta <- if(q1>0) {
      colSums( rR * as.vector(mT %*% a_t) )
    } else numeric(0)
    
    Uc <- sum(y_sust)
    
    # Agora concatenamos na ordem alpha, (betas), (phis), (thetas), sigma
    score <- c(Ualpha)
    if(k>0)  score <- c(score, Ubeta)
    if(p1>0) score <- c(score, Uphi)
    if(q1>0) score <- c(score, Utheta)
    score <- c(score, Uc)
    
    return(score)
  }
  
  #----------------------------------------------------------------------------
  # 5) Chamada ao 'optim'
  #----------------------------------------------------------------------------
  
  opt <- optim(
    par     = initial,
    fn      = loglik,
    gr      = escore.UGoarma,
    method  = "BFGS",
    hessian = TRUE,
    control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12)
  )
  
  if (opt$convergence != 0) {
    warning("FUNÇÃO NÃO CONVERGIU COM GRADIENTE ANALÍTICO!")
    opt <- optim(
      par     = initial,
      fn      = loglik, 
      method  = "BFGS",
      hessian = TRUE,
      control = list(fnscale = -1, maxit = maxit1, reltol = 1e-12)
    )
    if (opt$convergence != 0){
      warning("NÃO CONVERGIU NEM COM GRADIENTE NUMÉRICO!")
    } else {
      warning("CONVERGIU APENAS COM GRADIENTE NUMÉRICO!")
    }
  }
  
  # Armazenar resultados
  z$conv  <- opt$convergence
  z$par   <- opt$par
  z$loglik <- opt$value
  
  # Tamanho do Hessiano
  J_inv <- try(solve(-opt$hessian), silent=TRUE)
  if(inherits(J_inv, "try-error")){
    # se deu erro, não calcula EP
    z$stderror <- rep(NA, length(opt$par))
  } else {
    z$stderror <- sqrt(diag(J_inv))
  }
  
  z$zstat   <- z$par / z$stderror
  z$pvalues <- 2*(1 - pnorm(abs(z$zstat)))
  
  # Nomear parâmetros
  names_par <- c("alpha", 
                 if(k>0) names_beta,
                 if(p1>0) names_phi,
                 if(q1>0) names_theta,
                 "sigma")
  names(z$par)      <- names_par
  names(z$stderror) <- names_par
  names(z$zstat)    <- names_par
  names(z$pvalues)  <- names_par
  
  # Critérios de informação
  z$k <- length(z$par)
  z$aic <- -2 * z$loglik + 2 * z$k
  z$bic <- -2 * z$loglik + log(n) * z$k
  z$hq  <- -2 * z$loglik + log(log(n)) * z$k
  
  # Resumo do modelo
  z$model <- cbind(
    Estimate   = round(z$par, 4),
    `Std. Err` = round(z$stderror, 4),
    `z value`  = round(z$zstat, 4),
    `Pr(>|z|)` = round(z$pvalues, 4)
  )
  
  #----------------------------------------------------------------------------
  # 6) Fitted values e Forecast
  #----------------------------------------------------------------------------
  
  # Extrair parâmetros para predição
  alpha <- z$par[1]
  if(k>0)  beta  <- z$par[2:(1+k)] else beta <- numeric(0)
  i_phi <- (1 + k + 1) : (1 + k + p1)
  if(p1>0) phi <- z$par[i_phi] else phi <- numeric(0)
  i_theta <- (1 + k + p1 + 1) : (1 + k + p1 + q1)
  if(q1>0) theta <- z$par[i_theta] else theta <- numeric(0)
  sigma <- z$par[length(z$par)]
  
  # Criar vetores para armazenar valores
  errorhat <- rep(0, n)
  etahat   <- rep(NA, n)
  
  Xbeta <- if(k>0) as.vector(X %*% beta) else rep(0, n)
  if(p1>0) {
    Xbeta_ar <- suppressWarnings(matrix(Xbeta, (n-1), max(p,1,na.rm=TRUE)))
  } else {
    Xbeta_ar <- rep(0, (n-1))
  }
  
  for(i in (m+1):n){
    etahat[i] <- alpha
    if(k>0) etahat[i] <- etahat[i] + Xbeta[i]
    if(p1>0) {
      etahat[i] <- etahat[i] + (ynew_ar[(i-1), ar] - Xbeta_ar[(i-1), ar]) %*% phi
    }
    if(q1>0) {
      etahat[i] <- etahat[i] + t(theta) %*% errorhat[i-ma]
    }
    errorhat[i] <- ynew[i] - etahat[i]
  }
  
  # fitted no espaço original
  fitted_vals <- linkinv(etahat[(m+1):n])
  z$fitted    <- ts(c(rep(NA, m), fitted_vals), start=start(y), frequency=frequency(y))
  z$etahat    <- etahat
  z$errorhat  <- errorhat
  
  # Forecast (h passos)
  ynew_prev <- c(ynew, rep(NA, h))
  y_prev    <- rep(NA, n+h)
  y_prev[1:n] <- z$fitted
  
  # Se existirem covariáveis futuras
  if(k>0 && !any(is.na(X_hat))) {
    X_hat <- as.matrix(X_hat)
    X_prev <- rbind(X, X_hat)
  } else {
    X_prev <- X  # se não existir X_hat ou k=0
  }
  
  for(i in 1:h){
    idx <- n + i
    if(k>0){
      # alpha + X_hat
      etanext <- alpha + X_prev[idx, ] %*% beta
      if(p1>0){
        etanext <- etanext + phi %*% ( ynew_prev[idx - ar] - X_prev[idx - ar, ] %*% beta )
      }
    } else {
      # sem covariáveis
      etanext <- alpha
      if(p1>0){
        etanext <- etanext + phi %*% ynew_prev[idx - ar]
      }
    }
    if(q1>0){
      etanext <- etanext + theta %*% errorhat[idx - ma]
    }
    
    ynew_prev[idx] <- etanext
    y_prev[idx]    <- linkinv(etanext)
    
    # Se mantemos errorhat para previsões, assumindo E(error)=0
    errorhat[idx] <- 0
  }
  
  z$forecast <- y_prev[(n+1):(n+h)]
  
  #----------------------------------------------------------------------------
  # 7) Resíduos quantílicos (se tiver pUGo disponível)
  #----------------------------------------------------------------------------
  if(exists("pUGo")) {
    z$residuals <- as.vector(qnorm(pUGo(y[(m+1):n], z$fitted[(m+1):n], sigma)))
  } else {
    z$residuals <- NA
    warning("Função pUGo não encontrada, resíduos quantílicos não foram calculados.")
  }
  
  #----------------------------------------------------------------------------
  # 8) Diagnósticos (opcional)
  #----------------------------------------------------------------------------
  if(diag > 0){
    # Exibir resumo
    print(z$model)
    cat("Log-likelihood:", z$loglik, "\n")
    cat("AIC:", z$aic, " BIC:", z$bic, " HQ:", z$hq, "\n")
    cat("Gradiente convergiu? (0 = ok):", z$conv, "\n")
    
    if(diag > 1){
      # Exemplo de plot
      par(mfrow=c(1,1))
      plot(y, type="l", ylab="Série", xlab="Tempo", 
           main="Dados vs. Ajuste (Mediana)",
           ylim=range(y, na.rm=TRUE))
      lines(z$fitted, col="blue", lty=2)
      legend("topright", c("Observado","Ajuste (Mediana)"),
             col=c("black","blue"), lty=c(1,2), bty="n")
    }
  }
  
  return(z)
}

#------------------------------------------------------------------------------
# Exemplo de uso
#------------------------------------------------------------------------------
# Ajuste o caminho para suas funções:
source("simu.ugoarma.R")

# Exemplo com (p1=1, q1=1) e nenhuma covariável => total de 4 parâmetros:
# set.seed(2)
# y <- simu.ugoarma(100, phi=0.2, theta=0.4, alpha=1, sigma=6, tau=0.5, freq=12, link="logit")
# 
# 
# 
# fit <- uGoarma.fit(y, ma=1, ar=1, diag=1)
# 
# fit$model
