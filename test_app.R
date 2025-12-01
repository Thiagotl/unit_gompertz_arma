library(forecast)
library(BTSR)
library(tidyverse)
library(e1071)
source("ugo_fit.R")
library(lubridate)
library(lmtest)

# DATASETS ------

# INTEREST RATE
dados <- readxl::read_excel("dataset.xlsx", na = "-")
# SERIE
Y <- ts(dados[2] / 100, start = c(2011, 3), 
        end = c(2025, 3), frequency = 12)

m <- length(Y)
h1 <- 6
n <- m - h1

# TRAIN SERIE
y_train <- ts(Y[1:n], frequency = 12)

# TEST SERIE
y_test <- Y[(n + 1):m]

#### REGRESSORS MATRIX ----

# TREND AND SEASONAL

t      <- 1:length(y_train)
t_hat  <- (n + 1):(n + h1)

t_sc      <- t / max(t)         
t_hat_sc  <- t_hat / max(t)

C      <- cos(2 * pi * t / 12)
C_hat  <- cos(2 * pi * t_hat / 12)

X     <- as.matrix(cbind(trend = t_sc, seasonal = C))
X_hat <- as.matrix(cbind(trend = t_hat_sc, seasonal = C_hat))

#X <- as.matrix(cbind(seasonal = C, trend = t))
#X_hat <- as.matrix(cbind(seasonal = C_hat, trend = t_hat))


nX <- ncol(X) #nX <- dim(X)[2]
X0 <- rbind(X, X_hat)

#### ARIMA FIT ----
a01 <- auto.arima(y_train, seasonal = FALSE)
new1 <- Arima(y_test, model = a01) # one-step-ahead

a02 <- auto.arima(y_train, xreg = X, seasonal = FALSE)
new2 <- Arima(y_test, xreg = X_hat, model = a02) # one-step-ahead


#### BEST MODEL UGO ARMA----
source("best_ugoarma2.R")

pmax <- 3
qmax <- 3
y <- y_train

# FIT WITH REGRESSORS -----
best_ugoarmax <- best_ugo_2(y_train,
                            pmax = pmax, qmax = qmax,
                            nbest = 8, X = X, X_hat = X_hat
)


if (best_ugoarmax$q[1] == 0) {
  fit_ugoarmax <- uGoarma.fit(
    y_train, 
    ar = 1:best_ugoarmax$p[1], 
    ma = NA, X = X, X_hat = X_hat)
} else {
  fit_ugoarmax <- uGoarma.fit(y_train, ar = 1:best_ugoarmax$p[1], ma = 1:best_ugoarmax$q[1], X = X, X_hat = X_hat)
}


# FIT WITHOUT REGRESSORS -------

best_ugoarma <- best_ugo_2(y_train,
                           pmax = pmax, qmax = qmax,
                           nbest = 8
)

if (best_ugoarma$q[1] == 0) {
  fit_ugoarma <- uGoarma.fit(y_train, ar = 1:best_ugoarma$p[1], ma = NA)
} else {
  fit_ugoarma <- uGoarma.fit(y_train, ar = 1:best_ugoarma$p[1], ma = 1:best_ugoarma$q[1])
}


#### BETA E KW APPLICATION ----
quant <- .5 # quantil
ps <- 0:3
qs <- 0:3

order <- matrix(NA, nrow = length(ps) * length(qs), ncol = 6)
cont <- 1


colnames(order) <- c(
  "p", "q", "barma.AIC", "karma.AIC",
  "barmax.AIC", "karmax.AIC"
)



for (i in 0:3) {
  for (j in 0:3) {
    barma <- summary(BARFIMA.fit(
      y_train,
      p = i,
      d = F,
      q = j,
      info = T,
      start = list(
        phi = rep(0.3, i),
        theta = rep(0.3, j)
      ),
      control = list(method = "Nelder-Mead", stopcr = 1e-2),
      report = F
    ))
    
    
    barmax <- summary(BARFIMA.fit(
      y_train,
      p = i,
      d = F,
      q = j,
      info = T,
      start = list(
        phi = rep(0.3, i),
        theta = rep(0.3, j),
        beta = coef(lm(y_train ~ X + 0)),
        nu = 0.01
      ),
      control = list(
        method = "Nelder-Mead", 
        stopcr = 1e-2),
      xreg = X,
      report = F
    ))
    
    karma1 <- suppressWarnings(KARFIMA.fit(
      y_train,
      p = i,
      d = F,
      q = j,
      info = T,
      rho = quant,
      control = list(
        method = "Nelder-Mead", 
        stopcr = 1e-2),
      report = F
    ))
    
    karma <- summary(karma1)
    
    karmax1 <- suppressWarnings(KARFIMA.fit(
      y_train,
      p = i,
      d = F,
      q = j,
      info = T,
      xreg = X,
      rho = quant,
      control = list(method = "Nelder-Mead", stopcr = 1e-2),
                                            report = F
    ))
    karmax <- summary(karmax1)
    
    if (karma1$convergence == 1 || is.nan(karma$aic) == 1) karma$aic <- 0
    if (karmax1$convergence == 1 || is.nan(karmax$aic) == 1) karmax$aic <- 0
    
    order[cont, ] <- c(
      i, j, barma$aic, karma$aic,
      barmax$aic, karmax$aic
    )
    cont <- cont + 1
  }
}

order <- order[-1, ]
print(order)

orbarma <- order[which(order[, 3] == min(order[, 3])), c(1:3)]
orkarma <- order[which(order[, 4] == min(order[, 4])), c(1:2, 4)]
orbarmax <- order[which(order[, 5] == min(order[, 5])), c(1:2, 5)]
orkarmax <- order[which(order[, 6] == min(order[, 6])), c(1:2, 6)]

barma <- BARFIMA.fit(y_train,
                     p = orbarma[1],
                     d = F,
                     q = orbarma[2],
                     rho = quant,
                     info = T,
                     start = list(
                       alpha = 0,
                       phi = rep(0.3, orbarma[1]),
                       if (orbarma[2] == 0){
                         NULL
                       }else{
                         theta = rep(0.3, orbarma[2]) 
                       },
                       #theta = rep(0.3, 1),
                       nu = 0.01
                     ),
                     control = list(
                       method = "Nelder-Mead",
                       stopcr = 1e-2
                     ),
                     report = F
)

karma <- KARFIMA.fit(y_train,
                     p = orkarma[1],
                     d = F,
                     q = orkarma[2],
                     rho = quant,
                     control = list(
                       method = "Nelder-Mead",
                       stopcr = 1e-2
                     ),
                     info = T,
                     report = F
)


start_barmax <- list(
  alpha = 0,
  phi   = rep(0.3, orbarmax[1]),
  theta = if (orbarmax[2] == 0) NULL else rep(0.3, orbarmax[2]),
  beta  = coef(lm(y_train ~ X + 0)),
  nu    = 0.01
)


barmax <- BARFIMA.fit(y_train,
                      p = orbarmax[1],
                      d = F,
                      q = orbarmax[2],
                      xreg = X,
                      rho = quant,
                      start = start_barmax,
                      control = list(
                        method = "Nelder-Mead",
                        stopcr = 1e-2
                      ),
                      info = T,
                      report = F
)



karmax <- KARFIMA.fit(y_train,
                      p = orkarmax[1],
                      d = F,
                      q = orkarmax[2],
                      rho = quant,
                      control = list(
                        method = "Nelder-Mead",
                        stopcr = 1e-2
                      ),
                      xreg = X,
                      info = T,
                      report = F
)


results_insample <- rbind(
  forecast::accuracy(fit_ugoarmax$fitted, y_train),
  forecast::accuracy(barmax$fitted.values, y_train),
  forecast::accuracy(karmax$fitted.values, y_train),
  forecast::accuracy(a02$fitted, y_train),
  forecast::accuracy(barma$fitted.values, y_train),
  forecast::accuracy(karma$fitted.values, y_train),
  forecast::accuracy(a01$fitted, y_train)
)[, c(3, 2, 5)]

rownames(results_insample) <- c(
  "UGOARMA", "BARMAX", "KARMAX", "a02",
  "BARMA", "KARMA", "a01"
)

round(results_insample, 4)

### out-of-sample

barma_out <- BARFIMA.extract(
  yt = Y,
  coefs = list(
    alpha = barma$coefficients[1],
    phi = barma$coefficients[2:(orbarma[1] + 1)],
    theta = if (orbarma[2] == 0) {
      NULL
    } else {
      barma$coefficients[(orbarma[1] + 2):(orbarma[1] + 1 + orbarma[2])]
    },
    nu = barma$coefficients[(orbarma[1] + 2 + orbarma[2])]
  )
)


karma_out <- KARFIMA.extract(
  yt = Y,
  coefs = list(
    alpha = karma$coefficients[1],
    phi = karma$coefficients[2:(orkarma[1] + 1)],
    theta = if (orkarma[2] == 0) {
      NULL
    } else {
      karma$coefficients[(orkarma[1] + 2):(orkarma[1] + 1 + orkarma[2])]
    },
    nu = karma$coefficients[(orkarma[1] + 2 + orkarma[2])]
  )
)


barmax_out <- BARFIMA.extract(
  yt = Y, xreg = X0,
  coefs = list(
    alpha = barmax$coefficients[1],
    beta = barmax$coefficients[2:(nX + 1)],
    phi = barmax$coefficients[(nX + 2):(orbarmax[1] + nX + 1)],
    theta = if (orbarmax[2] == 0) {
      NULL
    } else {
      barmax$coefficients[(orbarmax[1] + (nX + 2)):(orbarmax[1] + nX + 1 + orbarmax[2])]
    },
    nu = barmax$coefficients[(orbarmax[1] + (nX + 2) + orbarmax[2])]
  )
)


karmax_out <- KARFIMA.extract(
  yt = Y, xreg = X0, rho = quant,
  coefs = list(
    alpha = karmax$coefficients[1],
    beta = karmax$coefficients[2:(nX + 1)],
    phi = karmax$coefficients[(nX + 2):(orkarmax[1] + nX + 1)],
    theta = if (orkarmax[2] == 0) {
      NULL
    } else {
      karmax$coefficients[(orkarmax[1] + nX + 2):(orkarmax[1] + nX + 1 + orkarmax[2])]
    },
    nu = karmax$coefficients[(orkarmax[1] + nX + 2 + orkarmax[2])]
  )
)


# OUT OF SAMPLE FIT WITHOUT REGRESSORS ------

ugoarma_out1 <- KARFIMA.extract(
  yt = Y, rho = quant,
  coefs = list(
    alpha = fit_ugoarma$coeff[1],
    phi = fit_ugoarma$coeff[2:(best_ugoarma$p[1] + 1)],
    theta = if (best_ugoarma$q[1] == 0) {
      NULL
    } else {
      fit_ugoarma$coeff[(best_ugoarma$p[1] + 2):(best_ugoarma$p[1] + 1 + best_ugoarma$q[1])]
    },
    nu = fit_ugoarma$coeff[(best_ugoarma$p[1] + 2 + best_ugoarma$q[1])]
  )
)


# OUT OF SAMPLE FIT WITH REGRESSORS  ------

ugoarma_out2 <- KARFIMA.extract(
  yt = Y, xreg = X0, rho = quant,
  coefs = list(
    alpha = fit_ugoarmax$coeff[1],
    beta = fit_ugoarmax$coeff[2:(nX + 1)],
    phi = fit_ugoarmax$coeff[(nX + 2):(best_ugoarmax$p[1] + 1 + nX)],
    theta = if (best_ugoarmax$q[1] == 0) {
      NULL
    } else {
      fit_ugoarmax$coeff[(best_ugoarmax$p[1] + 2 + nX):(best_ugoarmax$p[1] + 1 + best_ugoarmax$q[1] + nX)]
    },
    nu = fit_ugoarmax$coeff[(best_ugoarmax$p[1] + 2 + best_ugoarmax$q[1] + nX)]
  )
)


a <- n + 1:length(y_test)

results_outsample <- rbind(
  forecast::accuracy(ugoarma_out2$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(barmax_out$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(karmax_out$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(new2$fitted, y_test),
  forecast::accuracy(ugoarma_out1$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(barma_out$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(karma_out$mut[(n + 1)]:length(y_test), y_test),
  forecast::accuracy(new1$fitted, y_test)
)[, c(3, 2, 5)]
#


names_rows <- c(
  "UGOMAX", "BARMAX", "KARMAX", "ARIMAX",
  "UGO", "BARMA", "KARMA", "ARIMA"
)


row.names(results_outsample) <- names_rows

# xtable::xtable((results_outsample),digits=4)

print(round(results_outsample, 4))

round(results_outsample[, ], 4)