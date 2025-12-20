rm(list = ls())
library(BTSR)
source("ugo_fit.R")

# INTEREST RATE SERIES
dados <- readxl::read_excel("dataset.xlsx", na = "-")
Y <- ts(dados[2] / 100, frequency = 12, start = c(2011,03))

m <- length(Y)
h1 <- 6
n <- m - h1

# TRAIN SERIES
y_train <- ts(Y[1:n], frequency = 12)

# TEST SERIES
y_test <- Y[(n + 1):m]

#### REGRESSORS ----
t      <- 1:length(y_train)
t_hat  <- (n + 1):(n + h1)

C      <- cos(2 * pi * t / 12)
C_hat  <- cos(2 * pi * t_hat / 12)

#### TREND INVESTIGATION
### UGOARMA
fit_ugoarmax_trend <- uGoarma.fit(y_train, 
                                  ar = 1, 
                                  ma = NA, 
                                  X = cbind(C,t), 
                                  X_hat = cbind(C_hat, t_hat)
)

### BARMA
barmax_trend <- BARFIMA.fit(y_train,
                            p = 1,
                            d = F,
                            q = 0,
                            xreg = cbind(C,t),
                            start = list(
                              alpha = 0,
                              phi   = 0,
                              beta  = coef(lm(y_train ~ cbind(C,t) + 0)),
                              nu    = 1
                            ),
                            control = list(
                              method = "Nelder-Mead",
                              stopcr = 1e-2
                            ),
                            info = T,
                            report = F
) 

karmax_trend <- KARFIMA.fit(y_train,
                            p = 1,
                            d = F,
                            q = 0,
                            xreg = cbind(C,t),
                            control = list(
                              method = "Nelder-Mead",
                              stopcr = 1e-2
                            ),
                            info = T,
                            report = F
) 

### fitted values
fit_ugoarmax_trend$model
summary(barmax_trend)$coeff
summary(karmax_trend)$coeff 

ugoarmax_trend_fitted <- KARFIMA.extract(
  yt = Y, xreg = rbind(cbind(C,t),cbind(C_hat,t_hat)),
  coefs = list(
    alpha = fit_ugoarmax_trend$coeff[1],
    beta = fit_ugoarmax_trend$coeff[2:3],
    phi = fit_ugoarmax_trend$coeff[4],
    theta = NULL,
    nu = fit_ugoarmax_trend$coeff[5]
  )
)

barmax_trend_fitted <- BARFIMA.extract(
  yt = Y, xreg = rbind(cbind(C,t),cbind(C_hat,t_hat)),
  coefs = list(
    alpha = barmax_trend$coefficients[1],
    beta = barmax_trend$coefficients[2:3],
    phi = barmax_trend$coefficients[4],
    theta = NULL,
    nu = barmax_trend$coefficients[5]
  )
)

karmax_trend_fitted <- KARFIMA.extract(
  yt = Y, xreg = rbind(cbind(C,t),cbind(C_hat,t_hat)),
  coefs = list(
    alpha = karmax_trend$coefficients[1],
    beta = karmax_trend$coefficients[2:3],
    phi = karmax_trend$coefficients[4],
    theta = NULL,
    nu = karmax_trend$coefficients[5]
  )
)

results_trend_insample <- rbind(
  forecast::accuracy(ugoarmax_trend_fitted$mut[1:n], y_train),
  forecast::accuracy(barmax_trend_fitted$mut[1:n], y_train),
  forecast::accuracy(karmax_trend_fitted$mut[1:n], y_train)
)[, c(3, 2, 5)]


results_trend_outsample <- rbind(
  forecast::accuracy(ugoarmax_trend_fitted$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(barmax_trend_fitted$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(karmax_trend_fitted$mut[(n + 1):length(y_test)], y_test)
)[, c(3, 2, 5)]

results_trend_insample
xtable::xtable(results_trend_outsample, digits = 4)
