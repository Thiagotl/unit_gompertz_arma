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

#### REGRESSORS MATRIX ----
t      <- 1:length(y_train)
t_hat  <- (n + 1):(n + h1)

C      <- cos(2 * pi * t / 12)
C_hat  <- cos(2 * pi * t_hat / 12)

X     <- cbind(cos = C)
X_hat <- cbind(cos = C_hat)

nX <- ncol(X)
X0 <- rbind(X, X_hat)

#### FITTING THE UGO ARMA----
source("best_ugoarma.R")

pmax <- 3; qmax <- 3
y <- y_train

best_ugoarmax <- suppressWarnings(best_ugo_2(y_train,
                                             pmax = pmax, qmax = qmax,
                                             nbest = 8, X = X, X_hat = X_hat
))

if (best_ugoarmax$q[1] == 0) {
  fit_ugoarmax <- uGoarma.fit(
    y_train, 
    ar = 1:best_ugoarmax$p[1], 
    ma = NA, X = X, X_hat = X_hat)
} else {
  fit_ugoarmax <- uGoarma.fit(y_train, ar = 1:best_ugoarmax$p[1], ma = 1:best_ugoarmax$q[1], X = X, X_hat = X_hat)
}

#### FITTING THE BETA E KW ----
quant <- .5 # quantile
ps <- 0:3; qs <- 0:3

order <- matrix(NA, nrow = length(ps) * length(qs), ncol = 4)
cont <- 1

colnames(order) <- c(
  "p", "q", "barmax.AIC", "karmax.AIC"
)

start_btsr <- function(i,j,nu1=1){
  list(
    alpha = 0,
    phi   = rep(0, i),
    theta = rep(0, j),
    beta  = coef(lm(y_train ~ X + 0)),
    nu    = nu1
  )
}

for (i in 0:3) {
  for (j in 0:3) {
    barmax <- summary(BARFIMA.fit(
      y_train,
      p = i,
      d = F,
      q = j,
      info = T,
      start = start_btsr(i,j),
      control = list(
        method = "Nelder-Mead",
        stopcr = 1e-2),
      xreg = X,
      report = F
    ))
    
    karmax1 <- suppressWarnings(
      KARFIMA.fit(y_train,
                  p = i,
                  d = F,
                  q = j,
                  info = T,
                  control = list(
                    method = "Nelder-Mead",
                    stopcr = 1e-2),
                  xreg = X,
                  rho = quant,
                  report = F
      ))
    karmax <- summary(karmax1)
    
    if (karmax1$convergence == 1 || is.nan(karmax$aic) == 1) karmax$aic <- 0
    
    order[cont, ] <- c(
      i, j, 
      barmax$aic, karmax$aic
    )
    cont <- cont + 1
  }
}

order <- order[-1, ]

#### SELECTING THE FINAL MODELS ----

### selecting the BARMA best order
rank_barmax<-rank(order[,3])
## option 1
orbarmax1 <- order[which(rank_barmax==1), c(1:2, 3)]
barmax1 <- BARFIMA.fit(y_train,
                       p = orbarmax1[1],
                       d = F,
                       q = orbarmax1[2],
                       xreg = X,
                       rho = quant,
                       start = start_btsr(orbarmax1[1],orbarmax1[2]),
                       control = list(
                         method = "Nelder-Mead",
                         stopcr = 1e-2
                       ),
                       info = T,
                       report = F
) # several ARMA coefficients not significant

## option 2
orbarmax2 <- order[which(rank_barmax==2), c(1:2, 3)]
barmax2 <- BARFIMA.fit(y_train,
                       p = orbarmax2[1],
                       d = F,
                       q = orbarmax2[2],
                       xreg = X,
                       rho = quant,
                       start = start_btsr(orbarmax2[1],orbarmax2[2]),
                       control = list(
                         method = "Nelder-Mead",
                         stopcr = 1e-2
                       ),
                       info = T,
                       report = F
) # good candidate, all coeffs significants

## option 3
orbarmax3 <- order[which(rank_barmax==3), c(1:2, 3)]
barmax3 <- BARFIMA.fit(y_train,
                       p = orbarmax3[1],
                       d = F,
                       q = orbarmax3[2],
                       xreg = X,
                       rho = quant,
                       start = start_btsr(orbarmax3[1],orbarmax3[2]),
                       control = list(
                         method = "Nelder-Mead",
                         stopcr = 1e-2
                       ),
                       info = T,
                       report = F
) # good candidate, all coeffs significants, more parsimonious

barmax_insample<-rbind(
  forecast::accuracy(barmax2$fitted.values, y_train),
  forecast::accuracy(barmax3$fitted.values, y_train)
)[,c(3,2,5)] # barmax 3 outperforms for all measures and is more parsimonious

# barmax 3 is our final barmax model
barmax<-barmax3
orbarmax<-orbarmax3

### selecting the KARMA best order
rank_karmax<-rank(order[,4])

## option 1
orkarmax <- order[which(rank_karmax==1), c(1:2, 4)]
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

#### FINAL MODELS ----
###################################
### Table : fitted coefficients ###
###################################
fit_ugoarmax$model
round(summary(barmax)$coefficients,4)
round(summary(karmax)$coefficients,4)

### Table : in_sample predictions
names_rows <- c(
  "UGOMAX", "BARMAX", "KARMAX"
)

results_insample <- rbind(
  forecast::accuracy(fit_ugoarmax$fitted, y_train),
  forecast::accuracy(barmax$fitted.values, y_train),
  forecast::accuracy(karmax$fitted.values, y_train)
)[, c(3, 2, 5)]

rownames(results_insample) <- names_rows

#####################################
### Table : in_sample predictions ###
#####################################
round(results_insample,3)

#### OUT-OF-SAMPLE FORECASTING ----

### BARMA
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

### KARMA
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

### UGO-ARMA
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

results_outsample <- rbind(
  forecast::accuracy(ugoarma_out2$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(barmax_out$mut[(n + 1):length(y_test)], y_test),
  forecast::accuracy(karmax_out$mut[(n + 1):length(y_test)], y_test)
)[, c(3, 2, 5)]
#

row.names(results_outsample) <- names_rows

######################################
### Table : out_sample predictions ###
######################################
print(round(results_outsample, 4))

