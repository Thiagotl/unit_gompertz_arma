
source("ugo_fit.R")

y <- ts(c(.10, .12, .15, .13, .17, .16, .20, .18, .19, .23)) # Exemplo de série temporal

X <- matrix(rnorm(length(y) * 2), ncol = 2)  # Covariáveis fictícias
X_hat <- matrix(rnorm(5 * 2), ncol = 2)      # Covariáveis futuras

ugo_model <- uGoarma.fit(
  y = y,
  ar = 1, 
  ma = 1, 
  tau = 0.5, 
  link = "logit", 
  h = 5, 
  diag = 1, 
  X = X, 
  X_hat = X_hat
)

plot(y, type = "l", col = "black", main = "Observed vs Fitted")
lines(ugo_model$fitted, col = "blue", lty = 2)


plot(ugo_model$residuals, main = "Residuals")
abline(h = 0, col = "red", lty = 2)

cat("AIC:", ugo_model$aic, "BIC:", ugo_model$bic, "HQ:", ugo_model$hq)

print(ugo_model$forecast)
