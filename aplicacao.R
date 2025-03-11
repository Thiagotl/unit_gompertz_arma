
# APLICACAO 

#rm(list = ls()) 

library(tidyverse)
library(e1071)

source("ugo_fit.R")
source("functions.R")
source("best_ugoarma.R")


dados <- read_csv("EAR_southeast_may01-2000_aug31-2019.csv")
#View(dados)


y1 <- dados$val_eaconsimp4/100
y1

month = 5
year = 2000
tau = 0.5

# SERIE COMPLETA
Y <- ts(y1,
        start = c(year, month),
        frequency = 12)

plot(Y)

# TESTES DE TENDENCIA DETERMINISTICA 
tend_determ(Y)

# TESTE DE RAIZ UNITARIA 
raiz_unit(Y)

# TESTE DE SAZONALIDADE 

sazonalidade(Y)

# Tamanho da amostra

n_first <- length(Y)

# numeros de passos para previsão

h1 = 10

month_h = end(Y)[2]-h1+1

# remover as últimas obsevações de h1

n =  n_first-h1

# Dados de Treino 
y = ts(y1[1:n],
       start = c(year, month),
       frequency = 12)

# Dados de Teste

yh = ts(Y[(n+1):n_first],
        start = c(end(Y)[1],month_h),
        frequency = 12)

summary(yh)

plot(stl(y, s.window = 12))

plot(y)


monthplot(y)

acf(y)
pacf(y)



summary(y)  # resume
var(y)      # variance
skewness(y) # skewness
kurtosis(y) # kurtosis



# REGRESSORES - Variaveis Explicativas 


# Tendencia Temporal
t = 1:length(y) # sequencia de 1 até o o tamanho da serie de treinamento
t_hat = (n+1):(n+h1) # sequencia numerica para o periodo de testes - previsao

#Sazonalidade Deterministica 

C = cos(2*pi*t/12)  # Componente cosenoidal para treino
C_hat = cos(2*pi*t_hat/12)  # Componente cosenoidal para teste  

# Variavel Dummy para Crise Hidrica

S = sin(2*pi*t/12)  # Componente senoidal para treino
S_hat = sin(2*pi*t_hat/12)  # Componente senoidal para teste

# Variavel para periodo de crise hidrica (0 ou 1)
which(as.integer(time(y)) < 2002)

which(as.integer(time(y)) >= 2013)



l1 = length(which(as.integer(time(y)) < 2002))  # Número de observações até 2002
l2 = which(as.integer(time(y)) >= 2013)[1]  # Primeiro índice do período >= 2013

D = c(rep(1, l1),                             # Período de crise (2000-2001)
      rep(0, (length(y) - l1 - (length(y) - l2 + 1))), # Período sem crise (2002-2012)
      rep(1, length(y) - l2 + 1))             # Período de crise (2013 em diante)

D_hat = rep(1, h1)  # Supõe crise contínua nos próximos 10 meses de previsão




X = cbind(C, S, D)   # Matriz de regressoras para treino
X_hat = cbind(C_hat, S_hat, D_hat)  # Matriz de regressoras para teste


# escolha do melhor modelo ARMA para UGO-ARMA

pmax = 3
qmax = 3

ugo_best = best.ugo(y, sf = c(start = c(year,month),frequency = 12),
                        h = h1, pmax = pmax, qmax = qmax,
                        nbest = 10, tau = tau, link = "logit", X = X, X_hat = X_hat)






p_ugoarma = 1:2
q_ugoarma = 2:3

fit_ugoarma = uGoarma.fit(y,
                          ar = p_ugoarma,
                          ma = q_ugoarma,
                          tau = tau,
                          link  = "logit",
                          h = h1, 
                          diag = 1,
                          X = X, 
                          X_hat = X_hat)

fit_ugoarma$model

acf(fit_ugoarma$residuals)
pacf(fit_ugoarma$residuals)




# PREVISAO FORA DA AMOSTRA 

X_hat1 = cbind(C_hat,S_hat,D_hat)
# Test sets
test_data  = cbind(yh, X_hat)

# Criando matriz para armazenar previsões
forecasts_ugoarma = matrix(NA, nrow = 10, ncol = 10)
diag = 1

for (i in 1:nrow(test_data)) {
  
  # Atualizando os regressores para o período de previsão
  if(i == 1){
    X_hat = t(as.matrix(X_hat1[1:i,]))
  }else{
    X_hat = X_hat1[1:i,]
  }
  
  # Ajustando o modelo UGOARMA a cada iteração
  fit_ugoarma = uGoarma.fit(y,
                            ar    = 1:2,
                            ma    = 1:2,
                            tau   = tau,
                            link  = "logit",
                            h     = i, 
                            diag  = diag, 
                            X     = X,
                            X_hat = X_hat)
  
  # Armazenando previsões
  forecasts_ugoarma[i,1:i] = fit_ugoarma$forecast
}

# Exibindo previsões
print("Forecasted values (UGO-ARMA):")
round(forecasts_ugoarma, 4)



# CALCULO DAS METRICAS DE QUALIDADE DO AJUSTE
mse_ugoarma = numeric(h1)
mape_ugoarma = numeric(h1)

for (i in 1:h1) {
  actual_values = yh[1:i]
  forecasts_current_ugoarma = as.numeric(na.omit(forecasts_ugoarma[i,]))
  
  mse_ugoarma[i]  = mean((forecasts_current_ugoarma - actual_values)^2)
  mape_ugoarma[i] = mean(abs((forecasts_current_ugoarma - actual_values) / actual_values)) * 100
}

# Exibindo métricas de erro
results = rbind(mse_ugoarma, mape_ugoarma)
round(results, 4)



# GRAFICOS


plot(y,
     type = "l",
     ylab = "Stored hydroelectric energy",
     xlab = "Time",
     ylim = c(min(y), max(y)))
lines(fit_ugoarma$fitted, col = "blue", lty = 2)
legend("topright",
       c("Observed data","Predicted median"),
       pt.bg = "white", 
       lty = c(1,2), 
       bty = "n",
       col = c(1, "blue"))


out_of_sample_forecast <- forecasts_ugoarma[h1, 1:h1]

out_of_sample_forecast_ts <- ts(out_of_sample_forecast,
                                start = c(end(y)[1], end(y)[2]+1), 
                                frequency = 12)

# Plotar série completa (treino + teste)
plot(Y, main = "Previsão Fora da Amostra", xlab = "Tempo", ylab = "Energia")
lines(out_of_sample_forecast_ts, col = "red", lty = 2)
#points(yh, col = "blue", pch = 16)  # valores reais de teste
legend("topright", 
       legend = c("Observado", "Previsto"), 
       col = c("black", "red"), lty = c(1,2), bty = "n")





