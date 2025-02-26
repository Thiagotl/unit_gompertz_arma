# EAR Diário por Bacia
# EAR Diário por REE - Reservatório Equivalente de Energia
# EAR Diário por Reservatório
# EAR Diário por Subsistema
# ENA Diário por Bacia
# ENA Diário por REE - Reservatório Equivalente de Energia
# ENA Diário por Reservatório
# ENA Diário por Subsistema
# Energia Vertida Turbinável
# Equipamentos de Controle de Reativos da Rede de Operação
# Fator de Capacidade de Geração Eólica e Solar

rm(list = ls())
library(forecast)
library(dplyr)
library(BTSR)
library(tidyverse)
library(e1071)
source("ugo_fit.R")
source("functions.R")

dados <- read_csv("EAR_southeast_may01-2000_aug31-2019.csv")
View(dados)

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
y_train = ts(y1[1:n], #y 
             start = c(year, month),
             frequency = 12)

# Dados de Teste

y_test = ts(Y[(n+1):n_first], #yh
            start = c(end(Y)[1],month_h),
            frequency = 12)

summary(y_test)

plot(stl(y_train, s.window = 12))

plot(y_train)


monthplot(y_train)

acf(y_train)
pacf(y_train)

plot(y_test)

summary(y_train)  # resume
var(y_train)      # variance
skewness(y_train) # skewness
kurtosis(y_train) # kurtosis


# REGRESSORES - Variaveis Explicativas 


# Tendencia Temporal
t = 1:length(y_train) # sequencia de 1 até o o tamanho da serie de treinamento
t_test = (n+1):(n+h1) # sequencia numerica para o periodo de testes - previsao

#Sazonalidade Deterministica 

C = cos(2*pi*t/12)  # Componente cosenoidal para treino
C_test = cos(2*pi*t_test/12)  # Componente cosenoidal para teste  

# Variavel Dummy para Crise Hidrica

S = sin(2*pi*t/12)  # Componente senoidal para treino
S_test = sin(2*pi*t_test/12)  # Componente senoidal para teste

# Variavel para periodo de crise hidrica (0 ou 1)
which(as.integer(time(y_train)) < 2002)

which(as.integer(time(y_train)) >= 2013)

l1 = length(which(as.integer(time(y_train)) < 2002))  # Número de observações até 2002
l2 = which(as.integer(time(y_train)) >= 2013)[1]  # Primeiro índice do período >= 2013

D = c(rep(1, l1),                             # Período de crise (2000-2001)
      rep(0, (length(y_train) - l1 - (length(y_train) - l2 + 1))), # Período sem crise (2002-2012)
      rep(1, length(y_train) - l2 + 1))             # Período de crise (2013 em diante)

D_test = rep(1, h1)  # Supõe crise contínua nos próximos 10 meses de previsão


X = cbind(C, S, D)   # Matriz de regressoras para treino
X_test = cbind(C_test, S_test, D_test)  # Matriz de regressoras para teste

###### TESTES DE MODELO

# MODELO 1

mod1<-auto.arima(y_train)
mod1

mod1_new<-Arima(y_test, model = mod1)





