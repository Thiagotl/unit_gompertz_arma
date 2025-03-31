library(forecast)
library(BTSR)
library(tidyverse)
library(e1071)
source("ugo_fit.R")
source("functions.R")
library(readxl)


dados <- read_excel("ipeadata[31-03-2025-03-21].xls")

dados$`Taxa de desemprego`<-dados$`Taxa de desemprego`/100


hist(dados$`Taxa de desemprego`)

# SERIE COMPLETA
Y<-ts(dados$`Taxa de desemprego`, start = c(2012,3), frequency = 12)

plot(Y)

m<-length(Y)
h1 <- 12
n<-m-h1

# SERIE DE TREINO
y_train<-ts(Y[1:n], start=c(2012,3), frequency = 12)

# SERIE DE TESTE 
y_test<-Y[(n+1):m] 

tend_determ(Y)
raiz_unit(Y)
sazonalidade(Y)


# Tendencia Temporal
t = 1:length(y_train) # sequencia de 1 até o o tamanho da serie de treinamento
t_hat = (n+1):(n+h1) # sequencia numerica para o periodo de testes - previsao

#Sazonalidade Deterministica 

C = cos(2*pi*t/12)  # Componente cosenoidal para treino
C_hat = cos(2*pi*t_hat/12)  # Componente cosenoidal para teste  


S = sin(2*pi*t/12)  # Componente senoidal para treino
S_hat = sin(2*pi*t_hat/12)  # Componente senoidal para teste


X = cbind(C, S)   # Matriz de regressoras para treino
X_hat = cbind(C_hat, S_hat)  # Matriz de regressoras para teste


a01<-auto.arima(y_train)
new1<-Arima(y_test,model=a01) #one-step-ahead
#forecast(a01, h = length(hum_test))


a02<-auto.arima(y_train, xreg = X)
new2<-Arima(y_test,xreg = X_hat,model=a02) #one-step-ahead


#### APLICACAO MELHOR MODELO UGO ARMA----
source("best_ugoarma2.R")

pmax = 3
qmax = 3

y<-y_train
ugo_arma.model<-best_ugo_2(y_train, pmax = pmax, qmax = qmax,
                  nbest = 8, X=X, X_hat = X_hat) # 1 2 -1378.1261



#### APLICACAO BETA E KW---- 

quant<-.5 # quantil
# matrix de resultados
cont<-1

order <- matrix(NA, nrow = 16, ncol = 6)
colnames(order) <- c("p", "q", "barma.AIC", "karma.AIC", "barmax.AIC", "karmax.AIC")


for(i in 0:3){
  for(j in 0:3){
    
    
    barma<-summary(BARFIMA.fit(y_train,p=i,d=F,q=j,info=T,
                               report=F))
    karma1<-suppressWarnings(KARFIMA.fit(y_train,p=i,d=F,q=j,info=T,
                                         rho=quant,
                                         control = list(method="Nelder-Mead",stopcr=1e-2),
                                         report=F))
    karma<-summary(karma1)
   
    
    barmax<-summary(BARFIMA.fit(y_train,p=i,d=F,q=j,info=T,
                                xreg = X,
                                report=F))
    karmax1<-suppressWarnings(KARFIMA.fit(y_train,p=i,d=F,q=j,info=T,
                                          xreg = X,rho=quant,
                                          control = list(method="Nelder-Mead",stopcr=1e-2),
                                          report=F))
    karmax<-summary(karmax1)
    
    
  
    if(karma1$convergence==1 || is.nan(karma$aic)==1) karma$aic=0
    if(karmax1$convergence==1 || is.nan(karmax$aic)==1) karmax$aic=0
    
    
    order[cont, ] <- c(i, j, barma$aic, karma$aic, barmax$aic, karmax$aic)
    cont<-cont+1
    
  }
}

order<-order[-1,]
print(order)






