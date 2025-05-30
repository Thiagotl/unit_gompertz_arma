library(forecast)
library(BTSR)
library(tidyverse)
library(e1071)
source("ugo_fit.R")
#source("functions.R")
library(readxl)
library(lubridate)
library(lmtest)


rm(list = ls())
gc()

# DADOS ----
dados1 <- read_excel("STP-20250509171131887.xlsx", na = "-")
dados <- na.omit(dados1[4])/100

# Criar coluna de datas reais
datas <- seq(as.Date("2011-03-01"), by = "month", length.out = nrow(dados))
dados <- dados %>% mutate(data = datas)

# Dummies estruturais
dados <- dados %>%
  mutate(
    crise_politica = if_else(data >= ymd("2015-01-01") & data <= ymd("2017-03-01"), 1, 0),
    pandemia_covid = if_else(data >= ymd("2020-03-01") & data <= ymd("2022-03-01"), 1, 0),
    gov_bolsonaro = if_else(data >= ymd("2019-01-01") & data <= ymd("2022-12-31"), 1, 0)
  )

# SÉRIE TEMPORAL ----
Y <- ts(dados[[1]], start = c(2011, 3), frequency = 12)
m <- length(Y)
h1 <- 6
n <- m - h1



y_train <- ts(Y[1:n], frequency = 12)
y_test <- Y[(n+1):m]

t <- 1:length(y_train)
t_hat <- (n+1):(n+h1)

C <- cos(2*pi*t/12)
C_hat <- cos(2*pi*t_hat/12)

# Dummies
# crise <- dados$crise_politica[1:n]
# crise_hat <- dados$crise_politica[(n+1):m]

# bolsonaro <- dados$gov_bolsonaro[1:n]
# bolsonaro_hat <- dados$gov_bolsonaro[(n+1):m]

pandemia <- dados$pandemia_covid[1:n]
pandemia_hat <- dados$pandemia_covid[(n+1):m]

# Outliers identificados
# outliers <- matrix(0, nrow = length(y_train), ncol = 2)
# outliers[108, 1] <- 1
# outliers[110, 2] <- 1
# outliers_hat <- matrix(0, nrow = length(y_test), ncol = 2)

# Matriz de regressão
X <- cbind(C,  pandemia)
X_hat <- cbind(C_hat, pandemia_hat)

nX <- dim(X)[2]
X0 <- rbind(X, X_hat)



a01<-auto.arima(y_train, seasonal=F)
new1<-Arima(y_test,model=a01) #one-step-ahead

a02<-auto.arima(y_train, xreg = X, seasonal=F)
new2<-Arima(y_test,xreg = X_hat,model=a02) #one-step-ahead


#### APLICACAO MELHOR MODELO UGO ARMA----
source("best_ugoarma2.R")

pmax = 3
qmax = 3

y<-y_train

# AJUSTE COM REGRESSÃO -----
best_ugoarma<-best_ugo_2(y_train, pmax = pmax, qmax = qmax,
                         nbest = 8, X=X, X_hat = X_hat) # 1 2 


if(best_ugoarma$q[1]==0){
  fit_ugoarma<-uGoarma.fit(y_train, ar=1:best_ugoarma$p[1], ma=NA, X=X, X_hat = X_hat )
}else{
  fit_ugoarma<-uGoarma.fit(y_train, ar=1:best_ugoarma$p[1], ma=1:best_ugoarma$q[1], X=X, X_hat = X_hat)
}




# AJUSTE SEM REGRESSÃO -------

best_ugoarma_sr<-best_ugo_2(y_train, pmax = pmax, qmax = qmax,
                            nbest = 8) 

if(best_ugoarma_sr$q[1]==0){
  fit_ugoarma_sr<-uGoarma.fit(y_train, ar=1:best_ugoarma_sr$p[1], ma=NA )
}else{
  fit_ugoarma_sr<-uGoarma.fit(y_train, ar=1:best_ugoarma_sr$p[1], ma=1:best_ugoarma_sr$q[1] )
}





#### APLICACAO BETA E KW---- 

quant<-.5 # quantil
# matrix deresiduals# matrix de resultados
order<-matrix(NA,nrow = 16, ncol = 6) 
cont<-1

colnames(order)<-c("p", "q","barma.AIC", "karma.AIC", "barmax.AIC", "karmax.AIC")


for(i in 0:3){
  for(j in 0:3){
    
    barma<-summary(BARFIMA.fit(y_train,p=i,d=F,q=j,info=T,
                               start = list(phi = rep(0,i),
                                            theta = rep(0,j)),
                               report=F))
    
    karma1<-suppressWarnings(KARFIMA.fit(y_train,p=i,d=F,q=j,info=T,
                                         rho=quant,
                                         control = list(method="Nelder-Mead",stopcr=1e-2),
                                         report=F))
    karma<-summary(karma1)
    
    
    barmax<-summary(BARFIMA.fit(y_train,p=i,d=F,q=j,info=T,
                                start = list(phi = rep(0,i),
                                             theta = rep(0,j),
                                             beta = coef(lm(y_train~X+0))
                                ),
                                control = list(method="Nelder-Mead",stopcr=1e-2),
                                xreg = X,
                                report=F))
    
    karmax1<-suppressWarnings(KARFIMA.fit(y_train,p=i,d=F,q=j,info=T,
                                          xreg = X,rho=quant,
                                          control = list(method="Nelder-Mead",stopcr=1e-2),
                                          report=F))
    karmax<-summary(karmax1)
    
    if(karma1$convergence==1 || is.nan(karma$aic)==1) karma$aic=0
    if(karmax1$convergence==1 || is.nan(karmax$aic)==1) karmax$aic=0
    
    order[cont,]<-c(i,j,barma$aic,karma$aic,
                    barmax$aic,karmax$aic)
    cont<-cont+1
    
  }
}



order<-order[-1,]
print(order)

orbarma<-order[which(order[,3]==min(order[,3])),c(1:3)]
orkarma<-order[which(order[,4]==min(order[,4])),c(1:2,4)]
orbarmax<-order[which(order[,5]==min(order[,5])),c(1:2,5)]
orkarmax<-order[which(order[,6]==min(order[,6])),c(1:2,6)]

barma<-BARFIMA.fit(y_train,p=orbarma[1],d=F,q=orbarma[2],
                   info=T,report=F)
karma<-KARFIMA.fit(y_train,p=orkarma[1],d=F,q=orkarma[2],rho=quant,
                   control = list(method="Nelder-Mead",stopcr=1e-2),
                   info=T,report=F)

barmax<-BARFIMA.fit(y_train,p=orbarmax[1],d=F,q=orbarmax[2],
                    xreg=X,info=T,report=F)
karmax<-KARFIMA.fit(y_train,p=orkarmax[1],d=F,q=orkarmax[2],rho=quant,
                    control = list(method="Nelder-Mead",stopcr=1e-2),
                    xreg=X,info=T,report=F)


results_insample<-rbind(
  forecast::accuracy(fit_ugoarma$fitted, y_train),
  forecast::accuracy(barmax$fitted.values, y_train),
  forecast::accuracy(karmax$fitted.values, y_train),
  forecast::accuracy(a02$fitted, y_train),
  forecast::accuracy(barma$fitted.values, y_train),
  forecast::accuracy(karma$fitted.values, y_train),
  forecast::accuracy(a01$fitted, y_train)
)#[,c(3,2,5)]

rownames(results_insample)<-c("UGOARMA","BARMAX","KARMAX","a02",
                              "BARMA", "KARMA", "a01")



round(results_insample, 4)


### out-of-sample

barma_out<-BARFIMA.extract(yt=Y,
                           coefs = list(alpha = barma$coefficients[1], 
                                        phi= barma$coefficients[2:(orbarma[1]+1)], 
                                        theta = if(orbarma[2]==0) {NULL} else{
                                          barma$coefficients[(orbarma[1]+2):(orbarma[1]+1+orbarma[2])]},
                                        nu = barma$coefficients[(orbarma[1]+2+orbarma[2])])
)


karma_out<-KARFIMA.extract(yt=Y,
                           coefs = list(alpha = karma$coefficients[1], 
                                        phi= karma$coefficients[2:(orkarma[1]+1)],
                                        theta = if(orkarma[2]==0) {NULL} else{
                                          karma$coefficients[(orkarma[1]+2):(orkarma[1]+1+orkarma[2])]},
                                        nu = karma$coefficients[(orkarma[1]+2+orkarma[2])])
)


barmax_out<-BARFIMA.extract(yt=Y, xreg = X0,  
                            coefs = list(alpha = barmax$coefficients[1], 
                                         beta = barmax$coefficients[2:(nX+1)],
                                         phi= barmax$coefficients[(nX+2):(orbarmax[1]+nX+1)], 
                                         theta = if(orbarmax[2]==0) {NULL} else{
                                           barmax$coefficients[(orbarmax[1]+(nX+2)):(orbarmax[1]+nX+1+orbarmax[2])]},
                                         nu = barmax$coefficients[(orbarmax[1]+(nX+2)+orbarmax[2])]))


karmax_out<-KARFIMA.extract(yt=Y,xreg = X0,rho=quant,  
                            coefs = list(alpha = karmax$coefficients[1], 
                                         beta = karmax$coefficients[2:(nX+1)],
                                         phi= karmax$coefficients[(nX+2):(orkarmax[1]+nX+1)], 
                                         theta = if(orkarmax[2]==0) {NULL} else{
                                           karmax$coefficients[(orkarmax[1]+nX+2):(orkarmax[1]+nX+1+orkarmax[2])]},
                                         nu = karmax$coefficients[(orkarmax[1]+nX+2+orkarmax[2])]))





# FORA DA AMOSTRA SEM REGRESSÃO ------

ugoarma_out1<-KARFIMA.extract(yt=Y,rho=quant,
                              coefs = list(alpha=fit_ugoarma_sr$coeff[1],
                                           phi=fit_ugoarma_sr$coeff[2:(best_ugoarma_sr$p[1]+1)],
                                           theta=if(best_ugoarma_sr$q[1]==0){NULL}
                                           else{fit_ugoarma_sr$coeff[(best_ugoarma_sr$p[1]+2):(best_ugoarma_sr$p[1]+1+best_ugoarma_sr$q[1])]},
                                           nu=fit_ugoarma_sr$coeff[(best_ugoarma_sr$p[1]+2+best_ugoarma_sr$q[1])])
)


# FORA DA AMOSTRA COM REGRESSÃO ------

ugoarma_out2<-KARFIMA.extract(yt=Y,xreg = X0,rho=quant,
                              coefs = list(alpha=fit_ugoarma$coeff[1],
                                           beta=fit_ugoarma$coeff[2:(nX+1)],
                                           phi=fit_ugoarma$coeff[(nX+2):(best_ugoarma$p[1]+1+nX)],
                                           theta=if(best_ugoarma$q[1]==0){NULL}
                                           else{fit_ugoarma$coeff[(best_ugoarma$p[1]+2+nX):(best_ugoarma$p[1]+1+best_ugoarma$q[1]+nX)]},
                                           nu=fit_ugoarma$coeff[(best_ugoarma$p[1]+2+best_ugoarma$q[1]+nX)])
)







a<-n+1:length(y_test)





results_outsample<-rbind(
  forecast::accuracy(ugoarma_out2$mut[(n+1):length(y_test)],y_test),
  forecast::accuracy(barmax_out$mut[(n+1):length(y_test)],y_test ),
  forecast::accuracy(karmax_out$mut[(n+1):length(y_test)],y_test ),
  forecast::accuracy(new2$fitted,y_test),
  forecast::accuracy(ugoarma_out1$mut[(n+1):length(y_test)],y_test),
  forecast::accuracy(barma_out$mut[(n+1):length(y_test)],y_test ),
  forecast::accuracy(karma_out$mut[(n+1)]:length(y_test), y_test),
  forecast::accuracy(new1$fitted, y_test)
)[,c(3,2,5)]
# 
row.names(results_outsample)<-
  #row.names(results_insample)<-
  c("UGOMAX","BARMAX","KARMAX","ARIMAX",
    "UGO","BARMA","KARMA","ARIMA")

#xtable::xtable((results_outsample),digits=4)

print(round(results_outsample,4))



checkresiduals(fit_ugoarma$residuals)

hist(fit_ugoarma$residuals)

#Box.test(fit_ugoarma$residuals, lag = 1, type = "Ljung-Box")

shapiro.test(fit_ugoarma$residuals)

# library(nortest)
# ad.test(fit_ugoarma$residuals)


#which(fit_ugoarma$residuals < -3)

#time(y_train)[c(108, 110)]


