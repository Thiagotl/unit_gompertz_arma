library(forecast)
library(BTSR)
library(tidyverse)
library(e1071)
source("ugo_fit.R")
source("functions.R")
library(readxl)
library(lubridate)




dados <- read_excel("ipeadata[31-03-2025-03-21].xls")

dados$`Taxa de desemprego`<-dados$`Taxa de desemprego`/100

dados <- dados  |>
  mutate(
    taxa = as.numeric(gsub(",", ".", `Taxa de desemprego`)),
    data = parse_date_time(Data, orders = "ym")
  ) |>
  mutate(
    crise_politica = if_else(data >= ymd("2015-01-01") & data <= ymd("2017-03-01"), 1, 0),
    pandemia_covid = if_else(data >= ymd("2020-03-01") & data <= ymd("2022-03-01"), 1, 0),
    gov_bolsonaro = if_else(data >= ymd("2019-01-01") & data <= ymd("2022-12-31"), 1, 0)
  )

dados$data <- as.Date(parse_date_time(dados$data, orders = "ymd"))



ggplot(dados, aes(x = data, y = taxa)) +
  geom_line(color = "black") +
  
  annotate("rect", xmin = as.Date("2015-01-01"), xmax = as.Date("2017-03-01"), #2015-12-01,  2016-12-01
           ymin = -Inf, ymax = Inf, fill = "blue", alpha = 0.3) +
  
  annotate("rect", xmin = as.Date("2019-01-01"), xmax = as.Date("2022-12-31"),
           ymin = -Inf, ymax = Inf, fill = "green", alpha = 0.3) +
  
  annotate("rect", xmin = as.Date("2020-03-01"), xmax = as.Date("2022-03-01"),
           ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.3) +
  
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  
  labs(title = "Taxa de Desemprego no Brasil (2012–2025)",
       x = "Ano", y = "Taxa de Desemprego (%)") +
  theme_minimal()


#View(dados)

# SERIE COMPLETA
Y<-ts(dados$`Taxa de desemprego`, start = c(2012,3), frequency = 12)

m<-length(Y)
h1 <- 10
n<-m-h1

# SERIE DE TREINO
y_train<-ts(Y[1:n], start=c(2012,3), frequency = 12)

# SERIE DE TESTE 
y_test<-Y[(n+1):m] 

#tend_determ(Y)
raiz_unit(Y)
#sazonalidade(Y)

#### Matrizes de regressores ----

# TENDENCIA TEMPORAL 
t = 1:length(y_train)
t_hat = (n+1):(n+h1)

C     = cos(2*pi*t/12)
C_hat = cos(2*pi*t_hat/12)  

# CRISE POLITICA 

crise<-dados$crise_politica[1:length(y_train)]
crise_hat<-dados$crise_politica[length(y_train)+1:length(y_test)]

gov_bolsonaro <- dados$gov_bolsonaro[1:length(y_train)]
gov_bolsonaro_hat <- dados$gov_bolsonaro[length(y_train)+1:length(y_test)]

X = cbind(C, crise,gov_bolsonaro)   
X_hat = cbind(C_hat, crise_hat, gov_bolsonaro_hat)  

nX=dim(X)[2]
X0<-rbind(X,X_hat)

a01<-auto.arima(y_train)
new1<-Arima(y_test,model=a01) #one-step-ahead

a02<-auto.arima(y_train, xreg = X)
new2<-Arima(y_test,xreg = X_hat,model=a02) #one-step-ahead


#### APLICACAO MELHOR MODELO UGO ARMA----
source("best_ugoarma2.R")

pmax = 3
qmax = 3

y<-y_train
best_ugoarma<-best_ugo_2(y_train, pmax = pmax, qmax = qmax,
                nbest = 8, X=X, X_hat = X_hat) # 1 2 



fit_ugoarma<-uGoarma.fit(y_train, ar=1, ma=1:2, X=X, X_hat=X_hat)


checkresiduals(fit_ugoarma$residuals)

#fit_ugoarma<-uGoarma.fit(y_train, ar=c(1,2,4,5), ma=c(2,12), X=X, X_hat=X_hat)

#fit_ugoarma$model


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


checkresiduals(a02$residuals)

### out-of-sample

barma_out<-BARFIMA.extract(yt=Y,
                           coefs = list(alpha = barma$coefficients[1], 
                                        phi= barma$coefficients[2:(orbarma[1]+1)], 
                                        theta = if(orbarma[2]==0) {NULL} else{
                                          barma$coefficients[(orbarma[1]+2):(orbarma[1]+1+orbarma[2])]},
                                        nu = barma$coefficients[(orbarma[1]+2+orbarma[2])]))


karma_out<-KARFIMA.extract(yt=Y,
                           coefs = list(alpha = karma$coefficients[1], 
                                        phi= karma$coefficients[2:(orkarma[1]+1)],
                                        theta = if(orkarma[2]==0) {NULL} else{
                                          karma$coefficients[(orkarma[1]+2):(orkarma[1]+1+orkarma[2])]},
                                        nu = karma$coefficients[(orkarma[1]+2+orkarma[2])]))


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

ugoarma_out<-KARFIMA.extract(yt=Y,xreg = X0,rho=quant,  
                            coefs = list(alpha = karmax$coefficients[1], 
                                         beta = karmax$coefficients[2:(nX+1)],
                                         phi= karmax$coefficients[(nX+2):(orkarmax[1]+nX+1)], 
                                         theta = if(orkarmax[2]==0) {NULL} else{
                                           karmax$coefficients[(orkarmax[1]+nX+2):(orkarmax[1]+nX+1+orkarmax[2])]},
                                         nu = karmax$coefficients[(orkarmax[1]+nX+2+orkarmax[2])]))

ugoarma_out2<-KARFIMA.extract(yt=Y,xreg = X0,rho=quant,
                             coefs = list(alpha=fit_ugoarma$coeff[1],
                                          beta=fit_ugoarma$coeff[2:4],
                                          phi=fit_ugoarma$coeff[5],
                                          theta=fit_ugoarma$coeff[6:7],
                                          nu=fit_ugoarma$coeff[8]))






a<-n+1:length(y_test)


results_outsample<-rbind(
  forecast::accuracy(ugoarma_out2$mut[(n+1):length(y_test)],y_test),
  forecast::accuracy(barmax_out$mut[(n+1):length(y_test)],y_test ),
  forecast::accuracy(karmax_out$mut[(n+1):length(y_test)],y_test ),
  forecast::accuracy(new2$fitted,y_test),
  forecast::accuracy(barma_out$mut[(n+1):length(y_test)],y_test ),
  forecast::accuracy(karma_out$mut[(n+1)]:length(y_test), y_test),
  forecast::accuracy(new1$fitted, y_test)
)[,c(3,2,5)]
# 
row.names(results_outsample)<-
  row.names(results_insample)<-
  c("UGO","BARMAX","KARMAX","ARIMAX",
    "BARMA","KARMA","ARIMA")

xtable::xtable((results_outsample),digits=4)

round(results_outsample[,1:2],4)


