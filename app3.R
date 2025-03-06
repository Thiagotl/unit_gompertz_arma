#rm(list = ls())
library(forecast)
library(dplyr)
library(BTSR)
#source("ubxiiarma.fit.r")
source("ugo_fit.R")
######################
## Data preparation ##
######################
data <- readr::read_delim("combined_hourly_data.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE) %>% 
  mutate(timestamp=as.POSIXct(timestamp, tz="GMT",
                              origin="1970-01-01 00:00:00"),
         hum=hum/100)
#
# n<-round(dim(data)[1]*.8)-149
data<-data[150:1870,]
n<-round(dim(data)[1]*.8)
DM=NULL
#########################
## Train and test sets ##
#########################
datatrain<-cbind(data[1:n,])  
colnames(datatrain)<-paste0(colnames(datatrain),"_train")

datatest<-cbind(data[(n+1):(dim(data)[1]),])  
colnames(datatest)<-paste0(colnames(datatest),"_test")

datatrain$RSSI_03_train[is.na(datatrain$RSSI_03_train)]<-
  mean(na.omit(datatrain$RSSI_03_train))
datatrain$RSSI_04_train[is.na(datatrain$RSSI_04_train)]<-
  mean(na.omit(datatrain$RSSI_04_train))
datatrain$RSSI_05_train[is.na(datatrain$RSSI_05_train)]<-
  mean(na.omit(datatrain$RSSI_05_train))
datatrain$RSSI_06_train[is.na(datatrain$RSSI_06_train)]<-
  mean(na.omit(datatrain$RSSI_06_train))
datatrain$RSSI_08_train[is.na(datatrain$RSSI_08_train)]<-
  mean(na.omit(datatrain$RSSI_08_train))

suppressMessages(attach(datatrain))
suppressMessages(attach(datatest))
suppressMessages(attach(data))
# X<-as.matrix(datatrain[,6:13])
# Xtest<-as.matrix(datatest[,6:13])
X<-as.matrix(datatrain[,2])
Xtest<-as.matrix(datatest[,2])
X0<-rbind(X,Xtest)
nX<-dim(X)[2]
# X<-as.matrix(datatrain[,c(2,22)])
# Xtest<-as.matrix(datatest[,c(2,22)])
X0<-rbind(X,Xtest)
nX<-dim(X)[2]
######################
## stacionary tests ##
######################
truncation<-round(12*(n/100)^(.25)) # Schwert rule
adf.level<-tseries::adf.test(hum_train,k=truncation) # stacionary
adf.diff<-tseries::adf.test(diff(hum_train),k=truncation) # stacionary
kpss.level<-tseries::kpss.test(hum_train,lshort = F) # non-stacionary
kpss.diff<-tseries::kpss.test(diff(hum_train),lshort = F) # stacionary

table.stationarity<-data.frame(
  Series=c("In Level","1st difference"),
  ADF=c(adf.level$statistic,adf.diff$statistic),
  `p-value ADF`=c(adf.level$p.value,adf.diff$p.value),
  KPSS=c(kpss.level$statistic,kpss.diff$statistic),
  `p-value KPSS`=c(kpss.level$p.value,kpss.diff$p.value)
)

########################
## Fitting the models ##
########################
a01<-auto.arima(hum_train)
new1<-Arima(hum_test,model=a01) #one-step-ahead
forecast(a01, h = length(hum_test))

new1$fitted
a02<-auto.arima(hum_train, xreg = X)
new2<-Arima(hum_test,xreg = Xtest,model=a02) #one-step-ahead
# lmtest::coeftest(a02)
# xtable::xtable(lmtest::coeftest(a02)[,c(1,4)])
# xtable::xtable(summary(uwarmax)$coefficients[,c(1,4)])
quant<-.5
order<-matrix(NA,nrow = 16, ncol = 9) 
cont<-1
for(i in 0:3){
  for(j in 0:3){
    barma<-summary(BARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                               report=F))
    karma1<-suppressWarnings(KARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                                         rho=quant,
                                         control = list(method="Nelder-Mead",stopcr=1e-2),
                                         report=F))
    karma<-summary(karma1)
    uwarma1<-(UWARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,rho=quant,
                           report=F))
    uwarma<-summary(uwarma1)
    
    #ubxiiarma<-ubxiiarma.fit(ts(hum_test),ar=i,ma=i)
    
    ugoarma<-uGoarma.fit(hum_train,ar=i,ma=j)
    
    
    barmax<-summary(BARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                                xreg = X,
                                report=F))
    karmax1<-suppressWarnings(KARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                                          xreg = X,rho=quant,
                                          control = list(method="Nelder-Mead",stopcr=1e-2),
                                          report=F))
    karmax<-summary(karmax1)
    uwarmax1<-suppressWarnings(UWARFIMA.fit(hum_train,p=i,d=F,q=j,info=T,
                                            rho=quant, xreg = X,
                                            report=F))
    uwarmax<-summary(uwarmax1)
    if(karma1$convergence==1 || is.nan(karma$aic)==1) karma$aic=0
    if(karmax1$convergence==1 || is.nan(karmax$aic)==1) karmax$aic=0
    if(uwarmax1$convergence==1 || is.nan(uwarmax$aic)==1) uwarmax$aic=0
    #   print(c(karma1$convergence,karma$aic))
    order[cont,]<-c(i,j,barma$aic,karma$aic,uwarma$aic,
                    barmax$aic,karmax$aic,uwarmax$aic,ugoarma$aic)
    cont<-cont+1
  }
}
order<-order[-1,]
print(order)

orbarma<-order[which(order[,3]==min(order[,3])),c(1:3)]
orkarma<-order[which(order[,4]==min(order[,4])),c(1:2,4)]
oruwarma<-order[which(order[,5]==min(order[,5])),c(1:2,5)]
orbarmax<-order[which(order[,6]==min(order[,6])),c(1:2,6)]
orkarmax<-order[which(order[,7]==min(order[,7])),c(1:2,7)]
oruwarmax<-order[which(order[,8]==min(order[,8])),c(1:2,8)]

barma<-BARFIMA.fit(hum_train,p=orbarma[1],d=F,q=orbarma[2],
                   info=T,report=F)
karma<-KARFIMA.fit(hum_train,p=orkarma[1],d=F,q=orkarma[2],rho=quant,
                   control = list(method="Nelder-Mead",stopcr=1e-2),
                   info=T,report=F)
uwarma<-UWARFIMA.fit(hum_train,p=oruwarma[1],d=F,q=oruwarma[2],rho=quant,
                     info=T,report=F)
# ubxiiarma<-ubxiiarma.fit(ts(hum_train),ar=i,ma=i,tau = quant)
barmax<-BARFIMA.fit(hum_train,p=orbarmax[1],d=F,q=orbarmax[2],
                    xreg=X,info=T,report=F)
karmax<-KARFIMA.fit(hum_train,p=orkarmax[1],d=F,q=orkarmax[2],rho=quant,
                    control = list(method="Nelder-Mead",stopcr=1e-2),
                    xreg=X,info=T,report=F)
uwarmax<-UWARFIMA.fit(hum_train,p=oruwarmax[1],d=F,q=oruwarmax[2],rho=quant,
                      # control = list(method="Nelder-Mead"),
                      xreg=X,info=T,report=F)
results_insample<-rbind(
  forecast::accuracy(barmax$fitted.values, hum_train),
  forecast::accuracy(karmax$fitted.values, hum_train),
  forecast::accuracy(uwarmax$fitted.values, hum_train),
  forecast::accuracy(a02$fitted, hum_train),
  forecast::accuracy(barma$fitted.values, hum_train),
  forecast::accuracy(karma$fitted.values, hum_train),
  forecast::accuracy(uwarma$fitted.values, hum_train),
  forecast::accuracy(a01$fitted, hum_train)
)[,c(3,2,5)]

barma_out<-BARFIMA.extract(yt=hum,
                           coefs = list(alpha = barma$coefficients[1], 
                                        phi= barma$coefficients[2:(orbarma[1]+1)], 
                                        theta = if(orbarma[2]==0) {NULL} else{
                                          barma$coefficients[(orbarma[1]+2):(orbarma[1]+1+orbarma[2])]},
                                        nu = barma$coefficients[(orbarma[1]+2+orbarma[2])])
)

karma_out<-KARFIMA.extract(yt=hum,
                           coefs = list(alpha = karma$coefficients[1], 
                                        phi= karma$coefficients[2:(orkarma[1]+1)],
                                        theta = if(orkarma[2]==0) {NULL} else{
                                          karma$coefficients[(orkarma[1]+2):(orkarma[1]+1+orkarma[2])]},
                                        nu = karma$coefficients[(orkarma[1]+2+orkarma[2])])
)
uwarma_out<-UWARFIMA.extract(yt=hum,
                             coefs = list(alpha = uwarma$coefficients[1], 
                                          phi= uwarma$coefficients[2:(oruwarma[1]+1)],
                                          theta = if(oruwarma[2]==0) {NULL} else{
                                            uwarma$coefficients[(oruwarma[1]+2):(oruwarma[1]+1+oruwarma[2])]},
                                          nu = uwarma$coefficients[(oruwarma[1]+2+oruwarma[2])])
)
barmax_out<-BARFIMA.extract(yt=hum, xreg = X0,  
                            coefs = list(alpha = barmax$coefficients[1], 
                                         beta = barmax$coefficients[2:(nX+1)],
                                         phi= barmax$coefficients[(nX+2):(orbarmax[1]+nX+1)], 
                                         theta = if(orbarmax[2]==0) {NULL} else{
                                           barmax$coefficients[(orbarmax[1]+(nX+2)):(orbarmax[1]+nX+1+orbarmax[2])]},
                                         nu = barmax$coefficients[(orbarmax[1]+(nX+2)+orbarmax[2])])
)
karmax_out<-KARFIMA.extract(yt=hum,xreg = X0,rho=quant,  
                            coefs = list(alpha = karmax$coefficients[1], 
                                         beta = karmax$coefficients[2:(nX+1)],
                                         phi= karmax$coefficients[(nX+2):(orkarmax[1]+nX+1)], 
                                         theta = if(orkarmax[2]==0) {NULL} else{
                                           karmax$coefficients[(orkarmax[1]+nX+2):(orkarmax[1]+nX+1+orkarmax[2])]},
                                         nu = karmax$coefficients[(orkarmax[1]+nX+2+orkarmax[2])])
)
uwarmax_out<-UWARFIMA.extract(yt=hum,xreg = X0,rho=quant,  
                              coefs = list(alpha = uwarmax$coefficients[1],
                                           beta = uwarmax$coefficients[2:(nX+1)],
                                           phi= uwarmax$coefficients[(nX+2):(oruwarmax[1]+nX+1)],
                                           theta = if(oruwarmax[2]==0) {NULL} else{
                                             uwarmax$coefficients[(oruwarmax[1]+nX+2):(oruwarmax[1]+nX+1+oruwarmax[2])]},
                                           nu = uwarmax$coefficients[(oruwarmax[1]+nX+2+oruwarmax[2])])
)

results_outsample<-rbind(
  forecast::accuracy(barmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(karmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(uwarmax_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(new2$fitted, hum_test),
  forecast::accuracy(barma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(karma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(uwarma_out$mut[(n+1):(dim(data)[1])], hum_test),
  forecast::accuracy(new1$fitted, hum_test)
)[,c(3,2,5)]

row.names(results_outsample)<-
  row.names(results_insample)<-
  c("BARMAX","KARMAX","UWARMAX",
    "ARIMAX",
    "BARMA","KARMA","UWARMA","ARIMA")

xtable::xtable((results_outsample),digits=4)

round(results_outsample[,1:2],4)

