library(readxl)
umidade_atacama <- read_excel("umidade_atacama.xlsx")
View(umidade_atacama)

attach(umidade_atacama)

umidade_atacama$atacama<-umidade_atacama$atacama/100
umidade_atacama$chile<-umidade_atacama$chile/100

umid_mensal<-umidade_atacama |> 
  