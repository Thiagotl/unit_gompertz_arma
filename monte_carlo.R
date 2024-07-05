# # Checking the results
# Case 2: with regressors

phi=0.2;theta=0.4; alpha=1;sigma=6; tau=0.5
true_values=c(1,0.2,0.4,6)
n=50
R<-100
mu_result<-matrix(NA,R,4)
for (i in 1:R) {

  y<-simu.ugoarma(n,phi=phi,theta=theta, alpha=alpha,sigma=sigma, tau=tau,freq=12,
                    link="logit")
  fit1<-uGoarma.fit(y)
  mu_result[i,]<-fit1$coeff
}


mean_values<-c(apply(mu_result,2,mean))
b_values<-(true_values-mean_values)/true_values*100
eqm_values<-c(apply(mu_result,2,var))+(true_values-mean_values)^2
result1<- cbind(true_values,
                mean_values,
                b_values,
                eqm_values
)
colnames(result1)<-c("true value","mean","bias","eqm")
rownames(result1)<-c("alpha","phi","theta","sigma") 
print(round(result1,5))

