library(nnet)
library("MASS")

#function for generating simulated data, including both "true" and "error-prone" exposurs
data.generate<-function(sample_size=2000,seed=300,sd=20,phi=0.8,sd_rc=1,tau_par=0.6,beta_1=1,beta_par=1,correct_RC=1){
  
  options(digits=4) # only print 4 sig digits
  set.seed(seed)
  size<-sample_size
  #pre-treatment variables (confounders)
  cf1<-mvrnorm(n=size,mu=c(0,0,0),Sigma=matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1),ncol=3))
  cf2<-sample(c((-2):2),size,replace = T)
  cf3<-runif(size,min=-3,max=3)
  cf4<-rchisq(size,1)
  cf<-cbind(cf0=1,cf1,cf2,cf3,cf4)
  #coefficients of linear regression model
  tau<-tau_par*c(1,1,2,1.5,3,2,3) 
  
  # covariates in the RC model
  cova1<-cf1[,1]
  cova2<-rnorm(size,0,sd=4)
  cova3<-runif(size,min=-5,max=5)
  cova<-cbind(cova1,cova2,cova3)
  zeta<-c(2,1,3,phi) 
  
  treat.w<-as.numeric()
  treat<-as.numeric()
  for (i in 1:size){
    treat.w[i]<-sum(cf[i,]*tau)+rnorm(1,mean=0,sd=sd)
    if (correct_RC==1){
    treat[i]<-sum(c(cova[i,],treat.w[i])*zeta) +rnorm(1,mean=0,sd=sd_rc) 
    }else if (correct_RC==0){
    treat[i]<-sum(c(cova[i,],treat.w[i])*zeta)+treat.w[i]^2/20+rnorm(1,mean=0,sd=sd_rc)
    }
    #have other confounders in measurement error model 
  }
  #coefficients of linear model
  beta<-c(beta_1,c(3,2,1,4,2,1)*beta_par)
  #produce outcome Y
  Y<-as.numeric()
  for (i in 1:size){
    Y[i]<-sum(c(treat[i],cf[i,-1])*beta) +rnorm(1, mean = 0, sd = 1)
  }
  treat.cat<-as.numeric()
  for (i in 1:size){
    if (treat[i]>15){treat.cat[i]<-2}
    else if (treat[i]<=15 & treat[i]>-5){treat.cat[i]<-1}
    else if (treat[i]<=-5){treat.cat[i]<-0}
  }
  treat.w.cat<-as.numeric()
  for (i in 1:size){
    if (treat.w[i]>15){treat.w.cat[i]<-2}
    else if (treat.w[i]<=15 & treat.w[i]>-5){treat.w.cat[i]<-1}
    else if (treat.w[i]<=-5){treat.w.cat[i]<-0}
  }
  simulated.data<-data.frame(cbind(Y,treat,treat.cat,cf[,-1],cova,treat.w,treat.w.cat))
  colnames(simulated.data)[4:9]<-c("cf1","cf2","cf3","cf4","cf5","cf6")
  return(simulated.data)
}


simulated.data<-data.generate(sample_size=10000,phi=0.8,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=1)
var(simulated.data$treat)
var(simulated.data$treat.w)
cor(simulated.data$treat,simulated.data$treat.w)

summary(lm(treat~ treat.w  +cova1+cova2+cova3 , simulated.data))

#cor(simulated.data$treat.w,rowSums(simulated.data[,4:9]*0.6*c(1,1,2,1.5,3,2,3)))
table(simulated.data$treat.cat,simulated.data$treat.w.cat)


##################function for generating simulated data to pretend we have a misspecified terms in RC model.
data.generate.miss<-function(sample_size=2000,seed=300,sd=20,phi=1,sd_rc=5,gam_par=0.6,zeta_par=1,beta_par=1){
  
  options(digits=4) # only print 4 sig digits
  set.seed(seed)
  size<-sample_size
  #pre-treatment variables (confounders)
  cf1<-mvrnorm(n=size,mu=c(0,0,0),Sigma=matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1),ncol=3))
  cf2<-sample(c((-2):2),size,replace = T)
  cf3<-runif(size,min=-3,max=3)
  cf4<-rchisq(size,1)
  cf<-cbind(cf0=1,cf1,cf2,cf3,cf4)
  #coefficients of linear regression model
  gam<-gam_par*c(1,1,2,1.5,3,2,3) 
  
  # covariates in the RC model
  cova1<-cf1[,1]
  cova2<-rnorm(size,0,sd=4)
  cova3<-runif(size,min=-5,max=5)
  cova<-cbind(cova1,cova2,cova3)
  zeta<-zeta_par*c(2,1,3,phi) 
  
  treat.w<-as.numeric()
  treat<-as.numeric()
  for (i in 1:size){
    treat.w[i]<-sum(cf[i,]*gam)+rnorm(1,mean=0,sd=sd)
    treat[i]<-sum(c(cova[i,],treat.w[i])*zeta)-5*log(abs(cova2[i]))+5*exp(rnorm(1,0,1))-3*cova1[i]^2-20*sin(cova3[i])+rnorm(1,mean=0,sd=sd_rc) 
    #have other confounders in measurement error model 
  }
  #coefficients of linear model
  beta<-c(beta_1,c(3,2,1,4,2,1)*beta_par)
  #produce outcome Y
  Y<-as.numeric()
  for (i in 1:size){
    Y[i]<-sum(c(treat[i],cf[i,-1])*beta) +rnorm(1, mean = 0, sd = 1)
  }
  treat.cat<-as.numeric()
  for (i in 1:size){
    if (treat[i]>15){treat.cat[i]<-2}
    else if (treat[i]<=15 & treat[i]>-5){treat.cat[i]<-1}
    else if (treat[i]<=-5){treat.cat[i]<-0}
  }
  treat.w.cat<-as.numeric()
  for (i in 1:size){
    if (treat.w[i]>15){treat.w.cat[i]<-2}
    else if (treat.w[i]<=15 & treat.w[i]>-5){treat.w.cat[i]<-1}
    else if (treat.w[i]<=-5){treat.w.cat[i]<-0}
  }
  simulated.data<-data.frame(cbind(Y,treat,treat.cat,cf[,-1],cova,treat.w,treat.w.cat))
  colnames(simulated.data)[4:9]<-c("cf1","cf2","cf3","cf4","cf5","cf6")
  return(simulated.data)
}
