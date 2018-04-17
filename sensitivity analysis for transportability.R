#sensitivity analysis for if variance of regression calibration model coefficient is large

library(nnet)
library("MASS")

setwd(paste(getwd()))


source("data.generate.R")
source("RC_IPTW.R")
source("RC_subclassification.R")
source("RC_matching.R")


setwd(paste(getwd(),"/Simulation_result", sep=""))

gold_a = 22.56
gold_b = 21.50
#IPTW
regression.calibration.IPTW_trans<-function(model=RC.model,simulated.data=simulated.data,sd_ex=0){
  imputation.data<-simulated.data
  imputation.data$treat.estimate<-as.numeric((model$coefficients+c(0,rnorm(1,sd=sd_ex),0,0,0))%*%t(cbind(rep(1,2000),imputation.data[,c(13,10,11,12)])))
  treat.estimate.cat<-as.numeric()
  for (i in 1:nrow(imputation.data)){
    if (imputation.data$treat.estimate[i]>15){treat.estimate.cat[i]<-2}
    else if (imputation.data$treat.estimate[i]<=15 & imputation.data$treat.estimate[i]>(-5)){treat.estimate.cat[i]<-1}
    else if (imputation.data$treat.estimate[i]<=(-5)){treat.estimate.cat[i]<-0}
  }
  imputation.data$treat.estimate.cat<-treat.estimate.cat
  return(list(RC.IPTW.trun(treat.specfic=imputation.data$treat.estimate.cat,data=imputation.data),
              c(1)))
}

#transportablity for IPTW
IPTW.result.sd00<-matrix(NA,ncol=2,nrow=1000)
IPTW.result.sd01<-matrix(NA,ncol=2,nrow=1000)
IPTW.result.sd02<-matrix(NA,ncol=2,nrow=1000)
IPTW.result.sd03<-matrix(NA,ncol=2,nrow=1000)
IPTW.result.sd05<-matrix(NA,ncol=2,nrow=1000)

for (i in 1:1000){
  simulated.data<-data.generate(sample_size=2000,seed=i,sd=20,phi=0.8)
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  #gamma1[i]<-RC.model2$coefficients[2]
  
  IPTW.result.sd00[i,]<-regression.calibration.IPTW_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0)[[1]]
  IPTW.result.sd01[i,]<-regression.calibration.IPTW_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.1)[[1]]
  IPTW.result.sd02[i,]<-regression.calibration.IPTW_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.2)[[1]]
  IPTW.result.sd03[i,]<-regression.calibration.IPTW_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.3)[[1]]
  IPTW.result.sd05[i,]<-regression.calibration.IPTW_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.5)[[1]]
  
}

IPTW.sd.beta1<-cbind(IPTW.result.sd00[,1],IPTW.result.sd01[,1],IPTW.result.sd02[,1],IPTW.result.sd03[,1],IPTW.result.sd05[,1])
colnames(IPTW.sd.beta1)<-c("original","0.1 sigma","0.2 sigma","0.3 sigma","0.5 sigma")
IPTW.sd.beta2<-cbind(IPTW.result.sd00[,2],IPTW.result.sd01[,2],IPTW.result.sd02[,2],IPTW.result.sd03[,2],IPTW.result.sd05[,2])
colnames(IPTW.sd.beta2)<-c("original","0.1 sigma","0.2 sigma","0.3 sigma","0.5 sigma")


#########plot the box plots to assess the transportability assumption
pdf(paste0("beta_trans_IPTW.pdf"),width=14,height=8.5)
boxplot(IPTW.sd.beta1,main="RC, E(2)-E(1)",ylim=c(8,30), xaxt = "n",cex.axis=2.2,cex.main=2.5)
axis(side=1,at=1:5,labels=c("original",expression(paste("0.1+", sigma)),expression(paste("0.2+", sigma)),expression(paste("0.3+", sigma)),expression(paste("0.5+", sigma))),cex.axis=2.2)
abline(h=gold_a,lty=5,lwd=3,col="red")
boxplot(IPTW.sd.beta2,main="RC, E(3)-E(2)",ylim=c(8,30), xaxt = "n",cex.axis=2.2,cex.main=2.5)
axis(side=1,at=1:5,labels=c("original",expression(paste("0.1+", sigma)),expression(paste("0.2+", sigma)),expression(paste("0.3+", sigma)),expression(paste("0.5+", sigma))),cex.axis=2.2)
abline(h=gold_b,lty=5,lwd=3,col="red")
dev.off()
#sd(gamma1)



#transportablity for subclass
regression.calibration.10level_trans<-function(model=RC.model,simulated.data=simulated.data,sd_ex=0){
  imputation.data<-simulated.data
  imputation.data$treat.estimate<-as.numeric((model$coefficients+c(0,rnorm(1,sd=sd_ex),0,0,0))%*%t(cbind(rep(1,2000),imputation.data[,c(13,10,11,12)])))
  treat.estimate.cat<-as.numeric()
  for (i in 1:nrow(imputation.data)){
    if (imputation.data$treat.estimate[i]>15){treat.estimate.cat[i]<-2}
    else if (imputation.data$treat.estimate[i]<=15 & imputation.data$treat.estimate[i]>(-5)){treat.estimate.cat[i]<-1}
    else if (imputation.data$treat.estimate[i]<=(-5)){treat.estimate.cat[i]<-0}
  }
  imputation.data$treat.estimate.cat<-treat.estimate.cat
  return(list(RC.adjust.10level(treat.specfic=imputation.data$treat.estimate.cat,data=imputation.data),
              c(1)))
}
########
sub.result.sd00<-matrix(NA,ncol=2,nrow=1000)
sub.result.sd01<-matrix(NA,ncol=2,nrow=1000)
sub.result.sd02<-matrix(NA,ncol=2,nrow=1000)
sub.result.sd03<-matrix(NA,ncol=2,nrow=1000)
sub.result.sd05<-matrix(NA,ncol=2,nrow=1000)
for (i in 1:1000){
  simulated.data<-data.generate(sample_size=2000,seed=i,sd=20,phi=0.8)
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  #gamma1[i]<-RC.model2$coefficients[2]
  
  sub.result.sd00[i,]<-regression.calibration.10level_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0)[[1]]
  sub.result.sd01[i,]<-regression.calibration.10level_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.1)[[1]]
  sub.result.sd02[i,]<-regression.calibration.10level_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.2)[[1]]
  sub.result.sd03[i,]<-regression.calibration.10level_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.3)[[1]]
  sub.result.sd05[i,]<-regression.calibration.10level_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.5)[[1]]
  
}

sub.sd.beta1<-cbind(sub.result.sd00[,1],sub.result.sd01[,1],sub.result.sd02[,1],sub.result.sd03[,1],sub.result.sd05[,1])
colnames(sub.sd.beta1)<-c("original",expression(paste("0.1+", sigma)),"0.2 sigma","0.3 sigma","0.5 sigma")
sub.sd.beta2<-cbind(sub.result.sd00[,2],sub.result.sd01[,2],sub.result.sd02[,2],sub.result.sd03[,2],sub.result.sd05[,2])
colnames(sub.sd.beta2)<-c("original","0.1 sigma","0.2 sigma","0.3 sigma","0.5 sigma")

pdf(paste0("beta_trans_subclass.pdf"),width=14,height=8.5)
par(mfrow=c(1,2))
boxplot(sub.sd.beta1,main="RC, E(2)-E(1)",ylim=c(10,25), xaxt = "n",cex.axis=2.2,cex.main=2.5)
axis(side=1,at=1:5,labels=c("original",expression(paste("0.1+", sigma)),expression(paste("0.2+", sigma)),expression(paste("0.3+", sigma)),expression(paste("0.5+", sigma))),cex.axis=2.2)
abline(h=21.23,lty=5,lwd=3,col="red")
boxplot(sub.sd.beta2,main="RC, E(3)-E(2)",ylim=c(10,25), xaxt = "n",cex.axis=2.2,cex.main=2.5)
axis(side=1,at=1:5,labels=c("original",expression(paste("0.1+", sigma)),expression(paste("0.2+", sigma)),expression(paste("0.3+", sigma)),expression(paste("0.5+", sigma))),cex.axis=2.2)
abline(h=20.13,lty=5,lwd=3,col="red")
dev.off()
#sd(gamma1)

#transportablity for matching
regression.calibration.matching_trans<-function(model=RC.model,simulated.data=simulated.data,sd_ex=0){
  imputation.data<-simulated.data
  imputation.data$treat.estimate<-as.numeric((model$coefficients+c(0,rnorm(1,sd=sd_ex),0,0,0))%*%t(cbind(rep(1,2000),imputation.data[,c(13,10,11,12)])))
  treat.estimate.cat<-as.numeric()
  for (i in 1:nrow(imputation.data)){
    if (imputation.data$treat.estimate[i]>15){treat.estimate.cat[i]<-2}
    else if (imputation.data$treat.estimate[i]<=15 & imputation.data$treat.estimate[i]>(-5)){treat.estimate.cat[i]<-1}
    else if (imputation.data$treat.estimate[i]<=(-5)){treat.estimate.cat[i]<-0}
  }
  imputation.data$treat.estimate.cat<-treat.estimate.cat
  return(list(RC.matching(treat.specfic=imputation.data$treat.estimate.cat,data=imputation.data),
              c(1)))
}


matching.result.sd00<-matrix(NA,ncol=2,nrow=1000)
matching.result.sd01<-matrix(NA,ncol=2,nrow=1000)
matching.result.sd02<-matrix(NA,ncol=2,nrow=1000)
matching.result.sd03<-matrix(NA,ncol=2,nrow=1000)
matching.result.sd05<-matrix(NA,ncol=2,nrow=1000)
for (i in 1:1000){
  simulated.data<-data.generate(sample_size=2000,seed=i,sd=20,phi=0.8)
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  #gamma1[i]<-RC.model2$coefficients[2]
  
  matching.result.sd00[i,]<-regression.calibration.matching_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0)[[1]]
  matching.result.sd01[i,]<-regression.calibration.matching_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.1)[[1]]
  matching.result.sd02[i,]<-regression.calibration.matching_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.2)[[1]]
  matching.result.sd03[i,]<-regression.calibration.matching_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.3)[[1]]
  matching.result.sd05[i,]<-regression.calibration.matching_trans(model=RC.model2,simulated.data=simulated.data,sd_ex=0.5)[[1]]
}

matching.sd.beta1<-cbind(matching.result.sd00[,1],matching.result.sd01[,1],matching.result.sd02[,1],matching.result.sd03[,1],matching.result.sd05[,1])
#colnames(sub.sd.beta1)<-c("original",expression(paste("0.1+", sigma)),"0.2 sigma","0.3 sigma","0.5 sigma")
matching.sd.beta2<-cbind(matching.result.sd00[,2],matching.result.sd01[,2],matching.result.sd02[,2],matching.result.sd03[,2],matching.result.sd05[,2])
#colnames(sub.sd.beta2)<-c("original","0.1 sigma","0.2 sigma","0.3 sigma","0.5 sigma")

pdf(paste0("beta_trans_matching.pdf"),width=14,height=8.5)
par(mfrow=c(1,2))
boxplot(matching.sd.beta1,main="RC, E(2)-E(1)",ylim=c(10,25), xaxt = "n",cex.axis=2.2,cex.main=2.5)
axis(side=1,at=1:5,labels=c("original",expression(paste("0.1+", sigma)),expression(paste("0.2+", sigma)),expression(paste("0.3+", sigma)),expression(paste("0.5+", sigma))),cex.axis=2.2)
abline(h=21.23,lty=5,lwd=3,col="red")
boxplot(matching.sd.beta2,main="RC, E(3)-E(2)",ylim=c(10,25), xaxt = "n",cex.axis=2.2,cex.main=2.5)
axis(side=1,at=1:5,labels=c("original",expression(paste("0.1+", sigma)),expression(paste("0.2+", sigma)),expression(paste("0.3+", sigma)),expression(paste("0.5+", sigma))),cex.axis=2.2)
abline(h=20.13,lty=5,lwd=3,col="red")
dev.off()
#sd(gamma1)




