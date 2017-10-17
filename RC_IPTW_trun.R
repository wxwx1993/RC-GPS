#For GPS method using IPTW, truncated lowest 1% and highest 1% GPS

library(nnet)
library("MASS")

RC.IPTW.trun<-function(treat.specfic=treat.estimate,data=simulated.data){
  data$treat.use<-treat.specfic
  GPS_mod<-multinom(treat.use ~ cf1+cf2+cf3+cf4+cf5+cf6 , data)
  #GPS_mod
  data$GPS0<-GPS_mod$fitted.values[,1]
  data$GPS1<-GPS_mod$fitted.values[,2]
  data$GPS2<-GPS_mod$fitted.values[,3]
  
  GPS0.trun<-quantile(data$GPS0,c(0.01,0.99))
  GPS1.trun<-quantile(data$GPS1,c(0.01,0.99))
  GPS2.trun<-quantile(data$GPS2,c(0.01,0.99))

  data.X2<-data[which(data$treat.use==2 & data$GPS2>GPS2.trun[1] & data$GPS2<GPS2.trun[2]),]
  data.X1<-data[which(data$treat.use==1 & data$GPS1>GPS1.trun[1] & data$GPS1<GPS1.trun[2]),]
  data.X0<-data[which(data$treat.use==0 & data$GPS0>GPS0.trun[1] & data$GPS0<GPS0.trun[2]),]

  E.Y2<-sum(data.X2$Y*(1/data.X2$GPS2))/(nrow(data.X2)+nrow(data.X1)+nrow(data.X0))
  E.Y1<-sum(data.X1$Y*(1/data.X1$GPS1))/(nrow(data.X2)+nrow(data.X1)+nrow(data.X0))
  E.Y0<-sum(data.X0$Y*(1/data.X0$GPS0))/(nrow(data.X2)+nrow(data.X1)+nrow(data.X0))
  
  #E.Y2<-sum(data.X2$Y*(1/data.X2$GPS2))/nrow(data)
  #E.Y1<-sum(data.X1$Y*(1/data.X1$GPS1))/nrow(data)
  #E.Y0<-sum(data.X0$Y*(1/data.X0$GPS0))/nrow(data)
  
  return(c(E.Y1-E.Y0,E.Y2-E.Y1))
}

regression.calibration.IPTW.trun<-function(model=RC.model,simulated.data=simulated.data){
  imputation.data<-simulated.data
  imputation.data$treat.estimate<-predict(model,simulated.data)
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


IPTW.result1<-matrix(NA,ncol=2,nrow=1000)
IPTW.result2<-matrix(NA,ncol=2,nrow=1000)
IPTW.crude<-matrix(NA,ncol=2,nrow=1000)
IPTW.noerror<-matrix(NA,ncol=2,nrow=1000)

for (i in 1:1000){
  simulated.data<-data.generate(sample_size=2000,seed=i,sd=20,phi=0.7)
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  RC.model1<-lm(treat~ treat.w, validation)
  IPTW.result1[i,]<-regression.calibration.IPTW.trun(model=RC.model1,simulated.data=simulated.data)[[1]]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  gamma1[i]<-RC.model2$coefficients[2]
  IPTW.result2[i,]<-regression.calibration.IPTW.trun(model=RC.model2,simulated.data=simulated.data)[[1]]
  
  IPTW.noerror[i,]<-RC.IPTW.trun(treat.specfic=simulated.data$treat.cat,data=simulated.data)
  IPTW.crude[i,]<-RC.IPTW.trun(treat.specfic=simulated.data$treat.w.cat,data=simulated.data)
}


IPTW.beta1<-cbind(IPTW.noerror[,1],IPTW.crude[,1],IPTW.result1[,1],IPTW.result2[,1])
colnames(IPTW.beta1)<-c("error-free","error-prone","X~W","X~W+D")
IPTW.beta2<-cbind(IPTW.noerror[,2],IPTW.crude[,2],IPTW.result1[,2],IPTW.result2[,2])
colnames(IPTW.beta2)<-c("error-free","error-prone","X~W","X~W+D")

pdf("ATE.beta1_IPTW.pdf",width=11,height=8.5)
boxplot(IPTW.beta1,main="IPTW, E(Y2)-E(Y1)",ylim=c(10,25),cex.axis=2.2,cex.main=2.5)
abline(h=21.23,lty=5,lwd=3,col="red")
dev.off()
pdf("ATE.beta2_IPTW.pdf",width=11,height=8.5)
boxplot(IPTW.beta2,main="IPTW, E(Y3)-E(Y2)",ylim=c(10,25),cex.axis=2.2,cex.main=2.5)
abline(h=20.13,lty=5,lwd=3,col="red")
dev.off()
#boxplot(gamma1,main="coefficients of regression calibration, IPTW")



