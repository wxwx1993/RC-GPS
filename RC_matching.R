#For GPS method using matching


library(nnet)
library("MASS")

RC.matching<-function(treat.specfic=treat.estimate,data=simulated.data){
  data$treat.use<-treat.specfic
  GPS_mod<-multinom(treat.use ~ cf1+cf2+cf3+cf4+cf5+cf6 , data)
  #GPS_mod
  data$GPS0<-GPS_mod$fitted.values[,1]
  data$GPS1<-GPS_mod$fitted.values[,2]
  data$GPS2<-GPS_mod$fitted.values[,3]
  
  data.X2<-data[which(data$treat.use==2),]
  data.X1<-data[which(data$treat.use==1),]
  data.X0<-data[which(data$treat.use==0),]
  
  E.Y0<-NULL;E.Y1<-NULL;E.Y2<-NULL
  for (i in 1:2000){
    E.Y0[i]<-data.X0[which.min(abs(data.X0$GPS0-data$GPS0[i])),1]
    E.Y1[i]<-data.X1[which.min(abs(data.X1$GPS1-data$GPS1[i])),1]
    E.Y2[i]<-data.X2[which.min(abs(data.X2$GPS2-data$GPS2[i])),1]
  }

  return(c(mean(E.Y1)-mean(E.Y0),mean(E.Y2)-mean(E.Y1)))
}

regression.calibration.matching<-function(model=RC.model,simulated.data=simulated.data){
  imputation.data<-simulated.data
  imputation.data$treat.estimate<-predict(model,simulated.data)
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

match.result1<-matrix(NA,ncol=2,nrow=1000)
match.result2<-matrix(NA,ncol=2,nrow=1000)
match.crude<-matrix(NA,ncol=2,nrow=1000)
match.noerror<-matrix(NA,ncol=2,nrow=1000)

for (i in 1:1000){
  simulated.data<-data.generate(sample_size=2000,seed=i,sd=20,phi=0.7)
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  RC.model1<-lm(treat~ treat.w, validation)
  match.result1[i,]<-regression.calibration.matching(model=RC.model1,simulated.data=simulated.data)[[1]]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  gamma1[i]<-RC.model2$coefficients[2]
  match.result2[i,]<-regression.calibration.matching(model=RC.model2,simulated.data=simulated.data)[[1]]
  
  match.noerror[i,]<-RC.matching(treat.specfic=simulated.data$treat.cat,data=simulated.data)
  match.crude[i,]<-RC.matching(treat.specfic=simulated.data$treat.w.cat,data=simulated.data)
}

matching.beta1<-cbind(match.noerror[,1],match.crude[,1],match.result1[,1],match.result2[,1])
colnames(matching.beta1)<-c("error-free","error-prone","X~W","X~W+D")
matching.beta2<-cbind(match.noerror[,2],match.crude[,2],match.result1[,2],match.result2[,2])
colnames(matching.beta2)<-c("error-free","error-prone","X~W","X~W+D")

pdf("ATE.beta1_matching.pdf",width=11,height=8.5)
boxplot(matching.beta1,main="Matching, E(Y2)-E(Y1)",ylim=c(10,25),cex.axis=2.2,cex.main=2.5)
abline(h=21.23,lty=5,lwd=3,col="red")
dev.off()

pdf("ATE.beta2_matching.pdf",width=11,height=8.5)
boxplot(matching.beta2,main="Matching, E(Y3)-E(Y2)",ylim=c(10,25),cex.axis=2.2,cex.main=2.5)
abline(h=20.13,lty=5,lwd=3,col="red")
dev.off()


