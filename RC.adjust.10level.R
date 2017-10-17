#For GPS method using subclassfication (10 levels)

library(nnet)
library("MASS")


RC.adjust.10level<-function(treat.specfic=treat.estimate,data=imputation.data){
  data$treat.use<-treat.specfic
  GPS_mod<-multinom(treat.use ~ cf1+cf2+cf3+cf4+cf5+cf6 , data)
  
  GPS0<-GPS_mod$fitted.values[,1]
  GPS1<-GPS_mod$fitted.values[,2]
  GPS2<-GPS_mod$fitted.values[,3]
  
  GPS0.quantile<-quantile(GPS0,seq(0.1,0.9,0.1))
  GPS1.quantile<-quantile(GPS1,seq(0.1,0.9,0.1))
  GPS2.quantile<-quantile(GPS2,seq(0.1,0.9,0.1))
  
  Y.GPS0<-list(data[which(GPS0>GPS0.quantile[9]),],
               data[which(GPS0<GPS0.quantile[9] & GPS0>GPS0.quantile[8]),],
               data[which(GPS0<GPS0.quantile[8] & GPS0>GPS0.quantile[7]),],
               data[which(GPS0<GPS0.quantile[7] & GPS0>GPS0.quantile[6]),],
               data[which(GPS0<GPS0.quantile[6] & GPS0>GPS0.quantile[5]),],
               data[which(GPS0<GPS0.quantile[5] & GPS0>GPS0.quantile[4]),],
               data[which(GPS0<GPS0.quantile[4] & GPS0>GPS0.quantile[3]),],
               data[which(GPS0<GPS0.quantile[3] & GPS0>GPS0.quantile[2]),],
               data[which(GPS0<GPS0.quantile[2] & GPS0>GPS0.quantile[1]),],
               data[which(GPS0<GPS0.quantile[1]),])
  
  Y.GPS1<-list(data[which(GPS1>GPS1.quantile[9]),],
               data[which(GPS1<GPS1.quantile[9] & GPS1>GPS1.quantile[8]),],
               data[which(GPS1<GPS1.quantile[8] & GPS1>GPS1.quantile[7]),],
               data[which(GPS1<GPS1.quantile[7] & GPS1>GPS1.quantile[6]),],
               data[which(GPS1<GPS1.quantile[6] & GPS1>GPS1.quantile[5]),],
               data[which(GPS1<GPS1.quantile[5] & GPS1>GPS1.quantile[4]),],
               data[which(GPS1<GPS1.quantile[4] & GPS1>GPS1.quantile[3]),],
               data[which(GPS1<GPS1.quantile[3] & GPS1>GPS1.quantile[2]),],
               data[which(GPS1<GPS1.quantile[2] & GPS1>GPS1.quantile[1]),],
               data[which(GPS1<GPS1.quantile[1]),])
  
  Y.GPS2<-list(data[which(GPS2>GPS2.quantile[9]),],
               data[which(GPS2<GPS2.quantile[9] & GPS2>GPS2.quantile[8]),],
               data[which(GPS2<GPS2.quantile[8] & GPS2>GPS2.quantile[7]),],
               data[which(GPS2<GPS2.quantile[7] & GPS2>GPS2.quantile[6]),],
               data[which(GPS2<GPS2.quantile[6] & GPS2>GPS2.quantile[5]),],
               data[which(GPS2<GPS2.quantile[5] & GPS2>GPS2.quantile[4]),],
               data[which(GPS2<GPS2.quantile[4] & GPS2>GPS2.quantile[3]),],
               data[which(GPS2<GPS2.quantile[3] & GPS2>GPS2.quantile[2]),],
               data[which(GPS2<GPS2.quantile[2] & GPS2>GPS2.quantile[1]),],
               data[which(GPS2<GPS2.quantile[1]),])
  
  #the average effect when treatment is 0, E(Y(0))
  E.Y0<-mean(c(mean(Y.GPS0[[1]][which(Y.GPS0[[1]]$treat.use==0),1]),
               mean(Y.GPS0[[2]][which(Y.GPS0[[2]]$treat.use==0),1]),
               mean(Y.GPS0[[3]][which(Y.GPS0[[3]]$treat.use==0),1]),
               mean(Y.GPS0[[4]][which(Y.GPS0[[4]]$treat.use==0),1]),
               mean(Y.GPS0[[5]][which(Y.GPS0[[5]]$treat.use==0),1]),
               mean(Y.GPS0[[6]][which(Y.GPS0[[6]]$treat.use==0),1]),
               mean(Y.GPS0[[7]][which(Y.GPS0[[7]]$treat.use==0),1]),
               mean(Y.GPS0[[8]][which(Y.GPS0[[8]]$treat.use==0),1]),
               mean(Y.GPS0[[9]][which(Y.GPS0[[9]]$treat.use==0),1]),
               mean(Y.GPS0[[10]][which(Y.GPS0[[10]]$treat.use==0),1])),na.rm = T)
  
  #the average effect when treatment is 1, E(Y(1))
  E.Y1<-mean(c(mean(Y.GPS1[[1]][which(Y.GPS1[[1]]$treat.use==1),1]),
               mean(Y.GPS1[[2]][which(Y.GPS1[[2]]$treat.use==1),1]),
               mean(Y.GPS1[[3]][which(Y.GPS1[[3]]$treat.use==1),1]),
               mean(Y.GPS1[[4]][which(Y.GPS1[[4]]$treat.use==1),1]),
               mean(Y.GPS1[[5]][which(Y.GPS1[[5]]$treat.use==1),1]),
               mean(Y.GPS1[[6]][which(Y.GPS1[[6]]$treat.use==1),1]),
               mean(Y.GPS1[[7]][which(Y.GPS1[[7]]$treat.use==1),1]),
               mean(Y.GPS1[[8]][which(Y.GPS1[[8]]$treat.use==1),1]),
               mean(Y.GPS1[[9]][which(Y.GPS1[[9]]$treat.use==1),1]),
               mean(Y.GPS1[[10]][which(Y.GPS1[[10]]$treat.use==1),1])),na.rm = T)
  
  #the average effect when treatment is 2, E(Y(2))
  E.Y2<-mean(c(mean(Y.GPS2[[1]][which(Y.GPS2[[1]]$treat.use==2),1]),
               mean(Y.GPS2[[2]][which(Y.GPS2[[2]]$treat.use==2),1]),
               mean(Y.GPS2[[3]][which(Y.GPS2[[3]]$treat.use==2),1]),
               mean(Y.GPS2[[4]][which(Y.GPS2[[4]]$treat.use==2),1]),
               mean(Y.GPS2[[5]][which(Y.GPS2[[5]]$treat.use==2),1]),
               mean(Y.GPS2[[6]][which(Y.GPS2[[6]]$treat.use==2),1]),
               mean(Y.GPS2[[7]][which(Y.GPS2[[7]]$treat.use==2),1]),
               mean(Y.GPS2[[8]][which(Y.GPS2[[8]]$treat.use==2),1]),
               mean(Y.GPS2[[9]][which(Y.GPS2[[9]]$treat.use==2),1]),
               mean(Y.GPS2[[10]][which(Y.GPS2[[10]]$treat.use==2),1])),na.rm = T)
  
  
  return(c(E.Y1-E.Y0,E.Y2-E.Y1))
  
}

regression.calibration.10level<-function(model=RC.model,simulated.data=simulated.data){
  imputation.data<-simulated.data
  imputation.data$treat.estimate<-predict(model,simulated.data)
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

sub.result1<-matrix(NA,ncol=2,nrow=1000)
sub.result2<-matrix(NA,ncol=2,nrow=1000)
sub.crude<-matrix(NA,ncol=2,nrow=1000)
sub.noerror<-matrix(NA,ncol=2,nrow=1000)
gamma1<-NULL

for (i in 1:1000){
  simulated.data<-data.generate(sample_size=2000,seed=i,sd=20,phi=0.7)
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  RC.model1<-lm(treat~ treat.w, validation)
  sub.result1[i,]<-regression.calibration.10level(model=RC.model1,simulated.data=simulated.data)[[1]]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  gamma1[i]<-RC.model2$coefficients[2]
  sub.result2[i,]<-regression.calibration.10level(model=RC.model2,simulated.data=simulated.data)[[1]]
  
  sub.noerror[i,]<-RC.adjust.10level(treat.specfic=simulated.data$treat.cat,data=simulated.data)
  sub.crude[i,]<-RC.adjust.10level(treat.specfic=simulated.data$treat.w.cat,data=simulated.data)
}

sub.beta1<-cbind(sub.noerror[,1],sub.crude[,1],sub.result1[,1],sub.result2[,1])
colnames(sub.beta1)<-c("error-free","error-prone","X~W","X~W+D")
sub.beta2<-cbind(sub.noerror[,2],sub.crude[,2],sub.result1[,2],sub.result2[,2])
colnames(sub.beta2)<-c("error-free","error-prone","X~W","X~W+D")

pdf("ATE.beta1_subclass.pdf",width=11,height=8.5)
boxplot(sub.beta1,main="Subclass, E(Y2)-E(Y1)",ylim=c(12,25),cex.axis=2.2,cex.main=2.5)
abline(h=21.23,lty=5,lwd=3,col="red")
dev.off()

pdf("ATE.beta2_subclass.pdf",width=11,height=8.5)
boxplot(sub.beta2,main="Subclass, E(Y3)-E(Y2)",ylim=c(12,25),cex.axis=2.2,cex.main=2.5)
abline(h=20.13,lty=5,lwd=3,col="red")
dev.off()
#boxplot(gamma1,main="coefficients of regression calibration, Subclassification")



#bias and MSE
(mean(sub.beta1[,2])-21.23)/21.23


mean((sub.beta1[,2]-21.23)^2)

(mean(sub.beta1[,3])-21.23)/21.23

mean((sub.beta1[,3]-21.23)^2)

(mean(sub.beta2[,2])-20.13)/20.13

(mean(sub.beta2[,3])-20.13)/20.13





