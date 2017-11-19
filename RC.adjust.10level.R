#For GPS method using subclassfication (10 levels)

library(nnet)
library("MASS")

RC.adjust.10level<-function(treat.specfic=treat.estimate,data=imputation.data){
  data$treat.use<-treat.specfic
  if (length(table(data$treat.use))<3){return(c(NA,NA))}
      else{
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


subclass.10level.fun<-function(scenarios=1){
  
sub.result1<-matrix(NA,ncol=2,nrow=1000)
sub.result2<-matrix(NA,ncol=2,nrow=1000)
sub.result3<-matrix(NA,ncol=2,nrow=1000)
sub.crude<-matrix(NA,ncol=2,nrow=1000)
sub.noerror<-matrix(NA,ncol=2,nrow=1000)
#gamma1<-NULL
########default
if (scenarios==1){text<-"phi=0.8,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=1"}
########correlation low
if (scenarios==2){text<-"phi=0.2,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=1"}
######## poor model fit
if (scenarios==3){text<-"phi=0.8,sd_rc=10,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=1"}
######## outcome effect small
if (scenarios==4){text<-"phi=0.8,sd_rc=1,tau_par=0.8,beta_1=0.5,beta_par=1,correct_RC=1"}
####### confounder on outcome large
if (scenarios==5){text<-"phi=0.8,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=5,correct_RC=1"}
####### confounder on exposure large
if (scenarios==6){text<-"phi=0.8,sd_rc=1,tau_par=1.6,beta_1=1,beta_par=1,correct_RC=1"}
####### RC model misspecified
if (scenarios==7){text<-"phi=0.8,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=0"}

for (i in 1:1000){
  simulated.data<-eval(parse(text=paste0("data.generate(sample_size=2000,seed = i,",text,")")))
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  #main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  RC.model1<-lm(treat~ treat.w, validation)
  sub.result1[i,]<-regression.calibration.10level(model=RC.model1,simulated.data=simulated.data)[[1]]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  #gamma1[i]<-RC.model2$coefficients[2]
  sub.result2[i,]<-regression.calibration.10level(model=RC.model2,simulated.data=simulated.data)[[1]]
  
  if (scenarios==7){
  RC.model3<-lm(treat~ treat.w + I(treat.w^2) +cova1+cova2+cova3 , validation)
  sub.result3[i,]<-regression.calibration.10level(model=RC.model3,simulated.data=simulated.data)[[1]]
  }
  
  sub.noerror[i,]<-RC.adjust.10level(treat.specfic=simulated.data$treat.cat,data=simulated.data)
  sub.crude[i,]<-RC.adjust.10level(treat.specfic=simulated.data$treat.w.cat,data=simulated.data)
}

if (scenarios!=7){
sub.beta1<-cbind(sub.noerror[,1],sub.crude[,1],sub.result1[,1],sub.result2[,1])
colnames(sub.beta1)<-c("error-free","error-prone","X~W","X~W+D")
sub.beta2<-cbind(sub.noerror[,2],sub.crude[,2],sub.result1[,2],sub.result2[,2])
colnames(sub.beta2)<-c("error-free","error-prone","X~W","X~W+D")
}

if (scenarios==7){
sub.beta1<-cbind(sub.noerror[,1],sub.crude[,1],sub.result1[,1],sub.result2[,1],sub.result3[,1])
colnames(sub.beta1)<-c("error-free","error-prone","X~W","X~W+D","X~W2+D")
sub.beta2<-cbind(sub.noerror[,2],sub.crude[,2],sub.result1[,2],sub.result2[,2],sub.result3[,2])
colnames(sub.beta2)<-c("error-free","error-prone","X~W","X~W+D","X~W2+D")
}

true.data1<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 100,",text,")")))
true.data2<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 200,",text,")")))
true.data3<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 300,",text,")")))
true.data4<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 400,",text,")")))
true.data5<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 500,",text,")")))
true.data6<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 600,",text,")")))
true.data7<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 700,",text,")")))
true.data8<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 800,",text,")")))
true.data9<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 900,",text,")")))
true.data10<-eval(parse(text=paste0("data.generate(sample_size=10000,seed = 1000,",text,")")))

true.data<-rbind(true.data1,true.data2,true.data3,true.data4,true.data5,true.data6,true.data7,true.data8,true.data9,true.data10)

gold_model<-summary(lm(Y~I(treat.cat==1)+I(treat.cat==2)+cf1+cf2+cf3+cf4+cf5+cf6,data=true.data))

gold_a<-gold_model$coefficients[2]
gold_b<-gold_model$coefficients[3]-gold_model$coefficients[2]

if(scenarios %in% c(1,3,5,6)){
pdf(paste0("ATE.beta1_subclass_",scenarios,".pdf"),width=14,height=8.5)
boxplot(sub.beta1,main="Subclass, E(Y2)-E(Y1)",cex.axis=3,cex.main=4.5,ylim=c(16,25))
abline(h=gold_a,lty=5,lwd=3,col="red")
dev.off()

pdf(paste0("ATE.beta2_subclass_",scenarios,".pdf"),width=14,height=8.5)
boxplot(sub.beta2,main="Subclass, E(Y3)-E(Y2)",cex.axis=3,cex.main=4.5,ylim=c(16,25))
abline(h=gold_b,lty=5,lwd=3,col="red")
dev.off()
}else if(scenarios==2){
  
  pdf(paste0("ATE.beta1_subclass_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(sub.beta1,main="Subclass, E(Y2)-E(Y1)",cex.axis=3,cex.main=4.5,ylim=c(3,18))
  abline(h=gold_a,lty=5,lwd=3,col="red")
  dev.off()
  
  pdf(paste0("ATE.beta2_subclass_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(sub.beta2,main="Subclass, E(Y3)-E(Y2)",cex.axis=3,cex.main=4.5,ylim=c(-10,40))
  abline(h=gold_b,lty=5,lwd=3,col="red")
  dev.off()
}else if(scenarios==4){
  pdf(paste0("ATE.beta1_subclass_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(sub.beta1,main="Subclass, E(Y2)-E(Y1)",cex.axis=3,cex.main=4.5,ylim=c(7,13))
  abline(h=gold_a,lty=5,lwd=3,col="red")
  dev.off()
  
  pdf(paste0("ATE.beta2_subclass_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(sub.beta2,main="Subclass, E(Y3)-E(Y2)",cex.axis=3,cex.main=4.5,ylim=c(7,13))
  abline(h=gold_b,lty=5,lwd=3,col="red")
  dev.off()
}else if(scenarios==7){
  pdf(paste0("ATE.beta1_subclass_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(sub.beta1,main="Subclass, E(Y2)-E(Y1)",cex.axis=2.5,cex.main=4.5,ylim=c(-50,20))
  abline(h=gold_a,lty=5,lwd=3,col="red")
  dev.off()
  
  pdf(paste0("ATE.beta2_subclass_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(sub.beta2,main="Subclass, E(Y3)-E(Y2)",cex.axis=2.5,cex.main=4.5,ylim=c(10,65))
  abline(h=gold_b,lty=5,lwd=3,col="red")
  dev.off()
}

if (scenarios!=7){
MSE.beta1<-c(mean((sub.beta1[,1]-gold_a)^2,na.rm = T),
             mean((sub.beta1[,2]-gold_a)^2,na.rm = T),
             mean((sub.beta1[,3]-gold_a)^2,na.rm = T),
             mean((sub.beta1[,4]-gold_a)^2,na.rm = T))
MSE.beta2<-c(mean((sub.beta2[,1]-gold_b)^2,na.rm = T),
             mean((sub.beta2[,2]-gold_b)^2,na.rm = T),
             mean((sub.beta2[,3]-gold_b)^2,na.rm = T),
             mean((sub.beta2[,4]-gold_b)^2,na.rm = T))
bias.beta1<-c((mean(sub.beta1[,1],na.rm = T)-gold_a)/gold_a*100,
              (mean(sub.beta1[,2],na.rm = T)-gold_a)/gold_a*100,
              (mean(sub.beta1[,3],na.rm = T)-gold_a)/gold_a*100,
              (mean(sub.beta1[,4],na.rm = T)-gold_a)/gold_a*100)
bias.beta2<-c((mean(sub.beta2[,1],na.rm = T)-gold_b)/gold_b*100,
              (mean(sub.beta2[,2],na.rm = T)-gold_b)/gold_b*100,
              (mean(sub.beta2[,3],na.rm = T)-gold_b)/gold_b*100,
              (mean(sub.beta2[,4],na.rm = T)-gold_b)/gold_b*100)
}

if (scenarios==7){
  MSE.beta1<-c(mean((sub.beta1[,1]-gold_a)^2,na.rm = T),
               mean((sub.beta1[,2]-gold_a)^2,na.rm = T),
               mean((sub.beta1[,3]-gold_a)^2,na.rm = T),
               mean((sub.beta1[,4]-gold_a)^2,na.rm = T),
               mean((sub.beta1[,5]-gold_a)^2,na.rm = T))
  MSE.beta2<-c(mean((sub.beta2[,1]-gold_b)^2,na.rm = T),
               mean((sub.beta2[,2]-gold_b)^2,na.rm = T),
               mean((sub.beta2[,3]-gold_b)^2,na.rm = T),
               mean((sub.beta2[,4]-gold_b)^2,na.rm = T),
               mean((sub.beta2[,5]-gold_b)^2,na.rm = T))
  bias.beta1<-c((mean(sub.beta1[,1],na.rm = T)-gold_a)/gold_a*100,
                (mean(sub.beta1[,2],na.rm = T)-gold_a)/gold_a*100,
                (mean(sub.beta1[,3],na.rm = T)-gold_a)/gold_a*100,
                (mean(sub.beta1[,4],na.rm = T)-gold_a)/gold_a*100,
                (mean(sub.beta1[,5],na.rm = T)-gold_a)/gold_a*100)
  bias.beta2<-c((mean(sub.beta2[,1],na.rm = T)-gold_b)/gold_b*100,
                (mean(sub.beta2[,2],na.rm = T)-gold_b)/gold_b*100,
                (mean(sub.beta2[,3],na.rm = T)-gold_b)/gold_b*100,
                (mean(sub.beta2[,4],na.rm = T)-gold_b)/gold_b*100,
                (mean(sub.beta2[,5],na.rm = T)-gold_b)/gold_b*100)
}

return(list(MSE.beta1,MSE.beta2,bias.beta1,bias.beta2))
}





