#For GPS method using IPTW, truncated GPS if wieghts is >10
library(nnet)
library("MASS")

##################the function to implement IPTW approach
RC.IPTW.trun<-function(treat.specfic=treat.estimate,data=simulated.data){
  data$treat.use<-treat.specfic
  if (length(table(data$treat.use))<3){return(c(NA,NA))}
  else{
  GPS_mod<-multinom(treat.use ~ cf1+cf2+cf3+cf4+cf5+cf6 , data)
  #GPS_mod
  data$GPS0<-GPS_mod$fitted.values[,1]
  data$GPS1<-GPS_mod$fitted.values[,2]
  data$GPS2<-GPS_mod$fitted.values[,3]


  data.X2<-data[which(data$treat.use==2),]
  data.X1<-data[which(data$treat.use==1),]
  data.X0<-data[which(data$treat.use==0),]
  
  if (sum(data.X2$GPS2<0.1)>0){data.X2[which(data.X2$GPS2<0.1),]$GPS2<-0.1}
  if (sum(data.X1$GPS1<0.1)>0){data.X1[which(data.X1$GPS1<0.1),]$GPS1<-0.1}
  if (sum(data.X0$GPS0<0.1)>0){data.X0[which(data.X0$GPS0<0.1),]$GPS0<-0.1}
  
  E.Y2<-sum(data.X2$Y*(1/data.X2$GPS2))/nrow(data)
  E.Y1<-sum(data.X1$Y*(1/data.X1$GPS1))/nrow(data)
  E.Y0<-sum(data.X0$Y*(1/data.X0$GPS0))/nrow(data)
  }
  return(c(E.Y1-E.Y0,E.Y2-E.Y1))
}

################################use RC model to calibrate exposures, and implement RC.IPTW on calibrated exposure.
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

############function to simulate 7 scenarios to vary e.g. the strength of confounders, the variances, specification of RC model.
############output are bias and MSE to assess the performance of RC-GPS under different scenarios
IPTW.trun1.fun<-function(scenarios=1){

IPTW.result1<-matrix(NA,ncol=2,nrow=1000)
IPTW.result2<-matrix(NA,ncol=2,nrow=1000)
IPTW.result3<-matrix(NA,ncol=2,nrow=1000)
IPTW.crude<-matrix(NA,ncol=2,nrow=1000)
IPTW.noerror<-matrix(NA,ncol=2,nrow=1000)
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
  IPTW.result1[i,]<-regression.calibration.IPTW.trun(model=RC.model1,simulated.data=simulated.data)[[1]]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  IPTW.result2[i,]<-regression.calibration.IPTW.trun(model=RC.model2,simulated.data=simulated.data)[[1]]
  
  if (scenarios==7){
    RC.model3<-lm(treat~ treat.w + I(treat.w^2) +cova1+cova2+cova3 , validation)
    IPTW.result3[i,]<-regression.calibration.IPTW.trun(model=RC.model3,simulated.data=simulated.data)[[1]]
  }
  
  IPTW.noerror[i,]<-RC.IPTW.trun(treat.specfic=simulated.data$treat.cat,data=simulated.data)
  IPTW.crude[i,]<-RC.IPTW.trun(treat.specfic=simulated.data$treat.w.cat,data=simulated.data)
}

if (scenarios!=7){
IPTW.beta1<-cbind(IPTW.noerror[,1],IPTW.crude[,1],IPTW.result1[,1],IPTW.result2[,1])
colnames(IPTW.beta1)<-c("error-free","error-prone","X~W","X~W+D")
IPTW.beta2<-cbind(IPTW.noerror[,2],IPTW.crude[,2],IPTW.result1[,2],IPTW.result2[,2])
colnames(IPTW.beta2)<-c("error-free","error-prone","X~W","X~W+D")
}

if (scenarios==7){
  IPTW.beta1<-cbind(IPTW.noerror[,1],IPTW.crude[,1],IPTW.result1[,1],IPTW.result2[,1],IPTW.result3[,1])
  colnames(IPTW.beta1)<-c("error-free","error-prone","X~W","X~W+D","X~W2+D")
  IPTW.beta2<-cbind(IPTW.noerror[,2],IPTW.crude[,2],IPTW.result1[,2],IPTW.result2[,2],IPTW.result3[,2])
  colnames(IPTW.beta2)<-c("error-free","error-prone","X~W","X~W+D","X~W2+D")
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
  pdf(paste0("ATE.beta1_IPTW_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(IPTW.beta1,main="IPTW, E(Y2)-E(Y1)",cex.axis=3,cex.main=4.5,ylim=c(16,25))
  abline(h=gold_a,lty=5,lwd=3,col="red")
  dev.off()
  pdf(paste0("ATE.beta2_IPTW_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(IPTW.beta2,main="IPTW, E(Y3)-E(Y2)",cex.axis=3,cex.main=4.5,ylim=c(16,25))
  abline(h=gold_b,lty=5,lwd=3,col="red")
  dev.off()
}else if(scenarios==2){
  pdf(paste0("ATE.beta1_IPTW_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(IPTW.beta1,main="IPTW, E(Y2)-E(Y1)",cex.axis=3,cex.main=4.5,ylim=c(3,18))
  abline(h=gold_a,lty=5,lwd=3,col="red")
  dev.off()
  pdf(paste0("ATE.beta2_IPTW_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(IPTW.beta2,main="IPTW, E(Y3)-E(Y2)",cex.axis=3,cex.main=4.5,ylim=c(-10,40))
  abline(h=gold_b,lty=5,lwd=3,col="red")
  dev.off()
}else if(scenarios==4){
  pdf(paste0("ATE.beta1_IPTW_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(IPTW.beta1,main="IPTW, E(Y2)-E(Y1)",cex.axis=3,cex.main=4.5,ylim=c(7,13))
  abline(h=gold_a,lty=5,lwd=3,col="red")
  dev.off()
  pdf(paste0("ATE.beta2_IPTW_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(IPTW.beta2,main="IPTW, E(Y3)-E(Y2)",cex.axis=3,cex.main=4.5,ylim=c(7,13))
  abline(h=gold_b,lty=5,lwd=3,col="red")
  dev.off()
}else if(scenarios==7){
  pdf(paste0("ATE.beta1_IPTW_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(IPTW.beta1,main="IPTW, E(Y2)-E(Y1)",cex.axis=2.5,cex.main=4.5,ylim=c(-50,20))
  abline(h=gold_a,lty=5,lwd=3,col="red")
  dev.off()
  pdf(paste0("ATE.beta2_IPTW_",scenarios,".pdf"),width=14,height=8.5)
  boxplot(IPTW.beta2,main="IPTW, E(Y3)-E(Y2)",cex.axis=2.5,cex.main=4.5,ylim=c(10,65))
  abline(h=gold_b,lty=5,lwd=3,col="red")
  dev.off()
}

if (scenarios!=7){
  MSE.beta1<-c(mean((IPTW.beta1[,1]-gold_a)^2,na.rm = T),
               mean((IPTW.beta1[,2]-gold_a)^2,na.rm = T),
               mean((IPTW.beta1[,3]-gold_a)^2,na.rm = T),
               mean((IPTW.beta1[,4]-gold_a)^2,na.rm = T))
  MSE.beta2<-c(mean((IPTW.beta2[,1]-gold_b)^2,na.rm = T),
               mean((IPTW.beta2[,2]-gold_b)^2,na.rm = T),
               mean((IPTW.beta2[,3]-gold_b)^2,na.rm = T),
               mean((IPTW.beta2[,4]-gold_b)^2,na.rm = T))
  bias.beta1<-c((mean(IPTW.beta1[,1],na.rm = T)-gold_a)/gold_a*100,
                (mean(IPTW.beta1[,2],na.rm = T)-gold_a)/gold_a*100,
                (mean(IPTW.beta1[,3],na.rm = T)-gold_a)/gold_a*100,
                (mean(IPTW.beta1[,4],na.rm = T)-gold_a)/gold_a*100)
  bias.beta2<-c((mean(IPTW.beta2[,1],na.rm = T)-gold_b)/gold_b*100,
                (mean(IPTW.beta2[,2],na.rm = T)-gold_b)/gold_b*100,
                (mean(IPTW.beta2[,3],na.rm = T)-gold_b)/gold_b*100,
                (mean(IPTW.beta2[,4],na.rm = T)-gold_b)/gold_b*100)
}

if (scenarios==7){
  MSE.beta1<-c(mean((IPTW.beta1[,1]-gold_a)^2,na.rm = T),
               mean((IPTW.beta1[,2]-gold_a)^2,na.rm = T),
               mean((IPTW.beta1[,3]-gold_a)^2,na.rm = T),
               mean((IPTW.beta1[,4]-gold_a)^2,na.rm = T),
               mean((IPTW.beta1[,5]-gold_a)^2,na.rm = T))
  MSE.beta2<-c(mean((IPTW.beta2[,1]-gold_b)^2,na.rm = T),
               mean((IPTW.beta2[,2]-gold_b)^2,na.rm = T),
               mean((IPTW.beta2[,3]-gold_b)^2,na.rm = T),
               mean((IPTW.beta2[,4]-gold_b)^2,na.rm = T),
               mean((IPTW.beta2[,5]-gold_b)^2,na.rm = T))
  bias.beta1<-c((mean(IPTW.beta1[,1],na.rm = T)-gold_a)/gold_a*100,
                (mean(IPTW.beta1[,2],na.rm = T)-gold_a)/gold_a*100,
                (mean(IPTW.beta1[,3],na.rm = T)-gold_a)/gold_a*100,
                (mean(IPTW.beta1[,4],na.rm = T)-gold_a)/gold_a*100,
                (mean(IPTW.beta1[,5],na.rm = T)-gold_a)/gold_a*100)
  bias.beta2<-c((mean(IPTW.beta2[,1],na.rm = T)-gold_b)/gold_b*100,
                (mean(IPTW.beta2[,2],na.rm = T)-gold_b)/gold_b*100,
                (mean(IPTW.beta2[,3],na.rm = T)-gold_b)/gold_b*100,
                (mean(IPTW.beta2[,4],na.rm = T)-gold_b)/gold_b*100,
                (mean(IPTW.beta2[,5],na.rm = T)-gold_b)/gold_b*100)
}

return(list(MSE.beta1,MSE.beta2,bias.beta1,bias.beta2))

}

