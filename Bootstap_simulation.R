sim_num<-10
boot_num<-100
#sub.result1<-matrix(NA,ncol=2,nrow=1000)
sub.result2<-matrix(NA,ncol=2,nrow=sim_num*boot_num)
sub.crude<-matrix(NA,ncol=2,nrow=sim_num*boot_num)
sub.noerror<-matrix(NA,ncol=2,nrow=sim_num*boot_num)
#gamma1<-NULL

for (i in 1:sim_num){
  simulated.data<-data.generate(sample_size=2000,seed=i,sd_rc=15)
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  #main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  for (boots in 1:boot_num){
  validation_boot<-validation[sample.int(500,500,replace = T),]
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation_boot)
  sub.result2[(i-1)*boot_num+boots,]<-regression.calibration.10level(model=RC.model2,simulated.data=simulated.data)[[1]]
  }
  
  sub.noerror[c(((i-1)*boot_num+1):(i*boot_num)),]<-matrix(rep(RC.adjust.10level(treat.specfic=simulated.data$treat.cat,data=simulated.data),boot_num),ncol=2,byrow =T)
  sub.crude[c(((i-1)*boot_num+1):(i*boot_num)),]<-matrix(rep(RC.adjust.10level(treat.specfic=simulated.data$treat.w.cat,data=simulated.data),boot_num),ncol=2,byrow=T)
}

sub.beta1<-cbind(sub.noerror[,1],sub.crude[,1],sub.result2[,1])
colnames(sub.beta1)<-c("error-free","error-prone","X~W+D")
sub.beta2<-cbind(sub.noerror[,2],sub.crude[,2],sub.result2[,2])
colnames(sub.beta2)<-c("error-free","error-prone","X~W+D")

pdf("ATE.beta1_subclass_boot.pdf",width=11,height=8.5)
boxplot(sub.beta1,main="Subclass, E(Y2)-E(Y1)",cex.axis=2.2,cex.main=2.5)
abline(h=gold_a,lty=5,lwd=3,col="red")
dev.off()

pdf("ATE.beta2_subclass_boot.pdf",width=11,height=8.5)
boxplot(sub.beta2,main="Subclass, E(Y3)-E(Y2)",cex.axis=2.2,cex.main=2.5)
abline(h=gold_b,lty=5,lwd=3,col="red")
dev.off()
#boxplot(gamma1,main="coefficients of regression calibration, Subclassification")
linear.result2<-NULL
linear.noerror<-NULL
linear.crude<-NULL
for (i in 1:100){
  simulated.data<-data.generate(sample_size=2000,seed=i,sd_rc=20)
  validation<-simulated.data[sample.int(nrow(simulated.data),500),]
  #main.minusvali<-simulated.data[-as.numeric(row.names(validation)),]
  
  RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)
  simulated.data$treat.estimate<-predict(RC.model2,simulated.data)
  linear.result2[i]<-summary(lm(Y~treat.estimate+cf1+cf2+cf3+cf4+cf5+cf6,data=simulated.data))$coefficients[2]
  linear.noerror[i]<-summary(lm(Y~treat+cf1+cf2+cf3+cf4+cf5+cf6,data=simulated.data))$coefficients[2]
  linear.crude[i]<-summary(lm(Y~treat.w+cf1+cf2+cf3+cf4+cf5+cf6,data=simulated.data))$coefficients[2]
}
linear.beta<-cbind(linear.noerror,linear.crude,linear.result2)
colnames(linear.beta)<-c("error-free","error-prone","X~W+D")
boxplot(linear.beta,main="Linear",cex.axis=2.2,cex.main=2.5)
#abline(h=gold_linear,lty=5,lwd=3,col="red")


mean((linear.result2[1]-1)^2,na.rm = T)
mean((linear.noerror[1]-1)^2,na.rm = T)


