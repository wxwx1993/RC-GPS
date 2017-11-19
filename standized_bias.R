library("designmatch")
library("matrixStats")
library("SDMTools")

#########balance 
#### generate data
simulated.data<-data.generate(sample_size=2000,sd=20,phi=0.8,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=1)
validation<-simulated.data[sample.int(nrow(simulated.data),500),]
RC.model2<-lm(treat~ treat.w  +cova1+cova2+cova3 , validation)

imputation.data<-simulated.data
imputation.data$treat.estimate<-predict(RC.model2,simulated.data)
treat.estimate.cat<-as.numeric()
for (i in 1:nrow(imputation.data)){
  if (imputation.data$treat.estimate[i]>15){treat.estimate.cat[i]<-2}
  else if (imputation.data$treat.estimate[i]<=15 & imputation.data$treat.estimate[i]>(-5)){treat.estimate.cat[i]<-1}
  else if (imputation.data$treat.estimate[i]<=(-5)){treat.estimate.cat[i]<-0}
}
imputation.data$treat.estimate.cat<-treat.estimate.cat

data<-imputation.data
GPS_mod<-multinom(treat.estimate.cat ~ cf1+cf2+cf3+cf4+cf5+cf6 , data)
#GPS_mod
data$GPS0<-GPS_mod$fitted.values[,1]
data$GPS1<-GPS_mod$fitted.values[,2]
data$GPS2<-GPS_mod$fitted.values[,3]

########subclass balance
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

#subclass
subclass.asd<-function(GPSdataset=Y.GPS0,cf=1,cate=0){
  sum(c((mean(subset(GPSdataset[[1]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[1]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[1]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[1]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[1]])/nrow(data),
        (mean(subset(GPSdataset[[2]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[2]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[2]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[2]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[2]])/nrow(data),
        (mean(subset(GPSdataset[[3]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[3]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[3]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[3]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[3]])/nrow(data),
        (mean(subset(GPSdataset[[4]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[4]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[4]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[4]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[4]])/nrow(data),
        (mean(subset(GPSdataset[[5]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[5]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[5]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[5]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[5]])/nrow(data),
        (mean(subset(GPSdataset[[6]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[6]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[6]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[6]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[6]])/nrow(data),
        (mean(subset(GPSdataset[[7]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[7]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[7]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[7]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[7]])/nrow(data),
        (mean(subset(GPSdataset[[8]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[8]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[8]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[8]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[8]])/nrow(data),
        (mean(subset(GPSdataset[[9]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[9]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[9]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[9]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[9]])/nrow(data),
        (mean(subset(GPSdataset[[10]],treat.estimate.cat==cate)[,cf+3])-mean(subset(GPSdataset[[10]],treat.estimate.cat!=cate)[,cf+3]))/((sd(subset(GPSdataset[[10]],treat.estimate.cat==cate)[,cf+3])+sd(subset(GPSdataset[[10]],treat.estimate.cat!=cate)[,cf+3]))/2)*nrow(GPSdataset[[10]])/nrow(data)
  ),na.rm = T)
}


subclass.balance<-function(GPSdataset=Y.GPS0,cf=1,cate=0){
  subclass.asd(GPSdataset,cf,cate)
}


############IPTW balance
data.X2.trun<-data[which(data$treat.estimate.cat==2),]
data.X2.trun$weight<-1/data.X2.trun$GPS2
data.X1.trun<-data[which(data$treat.estimate.cat==1),]
data.X1.trun$weight<-1/data.X1.trun$GPS1
data.X0.trun<-data[which(data$treat.estimate.cat==0),]
data.X0.trun$weight<-1/data.X0.trun$GPS0

data_weight<-rbind(data.X0.trun,data.X1.trun,data.X2.trun)
data_weight[which(data_weight$weight>10),]$weight<-10

IPTW.balance<-function(cf=1,cate=0){
  (wt.mean(subset(data_weight,treat.estimate.cat==cate)[,cf+3],subset(data_weight,treat.estimate.cat==cate)$weight)
  -wt.mean(subset(data_weight,treat.estimate.cat!=cate)[,cf+3],subset(data_weight,treat.estimate.cat!=cate)$weight))/
    ((wt.sd(subset(data_weight,treat.estimate.cat==cate)[,cf+3],subset(data_weight,treat.estimate.cat==cate)$weight)
      +wt.sd(subset(data_weight,treat.estimate.cat!=cate)[,cf+3],subset(data_weight,treat.estimate.cat!=cate)$weight))/2)
}


#######matching bbalance
data.X2<-data[which(data$treat.estimate.cat==2),]
data.X1<-data[which(data$treat.estimate.cat==1),]
data.X0<-data[which(data$treat.estimate.cat==0),]

rank.Y0<-NULL;rank.Y1<-NULL;rank.Y2<-NULL
for (i in 1:nrow(data)){
  rank.Y0[i]<-which.min(abs(data.X0$GPS0-data$GPS0[i]))
  rank.Y1[i]<-which.min(abs(data.X1$GPS1-data$GPS1[i]))
  rank.Y2[i]<-which.min(abs(data.X2$GPS2-data$GPS2[i]))
  
}
E.Y0<-data.X0[rank.Y0,]
E.Y1<-data.X1[rank.Y1,]
E.Y2<-data.X2[rank.Y2,]
matching.data<-rbind(E.Y0,E.Y1,E.Y2)

matching.balance<-function(cf=1,cate=0){
  (mean(subset(matching.data,treat.estimate.cat==cate)[,cf+3])-mean(subset(matching.data,treat.estimate.cat!=cate)[,cf+3]))/
   ((sd(subset(matching.data,treat.estimate.cat==cate)[,cf+3])+sd(subset(matching.data,treat.estimate.cat!=cate)[,cf+3]))/2)
}

##########original data balance
init.balance<-function(cf=1,cate=0){
  (mean(subset(QD_data_complete,treat.estimate.cat==cate)[,cf+3])-mean(subset(QD_data_complete,treat.estimate.cat!=cate)[,cf+3]))/
    ((sd(subset(QD_data_complete,treat.estimate.cat==cate)[,cf+3])+sd(subset(QD_data_complete,treat.estimate.cat==cate)[,cf+3]))/2)
}

############# calculation
subclass.balance.0<-NULL
subclass.balance.1<-NULL
subclass.balance.2<-NULL
for (i in c(1:6)){
  subclass.balance.0<-c( subclass.balance.0,subclass.balance(Y.GPS0,cf=i,cate=0))
  subclass.balance.1<-c( subclass.balance.1,subclass.balance(Y.GPS1,cf=i,cate=1))
  subclass.balance.2<-c( subclass.balance.2,subclass.balance(Y.GPS2,cf=i,cate=2))
}

IPTW.balance.0<-NULL
IPTW.balance.1<-NULL
IPTW.balance.2<-NULL
for (i in c(1:6)){
  IPTW.balance.0<-c( IPTW.balance.0,IPTW.balance(cf=i,cate=0))
  IPTW.balance.1<-c( IPTW.balance.1,IPTW.balance(cf=i,cate=1))
  IPTW.balance.2<-c( IPTW.balance.2,IPTW.balance(cf=i,cate=2))
}


matching.balance.0<-NULL
matching.balance.1<-NULL
matching.balance.2<-NULL
for (i in c(1:6)){
  matching.balance.0<-c( matching.balance.0,matching.balance(cf=i,cate=0))
  matching.balance.1<-c( matching.balance.1,matching.balance(cf=i,cate=1))
  matching.balance.2<-c( matching.balance.2,matching.balance(cf=i,cate=2))
}

init.balance.0<-NULL
init.balance.1<-NULL
init.balance.2<-NULL
for (i in c(1:6)){
  init.balance.0<-c( init.balance.0,init.balance(cf=i,cate=0))
  init.balance.1<-c( init.balance.1,init.balance(cf=i,cate=1))
  init.balance.2<-c( init.balance.2,init.balance(cf=i,cate=2))
}

order.0<-sort(abs(init.balance.0), index.return =T)$ix
order.1<-sort(abs(init.balance.1), index.return =T)$ix
order.2<-sort(abs(init.balance.2), index.return =T)$ix

#############plot balance plot
par(oma = c(4, 1, 1, 1))
par(mfrow=c(1,3),las=1, mgp=c(3, 2, 0))
histcolors = c(rgb(0,0,1,.8), rgb(0,1,0,.8),rgb(1,0,0,.8),rgb(0,0,0,.8))
plot(abs(init.balance.0)[order.0],c(1,2,3,4,5,6),yaxt="n",xlab="",ylab="",type = "l",xlim=c(0,0.4), col = histcolors[4],main=expression(paste(X["c"]==1," vs ",X["c"]!=1)),lwd=3,cex.axis=3,cex.main=4.5)
lines(abs(subclass.balance.0)[order.0],c(1,2,3,4,5,6), col = histcolors[1],lwd=3)
lines(abs(IPTW.balance.0)[order.0],c(1,2,3,4,5,6), col = histcolors[2],lwd=3)
lines(abs(matching.balance.0)[order.0],c(1,2,3,4,5,6), col = histcolors[3],lwd=3)
axis(2, at=c(1,2,3,4,5,6), labels=c("cf1","cf2","cf3","cf4","cf5","cf6")[order.1],cex.axis=2.9)
#legend("topright",c("Subclass","IPTW","Matching","Original"),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=0.8)

plot(abs(init.balance.1)[order.1],c(1,2,3,4,5,6),yaxt="n",xlab="",ylab="",type = "l",xlim=c(0,0.4), col = histcolors[4],main=expression(paste(X["c"]==2," vs ",X["c"]!=2)),lwd=3,cex.axis=3,cex.main=4.5)
lines(abs(subclass.balance.1)[order.1],c(1,2,3,4,5,6), col = histcolors[1],lwd=3)
lines(abs(IPTW.balance.1)[order.1],c(1,2,3,4,5,6), col = histcolors[2],lwd=3)
lines(abs(matching.balance.1)[order.1],c(1,2,3,4,5,6), col = histcolors[3],lwd=3)
axis(2, at=c(1,2,3,4,5,6), labels=c("cf1","cf2","cf3","cf4","cf5","cf6")[order.1],cex.axis=2.9)
#legend("topright",c("Subclass","IPTW","Matching","Original"),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=0.8)


plot(abs(init.balance.2)[order.2],c(1,2,3,4,5,6),yaxt="n",xlab="",ylab="",type = "l",xlim=c(0,0.4), col = histcolors[4],main=expression(paste(X["c"]==3," vs ",X["c"]!=3)),lwd=3,cex.axis=3,cex.main=4.5)
lines(abs(subclass.balance.0)[order.2],c(1,2,3,4,5,6), col = histcolors[1],lwd=3)
lines(abs(IPTW.balance.2)[order.2],c(1,2,3,4,5,6), col = histcolors[2],lwd=3)
lines(abs(matching.balance.2)[order.2],c(1,2,3,4,5,6), col = histcolors[3],lwd=3)
axis(2, at=c(1,2,3,4,5,6), labels=c("cf1","cf2","cf3","cf4","cf5","cf6")[order.2],cex.axis=2.9)
#legend("topright",c("Subclass","IPTW","Matching","Original"),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=0.8)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0.1, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottomleft", c("Subclass","IPTW","Matching","Original"), xpd = TRUE, horiz = TRUE, inset = c(0, 
        0.0), bty = "n", col=histcolors,lwd=c(10,10,10,10),lty=c(1,1,1,1), cex = 2)


