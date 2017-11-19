
pdf("overlap_cor.pdf",width=14,height=8.5,paper='special')   
par(mfrow=c(1,3))
histcolors = c(rgb(0,0,1,.2), rgb(1,0,0,.2),rgb(0,1,0,.2))
with(subset(data,treat.cat==0), hist(GPS0, col = histcolors[1], breaks=seq(0,1,0.03), xlim = c(0,1), axes=T, ylab = "", xlab = "", main = "GPS1",cex.main=2.5,cex.axis=2.5))
par(new=TRUE)
with(subset(data,treat.cat==1), hist(GPS0, col = histcolors[2], breaks=seq(0,1,0.03), xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = ""))
par(new=TRUE)
with(subset(data,treat.cat==2), hist(GPS0, col = histcolors[3], breaks=seq(0,1,0.03), xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = ""))
legend("topright",c(expression(X["c"]==1),expression(X["c"]==2),expression(X["c"]==3)),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=2)

histcolors = c(rgb(0,0,1,.2), rgb(1,0,0,.2),rgb(0,1,0,.2))
with(subset(data,treat.cat==0), hist(GPS1, col = histcolors[1], breaks=seq(0,1,0.03), xlim = c(0,1), axes=T, ylab = "", xlab = "", main = "GPS2",cex.main=2.5,cex.axis=2.5))
par(new=TRUE)
with(subset(data,treat.cat==1), hist(GPS1, col = histcolors[2], breaks=seq(0,1,0.03), xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = ""))
par(new=TRUE)
with(subset(data,treat.cat==2), hist(GPS1, col = histcolors[3], breaks=seq(0,1,0.03), xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = ""))
legend("topright",c(expression(X["c"]==1),expression(X["c"]==2),expression(X["c"]==3)),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=2)

histcolors = c(rgb(0,0,1,.2), rgb(1,0,0,.2),rgb(0,1,0,.2))
with(subset(data,treat.cat==0), hist(GPS2, col = histcolors[1], breaks=seq(0,1,0.03), xlim = c(0,1), axes=T, ylab = "", xlab = "", main = "GPS3",cex.main=2.5,cex.axis=2.5))
par(new=TRUE)
with(subset(data,treat.cat==1), hist(GPS2, col = histcolors[2], breaks=seq(0,1,0.03), xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = ""))
par(new=TRUE)
with(subset(data,treat.cat==2), hist(GPS2, col = histcolors[3], breaks=seq(0,1,0.03), xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = ""))
legend("topright",c(expression(X["c"]==1),expression(X["c"]==2),expression(X["c"]==3)),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=2)

dev.off()



pdf("balance_cor.pdf",width=14,height=8.5,paper='special')   
par(mfrow=c(1,3))
histcolors = c(rgb(0,0,1,.2), rgb(1,1,0,.2))
with(subset(data,treat.cat==0), hist(GPS0, col = histcolors[1], breaks=25, xlim = c(0,1), axes=T, ylab = "", xlab = "", main = "GPS0",cex.main=2.5,cex.axis=2.5,freq=F))
par(new=TRUE)
with(subset(data,treat.cat!=0), hist(GPS0, col = histcolors[2], breaks=25, xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = "",freq=F))
#par(new=TRUE)
#with(subset(data,treat.cat==2), hist(GPS0, col = histcolors[3], breaks=25, xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = "",freq=F))
legend("topright",c("T=0","T=1 or 2"),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=2)

histcolors = c(rgb(1,0,0,.2),rgb(0,1,1,.2))
with(subset(data,treat.cat==1), hist(GPS1, col = histcolors[1], breaks=15, xlim = c(0,1), axes=T, ylab = "", xlab = "", main = "GPS1",cex.main=2.5,cex.axis=2.5,freq=F))
par(new=TRUE)
with(subset(data,treat.cat!=1), hist(GPS1, col = histcolors[2], breaks=15, xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = "",freq=F))
#par(new=TRUE)
#with(subset(data,treat.cat==2), hist(GPS1, col = histcolors[3], breaks=25, xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = "",freq=F))
legend("topright",c("T=1","T=0 or 2"),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=2)


histcolors = c(rgb(0,1,0,.2),rgb(1,0,1,.2))
with(subset(data,treat.cat==2), hist(GPS2, col = histcolors[1], breaks=25, xlim = c(0,1), axes=T, ylab = "", xlab = "", main = "GPS2",cex.main=2.5,cex.axis=2.5,freq=F))
par(new=TRUE)
with(subset(data,treat.cat!=2), hist(GPS2, col = histcolors[2], breaks=25, xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = "",freq=F))
#par(new=TRUE)
#with(subset(data,treat.cat==2), hist(GPS2, col = histcolors[3], breaks=25, xlim = c(0,1), axes=FALSE, ylab = "", xlab = "", main = "",freq=F))
legend("topright",c("T=2","T=1 or 2"),col=histcolors,lwd=c(10,10,10),lty=c(1,1,1),cex=2)

dev.off()