
## Generalized Synthetic Control Method
## Replicatino Material: An Simulated Example

## Author: Yiqing Xu

## The codes below replicate Figure 1 in the paper and Figures A2 and
## A3 in the Online Appendix


rm(list=ls(all=TRUE)) ## clear memeory; be cautious!

## load package and data
library(gsynth)
data(gsynth)
require(parallel)
require(foreach)
require(ggplot2)
require(GGally)

## run the algorithm
system.time(
    out <- gsynth(Y ~ D + X1 + X2, data = simdata,
                  index=c("id","time"), inference="parametric",
                  se = TRUE, nboots = 1000, r = c(0, 5), CV = TRUE,
                  force = "two-way", parallel = TRUE, cores = 4)
)

## save the results
Y<-out$Y.dat
tb<-out$est.att
Yb<-out$Y.bar[,1:2] ## treated average and counterfactual average
tr<-out$tr
pre<-out$pre
T<-out$T
T0<-out$T0
p<-out$p
m<-out$m
Ntr<-out$Ntr
F.hat<-out$factor
L.tr<-out$lambda.tr
L.co<-out$lambda.co
time<-out$time
time.bf<-time[unique(T0)]
show <- 1:30

##############################
## Figure 1: Raw data and ATT
##############################

pdf("sim_att.pdf",width=10,height=7)
true.effect<-matrix(simdata$eff,30,50)[,1:5]
par(mfcol=c(2,1),mar=c(2,3,1,1),lend=1)
# raw plot
plot(time[show],Yb[show,1],type="n",ylim=c(-10, 50),main="",ylab="",xlab="")
for (j in which(tr==0)) lines(time[show],Y[show,j], col="#AAAAAA30")
for (j in which(tr==1)) {lines(time[show],Y[show,j],col="gray60",lwd=0.8)}
lines(time[show],Yb[show,2],col="gray20",lty=5,lwd=2.5)
lines(time[show],Yb[show,1],col=1,lty=1,lwd=2.5)
abline(v=time.bf,col="gray50",lty=1,lwd=1)
legend("topleft",legend=c("Treated Average",
                          "Estimated Y(0) Average for the Treated",
                          "Treated",
                          "Control"),
       cex=1.3, seg.len=2, col=c("1","gray20","gray60","#AAAAAA60"), 
       fill=NA,border=NA, lty=c(1,5,1,1), lwd=c(2.5,2.5,2,2), merge=T,bty="n")
## gap plot
par(lend=1)
plot(time[show],tb[show,1],type="n",ylim=c(-2.2,12),main="",ylab="",xlab="")
polygon(c(rev(time[show]),
          time[show]),
        c(rev(tb[show,4]), tb[show,3]), col = '#80808050', border = NA)
abline(h=0,col="gray50",lty=1)
abline(v=time.bf,col="gray50",lty=1,lwd=1)
lines(time[show],tb[show,1],lwd=2)
lines(1:T,rowMeans(true.effect),col="gray20",lty=5,lwd=2.5)
legend("topleft",
       legend=c("Estimated ATT","True ATT","95% Confidence Intervals"),
       seg.len=2, cex=1.3, col=c(1,"gray20","#80808050"),
       lty=c(1,5, 1),lwd=c(2.5,2.5,15), bty="n",border=NA)
graphics.off()

###################################
## Figure A2: Factors and loadings
###################################

pdf("sim_factor.pdf",width=7,height=7)
par(mar=c(2,2,1,1))
norm<-sqrt(diag(t(L.co)%*%L.co)/(out$N-out$Ntr))
ymax<-max(abs(F.hat))*max(norm)*1.5
ylim<-c(-ymax,ymax)
xlim<-range(time[show])
plot(time[show],Yb[show,2],type="n",ylim=ylim,xlim=xlim, ylab="",xlab="")
abline(h=0,col=1,lty=2)
lines(time[show],F.hat[,1]*norm[1],col=1, lwd = 2)
lines(time[show],F.hat[,2]*norm[2],col="gray50",lty=5, lwd = 2)
legend("topleft",legend=c("Factor 1","Factor 2"), cex = 1.8,
       seg.len=2, col=c(1,"gray60"), fill=NA, border=NA,
       lty=c(1,5),lwd=1.5,merge=T, bty="n")  
graphics.off()

pdf("sim_loading.pdf",width=7,height=7)
par(mfcol=c(2,2),mar=c(0,0,0,0), oma=c(1,2,3,1))
for (i in 1:2) {
    for (j in 1:2) {
        if (i==j) {          
            ylimd<-range(c(density(out$lambda.co[,i])$y))
            plot(density(out$lambda.co[,i]),main="",xlab="",ylab="",axes=FALSE,
                 ,ylim=ylimd)
            if (out$Ntr>=10) lines(density(out$lambda.tr[,i]),
                                   col=1,lty=2,lwd=2)  
            else abline(v=out$lambda.tr[,i],col=1,lty=2)
            box() 
        } else {
            plot(out$lambda.co[,i],out$lambda.co[,j],cex=1,main="",
                 xlab="",ylab="",axes=FALSE,col=1,
                 xlim=range(c(out$lambda.co[,i],out$lambda.tr[,i])),
                 ylim=range(c(out$lambda.co[,j],out$lambda.tr[,j])))
            points(out$lambda.tr[,i],out$lambda.tr[,j],col=1,cex=1.2, pch=16)
            box()
            legend("topleft", c(paste("L",i," vs ","L",j,sep="")),
                       cex=1.8,bty="n")
             } 
    }
}
graphics.off()

###################################
## Figure A3: Individual effects
###################################

pdf("sim_iTE.pdf",width=7,height=10)
Ntr<-out$Ntr
Y.tr<-out$Y.tr
Y.ct<-out$Y.ct
eff<-out$eff
true.effect<-matrix(panel$eff,T,N)[,1:5]
par(mfcol=c(Ntr,2))
par(mar=c(3,2,2,1))
for (i in 1:Ntr) {
  mymain<-paste("id = ",i,sep="")
  plot(Y.tr[,i],type="l",ylim=range(c(Y.tr,Y.ct)),col=1,
       main=mymain,ylab="",xlab="")
  lines(Y.ct[,i],lty=2,col="gray30")
  if(i==1){legend("topleft",legend=c("Treated Outcomes","Counterfactuals"),
                  col=c(1,"gray30"),lty=c(1,2),lwd=c(1,1.2),merge=T,border=NA,bty="n")}
}
for (i in 1:Ntr) {
  plot(eff[,i],type="l",ylim=c(-5,15),main="",ylab="",xlab="")
  lines(true.effect[,i],col="gray50",lty=2,lwd=1.2)
  if(i==1){legend("topleft",legend=c("Estimated Treatment Effect",
                                     "Actual Treatment Effect"),
                  col=c(1,"gray50"),lty=c(1,1.2),
                  lwd=c(1,1.5),merge=T,border=NA,bty="n")}
}
graphics.off()







