
## Generalized Synthetic Control Method
## Replicatino Material: EDR on Voter Turnout

## Author: Yiqing Xu

## The codes below replicate Tables 2 and 3, Figures 2 and 3 in the
## paper, as well as Figures A5 and A6 in the Online Appendix

rm(list=ls(all=TRUE)) ## eliminating everything in the memory; be cautious

library(gsynth)
data(gsynth)


##########################
## Table 2
##########################

## interFE is a function in the gsynth package that implements
## interactive fixed effects models, of which additive fixed effects
## models are special case (when r = 0).

## Column 1 (DID)
out.did1<-interFE(turnout ~ policy_edr,
                  data = turnout, index = c("abb", "year"),
                  r = 0, force="two-way", SE = TRUE, nboots = 2000,
                  seed = 2139)

## Column 2 (DID)
out.did2<-interFE(turnout ~ policy_edr + policy_mail_in + policy_motor,
                  data = turnout, index = c("abb", "year"),
                  r = 0, force="two-way", SE = TRUE, nboots = 2000,
                  seed = 2139)

## Column 3 (GSC)
out.syn1<-gsynth(turnout ~ policy_edr,
                 data = turnout, index =c("abb", "year"),
                 force = "two-way", CV = TRUE, r=c(0,5), se = TRUE,
                 parallel = FALSE, nboots=2000, seed = 2139)

## Column 4 (GSC)
out.syn2<-gsynth(turnout ~ policy_edr + policy_mail_in + policy_motor,
                 data = turnout, index =c("abb", "year"),
                 force = "two-way", CV = TRUE, r=c(0,5), se = TRUE,
                 parallel = FALSE, nboots=2000, seed = 2139)


## first look (not in the paper)
plot(out.syn2, type = "gap", xlim = c(-13,4))
plot(out.syn2, type = "raw")
plot(out.syn2, type = "counterfactual")
plot(out.syn2, type = "factors")
plot(out.syn2, type = "loadings")


## DID (dynamic effect)
out.did3<-gsynth(turnout ~ policy_edr + policy_mail_in + policy_motor,
                 data = turnout, index =c("abb", "year"),
                 force = "two-way", CV = FALSE, r = 0, se = TRUE,
                 parallel = FALSE, nboots=2000,
                 inference = "nonparametric", seed = 2139)
plot(out.did3, type = "gap", xlim = c(-13,4))
plot(out.did3, type = "counterfactual")

##########################
## Table 3
##########################

## ME, MN, WI
sub1 <- gsynth(turnout ~ policy_edr + policy_mail_in + policy_motor,
               data = turnout[which(!turnout$abb%in%c("ID","NH","WY","MT","IA","CT")),],
               index =c("abb", "year"), force = "two-way", CV = FALSE,
               r = 2, se = TRUE, parallel = FALSE, nboots=2000, seed = 2139)

## ID, NH, WY
sub2 <- gsynth(turnout ~ policy_edr + policy_mail_in + policy_motor,
               data = turnout[which(!turnout$abb%in%c("ME","MN","WI","MT","IA","CT")),],
               index =c("abb", "year"), force = "two-way", CV = FALSE,
               r = 2, se = TRUE, parallel = FALSE, nboots=2000, seed = 2139)

## MT, IA, CT
sub3 <- gsynth(turnout ~ policy_edr + policy_mail_in + policy_motor,
               data = turnout[which(!turnout$abb%in%c("ME","MN","WI","ID","NH","WY")),],
               index =c("abb", "year"), force = "two-way", CV = FALSE,
               r=2, se = TRUE, parallel = FALSE, nboots=2000, seed = 2139)

##############################
## Figiure 2
##############################

## DID ##
pdf("fg_edr_main_did.pdf",width=14,height=5)
## counterfactual
par(mfcol=c(1,2),mar=c(4,4,1,1),lend=1)
out<-out.did3
time<-c(-13:4)
plot(1,type="n",xlab="",ylab='',axes=F,xlim=range(time),ylim=c(55,75))
box()
axis(1,at=seq(-12,4,2));mtext("Term Relative to Reform",1,2.5,cex=1.5)
axis(2);mtext("Turnout %",2,2.5,cex=1.5)
abline(v=0,col="gray",lwd=2,lty=2)
lines(time,out$Y.tr.cnt[1:18],lwd=2)
lines(time,out$Y.ct.cnt[1:18],col=1,lty=5,lwd=2)
legend("topleft",legend=c("Treated Average","Estimated Y(0) Average for the Treated"),cex=1.5,
       seg.len=2, col=c(1,1),lty=c(1,5),lwd=2,bty="n")
## gap
newx<-c(-13:4)
plot(1,type="n",xlab="",ylab='',axes=F,xlim=c(-13,4),ylim=c(-8,8))
box()
axis(1,at=seq(-12,4,2));mtext("Term relative to reform",1,2.5,cex=1.5)
axis(2,at=seq(-8,8,4));mtext("Turnout %",2,2.5,cex=1.5)
abline(v=0,col="gray",lty=2,lwd=2)
abline(h=0,col="gray20",lty=2,lwd=1)
polygon(c(rev(newx), newx), c(rev(out$est.att[1:18,3]), out$est.att[1:18, 4]),
        col = "#55555530", border = NA)
lines(newx,out$est.att[1:18,1],col=1,lty=1,lwd=2)
legend("topleft",legend=c("Estimated ATT","95% Confidence Intervals"), cex=1.5, seg.len=2,
       col=c(1,"#55555530"),lty=c(1,5),lwd=c(2,20),bty="n")
graphics.off()


## GSC ##
pdf("fg_edr_main_syn.pdf",width=14,height=5)
## counterfactual
par(mfcol=c(1,2),mar=c(4,4,1,1),lend=1)
out<-out.syn1
time<-c(-13:4)
plot(1,type="n",xlab="",ylab='',axes=F,xlim=range(time),ylim=c(55,75))
box()
axis(1,at=seq(-12,4,2));mtext("Term Relative to Reform",1,2.5,cex=1.5)
axis(2);mtext("Turnout %",2,2.5,cex=1.5)
abline(v=0,col="gray",lwd=2,lty=2)
lines(time,out$Y.tr.cnt[1:18],lwd=2)
lines(time,out$Y.ct.cnt[1:18],col=1,lty=5,lwd=2)
legend("topleft",legend=c("Treated Average","Estimated Y(0) Average for the Treated"),cex=1.5,
       seg.len=2, col=c(1,1),lty=c(1,5),lwd=2,bty="n") 
## gap
newx<-c(-13:4)
plot(1,type="n",xlab="",ylab='',axes=F,xlim=c(-13,4),ylim=c(-8,8))
box()
axis(1,at=seq(-12,4,2));mtext("Term relative to reform",1,2.5,cex=1.5)
axis(2,at=seq(-8,8,4));mtext("Turnout %",2,2.5,cex=1.5)
abline(v=0,col="gray",lty=2,lwd=2)
abline(h=0,col="gray20",lty=2,lwd=1)
polygon(c(rev(newx), newx), c(rev(out$est.att[1:18,3]), out$est.att[1:18, 4]),
        col = "#55555530", border = NA)
lines(newx,out$est.att[1:18,1],col=1,lty=1,lwd=2)
legend("topleft",legend=c("Estimated ATT","95% Confidence Intervals"), cex=1.5, seg.len=2,
       col=c(1,"#55555530"),lty=c(1,5),lwd=c(2,20),bty="n")
graphics.off()

##################################
## Figiure 3: Factors and Loadings
##################################

## factors

pdf("fg_edr_factors.pdf")
out<-out.syn2
L.co<-out$lambda.co
norm<-sqrt(diag(t(L.co)%*%L.co)/(out$N-out$Ntr))
time<-out$time
par(mar=c(3,3,1,1))
plot(time,out$Yb[,2],type="n",ylim=c(-15,15),xlim=c(1920,2016),axes=0,ylab="",xlab="")
box()
axis(1,at=seq(1920,2016,8));mtext("Year",1,2,cex=1.2)
axis(2);mtext("Turnout %",2,2,cex=1.2)
abline(h=0,col=1,lty=2)
lines(time,out$factor[,1]*norm[1],col=1,lwd=2)
lines(time,out$factor[,2]*norm[2],col=1,lty=5,lwd=2)
legend("topleft",legend=c("Factor 1","Factor 2"), cex=1.5, seg.len=2,
       col=c(1,1),lty=c(1,5),lwd=2,bty="n")
graphics.off()

## loadings

pdf("fg_edr_loadings2.pdf")
out<-out.syn2
par(mar=c(3,3,1,1),lend=1)
plot(1,main="",type="n",xlab="",ylab="",axes=FALSE,xlim=c(-20,20),ylim=c(-5,6))
text(out$lambda.co[,1],out$lambda.co[,2], labels=out$id.co, cex= 1.3,col="gray50")
text(out$lambda.tr[,1],out$lambda.tr[,2], labels=out$id.tr, cex= 1.3,col=1, font=2)
box()
axis(1);mtext("Loadings for factor 1",1,2,cex=1.2)
axis(2);mtext("Loadings for factor 2",2,2,cex=1.2)
# mark south
south<-order(out$lambda.co[,1],decreasing=1)[1:11]
x0<-out$lambda.co[south,1]
y0<-out$lambda.co[south,2]-0.2
segments(x0-0.95,y0,x0+0.95,y0,col="gray50",lwd=1.5)
text(-17.3,5.9, "Treated", cex= 1.5,col=1, font=2)
text(-12.5,5.3, "Control (non-South)", cex= 1.5,col="gray50")
text(-14.5,4.7, "Control (South)", cex= 1.5,col="gray50")
segments(-20.5,4.4,-8.5,4.4,col="gray50",lwd=1.5)
graphics.off()

##############################
## Figiure A5: Raw Data
##############################

pdf("fg_edr_raw.pdf",width=10)
out<-out.syn2
time<-out$time
Y<-out$Y.dat
Y.tr<-out$Y.dat[,which(out$tr==1)]
tr<-out$tr
par(mar=c(3,3,1,1))
plot(1,type="n",ylim=c(0,100),xlim=c(1920,2016),axes=F,main="",ylab="",xlab="")
box()
axis(1,at=seq(1920,2016,8));mtext("Year",1,2,cex=1.2)
axis(2);mtext("Turnout %",2,2,cex=1.2)
for (j in which(tr==0)) lines(time,Y[,j],col="#99999980",lwd=1.3)
for (j in 1:out$Ntr) lines(time[1:(max(which(out$pre[,j]==1))+1)],
                           Y.tr[1:(max(which(out$pre[,j]==1))+1),j],col= 1,lty=1,lwd=1.3)
for (j in 1:out$Ntr) lines(time[which(out$pre[,j]==0)],
                           Y.tr[which(out$pre[,j]==0),j],col=1,lty=5,lwd=1.5)
legend("bottomright",legend=c("Controls","Treated (pre)","Treated (post)"),
       cex=1.5, seg.len=2, col=c("#99999980",1,1), 
       fill=NA,border=NA, lty=c(1,1,5),lwd=c(1.5,1.5,1.5), merge=T,bty="n")
graphics.off()

##############################
## Figiure A6: Each State
##############################

out<-out.syn2
names(out)
time<-out$time
T0<-out$T0
T<-out$T
states<-out$id.tr

for (abb in states) { # treated states
  
    id<-which(out$id.tr==abb)
    t0<-T0[id]
    
    pdf(paste("fg_edr_case_",abb,".pdf",sep=""),width=14,height=5)
    par(mfcol=c(1,2),mar=c(4,4,1,1),lend=1)
    
    ## counterfactual
    yrg<-range(c(out$Y.tr[,id],out$Y.ct[,id]))
    ylim=c(yrg[1]-10,yrg[2]+5)
    plot(1,type="n",xlab="",ylab='',axes=F,xlim=c(1920,2016),ylim=ylim)
    box()
    axis(1,at=seq(1920,2016,8));mtext("Year",1,2.5,cex=1.5)
    axis(2);mtext("Turnout %",2,2,cex=1.5)
    abline(v=time[t0],col="gray",lwd=2)
    lines(time,out$xi+mean(out$Y.tr[,id]),col="gray50",lwd=2,lty=5)
    lines(time,out$Y.tr[,id],lwd=3)  
    lines(time,out$Y.ct[,id],col="gray50",lty=1,lwd=2)
    text(1925,yrg[2]+2,abb,cex=2)
    legend("bottomright",
           legend=c("Actual Outcome","Estimated Y(0) (DID)","Estimated Y(0) (GSC)"),
           cex=1.3, seg.len=2, col=c(1,"gray50","gray50"),lty=c(1,5,1),
           lwd=c(3,2,2),bty="n")  
    
    yrg<-range(c(out$est.ind[,,id],out$Y.tr[,id]-out$xi-mean(out$Y.tr[,id])))
    ylim<-c(yrg[1]-15,yrg[2]+5)
    newx<-seq(1920,2012,4)
    plot(1,type="n",xlab="",ylab='',axes=F,xlim=c(1920,2016),ylim=ylim)
    box()
    axis(1,at=seq(1920,2016,8));mtext("Year",1,2.5,cex=1.5)
    axis(2);mtext("Turnout %",2,2,cex=1.5)
    abline(v=time[t0],col="gray",lwd=2)
    abline(h=0,col="gray20",lty=2,lwd=1)
    lines(time,out$Y.tr[,id]-out$xi-mean(out$Y.tr[,id]),col=1,lwd=2,lty=5)
    polygon(c(rev(newx), newx), c(rev(out$est.ind[ ,3,id]), out$est.ind[ ,4,id]),
            col = "#55555530", border = NA)
    lines(time,out$est.ind[,1,id],col=1,lty=1,lwd=2)  
    legend("bottomright",
           legend=c("Treatment Effect (DID)","Treatment Effect (GSC)",
                    "95% Confidence Interval (GSC)"),
           cex=1.3, seg.len=2, col=c(1,1,"#55555550"),lty=c(5,1,1),lwd=c(2,2,15),
           bty="n")  

  graphics.off()   
}

