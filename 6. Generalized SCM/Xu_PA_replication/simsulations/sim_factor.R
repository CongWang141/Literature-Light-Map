## Generalized Synthetic Control Method
## Replication Materials
## Correct Choice of the number of factors

## Author: Yiqing Xu

## The codes below replicate Table A5 in the Online Appendix

## Varying T_0, N_co and N_tr  (warning: takes a long time!)


rm(list=ls(all=TRUE))
library(doParallel)
library(parallel)
library(foreach)
library(xtable)
library(abind)
library(gsynth)
source("sim_sampling.R")

sims<-5000
set.seed(123)

# register multiple cores
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)

# for each NNtr, consider 8 cases
TT0 <-rep(c(10,30,50,15,15,15),3)
NNco<-rep(c(40,40,40,40,80,120),3)
NNtr<-rep(c(5,20,40),each=6)
ncases<-length(TT0)

# storage
result<-array(NA,dim=c(sims,8,length(TT0))) # 18 case; sims*8*18

# loop
begin.time<-Sys.time()
cat("Cores: ",cores,"; cases: ",ncases,"; Sims: ",sims,sep="")
for (case in 1:ncases) {
    
    # annouce case
    cat("\nCase ",case,": T0 =",TT0[case]," Nco =",NNco[case]," Ntr =",NNtr[case],"\n",sep="")
    
    onecase<-foreach (j=1:sims,.combine="rbind",.inorder=FALSE,.export=c("gsynth")) %dopar% { # parallel computing
      
      out<-rep(NA,8) # MSPE OLS, MSPE: 0-5, r.cv 
      
      # general random sample: 
      panel<-simulate(Ntr=NNtr[case],Nco=NNco[case],T0=TT0[case],p=2,r=2,m=0,w=0.5,D.sd=1,beta=c(1,3),mu=5,att=c(1:10),fsize=1,FE=1,fixF=FALSE)
      
      # OLS
        out[1]<-gsynth(Y ~ D + X1 + X2, data=panel,index=c("id","time"),
                       se=0,r=c(0,0),CV=1,force="none")$CV.out[1,"MSPE"]
      
      # run the model: the number of factors are known to be 2
        synth.out<-gsynth(Y ~ D + X1 + X2, data=panel,index=c("id","time"),
                          se=0,r=c(0,5),CV=1,force="two-way")
      out[2:7]<-synth.out$CV.out[,"MSPE"]
      out[8]<-synth.out$r.cv
      
      return(out)
    }
        
    # Save after each case
    result[,,case]<-onecase
    save(result,file="sim_factor_result.RData")

} # end of all cases

stopCluster(cl) # stop parallel computing
print(Sys.time()-begin.time)

######################
# Output
######################

load("sim_factor_result.RData")
out1<-apply(result[,8,],2,function(vec){mean(ifelse(vec==2,1,0))})
out2<-apply(result[,8,],2,function(vec){mean(ifelse(vec>2,1,0))})
out3<-apply(result[,8,],2,function(vec){mean(ifelse(vec<2,1,0))})

out1 # correct cases
out2 # case of overfitting
out2 # case of choosing less than 2 factors

output<-as.data.frame(cbind(TT0[1:6],NNco[1:6],rep(NA,length(TT0)/3),out1[1:6],out1[7:12],out1[13:18]))
colnames(output)<-c("T0","Nco","","5","20","40")
round(output,3)

library(xtable)
sink("output_sim_factor.txt")
print(xtable(output,digits=c(0,0,0,rep(3,4))),include.rownames=FALSE)
sink()



pdf("sim_factor_all.pdf",height=6)
MSPE.out<-apply(result[,1:7,],c(2,3),mean)
par(mar=c(2,2,0,0))
plot(1,type="n",xlim=c(1,7),ylim=c(0,15),xlab="",ylab="",axes=F)
box(); axis(2); mtext("MSPE",2,2)
#abline(v=4,lwd=3,col="grey")
for (i in 1:18) lines(MSPE.out[,i])
axis(1,at=c(1:7),c("OLS","Two-way FEs","Plus 1","2","3","4","5"))
title("All 18 Cases",line=-1.5)
graphics.off()

for (case in 1:dim(result)[3]) {
  png(paste("sim_factor_case",case,".png",sep=""),height=300,width=1000)
  par(mar=c(2,3,0,1))
  plot(1,type="n",xlim=c(1,7),ylim=c(0,20),xlab="",ylab="",axes=F)
  box(); axis(2); mtext("MSPE",2,2)
  for (j in 1:dim(result)[1]) lines(1:7, result[j,1:7,case],col="#AAAAAA10")
  #abline(v=4,lty=2,col=2)
  lines(1:7, apply(result[,1:7,case],2,mean),lwd=2)
  points(1:7, apply(result[,1:7,case],2,mean),lwd=2,pch=16)
  axis(1,at=c(1:7),c("OLS","Two-way FEs","Plus 1","2","3","4","5"))
  title(main=list(paste("T0=",TT0[case],", Ntr=",NNtr[case],", Nco=",NNco[case], sep=""),cex=2),line=-1.5)
  graphics.off()
}
