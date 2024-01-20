## Generalized Synthetic Control Method
## Replication Materials
## Comparison with Diff-in-Diffs

## Author: Yiqing Xu

## The codes below replicate Table A2 in the Online Appendix

## Varying w and N_tr (warning: takes a long time!)

# test Bais for SATT 15 (treatment period 5)

rm(list=ls(all=TRUE))
library(doParallel)
library(parallel)
library(foreach)
library(xtable)
library(abind)
library(plm)
library(gsynth)
source("sim_sampling.R")

sims<-5000
set.seed(123)

# register multiple cores
#cores<-detectCores()
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)


# for each model, consider 16 cases
WW<-rep(c(1,0.8,0.6),each=4)
NNtr<-rep(20,12)
NNco<-rep(c(40,80,120,200),3)
T0<-15
T<-T0+10
ncases<-length(WW)

# storage
result<-array(NA,dim=c(length(WW),3,2)) # bias, sd, RMSE; synth and did

# bind funciton
f <- function(...){
  abind(...,along=3)  # bind as an array
}

# loop
begin.time<-Sys.time()
cat("Cores: ",cores,"; cases: ",ncases,"; Sims: ",sims,sep="")
for (case in 1:length(WW)) {
    
    ## annouce case
    w<-WW[case]; Ntr<-NNtr[case]; Nco<-NNco[case]; N<-Ntr+Nco;
    cat("\nCase " ,case,": Ntr = ",Ntr,": Nco = ",Nco,"; w = ",w,"\n",sep="")
 
    panel<-simulate(Ntr=Ntr,Nco=Nco,T0=T0,p=2,r=2,m=0,w=w,D.sd=1,beta=c(1,3),
                    mu=5,att=c(1:10),fsize=1,FE=1,fixF=TRUE, fixL=TRUE)
    
    onecase<-foreach (i=1:sims,.combine="f",.inorder=FALSE,
                      .packages=c("plm","gsynth")) %dopar% { # parallel computing
         

        ## general random sample:
        panel$Y <- panel$Ybar + rnorm(N*T) 
        effect<-apply(as.matrix(matrix(panel$eff,T,N)[,1:Ntr]),1,sum)/Ntr # ATT
        
        ## GSC: the number of factors are known to be 2
        out.syn<-gsynth(Y ~ D + X1 + X2, data=panel,index=c("id","time"),
                        se=FALSE, r = 2, CV = FALSE, force="two-way")
        att.syn<-out.syn$att[(T0+1):T]
        bias.syn<-att.syn-effect[(T0+1):T] # bias
        
        ## DID
        dummies<-model.matrix(~factor(panel$time)-1)[,(1+T0):T]
        dummies[which(panel$D==0),]<-0
        colnames(dummies)<-c("d1","d2","d3","d4","d5","d6","d7","d8","d9","d10")
        panel.ex<-cbind(panel,dummies)
        out.did<-plm(panel$Y~X1+X2+d1+d2+d3+d4+d5+d6+d7+d8+d9+d10,data=panel.ex,index=c("id","time"),effect="twoway",model="within")
        att.did<-out.did$coefficients[3:12]
        bias.did<-att.did-effect[(T0+1):T]
        
        out<-cbind(att.syn,bias.syn,att.did,bias.did)
        return(out)      
  }
  # gsynth
  result[case,1,1]<-mean(onecase[5,2,]) # bias, focus on ATT_15
  result[case,2,1]<-sd(onecase[5,1,])   # sd of ATT_15 
  result[case,3,1]<-sqrt(mean((onecase[5,2,])^2)) # RMSE
  # did
  result[case,1,2]<-mean(onecase[5,4,]) # bias, focus on ATT_15
  result[case,2,2]<-sd(onecase[5,3,])   # sd of ATT_1
  result[case,3,2]<-sqrt(mean((onecase[5,4,])^2)) # RMSE
  
  # Save after each case
  save(result,file="sim_did_result.RData")
  
} # end of all cases
stopCluster(cl) # stop parallel computing
cat("\n");print(Sys.time()-begin.time)

####################
# Output
####################

load("sim_did_result.RData")
output<-as.data.frame(cbind(rep(15,ncases),NNtr,NNco,WW,
                            rep(NA,ncases),result[,,1],
                            rep(NA,ncases),result[,,2]))
colnames(output)<-c("T0","Ntr","Nco","w",rep(c("sp","Bias","SD","RMSE"),2))
round(output,3)


library(xtable)
sink("output_sim_did.txt")
print(xtable(as.matrix(output),digits=c(rep(0,4),2,rep(3,(dim(output)[2]-5+1)))),include.rownames=0)
sink()

