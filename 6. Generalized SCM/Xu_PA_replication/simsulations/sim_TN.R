## Generalized Synthetic Control Method
## Replication Materials
## Finite Sample Properties

## Author: Yiqing Xu

## The codes below replicate Table A1 in the Online Appendix

## Varying T_0, N_co and N_tr (warning: takes a long time!)


## test Bais for ATT 15 (treatment period 5)

## Outputs:
## LaTeX: output_sim_TN.txt
## R file: sim_TN_result.RData

#install.packages("foreach")
#install.packages("xtable")
#install.packages("doParallel")
#install.packages("abind")
#install.packages("doRedis")

rm(list=ls(all=TRUE))
Sys.setenv(LANGUAGE='en')
library(doParallel)
library(parallel)
library(foreach)
library(xtable)
library(abind)
library(gsynth)
source("sim_sampling.R")

# register multiple cores
#cores<-detectCores()
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)
cat("Nodes:",cores)

sims<-5000
set.seed(123)

# for each NNtr, consider 36 cases in total
r=2
w=0.8 # 80% overlap
TT0 <-rep(c(15,30,50),each=4)
NNco<-rep(c(40,80,120,200),3)
NNtr<-c(1,5,20)

# set progress bar and combind function
f <- function(){function(...) {abind(...,along=3)}}

# storage
result<-array(NA,dim=c(length(TT0),3,length(NNtr))) # bias and sd and RMSE
betas<-matrix(NA,sims,length(TT0)*length(NNtr))

# loop
k=1
begin.time<-Sys.time()
for (i.ntr in 1:length(NNtr)) {
  Ntr<-NNtr[i.ntr]
  for (case in 1:length(NNco)) {
    T0 <-TT0[case]
    Nco<-NNco[case]
    T<-T0+10
    N<-Ntr+Nco
    
    ## annouce case
    cat("\nCase: ",k," T0 =",T0," Nco =",Nco," Ntr =",Ntr," w=",w,"\n",sep="")

    panel<-simulate(Ntr=Ntr,Nco=Nco,T0=T0,p=2,r=2,m=0,w=w,D.sd=1,beta=c(1,3),
                     mu=5,att=c(1:10),fsize=1,FE=1,fixF=TRUE, fixL=TRUE)
    
    onecase<-foreach(i=1:sims,
                     .combine=f(),
                     .packages=c("gsynth"),
                     .inorder=FALSE) %dopar% { # parallel computing
                         
                         ## general random sample: (5tr+45co)*(10pre+10post)
                         panel$Y <- panel$Ybar + rnorm(N*T) 
                         effect<-apply(as.matrix(matrix(panel$eff,T,N)[,1:Ntr]),1,sum)/Ntr # SATT
                         
                         ## run the model: the number of factors are known to be 2
                         out<-gsynth(Y ~ D + X1 + X2, data=panel,
                                     index = c("id","time"),force = "two-way",
                                     se=FALSE, r = 2, CV = FALSE)
                         att<-out$att
                         bias<-out$att-effect
                         result<-cbind(att,bias,out$beta[2])
                         return(result)         
    }
    
    # Save after each case    
    result[case,1,i.ntr]<-mean(onecase[(T0+5),2,]) # bias, focus on ATT_15
    result[case,2,i.ntr]<-sd(onecase[(T0+5),1,])   # sd of ATT_15 
    result[case,3,i.ntr]<-sqrt(mean((onecase[(T0+5),2,])^2)) # RMSE
    betas[,k]<-onecase[1,3,]
    
    save(result,betas=betas,file="sim_TN_result.RData")
    k=k+1
    
  }  # end of T0, N.co combination
} # end of N.tr case
stopCluster(cl) # stop parallel computing
print(Sys.time()-begin.time)

####################
# Output
####################

load("sim_TN_result.RData")
output<-cbind(TT0,NNco)
for (i in 1:length(NNtr)) {
  output<-  cbind(output,rep(NA,length(TT0)),result[,,i])
}
colnames(output)<-c("T0","Nco",rep(c("sp","Bias","SD","RMSE"),length(NNtr)))
round(output,3)
sink("output_sim_TN.txt")
print(xtable(output,digits=c(0,0,0,rep(3,4*length(NNtr)))),include.rownames=FALSE)
sink()
