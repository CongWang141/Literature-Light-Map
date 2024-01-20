## Generalized Synthetic Control Method
## Replication Materials
## 95% Coverage

## Author: Yiqing Xu

## Coverage

# test CI coverage for period 0, 5 and 10
# Warning: this gonna take a long time!

rm(list=ls(all=TRUE))
library(doParallel)
library(foreach)
library(xtable)
library(abind)
library(MASS)
library(gsynth)
source("sim_sampling.R")

inference="parametric"
sims<-5000
nboots<-2000
set.seed(123)

# register multiple cores
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)

# set progress bar and combind function
f <- function(){
  function(...) abind(...,along=3)
}

TT0<-c(30)
NNco<-c(120)
Tmax<-max(TT0)+10
ncases<-length(TT0)


Ntr<-40
D.sd<-1   # variance of treatment effect
r<-2
p<-2
w<-1

cat("Cores: ",cores,"; cases: ",ncases,"\n",sep="")
begin.time<-Sys.time()
result<-matrix(NA,Tmax,ncases)
for (case in 1:ncases) {
  
    Nco<-NNco[case];N<-Ntr+Nco;T0<-TT0[case];T<-T0+10;
    cat(inference,": sims = ",sims,", nboots = ",nboots, ", Ntr = ",Ntr,
        ", Nco = ",Nco,", T0 = ",T0,"\n",sep="")
  
    onecase<-foreach (i=1:sims,
                      .combine=f(),
                      .inorder=FALSE,
                      .packages=c("MASS","gsynth")
                      ) %dopar% { # parallel computing
  
                          output<-matrix(NA,T,5)
                          
                                        # generate a random sample: 
                          panel<-simulate(Ntr=Ntr,Nco=Nco,T0=T0,p=p,r=r,
                                          m=0,w=w,D.sd=D.sd,
                                          beta=c(1,3),
                                          mu=5,att=c(1:10),
                                          fsize=1,FE=1,fixF=FALSE)
                          if (inference=="nonparametric") {
                              effect<-c(rep(0,T0),1:10) 
                          } else {
                              effect<-apply(as.matrix(matrix(panel$eff,T,N)[,1:Ntr]),1,sum)/Ntr  
                          }

                          begin.time<-Sys.time()
                          ## run the model: the number of factors are known to be 2
                          out<-gsynth(Y ~ D + X1 + X2, data = panel,
                                      index = c("id","time"),
                                      se = 1, r=r, CV=0, force="two-way",nboots=nboots,
                                      inference=inference,
                                      parallel = FALSE)
                          time<-Sys.time()-begin.time
                          
                                        # storage
                          ATT<-out$est.att[,1]
                          output[,1]<-ATT
                          output[,2]<-ATT-effect #Bias
                          output[,3]<-out$est.att[,3] #CI lower
                          output[,4]<-out$est.att[,4] #CI upper
                          output[,5]<-ifelse(effect>=out$est.att[,3] & effect<=out$est.att[,4],1,0) #good CI
                          return(output)
  }
  # coverage
  result[(Tmax-T+1):Tmax,case]<-apply(onecase[,5,],1,mean)
  
  save(result,file="sim_cover_result.RData")
}
stopCluster(cl) # stop parallel computing
time<-Sys.time()-begin.time
print(time)

save(time,result,file="sim_cover_result.RData")

load("sim_cover_result.RData")
apply(result,2,mean,na.rm=TRUE)

plot(onecase[,1,1],ylim=c(-5,15))
lines(onecase[,3,1])
lines(onecase[,4,1])






