## Generalized Synthetic Control Method
## Replication Materials
## Comparison with Interative Fixed Effects

## Author: Yiqing Xu

## The codes below replicate Table A3 in the Online Appendix

## Varying var(delta), and N_tr (warning: takes a long time!)

# test Bais for ATT 15 (treatment period 5)

rm(list=ls(all=TRUE))
require(doParallel)
require(parallel)
require(foreach)
require(xtable)
require(abind)
require(gsynth)
source("sim_sampling.R")

sims<-5000
set.seed(123)

# register multiple cores
cores<-8
cl<-makeCluster(cores)
registerDoParallel(cl)
cat("Cores:",cores)

# set progress bar and combind function
f <- function(){
  function(...) {
    abind(...,along=3)  # bind as an array
  }
}

# for each model, consider 8 cases
DDsd<-rep(c(0,5),each=4)
NNtr<-rep(20,8)
NNco<-rep(c(40,80,120,200),2)
T0<-15
w<-1

# storage
result<-array(NA,dim=c(length(DDsd),3,2)) # bias, sd, RMSE; synth and inter

# loop
begin.time<-Sys.time()
for (case in 1:length(DDsd)) {
    
  ## annouce case
    cat("\nCase " ,case,": D.sd = ",DDsd[case],"; Nco = ",NNco[case],"\n",sep="")
    Ntr=NNtr[case]; Nco=NNco[case]; N=Ntr+Nco; T=T0+10; Dsd<-DDsd[case]

    panel<-simulate(Ntr=Ntr,Nco=Nco,T0=T0,p=2,r=2,m=0,w=w,D.sd=Dsd,beta=c(1,3),
                    mu=5,att=c(1:10),fsize=1,FE=1,fixF=TRUE, fixL=TRUE)
    
    onecase<-foreach (i=1:sims,
                      .combine=f(),
                      .inorder=FALSE,
                      .packages=c("gsynth")) %dopar% { # parallel computing
                           
        ## general random sample:
        panel$Y <- panel$Ybar + rnorm(N*T) 
        effect<-apply(matrix(panel$eff,T,N)[,1:Ntr],1,sum)/Ntr # ATT
        
        ## GSC: the number of factors are known to be 2
        out.syn<-gsynth(Y ~ D + X1 + X2, data=panel, index=c("id","time"),
                         se = FALSE, r = 2,CV = FALSE, force="two-way")
        att.syn<-out.syn$att[(T0+1):T]
        bias.syn<-att.syn-effect[(T0+1):T] # bias
        
        ## InterFE
        Y<-matrix(panel$Y,T,N)
        X<-array(0,dim=c(T,N,(2+10)))
        X[,,1]<-matrix(panel$X1,T,N)
        X[,,2]<-matrix(panel$X2,T,N)
        for (t in 1:10) { # only post treatment periods
            X[(T0+t),1:Ntr,(2+t)]<-rep(1,Ntr)
        }
        ## inter_fe is an interval functions that is now hidden in the package                  
        out.int<-inter_fe(Y,X,r=2,beta0=as.matrix(rep(0,12)),
                          force=3, trends=0)
        att.int<-out.int$beta[3:12] # att
        bias.int<-att.int-effect[(T0+1):T] # bias
        
        out<-cbind(att.syn,bias.syn,att.int,bias.int)
        return(out)      
    }
    ## gsynth
    result[case,1,1]<-mean(onecase[5,2,]) # bias, focus on ATT_15
    result[case,2,1]<-sd(onecase[5,1,])   # sd of ATT_15 
    result[case,3,1]<-sqrt(mean((onecase[5,2,])^2)) # RMSE
    
    ## inter
    result[case,1,2]<-mean(onecase[5,4,]) # bias, focus on ATT_15
    result[case,2,2]<-sd(onecase[5,3,])   # sd of ATT_1
    result[case,3,2]<-sqrt(mean((onecase[5,4,])^2)) # RMSE
    
    ## Save after each case
    save(DDsd,NNtr,result,file="sim_inter_result.RData")
  
} # end of all cases
stopCluster(cl) # stop parallel computing
print(Sys.time()-begin.time)


####################
# Output
####################

load("sim_inter_result.RData")

# w=1 only
output<-as.data.frame(cbind(rep(15,8),NNtr,NNco,DDsd[1:12],
                             rep(NA,8),result[1:8,,1],rep(NA,8),result[1:8,,2]))
colnames(output)<-c("T0","NNtr","NNco","DDsd",rep(c("sp","Bias","SD","RMSE"),2))
round(output,3)

sink("output_sim_inter.txt")
print(xtable(as.matrix(output),digits=c(rep(0,5),rep(3,(dim(output)[2]-5+1)))),include.rownames=0)
sink()

