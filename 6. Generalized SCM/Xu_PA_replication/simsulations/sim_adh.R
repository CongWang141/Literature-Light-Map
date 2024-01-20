## Generalized Synthetic Control Method
## Replication Materials
## Comparison with Synthetic Control -- ADH (2010)

## Author: Yiqing Xu

## The codes below replicate Table A3 in the Online Appendix

## Varying r and w (warning: takes a long time!)

## test Bais for ATT 15 (treatment period 5)


suppressWarnings(sink())
rm(list=ls(all=TRUE))
Sys.setenv(LANGUAGE='en')
library(Synth)
library(doParallel)
library(parallel)
library(foreach)
library(abind)
library(xtable)
library(gsynth)
source("sim_sampling.R")

sims<-5000
set.seed(1234)

# register multiple cores
cores<-8
Sys.setenv(GOTO_NUM_THREADS=cores)
cl<-makeCluster(cores)
registerDoParallel(cl)

# set progress bar and combind function
f <- function(){
  function(...) {
    abind(...,along=3)  # bind as an array
  }
}

# fixed
T0<-15
T<-T0+10
Nco<-40
Ntr<-1
N<-Ntr+Nco
p<-0     # no covariates, constant not included
D.sd<-1

# for each model, consider 8 cases
RR<-c(1,2,3,4,2,2,2,2)
WW<-c(1,1,1,1,0.75,0.50,0.25,0)
ncases<-length(RR)
# storage
result<-matrix(NA,length(RR),7) # bias sd rmse; gsynth and synth(adh) + synth cannot find solutions
                     
# loop
begin.time<-Sys.time()
cat("Cores: ",cores,"; cases: ",ncases,"; Sims: ",sims,sep="")
for (case in 1:8) {
    
  
      
    ## annouce case
    r<-RR[case]; w<-WW[case]
    cat("\nCase ",case,": r = ",r,", w = ",w,"\n",sep="")

    
    panel<-simulate(Ntr=Ntr,Nco=Nco,T0=T0,p=2,r=2,m=0,w=w,D.sd=D.sd,beta=c(1,3),
                    mu=5,att=c(1:10),fsize=1,FE=1,fixF=TRUE, fixL=TRUE)
    
    ## record the number of no solution cases for Synth
    onecase<-foreach (i=1:sims,.combine=f(),.inorder=FALSE,
                      .packages=c("Synth","gsynth")) %dopar% { # parallel computing
      
      output<-matrix(NA,10,5)
      output[,5]<-0
      
      # general random sample: (5tr+45co)*(10pre+10post)
      panel$Y <- panel$Ybar + rnorm(N*T) 
      effect<-apply(as.matrix(matrix(panel$eff,T,N)[,1:Ntr]),1,sum)/Ntr # SATT
    
      # Gsynth: correct specification
      out.gsyn<-gsynth(Y ~ D, data=panel,index = c("id","time"),
                                           se = FALSE,r=r, CV=0, force="two-way")
      output[,1]<-out.gsyn$att[(T0+1):T] # att
      output[,2]<-out.gsyn$att[(T0+1):T]-effect[(T0+1):T] # bias
      #cat("-")
      
      # Synth (ADH2010)
      tr.id<-unique(panel$id)[1:Ntr]
      co.id<-unique(panel$id)[(Ntr+1):N]
      time<-unique(panel$time)
      pre<-time[1:T0]
      
      # Run Synth: see if there's an error
      synth.ite<-matrix(NA,T,Ntr) # to store individual treatment effect
      
      possibleError <- tryCatch({
        
      for (j in 1:Ntr) { # synth loop for each treated unit
            dataprep.out<-
              dataprep(
                foo = panel,
                dependent = "Y",
                unit.variable = "id",
                time.variable = "time",
                special.predictors = list(
                  list("Y", 1, "mean"),
                  list("Y", 8, "mean"),
                  list("Y", 15, "mean")
                ),
                treatment.identifier = tr.id[j],
                controls.identifier = co.id,
                time.predictors.prior = pre,
                time.optimize.ssr = c(2:6,9:14),
                time.plot = time
              )
            out.synth <- synth(dataprep.out)
            synth.ite[,j]<- dataprep.out$Y1plot-(dataprep.out$Y0plot%*%out.synth$solution.w)
            #cat("+")                  
      } # end of synth loop 
      }, error=function(e) e
      )
      
      if(inherits(possibleError, "error")) { # if wrong, skip
        output[,5]<-1
      } else { # if correct
        output[,3]<-apply(synth.ite,1,mean)[(T0+1):T]  # att
        output[,4]<-output[,3]-effect[(T0+1):T] # bias
      }
      
     return(output)
            
  } # end of one simulation
  
  # gsynth
  result[case,1]<-mean(onecase[5,2,]) # bias, focus on SATT_15
  result[case,2]<-sd(onecase[5,1,])   # sd of SATT_15 
  result[case,3]<-sqrt(mean((onecase[5,2,])^2))
  # synth
  result[case,4]<-mean(onecase[5,4,],na.rm=TRUE) # bias, focus on SATT_15
  result[case,5]<-sd(onecase[5,3,],na.rm=TRUE)  # sd of SATT_1
  result[case,6]<-sqrt(mean((onecase[5,4,])^2,na.rm=TRUE))
  result[case,7]<-mean(apply(onecase[,5,],2,mean))
  
    
  # Save after each case
  save(result,file="sim_adh_result.RData")
    
  
} # end of 8 case
stopCluster(cl) # stop parallel computing
print(Sys.time()-begin.time)

####################
# Output
####################

load("sim_adh_result.RData")
output<-as.data.frame(cbind(rep(15,ncases),rep(Ntr,ncases),rep(Nco,ncases),RR,WW,
                            rep(NA,ncases),result[,1:3],rep(0,ncases),rep(NA,ncases),result[,4:7]))

colnames(output)<-c("T0","Ntr","Nco","r","w",rep(c("sp","Bias","SD","RMSE","Fail"),2))
round(output,3)

library(xtable)
sink("output_sim_adh.txt")
print(xtable(output,digits=c(rep(0,5),2,rep(c(1,3,3,3,3),2))),include.rownames=FALSE)
sink()


