Replication Files for:

Xu, Yiqing, "Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models," Political Analysis

* Program Files *

gsynth_1.0.tar.gz	

Please install the package in R (in RStudio, Tools ¡ª Install Packages ¡ª Install from File). The program files include sources files written in R (i.e., gsynth.R, interFE.R) and C++ (i.e., xtinter.cpp) and two datasets. (In Windows Operating System, you may need to install the "Rtools¡±, ¡°Rcpp" and "RcppArmadillo" packages first before installing "gsynth".)

gsynth is dependent on the following packages
(1) doParallel_1.0.10
(2) parallel_1.4.0
(3) ggplot2_2.1.0
(4) GGally_1.0.1

* Replication Files for the Two Examples *

(1) gsynth.RData	storing two datasets: "simdata" and "turnout" (redundant since they are already included in the gsynth package)
(2) rep_simdata.R 	replicating results using the simulated data
(3) rep_turnout.R 	replicating results of the EDR example


* Replication Files for Simulations * 

(1) sim_sampling.R	Generating simulated samples
(2) sim_TN.R		Finite sample properties
(3) sim_DID.R		Comparison with diff-in-diffs
(4) sim_inter.R		Comparison with interactive fixed effects
(5) sim_adh.R		Comparison with synthetic matching
(6) sim_factor.R	Choosing the correct number of factors
(7) sim_coverage.R	Coverage using the parametric bootstrap procedure
(8) FLSource.RData	Storing fixed factors and factor loadings

The above simulation programs are dependent on the following packages:
(1) Synth_1.1.5
(2) doParallel_1.0.10
(3) parallel_1.4.0
(4) foreach_1.4.3
(5) abind_1.4.3
(6) xtable_1.8.2
(7) MASS_7.3.45
(8) plm_1.4.0

