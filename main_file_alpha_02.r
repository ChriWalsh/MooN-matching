# Clear memory
rm(list=ls())

# This program forms the basis of the simulation study in 
#
# Walsh, C. and Jentsch, C. (2022+), Nearest neighbor matching: Does 
# direct M-out-of-N bootstrapping without bias correction work when the naive 
# naive bootstrap fails?
#
# The program runs a simulation study for a fixed ratio of treated to control
# given by alpha = 0.2.
#
# Inputs:
# -------
# 
# ALPHA         - ratio of treated to controls 
#
# SIM           - number of simulation runs
# N_Grid        - vector containing values for sample size
# B             - integer specifying the number of replications when resampling
# TAU           - the ATET for the simulation
# EXPONENT_LIST - vector with gammas to set the M for the M-out-of-N
#
# Outputs:
# --------
# tau_hat       - vector of estimated ATET for the different sample sizes  
# var_MooN      - the resampling based MooN-type variance estimates for the 
#                 different sample sizes and (for the different values of M 
#                 corresponding to the supplied gammas)  


#=================================================
#
# Step 0: Load necessary packages and functions and set 
#         parameters of simulation study.
#
#==================================================

library("parallel")

source("moon_matching.r")
# Contains the function 
#  
#moon_matching(X0,X1,Y0,Y1,m0,m1,...)
#
#, which calculates the M-out-of-N bootstrap estimator of the ATET. 

source("simulation_function.r")
# Contains the function 
#
# simulation_function(N_Grid,ALPHA,B,TAU,EXPONENT_LIST,...)
#
#, which does one entire simulation run.

# Set parameters
#-------------------

# Set the ratio of treated to controls
ALPHA <- 0.2

# Number of simulation runs
SIM <- 1000
# Number of bootstrap replications
B <- 1000
# Treatment effect
TAU <- 1
# Number or observations (N=N_1+N_0)
N_Grid <-c(100,250,500,1000,2000)
# Set the exponents for the MooN
EXPONENT_LIST <- c(0.5,0.6,0.7,0.8,0.9,1)

#=======================================================
#
# Step 1: Run the simulation study and store the results
#
#=======================================================

# Set-up the computing cluster for the simulation
cores.avail <- detectCores()
mc <- cores.avail-1
cl <- makeCluster(mc)

# Set the seed 
clusterSetRNGStream(cl,123)
# Pass the required functions to the cluster
clusterExport(cl=cl, varlist=c("moon_matching","simulation_function"), envir=.GlobalEnv)

# Run the simulation study
#----------------------------
estimation_results <- parLapply(cl,1:SIM,simulation_function,N_Grid=N_Grid,B=B,ALPHA=ALPHA,TAU=TAU,EXPONENT_LIST=EXPONENT_LIST) 
stopCluster(cl)


# Save the results from the simulation study
file_name <- paste("ATET_MooN_alpha_",ALPHA,".RData",sep="")
save(estimation_results,file=file_name)


