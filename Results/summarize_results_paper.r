# Clear memory
rm(list=ls())

# This program collects the estimation results and constructs:
#
#   - Figure 1 of the paper illustrating the performance of the MooN
#     based variance estimators
#   - Table 1 of the paper containing the coverage probabilities of 
#     of the confidence intervals.
#
# The program requires loading the estimation results. Paths to these need 
# adjusting accordingly in lines 40, 43, 46 and 49.


# Load xtable package to construct tables for paper
library(xtable)

# BY HAND: Copy in the settings of the simulation study
#-----------------------------------------------
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

# Construct a list with all the ALPHA settings considered
ALPHA_LIST <- c(0.2,1,2.84,5)


# Read in and combine estimation results for all the ALPHA settings
#------------------------------------------------------------------
results <- list()

# Load data set
load("ATET_MooN_alpha_0.2.RData")
results[[1]] <- estimation_results
rm(estimation_results)
load("ATET_MooN_alpha_1.RData")
results[[2]] <- estimation_results
rm(estimation_results)
load("ATET_MooN_alpha_2.84.RData")
results[[3]] <- estimation_results
rm(estimation_results)
load("ATET_MooN_alpha_5.RData")
results[[4]] <- estimation_results
rm(estimation_results)


#=========================================================================
#
# (1) Construct plots to illustrate point estimates of asymptotic variance
#     in Figure 1 of the paper.
#=========================================================================

len_N  <- length(N_Grid)
len_exponents <- length(EXPONENT_LIST)

# Construct the plots
file_name <- paste("paper_plot_point_estimates.pdf",sep="")

pdf(file=file_name,width=10,height=8,paper="special")
par(mfrow=c(1,4),mar=c(2,2,2,0))
count_alpha <- 0
# Loop over the alphas considered
for (ALPHA in ALPHA_LIST){
  count_alpha <- count_alpha + 1
  
  # Initialize storage of variance estimates
  variance_MooN <- list()
  
  # Extract the results MooN based variance estimates for a fixed alphs
  for(exponent in 1:len_exponents){
    variance_MooN[[exponent]] <- matrix(NA,nrow=SIM,ncol=len_N)
    for(sim in 1:SIM){
      variance_MooN[[exponent]][sim,] <- results[[count_alpha]][[sim]]$var_MooN[,exponent]
    }
  }
 
  # Construct a data matrix of the results, which is amenable for plotting
  var_MooN_plot <- variance_MooN[[1]][,len_N:1]
  for(exponent in 2:len_exponents)
    var_MooN_plot <- cbind(var_MooN_plot,variance_MooN[[exponent]][,len_N:1])
  
  # Construct the boxplots for given ALPHA
  header <- paste0("alpha = ",ALPHA)
  boxplot(var_MooN_plot,
          at=c(1:5,8:12,15:19,22:26,29:33,36:40),horizontal=TRUE,
          main=header,cex.main=1.33,cex.axis=1.33,cex.lab=1,yaxt='n',
          xlim=c(1,40))#,ylim=x_lims_plot[count_alpha,])
  if(count_alpha==1)axis(side=2,at=c(3,10,17,24,31,38),labels=c("MooN(0.5)","MooN(0.6)","MooN(0.7)","MooN(0.8)","MooN(0.9)","Efron"),tck=0,cex.axis=1.33)
  abline(v=1+1.5*ALPHA,lwd=2)
  abline(h=6.5)
  abline(h=13.5)
  abline(h=20.5)
  abline(h=27.5)
  abline(h=34.5)
  # Close alpha loop
}  
dev.off()



#=================================================================
#
# (2) Calculate coverage probabilities of confidence intervals as 
#     in Table 1 of the paper 
#=================================================================



# (2)a) Auxiliary functions
#--------------------------

#Function to calculate the boundaries of the confidence interval
conf_int_gauss_boundary <- function(N,alpha,var,prob,tau_hat){
  # Calculate the scaling
  N1 <- round(alpha/(1+alpha)*N)
  conf_boundary <- tau_hat + qnorm(prob)*sqrt(var/N1)
}

# Function to check whether the confidence interval covers the ATET
check_coverage <- function(low,upp,value){
  covereage_indicator <- (low < value & upp > value)
}


# (2)b) For each alpha calculate the confidence intervals and whether they cover
#------------------------------------------------------------------------------

# Initialize storage keeping track of how many CIs cover
cov_prob <- list()
# Initialize storage to store length of the confidence intervals
CI_length <- list()

count_alpha <- 0
# Loop over the alphas
for (ALPHA in ALPHA_LIST){
  count_alpha <- count_alpha + 1

  # Initialize storage keeping track of how many CIs cover with zeroes
  cov_prob[[count_alpha]] <- matrix(0,nrow=len_N,ncol=len_exponents) 
  
  # Initialize storage to store length of the confidence intervals
  CI_length[[count_alpha]] <- list()

  # Retrieve the ATET estimates 
  estimates_ATET <- matrix(NA,nrow=SIM,ncol=len_N) 
  
  for(sim in 1:SIM){
    estimates_ATET[sim,] <- results[[count_alpha]][[sim]]$tau_hat
  }

  for(exponent in 1:len_exponents){
    CI_length[[count_alpha]][[exponent]] <- matrix(NA,nrow=SIM,ncol=len_N)
    for(sim in 1:SIM){
      for(n_count in 1:len_N){
        # Calculate the confidence intervals
        conf_gauss_low <- conf_int_gauss_boundary(N=N_Grid[n_count],alpha=ALPHA,
                             var=results[[count_alpha]][[sim]]$var_MooN[n_count,exponent],
                             prob=0.025,tau_hat=estimates_ATET[sim,n_count])
        conf_gauss_upp <- conf_int_gauss_boundary(N=N_Grid[n_count],alpha=ALPHA,
                             var=results[[count_alpha]][[sim]]$var_MooN[n_count,exponent],
                             prob=0.975,tau_hat=estimates_ATET[sim,n_count])
        
        cov_prob[[count_alpha]][n_count,exponent] <- cov_prob[[count_alpha]][n_count,exponent] + 
                            check_coverage(low=conf_gauss_low,upp=conf_gauss_upp,value=TAU)  
        CI_length[[count_alpha]][[exponent]][sim,n_count] <- conf_gauss_upp - conf_gauss_low
      # Close the n_count loop
      }
   # Close the sim loop 
   }
  # Close the exponent loop 
  }
  # Compute the coverage probabilities for each fixed alpha
  cov_prob[[count_alpha]] <- 1/SIM*cov_prob[[count_alpha]]
# Close alpha loop
}  

# (2)c) Put all the coverage probabilities in a large table 
#------------------------------------------------------------------------------

cov_prob_table <- cov_prob[[1]]
for(alpha_count in 2:length(ALPHA_LIST))
  cov_prob_table <- rbind(cov_prob_table,cov_prob[[alpha_count]])

# Add the sample sizes as an additional column
cov_prob_table <- cbind(rep(N_Grid,length(ALPHA_LIST)), cov_prob_table)
colnames(cov_prob_table) <- c("N","MooN(0.5)","MooN(0.6)","MooN(0.7)","MooN(0.8)","MooN(0.9)","Efron")

cov_prob_table <- xtable(cov_prob_table)
digits(cov_prob_table) <- c(0,0,rep(3,len_exponents))

print(cov_prob_table,file="Coverage_Table.tex",include.rownames = FALSE)





