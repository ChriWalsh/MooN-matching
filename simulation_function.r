simulation_function <- function(N_Grid,ALPHA,B,TAU,EXPONENT_LIST,...){
  
  # This function does one simulation run for a fixed treated to control ratio 
  # and various sample sizes and choices of bootstrap sample sizes.   
  #
  # Inputs:
  # -------
  #
  # N_Grid        - vector containing values for sample size
  # ALPHA         - ratio of treated to controls  
  # B             - integer specifying the number of replications when resampling
  # TAU           - the ATET for the simulation
  # EXPONENT_LIST - vector with gammas to set the M for the M-out-of-N
  #
  # Outputs:
  # --------
  # tau_hat       - the estimataed ATET  
  # var_MooN      - the resampling based MooN-type variance estimates (for 
  #                 different values of M)
 
 
  # Initialize storage variables for the estimated quantities for
  # all data sample sizes and bootstrap resample sizes.
  #--------------------------------------------------------
  
  # The nearest-neighbor matching estimator of the ATET
  tau_hat   <- rep(NA,length=length(N_Grid))
  
  # MooN-type variance estimator for all exponents
  var_MooN <- matrix(NA,nrow = length(N_Grid),ncol = length(EXPONENT_LIST))
  
  #==========================================================
  #
  # Step 0: Set-up parameters for the simulations
  #
  #======================================================= 
  
  # Determine the largest sample size to be considered
  max_N <- max(N_Grid)
  
  # ATET
  TAU <- TAU
  
  # Ratio of treated to control
  alpha <- ALPHA
  
  # Simulate the covariate
  X <- runif(max_N)
  
  # Simulate the potential outcomes 
  Y0 <- rnorm(max_N) 
  Y1 <- rep(TAU,max_N)
  
  
  #===================================================
  #
  # Start the actual simulation loop
  #
  #==================================================== 
  
  # Set a counter to keep track of the simulation design
  counter_design <- 0
  
  # Loop over the possible sample sizes
  for(n in N_Grid){
  
    # Increment counter tracking simulation design
    counter_design <- counter_design + 1
      
    #===========================================================
    #
    # Step 1: Simulate the data given number of observations and 
    #         ratio of treated to controls
    #
    #=========================================================== 
      
      
    # Calculate the number of treated and controls according to the simulation design
    #========================================================================
    n1 <- round(alpha/(1+alpha)*n)
    n0 <- n - n1
      
    # Construct the data for the simulation run given the design parameters
    #------------------------------------------------------------------
      
    # Set the treatment indicator 
    W_sim <- rep(0,n)
    W_sim[1:n1] <- 1
  
    # Choose the regressors for the treated and control
    X1_sim <- X[1:n1]
    X0_sim <- X[max_N:(max_N-n0+1)]
    # Combine to get vector with regressor values
    X_sim  <- c(X1_sim,X0_sim)
  
    # Set the potential outcomes and the observed outcome
    Y0_sim <- Y0[max_N:(max_N-n0+1)] 
    Y1_sim <- Y1[1:n1]
    Y_sim  <- c(Y1_sim,Y0_sim)
      
    #===========================================================
    #
    # Step 2: Calculate the nearest neighbor matching estimator
    #
    #=========================================================== 
      
    # Index for the nearest control unit
    min_index <- rep(NA,n1)
    # For each treated find the closest control.  
    # Note: In the original data there will be no ties.
    for(i in 1:n1){
      temp <- (X0_sim - X1_sim[i])^2
      min_index[i] <- which.min(temp)  
    }
      
    # Calculate the K_i term
    K_i <- tabulate(min_index,nbins=n0)
      
    # Calculate the individual contribution to the matching estimator for all the sample units.
    ATET_hat_i           <-  rep(0,length=n)
    ATET_hat_i[W_sim==1] <-  Y_sim[W_sim==1]
    ATET_hat_i[W_sim==0] <- -K_i*Y_sim[W_sim==0]
      
    # Calculate the ATE estimator
    ATET_hat <- 1/n1*sum(ATET_hat_i)
      
    # Store the ATE estimator
    tau_hat[counter_design] <- ATET_hat 
      
   
    #-------------------------------------------------------------------------
    #
    # Step 3: Calculate the M-out-of-N bootstrap estimator
    #
    #-------------------------------------------------------------------------
      
    # Loop over the exponent
    counter_exponent <- 0
      
    for(exponent in EXPONENT_LIST){
     counter_exponent <- counter_exponent + 1
        
     # Initialize storage for the Resampling based estimator of the ATET
     tau_MooN <- rep(NA,length=B) 
        
     # Set the parameters for the bootstrap 
     m  <- round(n^exponent)
      
     # Actual alpha in the data
     alpha <- n1/n0
        
     m1 <- round(alpha/(1+alpha)*m)
     m0 <- m-m1 
       
     if(m1==0) stop('No units in MooN treated sample.')
     if(m0==0) stop('No units in MooN control sample.')
       
     # Do the bootstrap   
     for(b in 1:B){
      tau_MooN[b] <- moon_matching(X0_sim,X1_sim,Y0_sim,Y1_sim,m0,m1)
      # Close the bootstrap loop   
      }  
        
      # Store the MooN-type bootstrap variance estimator
      var_MooN[counter_design,counter_exponent] <- m1 * mean((tau_MooN-rep(ATET_hat,length=B))^2)
      
    # Close the exponent loop
    }  
      
  # Close the loop over the sample size
  }    
  
  return(list(tau_hat=tau_hat,var_MooN=var_MooN))
} 
        







