moon_matching <- function(X0,X1,Y0,Y1,m0,m1,...){
    
    # This function calculates the M-out-of-N type bootstrap ATET estimator 
    #
    # Inputs: 
    # - X0, X1: Vector of regressors for controls and treated
    # - Y0, Y1: Vector of outcomes for controls and treated
    # - m0, m1: Size of control and treated groups in bootstrap sample
    #
    # Outputs:
    # - tau_MooN: ATET nearest neighbor estimate based on bootstrap sample 
  
  
    # (1) Extract size of control and treated group
    #------------------------------------------------
    n0 <- length(Y0)
    if(n0 != length(X0)) 
      stop('Number of regressors and outcomes do not match in control group.')
    n1 <- length(Y1)
    if(n1 != length(X1)) 
      stop('Number of regressors and outcomes do not match in control group.')
    
    # Check to make sure that MooN bootstrap group sizes are less than original.
    if(n0 < m0 | n1 < m1)
      stop('The group sizes in the bootstrap are larger than in original data.')
    
    
    # (2) Construct the MooN bootstrap sample
    #----------------------------------------
    # Draw indices for the control and treated unites respectively
    boot_index_control   <- sample(1:n0,m0,replace=TRUE)
    boot_index_treatment <- sample(1:n1,m1,replace=TRUE)
    
    # Construct the regressors and outcomes for both groups
    X0_boot <- X0[boot_index_control]  
    Y0_boot <- Y0[boot_index_control]
    
    X1_boot <- X1[boot_index_treatment]
    Y1_boot <- Y1[boot_index_treatment]
    
    # (3) Do the matching (allowing for ties) in the bootstrap sample
    #----------------------------------------------------------------
    
    # (3)a) Calculate the closest controls to each treated
    
    # Initialize list that contains matched indices for each treated
    matching_indices <- vector("list",length=m1)
    # Initialize vector for number of matched units for eacht treated
    num_index <- rep(0,length = m1)
    
    # Loop over the treated units
    for(i in 1:m1){
      # Determine (squared) distance of the treated to all the controls 
      temp <- (X0_boot - X1_boot[i])^2
      # Find the distance to closest controls
      min_index_value <- temp[which.min(temp)]
      # Store indices of set of closest controls
      matching_indices[[i]] <- c(1:m0)[temp==min_index_value] 
      # Store cardinality of matching set
      num_index[i] <- length(matching_indices[[i]])
    }
    
    # (3)b) Calculate K^*_i - the bootstrap counterpart of the number of times
    # each treated is a match to a treated
    
    # Initialize K^*_i
    K_i_star     <- rep(0,m0) 
    
    # Calculate K^*_i
    for(i in 1:m0){
      for(l in 1:m1){
        K_i_star[i] <- K_i_star[i] + 1/num_index[l]*(i %in% matching_indices[[l]])
      }
    }
    
    # (3)c) Calculate the MooN-type estimator 
    #----------------------------------------------------------------
    
    tau_MooN <- 1/m1*sum(Y1_boot) - 1/m1*sum(K_i_star*Y0_boot)
    
    return(tau_MooN)
}

