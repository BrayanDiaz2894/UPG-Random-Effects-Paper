
# Gibbs Sampling

 start_time <- Sys.time()
for (iter in 1:iterations) {
  print(iter*100/iterations)
  print("Basic")
  omega_prev <- if (iter == 1) matrix(1, nrow = n, ncol = t_obs) else omega_store[iter - 1, , ]
  
  # Step 1: Sample alpha_i (random effects)
  for (i in 1:n) {
    
    # 1. Posterior precision
    B_i <- V_alpha_inv + sum(omega_prev[i, ]) * ZZ_precomp[i,,]
    
    # 2. Posterior mean component
    temp <- Y[i, ] - 0.5 - omega_prev[i, ] * as.vector(X[i, , ] %*% Beta)
    A_i <- sum(temp) * Z[i, ]
    
    # 3. Posterior covariance and mean
    cov_alpha <- chol_inverse(B_i)  # Σ = B⁻¹
    mu_i <- cov_alpha %*% A_i           # μ = Σ * A_i
    
    # 4. Cholesky-based sampling
    R <- chol(cov_alpha)                # chol() returns R such that Σ = R'R
    z <- rnorm(q)                       # standard normal vector
    Alpha[i, ] <- as.numeric(mu_i + t(R) %*% z)  # αᵢ = μ + Rᵗ z
  }
  
  Alpha_save[iter, , ] <- Alpha
  

  # Step 2: Sample beta (fixed effects)

  # Vectorized computation for Sigma_beta_inv:
  omega_vec <- as.vector(omega_prev)
  # This first trick comes from this identity: 
  # \left(\sqrt{\boldsymbol{\omega}} \cdot \boldsymbol{X}\right)^\top\left(\sqrt{\boldsymbol{\omega}} \cdot \boldsymbol{X}\right)\;=\;\sum_{i,t} \omega_{it} \, \boldsymbol{X}_{it} \boldsymbol{X}_{it}^\top
  Sigma_beta_inv_data <- t(X_flat * sqrt(omega_vec)) %*% (X_flat * sqrt(omega_vec))
  Sigma_beta_inv <- Sigma0_inv + Sigma_beta_inv_data
  
  # For mu_beta, first compute the contribution from the data.
  #We aim to \mu_\beta = \Sigma_\beta \left( \sum_{i=1}^n \sum_{j=1}^t \left( y_{ij} - \frac{1}{2} - \omega_{ij} \cdot (Z_i^\top \alpha_i) \right) X_{ij} + \Sigma_0^{-1} \mu0 \right)

  # For each individual, compute the scalar sum_{k} (Z[i,] * Alpha[i,]) #Note that here only the dimmension j is being summed,
  # So in the double summation we can take it out. 
  dot_Z_alpha <- rowSums(Z * Alpha)  # length n
  
  # Replicate dot_Z_alpha across t columns to match Y and omega_prev.
  # Each row i is repeated t times.
  dot_Z_alpha_mat <- matrix(rep(dot_Z_alpha, each = t_obs), nrow = n, ncol = t_obs, byrow = TRUE)
  #This cretes a matrix where dot_Z_apha is the same for all individuals across t times. We will use it in a moment. 
  
  # Compute the t-by-n matrix of scalars for each (i,j):
  #Here y, omega_prev and dot_Z_alpha_mat are nxt, but dot_Z_alpha_mat has the same values across all columns (time invariant)
  #Then the next expression is 0.5*(2*Y_ij-1)*omega_i,j*Z_ij this is, we have a matrix where the rows is the individual and columns time, for the expression mentioned before.
  temp_mat <- 0.5 * (2 * Y - 1) - omega_prev * dot_Z_alpha_mat
  
  # Flatten temp_mat to a vector of length n*t:
  #We vectorize it, so we dont have 5 collumns any more.
  temp_vec <- as.vector(temp_mat)
  
  # Data contribution to mu_beta: sum_{i,j} (temp_scalar[i,j] * X[i,j,])
  mu_beta_data <- t(X_flat) %*% temp_vec
  
  # Finally, combine with the prior contribution:
  mu_beta <- chol_inverse(Sigma_beta_inv) %*% (mu_beta_data + Sigma0_inv %*% mu0)
  Sigma_beta <- chol_inverse(Sigma_beta_inv)  # posterior covariance
  R <- chol(Sigma_beta)                            # upper-triangular: Sigma_beta = R'R
  z <- rnorm(p)                                    # standard normal vector
  Beta <- as.numeric(mu_beta + t(R) %*% z) # correct draw
  Beta_save[iter, ] <- Beta
  
  # Step 3: Sample omega_ij
  Zalpha_vec <- dot_Z_alpha[row_id]             # ✓ FIX  use 'times = t', not 'each = t'
  XB_vec <- X_flat %*% Beta
  psi_vec      <- XB_vec + Zalpha_vec
  omega_vec    <- rpg(n*t_obs, h = 1, z = psi_vec)     # BayesLogit::rpg is vectorised
  omega_matrix <- matrix(omega_vec, nrow = n, ncol = t_obs, byrow = FALSE)
  omega_store[iter, , ] <- omega_matrix
  omega_prev            <- omega_matrix
  
  # Step 4: Sample V_alpha
  V_Alpha <- riwish(nu0 + n, Lambda0 + t(Alpha) %*% Alpha)
  V_alpha_save[iter, , ] <- V_Alpha
  V_alpha_inv <- chol_inverse(V_Alpha)
}
  end_time <- Sys.time()
  print(end_time - start_time)
  timesamples_naug <- time_to_seconds(end_time - start_time) 
  

    
 
#4. Sample V_\alpha

   # Extract samples for each element of V_alpha
   V_11_samples <- V_alpha_save[(burnin + 1):iterations, 1, 1]
   V_12_samples <- V_alpha_save[(burnin + 1):iterations, 1, 2]
   V_21_samples <- V_alpha_save[(burnin + 1):iterations, 2, 1]
   V_22_samples <- V_alpha_save[(burnin + 1):iterations, 2, 2]

  

#5. Plot Beta


   # Extract beta samples for each parameter
   beta_1_samples <- Beta_save[(burnin + 1):iterations, 1]
   beta_2_samples <- Beta_save[(burnin + 1):iterations, 2]
   beta_3_samples <- Beta_save[(burnin + 1):iterations, 3]



   

   a1_nag <- Alpha_save[(burnin + 1):iterations, 1, ]
   a2_nag <- Alpha_save[(burnin + 1):iterations, 2, ]
   beta_1_samplesnaug <- beta_1_samples
   beta_2_samplesnaug <- beta_2_samples
   beta_3_samplesnaug <- beta_3_samples
   V1_naug <- V_11_samples
   V2_naug <- V_22_samples
   V21_naug <- V_21_samples
   alpha_samples_naug <- Alpha_save
   
   unbalance <- abs(beta_true[1])
   
   filename <- paste0(
     "Output/V_samplesnaug_",
     n, "_", t_obs, "_", iterations, "_", unbalance, ".RData"
   )
   
   save(a1_nag, a2_nag, V21_naug, V2_naug, V1_naug, beta_1_samplesnaug, beta_2_samplesnaug, beta_3_samplesnaug,timesamples_naug, alpha_samples_naug, file = filename)
   rm(a1_nag, a2_nag, V21_naug, V2_naug, V1_naug, beta_1_samplesnaug, beta_2_samplesnaug, beta_3_samplesnaug,timesamples_naug, alpha_samples_naug)
 