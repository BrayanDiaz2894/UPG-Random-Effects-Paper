
start_time <- Sys.time()
# VI. Gibbs Sampler Loop
for (iter in 1:iterations) {
  print(iter*100/iterations)
  print("Hybrid")
  ## Step 1. Sample h given Y, beta, alpha, omega
  # --- Vectorized draw of all h[i,j] ---------------------
  
  ## Step 1. Sample h given Y, beta, alpha, omega via inverse‐CDF
  
  # 1) (X_mat built outside)
  # 2) fixed‐effects part
  xb_vec    <- X_mat %*% Beta        # length N = n*t_obs
  
  # 3) random‐effects part
  alpha_eff <- rowSums(Z * Alpha)    # length n
  alpha_vec <- rep(alpha_eff, times = t_obs) # length N
  
  # 4) full mean
  mu_vec    <- xb_vec + alpha_vec            # length N
  
  # 5) vector of s.d.’s
  sd_vec    <- as.vector(sqrt(1/omega))      # length N (or length n, recycled as in your original)
  # if omega is length n, do:
  # sd_vec <- rep(sqrt(1/omega), times = t_obs)
  
  # 6) flatten Y
  y_vec     <- as.vector(Y)                  # length N
  
  # 7) inverse‐CDF sampling
  # 7a) standardised truncation point a = (0 - mu)/sd
  a_vec     <- - mu_vec / sd_vec
  
  # 7b) lower‐tail CDF at zero
  F_vec     <- pnorm(a_vec)
  
  # 7c) draw uniforms on the correct intervals
  u_vec           <- numeric(length = length(mu_vec))
  idx1            <- which(y_vec == 1)
  idx0            <- which(y_vec == 0)
  u_vec[idx1]     <- F_vec[idx1] + runif(length(idx1)) * (1 - F_vec[idx1])
  u_vec[idx0]     <-                runif(length(idx0)) *   F_vec[idx0]
  eps <- .Machine$double.eps          # o usa .Machine$double.eps si quieres el mínimo representable
  u_vec[u_vec == 1] <- u_vec[u_vec == 1] - eps
  
  # 7d) invert to get truncated‐Normal draws
  h_vec <- mu_vec + sd_vec * qnorm(u_vec)
  
  # 8) reshape back into an n×t_obs matrix
  h     <- matrix(h_vec, nrow = n, ncol = t_obs)
  
  
  ## Step 2. Sample omega from PG(1, h - (X*beta + Z*alpha))
  
  # 1) Flatten X into a (n*t_obs) × p matrix, so each row is X[i,j,]
  #done outside the loop.
  
  # 2) Fixed-effect part: X β
  xb_vec    <- X_mat %*% Beta           # length = n*t_obs
  
  # 3) Random-effect part: αᵢ effect replicated for each j
  alpha_eff <- rowSums(Z * Alpha)       # length n
  alpha_vec <- rep(alpha_eff, times = t_obs)    # length = n*t_obs
  
  # 4) Compute the vector of “z_val”’s
  h_vec     <- as.vector(h)                     # column-major flatten of h[i,j]
  lp_vec    <- xb_vec + alpha_vec               # linear predictor at each (i,j)
  z_vec     <- h_vec - lp_vec
  
  # 5) One call to rpg for all n*t_obs draws
  omega_vec <- rpg(
    n = length(z_vec), 
    h = 2, 
    z = abs(z_vec)
  )
  
  # 6) Reshape back into your n × t_obs matrix
  omega     <- matrix(omega_vec, nrow = n, ncol = t_obs)
  
  ## Step 3. Sample alpha_i for each individual i
  # 1) linear predictor for fixed effects
  xb_vec   <- X_mat %*% Beta           # length = n*t_obs
  xb_mat   <- matrix(xb_vec, nrow = n, ncol = t_obs)
  
  # 2) residual matrix
  resid_mat   <- h - xb_mat                    # n × t_obs
  
  # 3) row-sums for omega and for omega * residual
  sum_omega   <- rowSums(omega)                # length = n
  sum_resid   <- rowSums(omega * resid_mat)    # length = n
  
  # 4) build A_mat where row i = A_i = Z[i,] * sum_resid[i]
  A_mat       <- Z * sum_resid                 # n × q
  
  # 5) posterior precision constant
  V_alpha_inv <- chol_inverse(V_alpha)        # q × q
  for(i in 1:n) {
    # Posterior precision and covariance
    Bi      <- V_alpha_inv + sum_omega[i] * ZZ_precomp[i,,]
    cov_i   <- chol_inverse(Bi)
    
    # Posterior mean
    mu_i <- cov_i %*% A_mat[i, ]
    
    # Cholesky-based draw (correct)
    R <- chol(cov_i)                    # R such that cov = R'R
    z <- rnorm(q)
    Alpha[i, ] <- as.numeric(mu_i + t(R) %*% z)
  }
  
  
  
  ## Step 4. Sample V_alpha from its full conditional:
  # V_alpha | ... ~ IW(nu0 + n, Lambda0 + sum_{i=1}^n alpha_i alpha_i^T)
  # Sum of outer­products of each row of Alpha:
  S_alpha <- crossprod(Alpha)  
  # (this is t(Alpha) %*% Alpha)
  
  # Then draw V_alpha as before:
  V_alpha <- riwish(nu0 + n, Lambda0 + S_alpha)
  
  
  ## Step 5. Sample beta from its full conditional:
  # beta | ... ~ N(mu_beta, Sigma_beta)
  
  omega_vec<- as.vector(omega)                        # length n*t_obs
  h_vec    <- as.vector(h)                            # length n*t_obs
  
  # 1) Flattened alpha-effect and residual
  alpha_eff <- rowSums(Z * Alpha)              # length n
  alpha_vec <- rep(alpha_eff, times = t_obs)          # length n*t_obs
  resid_vec <- h_vec - alpha_vec                      # length n*t_obs
  
  # 2) Precision update: Σβ⁻¹ = Sigma0_inv + ∑ ω_{ij} x_{ij} x_{ij}ᵀ
  Sigma_beta_inv <- Sigma0_inv +
    t(X_mat) %*% (X_mat * omega_vec)
  
  # 3) Sum term: ∑ ω_{ij} x_{ij} resid_{ij}
  sum_term <- as.numeric(
    t(X_mat) %*% (omega_vec * resid_vec)
  )
  
  # 4) Posterior covariance and mean
  Sigma_beta <- chol_inverse(Sigma_beta_inv)
  mu_beta    <- Sigma_beta %*% (Sigma0_inv %*% mu0 + sum_term)
  
  # 5) One MVN draw
  R <- chol(Sigma_beta)             # R is upper triangular such that Sigma_beta = R'R
  z <- rnorm(p)                     # Standard normal vector of length p
  Beta <- as.numeric(mu_beta + t(R) %*% z)
  
  
  
  # Save samples:
  Beta_save[iter, ] <- Beta
  Alpha_save[iter, , ] <- Alpha
  V_alpha_save[[iter]] <- V_alpha
}

end_time <- Sys.time()
print(end_time - start_time)
timesamples_aug <- end_time - start_time

save.image(file = paste0("workspace_augmented", t_obs, "obs", n, "_iter", iterations, ".RData"))

#load("V_samples.RData")


# --- Plot for Beta ---
# Extract beta samples for each parameter (post burn-in)
beta_1_samples <- Beta_save[(burnin + 1):iterations, 1]
beta_2_samples <- Beta_save[(burnin + 1):iterations, 2]
beta_3_samples <- Beta_save[(burnin + 1):iterations, 3]
mean(beta_1_samples)



# Assume V_alpha_save is your list of matrices
# from the Gibbs sampler and V_alpha_true is the true 2x2 matrix.

# Extract the number of iterations from your list of samples

# Extract components: note that V_alpha is symmetric, so [1,2] and [2,1] are the same.
v11 <- sapply(V_alpha_save, function(mat) mat[1, 1])
v12 <- sapply(V_alpha_save, function(mat) mat[1, 2])
v22 <- sapply(V_alpha_save, function(mat) mat[2, 2])



beta_1_samplesaug <- beta_1_samples
beta_2_samplesaug <- beta_2_samples
beta_3_samplesaug <- beta_3_samples
V1_aug <- v11
V2_aug <- v22
V21_aug <- v12
a1_ag <- Alpha_save[(burnin + 1):iterations, 1, ]
a2_ag <- Alpha_save[(burnin + 1):iterations, 2, ]
alpha_samples_aug <- Alpha_save

unbalance <- abs(beta_true[1])

filename <- paste0(
  "Output/V_samplesaug_",
  n, "_", t_obs, "_", iterations, "_", unbalance, ".RData"
)

save(a1_ag,a2_ag,V21_aug, V2_aug, V1_aug, beta_1_samplesaug, beta_2_samplesaug, beta_3_samplesaug, alpha_samples_aug, timesamples_aug, file = filename)





