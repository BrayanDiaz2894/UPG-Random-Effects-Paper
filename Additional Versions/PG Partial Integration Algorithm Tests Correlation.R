#-------------------------------------------------------------
# Gibbs Sampler with Random Effects for Logistic Posterior Estimation
#-------------------------------------------------------------

# I. Load required libraries
library(BayesLogit)   # For Polya-Gamma sampling (rpg)
library(MASS)         # For multivariate normal (mvrnorm)
library(MCMCpack)     # For inverse-Wishart sampling (riwish)
library(truncnorm)    # For truncated normal sampling
library(coda)         # For MCMC diagnostics (if needed)
library(ggplot2)      # For plotting (if needed)

# A helper function for matrix inversion using Cholesky (if needed)
inversechol <- function(mat, reg = 0) {
  mat_reg <- mat + reg * diag(nrow(mat))
  L <- chol(mat_reg)  # L such that t(L) %*% L = mat_reg
  chol2inv(L)
}

sample_truncnorm_inverse <- function(a, b, mean = 0, sd = 1, eps = 1e-10) {
  # Step 1: Standardize the bounds
  a_star <- (a - mean) / sd
  b_star <- (b - mean) / sd
  
  # Step 2: CDF of standardized bounds
  Fa <- pnorm(a_star)
  Fb <- pnorm(b_star)
  
  # Defensive check: if interval is too narrow, return midpoint
  if (abs(Fb - Fa) < eps) {
    warning("CDF bounds too close. Returning midpoint.")
    z <- (a_star + b_star) / 2
  } else {
    # Step 3: Sample from Uniform(Fa, Fb)
    u <- runif(1, min = Fa, max = Fb)
    
    # Step 4: Clip u to avoid Inf in qnorm
    u <- min(max(u, eps), 1 - eps)
    
    # Step 5: Invert
    z <- qnorm(u)
  }
  
  # Step 6: Transform back
  x <- mean + sd * z
  return(x)
}

#' Invert a symmetric positive-definite matrix via Cholesky

chol_inverse <- function(M) {
  # Cholesky factor R:  M = Rᵀ R   (R is upper-triangular)
  R <- chol(M)
  
  # Fast inversion via chol2inv:  (Rᵀ R)^{-1} = R^{-1} R^{-ᵀ}
  M_inv <- chol2inv(R)
  
  return(M_inv)
}

#' Fast multivariate-normal sampler via Cholesky


rmvnorm_chol <- function(mu, Sigma, n_draw = 1)
{
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma))
    stop("'Sigma' must be a square matrix")
  d <- length(mu)
  if (d != nrow(Sigma))
    stop("length(mu) must match nrow(Sigma)")
  
  R <- chol(Sigma)                       # Σ = Rᵀ R  (R upper-triangular)
  Z <- matrix(rnorm(n_draw * d), nrow = n_draw)   # N(0, I) draws
  X <- Z %*% t(R)                        # N(0, Σ)
  X <- sweep(X, 2, mu, "+")              # add mean
  
  return(X)
}


# II. Set seed for reproducibility
set.seed(123)

# III. Set simulation parameters
n      <- 100  # Number of individuals
t_obs  <- 30    # Observations per individual
p      <- 3     # Number of fixed-effect predictors (including intercept)
q      <- 2     # Number of random-effect covariates (random effects dimension)

# IV. Simulate data
# 1. Simulate predictor matrix X (n x t_obs x p) and ensure first column is intercept = 1
X <- array(rnorm(n * t_obs * p), dim = c(n, t_obs, p))
X[,,1] <- 1  # first predictor is intercept (set to 1 for all observations)

# 2. Simulate individual-level covariates Z for random effects (n x q)
Z <- matrix(rnorm(n * q), nrow = n, ncol = q)

# 3. Set "true" covariance for random effects and simulate true random effects alpha_true

V_alpha_true <- matrix(c(1, 0.5,   # Covariance matrix (q x q) for random effects
                         0.5, 3), nrow = q, ncol = q)
alpha_true <- mvrnorm(n = n, mu = rep(0, q), Sigma = V_alpha_true)  # random effects for each individual (n x q)

# 4. Set "true" fixed effects beta_true
beta_true <- c(-7, 2, -1.5)  # length p = 3 (intercept, X2, X3)

# 5. Generate latent linear predictor eta and binary response Y
eta <- matrix(0, nrow = n, ncol = t_obs)  # linear predictor for each observation
for (i in 1:n) {
  # Calculate eta[i,] = X[i,,] %*% beta_true + Z[i,] %*% alpha_true[i,]
  # (Note: Z[i,] %*% alpha_true[i,] gives a scalar, apply to all obs of individual i)
  eta[i, ] <- X[i,, ] %*% beta_true + as.numeric(Z[i, ] %*% alpha_true[i, ])
}
# Generate binary outcomes Y from logistic model: P(Y=1) = logistic(eta)
prob <- 1 / (1 + exp(-eta))
Y <- matrix(rbinom(n * t_obs, size = 1, prob = prob), nrow = n, ncol = t_obs)

prop_ones <- mean(Y)               # proporción en [0,1]
pct_ones  <- 100 * prop_ones       # porcentaje

cat(sprintf("Proporción de 1's: %.3f  (%.2f%%)\n",
            prop_ones, pct_ones))


# V. Set up priors and initial values for the Gibbs sampler
# 1. Priors for fixed effects beta ~ N(mu0, Sigma0)
mu0     <- rep(0, p)
Sigma0  <- diag(p)          # relatively diffuse prior (variance 10 on diagonal)
Sigma0_inv <- solve(Sigma0)

# 2. Prior for gamma (global location parameter) ~ N(0, G0)
G0    <- 10                      # variance of gamma's prior (could be tuned as needed)

# 3. Prior for random effects covariance V_alpha ~ Inverse-Wishart(nu0, Lambda0)
nu0     <- q + 3                 # degrees of freedom for IW prior (q+2 ensures finite mean)
Lambda0 <- diag(2, q)               # scale matrix for IW prior (identity for simplicity)

# 4. Initial values for parameters
Beta    <- rep(0, p)             # start fixed effects at 0
Alpha   <- matrix(0, n, q)       # start random effects at 0 for all individuals
V_alpha <- diag(q)               # start covariance as identity matrix
gamma   <- 0                     # start gamma (location expansion parameter) at 0

# 5. Pre-compute X_all (flattened design matrix) for efficiency in Beta updates
N_total <- n * t_obs
X_all   <- matrix(NA, nrow = N_total, ncol = p)
for (i in 1:n) {
  # Fill rows for individual i
  idx_start <- (i-1)*t_obs + 1
  idx_end   <- i * t_obs
  X_all[idx_start:idx_end, ] <- X[i,, ]
}

# VI. Prepare storage for MCMC samples
iterations <- 1500
burnin     <- floor(iterations * 0.1)    # 10% burn-in
Beta_save    <- matrix(NA, nrow = iterations, ncol = p)    # store beta draws
Alpha_save   <- array(NA, dim = c(iterations, n, q))       # store alpha draws
V_alpha_save <- vector("list", length = iterations)        # store V_alpha draws (as matrices)
gamma_save   <- numeric(iterations) 
# --- NUEVO: almacenamiento para h y omega --------------------
h_save     <- array(NA_real_, dim = c(iterations, n, t_obs))   # iter × n × t_obs
omega_save <- array(NA_real_, dim = c(iterations, n, t_obs))   # iter × n × t_obs

# store gamma draws

# VII. Gibbs Sampler loop
X_flat <- matrix(X, nrow = n * t_obs, ncol = p)


matriz <- matrix(data = 1:(3 * iterations), nrow = iterations, ncol = 3)

start_time <- Sys.time()
#iter <- 2
for (iter in 1:iterations) {
  
  print(iter)
  # 1. linear predictor η_ij
  ## ----------  STEP (Z)  : sample latent utilities  z_ij  -----------------
  ## Linear predictor η_ij  (same as before)
  
  xbeta_vec <- X_flat %*% Beta # 2. Fixed-effects contribution:  Xβ → vector of length n·t_obs             
  eta_current <- matrix(xbeta_vec, nrow = n, ncol = t_obs) # 3. Reshape back to an n × t_obs matrix
  alpha_eff <- rowSums(Z * Alpha) # 4. Row-wise dot-product of Z and Alpha (one scalar per individual)
  eta_current <- sweep(eta_current, 1, alpha_eff, "+") # 5. Broadcast α_i across all observations of individual i
  
  ## λ_ij   and   π_ij  = F(η_ij)
  lambda_mat <- exp(eta_current)                     # λ_ij = exp(η_ij)
  pi_mat     <- lambda_mat / (1 + lambda_mat)        # π_ij = λ /(1+λ) = plogis(η)
  
  ## Draw U_ij  ∼  U(0,1)   (avoid the exact endpoints)
  epsU  <- .Machine$double.eps
  U_mat <- matrix(runif(n * t_obs, epsU, 1 - epsU), n, t_obs)
  
  ## Construct V_ij  as in the paper:
  ##   V_ij = y_ij + U_ij * ( 1 - y_ij - π_ij )
  V_mat <- Y + U_mat * (1 - Y - pi_mat)
  
  ## In words:
  ##   y = 0  →  V ∼ U(0 , 1-π)
  ##   y = 1  →  V ∼ U(1-π , 1)
  
  ## Logistic inverse-CDF  ε_ij  and latent utility  h_ij
  eps_mat <- log( V_mat / (1 - V_mat) )          # F^{-1}_ε (V)
  h_mat   <- eta_current + eps_mat               # h_ij = η_ij + ε_ij
  # -----------------------------------------------
  
  
  
  # 2. **Sample Polya-Gamma ω_{ij} | h, Beta, Alpha** for each observation
  #    ω_{ij} ~ PG(2, |h_{ij} - X_{ij}^T Beta - Z_{ij}^T Alpha_i|).
  #    Note: X_{ij}^T Beta + Z_{ij}^T Alpha_i = eta_current, so use |h - eta|.
  diff_mat <- abs(h_mat - eta_current)
  # Flatten diff_mat to vector and draw PG in one call for efficiency
  diff_vec <- as.vector(diff_mat)
  omega_vec <- rpg(length(diff_vec), 2, diff_vec)   #  h = 2 is correct here
  omega_mat <- matrix(omega_vec, nrow = n, ncol = t_obs)
  
  # 3. **Sample auxiliary location parameter gamma_tilde ~ N(0, G0)**
  gamma_tilde <- rnorm(1, mean = 0, sd = sqrt(G0))
  
  # 4. **Propose shifted utilities:** z_tilde (or h_tilde) = h + gamma_tilde
  h_tilde <- h_mat + gamma_tilde
  
  #    Compute truncation bounds for gamma_new:
  #    L = max_{i,j: Y_{ij}=0} h_tilde_{ij},  U = min_{i,j: Y_{ij}=1} h_tilde_{ij}.
  L_bound <- max(h_tilde[Y == 0])   # if no Y==0, this will be -Inf
  U_bound <- min(h_tilde[Y == 1])   # if no Y==1, this will be Inf
  if (is.infinite(L_bound)) L_bound <- -Inf
  if (is.infinite(U_bound)) U_bound <- Inf
  
  # 5. **Sample gamma_new | ω, h_tilde, Y** from N(g_N, G_N) truncated to [L, U].
  
  # mu_beta1_star = sum_ij omega_ij * h_tilde_ij + Sigma0_inv %*% mu0
  omega_vec <- as.vector(omega_mat)              # vectorizar omega
  mu_beta1_star <- sum(omega_vec * h_tilde) + chol_inverse(Sigma0) %*% mu0
  # mu_beta2_star = sum_ij omega_ij * X_ij^T  (X_all' * omega)
  mu_beta2_star <- t(X_all) %*% omega_vec        # X_all es (n*t_obs x p)
  # Sigma_beta_star = sum_ij omega_ij * X_ij X_ij^T + Sigma0_inv
  # Equivalent to: X_all' diag(omega) X_all + Sigma0_inv
  Sigma_beta_star <- t(X_all) %*% (X_all * omega_vec) + chol_inverse(Sigma0)
  
  # Linear Term
  A <- sum(omega_vec * h_tilde) +
    t(mu_beta2_star) %*% chol_inverse(Sigma_beta_star) %*% mu_beta1_star
  
  # Squared Term
  B <- sum(omega_vec) -
    t(mu_beta2_star) %*% chol_inverse(Sigma_beta_star) %*% mu_beta2_star +
    1 / G0
  
  matriz[iter,1] <- sum(omega_vec)
  matriz[iter,2] <-  t(mu_beta2_star) %*% chol_inverse(Sigma_beta_star) %*% mu_beta2_star
  matriz[iter,3] <- 1 / G0
  
  # 3. Media y varianza
  mu_gamma     <- as.numeric(A / B)
  var_gamma    <- 1 / as.numeric(B)
  
  gamma_new <- sample_truncnorm_inverse(
    a = L_bound,
    b = U_bound,
    mean = mu_gamma,
    sd = sqrt(var_gamma)
  )
  
  # 6. **Shift utilities:** define new h = h_tilde - gamma_new (to center latent utilities back around 0 threshold)
  h_mat <- h_tilde - gamma_new
  #    Update gamma to the newly sampled value (store for output)
  gamma <- gamma_new
  
  # 7. **Sample Beta | h, ω, Alpha** from its full conditional (Gaussian).
  ## ---------------------------------------------------------------
  ## Step to sample β     (clear version, non-vectorized)
  ## ---------------------------------------------------------------
  #
  #   We want:    Q_beta  =  Σ_{ij} ω_ij X_ij X_ijᵀ  +  Σ0^{-1}
  #               s       =  Σ_{ij} ω_ij X_ij (h_ij - Z_ijᵀ α_i)  +  Σ0^{-1} μ0
  #               m_beta  =  Q_beta^{-1} s
  #               β ~ N( m_beta ,  Q_beta^{-1} )
  #
  
  ## --- 0.  Helpers / dimensions ---------------------------------------------
  N_tot <- n * t_obs                     # total number of observations
  # X_flat  : (N_tot × p)   already built earlier
  # omega_vec : length N_tot
  # h_mat     : (n × t_obs)
  # Z, Alpha  : (n × q)
  
  ## --- 1.  Scalar random effect per individual  -----------------------------
  # zα_i  = Σ_k  Z_{ik} * Alpha_{ik}
  zalpha_eff <- rowSums(Z * Alpha)       # length-n vector
  
  ## --- 2.  Residuals r_ij  (= h_ij − zα_i)  -------------------------------
  r_mat <- sweep(h_mat, 1, zalpha_eff, "-")   # n × t_obs
  r_vec <- as.vector(r_mat)                   # flatten to length N_tot
  
  ## --- 3.  Accumulate  Σ ω x xᵀ   (Q_beta)  -------------------------------
  # Element-wise multiply each row of X_flat by its ω_ij
  weighted_X <- X_flat * omega_vec           # (N_tot × p)
  
  Q_beta <- Sigma0_inv +                      # prior precision
    t(X_flat) %*% weighted_X          # cross-product Σ ω x xᵀ   (p × p)
  
  ## --- 4.  Accumulate  Σ ω x r   (s_vec)  ----------------------------------
  s_vec <- as.numeric(
    Sigma0_inv %*% mu0 +             # prior term
      t(X_flat) %*% (omega_vec * r_vec)# Σ ω x r    (p × 1)
  )
  
  ## --- 5.  Posterior moments  ----------------------------------------------
  Q_inv  <- chol_inverse(Q_beta)                     # Σ_β  (p × p)
  m_beta <- Q_inv %*% s_vec                   # μ_β  (p × 1)
  
  ## --- 6.  One draw  β ∼ N(m_beta, Σ_β)  -----------------------------------
  Beta <- as.numeric(mvrnorm(1, mu = m_beta, Sigma = Q_inv))
  
  
  # 8. **Sample Alpha_i | Beta, h, ω for each i** from its Gaussian full conditional.
  #    For each individual i:
  for (i in 1:n) {
    # Extract vectors for this individual's observations
    omega_i <- omega_mat[i, ]            # length t_obs
    h_i     <- h_mat[i, ]               # length t_obs
    X_i     <- X[i,, ]                  # (t_obs x p) matrix
    # Compute precision: Q_alpha_i = V_alpha^{-1} + ∑_j ω_ij Z_ij Z_ij^T.
    # (If Z is constant per individual, Z_ij = Z[i,] for all j.)
    Z_i <- matrix(Z[i, ], nrow = q, ncol = 1)    # (q x 1) vector for this individual
    # Sum over j: since Z_i is constant, ∑_j ω_ij * Z_i Z_i^T = (∑_j ω_ij) * (Z_i %*% t(Z_i))
    sum_omega_i <- sum(omega_i)
    Q_alpha_i <- solve(V_alpha) + sum_omega_i * (Z_i %*% t(Z_i))
    # Compute mean: m_alpha_i = Q_alpha_i^{-1} [ ∑_j ω_ij * Z_ij * (h_ij - X_ij^T Beta) ]
    # (Note: h_ij - X_ij^T Beta is the part of latent utility not explained by fixed effects.)
    resid_i <- h_i - (X_i %*% Beta)     # vector length t_obs
    # Since Z_i is constant, ∑_j ω_ij * (h_ij - X_ij β) * Z_i = (Z_i) * ∑_j ω_ij * resid_i
    S_i <- Z_i * sum(omega_i * resid_i)   # (q x 1) vector
    m_alpha_i <- solve(Q_alpha_i, S_i)
    # Draw Alpha_i ~ MVN(m_alpha_i, Q_alpha_i^{-1})
    cov_alpha_i <- solve(Q_alpha_i)
    Alpha[i, ] <- as.numeric(mvrnorm(1, mu = m_alpha_i, Sigma = cov_alpha_i))
  }
  
  # 9. **Sample V_alpha | Alpha** from Inverse-Wishart(ν0 + n, Λ0 + ∑_{i=1}^n α_i α_i^T).
  S_post <- Lambda0 + t(Alpha) %*% Alpha   # q x q scale matrix (sum of outer products of Alpha_i plus Lambda0)
  V_alpha <- riwish(nu0 + n, S_post)
  
  # -- Save current draws
  # -- Save current draws --------------------------------------
  Beta_save[iter, ]    <- Beta
  Alpha_save[iter, , ] <- Alpha
  V_alpha_save[[iter]] <- V_alpha
  gamma_save[iter]     <- as.numeric(gamma)
  
  h_save[iter,     , ] <- h_mat       # ← NUEVO
  omega_save[iter, , ] <- omega_mat   # ← NUEVO


}
end_time <- Sys.time()
print(end_time - start_time)
timesamples_pupg <- end_time - start_time

# VIII. Discard burn-in samples
Beta_samples    <- Beta_save[(burnin+1):iterations, , drop=FALSE]
Alpha_samples   <- Alpha_save[(burnin+1):iterations, , , drop=FALSE]
V_alpha_samples <- V_alpha_save[(burnin+1):iterations]
gamma_samples   <- gamma_save[(burnin+1):iterations]
h_samples     <- h_save[(burnin + 1):iterations,     , , drop = FALSE]
omega_samples <- omega_save[(burnin + 1):iterations, , , drop = FALSE]

# The variables Beta_samples, Alpha_samples, V_alpha_samples, and gamma_samples 
# now contain the posterior samples after burn-in.
# For example, we can examine the posterior means:
colMeans(Beta_samples)        # posterior mean of beta
apply(Alpha_samples, 3, mean) # posterior mean of alpha (averaged over individuals, for each component)
Reduce("+", V_alpha_samples) / length(V_alpha_samples)  # posterior mean of V_alpha

# (Optional) Use coda to summarize or plot the MCMC:
beta_mcmc <- mcmc(Beta_samples)
summary(beta_mcmc)
# Example trace plot for beta:
plot(beta_mcmc) 





####################### Heat map

#───────────────────────────────────────────────────────────────
# 1.  Resúmenes rápidos que faltan
#───────────────────────────────────────────────────────────────

## Vα  – elementos únicos (q = 2 ⇒ 3 entradas)
v11_chain <- sapply(V_alpha_samples, `[[`, 1)          # V[1,1]
v12_chain <- sapply(V_alpha_samples, function(M) M[1,2])
v22_chain <- sapply(V_alpha_samples, function(M) M[2,2])

## α  – media y sd sobre los 100×2 componentes
alpha_mat  <- cbind(
  Alpha_samples[ , , 1],   # iter × 100  (α₁)
  Alpha_samples[ , , 2]    # iter × 100  (α₂)
)                         # iter × 200
alpha_mean <- rowMeans(alpha_mat)
alpha_sd   <- apply(alpha_mat, 1, sd)

## ω  y  h  – media y sd sobre i,j
omega_mean <- apply(omega_samples, 1, mean)
omega_sd   <- apply(omega_samples, 1,  sd)
h_mean     <- apply(h_samples,     1, mean)
h_sd       <- apply(h_samples,     1,  sd)

#───────────────────────────────────────────────────────────────
# 2.  Matriz compacta con todos los draws
#───────────────────────────────────────────────────────────────
post_draws <- cbind(
  Beta_samples,                 # β1, β2, β3  (ya es iter × 3)
  v11 = v11_chain,
  v12 = v12_chain,
  v22 = v22_chain,
  alpha_mean = alpha_mean,
  alpha_sd   = alpha_sd,
  omega_mean = omega_mean,
  omega_sd   = omega_sd,
  h_mean     = h_mean,
  h_sd       = h_sd
)
colnames(post_draws)[1:3] <- c("beta1", "beta2", "beta3")
#───────────────────────────────────────────────────────────────
# 3.  Correlaciones y heat-map
#───────────────────────────────────────────────────────────────
if (!requireNamespace("corrplot", quietly = TRUE))
  install.packages("corrplot")
library(corrplot)

R <- cor(post_draws, use = "pairwise.complete.obs")

corrplot(
  R,
  method      = "color",
  type        = "upper",
  addCoef.col = "black",
  number.cex  = 0.55,
  tl.col      = "black",
  tl.cex      = 0.75,
  tl.srt      = 45,
  mar         = c(0,0,2,0),
  title       = "Posterior correlations: β, Vα, α, ω, h"
)

################################ ESS

#───────────────────────────────────────────────────────────────
# 4.  ESS para cada columna de post_draws
#───────────────────────────────────────────────────────────────
if (!requireNamespace("coda", quietly = TRUE))
  install.packages("coda")
library(coda)

ess_vec <- sapply(as.data.frame(post_draws), effectiveSize)

print(round(ess_vec, 1))   # muestra el vector ordenado

## (opcional) guarda para comparar corridas
# saveRDS(ess_vec, file = "Output/ESS_vector_current_run.rds")

