#------------------------------------------------------------------
# Gibbs Sampler with Random Effects for Logistic Posterior Estimation
#------------------------------------------------------------------

rm(list = ls())  # Remove all variables
gc()             # Free memory


# I. Load required libraries
library(BayesLogit)  # For Polya-Gamma sampling (rpg)
library(MASS)        # For mvrnorm
library(ggplot2)
library(MCMCpack)    # For riwish()
library(coda)
library(truncnorm)   # For truncated normal sampling
library(coda)
library(ggplot2)
#LAST VERSION

# II. Set seed for reproducibility
set.seed(123)

# III. Set parameters
n <- 100       # Number of individuals
t_obs <- 30     # Observations per individual
p <- 3         # Number of predictors (fixed effects)
q <- 2         # Number of individual-level covariates (random effects)
reg <- 0       # Regularization parameter (if needed)

# IV. Simulate data
# i. Basic data:
X <- array(rnorm(n * t_obs * p), dim = c(n, t_obs, p))  # Predictor array (n x t_obs x p)
X[, , 1] <- 1   # Force intercept = 1

Z <- matrix(rnorm(n * q), nrow = n, ncol = q)          # Individual-level covariates (n x q)

# ii. Additional elements:
V_alpha_true <- matrix(c(1, 0.5,
                         0.5, 3), nrow = q, ncol = q)    # True covariance for random effects
alpha_true <- mvrnorm(n, mu = rep(0, q), Sigma = V_alpha_true)   # True random effects (n x q)
beta_true <- c(-7, 2, -1.5)                              # True fixed effects (length p)

# iii. Generate response Y deterministically using thresholding
eta <- matrix(0, nrow = n, ncol = t_obs)
for (i in 1:n) {
  # t(X[i, , ]) %*% beta_true yields a vector of length t_obs
  eta[i, ] <- t(X[i, , ] %*% beta_true) + as.vector(Z[i, ] %*% alpha_true[i, ])
}
# Use thresholding: Y = 1 if eta > 0, else 0.
prob <- 1 / (1 + exp(-eta))
Y <- matrix(rbinom(n * t_obs, size = 1, prob = prob), nrow = n, ncol = t_obs)

prop_ones <- mean(Y)               # proporción en [0,1]
pct_ones  <- 100 * prop_ones       # porcentaje

cat(sprintf("Proporción de 1's: %.3f  (%.2f%%)\n",
            prop_ones, pct_ones))

# # Simulate the logistic error for each observation using rlogis (location 0, scale 1)
# epsilon <- matrix(rlogis(n * t_obs, location = 0, scale = 1), nrow = n, ncol = t_obs)
# h <- eta + epsilon
# Y <- ifelse(h > 0, 1, 0)

# A helper function for matrix inversion using Cholesky (if needed)
inversechol <- function(mat, reg = 0) {
  mat_reg <- mat + reg * diag(nrow(mat))
  L <- chol(mat_reg)  # L such that t(L) %*% L = mat_reg
  chol2inv(L)
}

# V. Gibbs Sampler Settings & Initialization
n_iter <- 1500   # total iterations
burn   <-  n_iter*0.1 # burn-in iterations

# Storage for parameters:
beta_samples <- matrix(NA, n_iter, p)           # Fixed effects (p-dimensional)
alpha_samples <- array(NA, dim = c(n_iter, n, q))  # Random effects: each individual gets a q-vector
V_alpha_samples <- vector("list", n_iter)       # To store the q x q V_alpha
# Storage for latent variables
h_samples     <- array(NA_real_, dim = c(n_iter, n, t_obs))   # h draws
omega_samples <- array(NA_real_, dim = c(n_iter, n, t_obs))   # omega draws


# Initialize parameters:
beta_current <- rep(0, p)
alpha_current <- mvrnorm(n, mu = rep(0, q), Sigma = diag(1, q))
V_alpha_current <- diag(1, q)

# Initialize latent variables:
h <- matrix(0, nrow = n, ncol = t_obs)        # latent h (n x t_obs)
omega <- matrix(1, nrow = n, ncol = t_obs)      # Polya-Gamma latent vars (n x t_obs)

# Priors:
# For beta ~ N(mu0, Sigma0)
mu0 <- rep(0, p)
Sigma0 <- diag(1, p)
invSigma0 <- solve(Sigma0)

# For V_alpha ~ Inverse-Wishart(nu0, Lambda0)
nu0 <- 5         # degrees of freedom (at minimum q+1)
Lambda0 <- diag(2, q)

#Fix stuff outside the loop
X_mat      <- matrix(X, nrow = n * t_obs, ncol = p)
ZZ_precomp <- array(NA, dim = c(n, q, q))
for(i in 1:n) {
  ZZ_precomp[i,,] <- tcrossprod(Z[i,])
}
start_time <- Sys.time()
# VI. Gibbs Sampler Loop
for (iter in 1:n_iter) {
  print(iter)
  ## Step 1. Sample h given Y, beta, alpha, omega
  # --- Vectorized draw of all h[i,j] ---------------------
  
  ## Step 1. Sample h given Y, beta, alpha, omega via inverse‐CDF
  
  # 1) (X_mat built outside)
  # 2) fixed‐effects part
  xb_vec    <- X_mat %*% beta_current        # length N = n*t_obs
  
  # 3) random‐effects part
  alpha_eff <- rowSums(Z * alpha_current)    # length n
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
  
  # 7d) invert to get truncated‐Normal draws
  h_vec <- mu_vec + sd_vec * qnorm(u_vec)
  
  # 8) reshape back into an n×t_obs matrix
  h     <- matrix(h_vec, nrow = n, ncol = t_obs)
  
  
  ## Step 2. Sample omega from PG(1, h - (X*beta + Z*alpha))
  
  # 1) Flatten X into a (n*t_obs) × p matrix, so each row is X[i,j,]
  #done outside the loop.
  
  # 2) Fixed-effect part: X β
  xb_vec    <- X_mat %*% beta_current           # length = n*t_obs
  
  # 3) Random-effect part: αᵢ effect replicated for each j
  alpha_eff <- rowSums(Z * alpha_current)       # length n
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
  xb_vec   <- X_mat %*% beta_current           # length = n*t_obs
  xb_mat   <- matrix(xb_vec, nrow = n, ncol = t_obs)
  
  # 2) residual matrix
  resid_mat   <- h - xb_mat                    # n × t_obs
  
  # 3) row-sums for omega and for omega * residual
  sum_omega   <- rowSums(omega)                # length = n
  sum_resid   <- rowSums(omega * resid_mat)    # length = n
  
  # 4) build A_mat where row i = A_i = Z[i,] * sum_resid[i]
  A_mat       <- Z * sum_resid                 # n × q
  
  # 5) posterior precision constant
  V_alpha_inv <- solve(V_alpha_current)        # q × q
  for(i in 1:n) {
    # Posterior precision and covariance
    Bi      <- V_alpha_inv + sum_omega[i] * ZZ_precomp[i,,]
    cov_i   <- inversechol(Bi, reg = 1e-6)
    
    # Posterior mean
    mu_i <- cov_i %*% A_mat[i, ]
    
    # Cholesky-based draw (correct)
    R <- chol(cov_i)                    # R such that cov = R'R
    z <- rnorm(q)
    alpha_current[i, ] <- as.numeric(mu_i + t(R) %*% z)
  }
  
  
  
  ## Step 4. Sample V_alpha from its full conditional:
  # V_alpha | ... ~ IW(nu0 + n, Lambda0 + sum_{i=1}^n alpha_i alpha_i^T)
  # Sum of outer­products of each row of alpha_current:
  S_alpha <- crossprod(alpha_current)  
  # (this is t(alpha_current) %*% alpha_current)
  
  # Then draw V_alpha as before:
  V_alpha_current <- riwish(nu0 + n, Lambda0 + S_alpha)
  
  
  ## Step 5. Sample beta from its full conditional:
  # beta | ... ~ N(mu_beta, Sigma_beta)
  
  omega_vec<- as.vector(omega)                        # length n*t_obs
  h_vec    <- as.vector(h)                            # length n*t_obs
  
  # 1) Flattened alpha-effect and residual
  alpha_eff <- rowSums(Z * alpha_current)              # length n
  alpha_vec <- rep(alpha_eff, times = t_obs)          # length n*t_obs
  resid_vec <- h_vec - alpha_vec                      # length n*t_obs
  
  # 2) Precision update: Σβ⁻¹ = invSigma0 + ∑ ω_{ij} x_{ij} x_{ij}ᵀ
  Sigma_beta_inv <- invSigma0 +
    t(X_mat) %*% (X_mat * omega_vec)
  
  # 3) Sum term: ∑ ω_{ij} x_{ij} resid_{ij}
  sum_term <- as.numeric(
    t(X_mat) %*% (omega_vec * resid_vec)
  )
  
  # 4) Posterior covariance and mean
  Sigma_beta <- solve(Sigma_beta_inv)
  mu_beta    <- Sigma_beta %*% (invSigma0 %*% mu0 + sum_term)
  
  # 5) One MVN draw
  R <- chol(Sigma_beta)             # R is upper triangular such that Sigma_beta = R'R
  z <- rnorm(p)                     # Standard normal vector of length p
  beta_current <- as.numeric(mu_beta + t(R) %*% z)
  
  
  
  # Save samples:
  beta_samples[iter, ] <- beta_current
  alpha_samples[iter, , ] <- alpha_current
  V_alpha_samples[[iter]] <- V_alpha_current
  h_samples[iter,     , ] <- h          # n × t_obs
  omega_samples[iter, , ] <- omega      # n × t_obs
  
}

end_time <- Sys.time()
print(end_time - start_time)
timesamples_aug <- end_time - start_time

# Extrae las 3 entradas únicas de la matriz 2×2
v11 <- sapply(V_alpha_samples, function(mat) mat[1, 1])  # componente (1,1)
v12 <- sapply(V_alpha_samples, function(mat) mat[1, 2])  # componente (1,2) ≡ (2,1)
v22 <- sapply(V_alpha_samples, function(mat) mat[2, 2])  # componente (2,2)



# --- recorta burn-in ---
# --- recorta burn-in ---
post   <- (burn + 1):n_iter
beta_1_samplesaug <- beta_samples[post, 1]
beta_2_samplesaug <- beta_samples[post, 2]
beta_3_samplesaug <- beta_samples[post, 3]
V1_aug            <- v11[post]
V2_aug            <- v22[post]
V21_aug           <- v12[post]
a1_ag             <- alpha_samples[post, , 1]
a2_ag             <- alpha_samples[post, , 2]

# nuevos objetos
h_samples_aug     <- h_samples[post,     , ]    # iter × n × t_obs
omega_samples_aug <- omega_samples[post, , ]    # iter × n × t_obs




####################################  Correlaciones. 

# ============================================================
#    ==  HEAT-MAP DE CORRELACIONES  ==========================
# ============================================================

## 1) Construye las columnas post-burn -----------------------
Beta_samples <- cbind(
  beta1 = beta_1_samplesaug,
  beta2 = beta_2_samplesaug,
  beta3 = beta_3_samplesaug
)

# α  (media y sd sobre los 100 individuos × 2 componentes)
alpha_mat  <- cbind(a1_ag, a2_ag)      # iter × 200
alpha_mean <- rowMeans(alpha_mat)
alpha_sd   <- apply(alpha_mat, 1, sd)

# ω y h  (media y sd sobre i, j)
omega_mean <- apply(omega_samples_aug, 1, mean)
omega_sd   <- apply(omega_samples_aug, 1, sd)
h_mean     <- apply(h_samples_aug,     1, mean)
h_sd       <- apply(h_samples_aug,     1, sd)

# Vα elementos
v11_chain <- V1_aug
v12_chain <- V21_aug
v22_chain <- V2_aug

## 2) Matriz de draws compacta --------------------------------
post_draws <- cbind(
  Beta_samples,
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

## 3) Correlaciones y heat-map -------------------------------
if (!requireNamespace("corrplot", quietly = TRUE))
  install.packages("corrplot")
library(corrplot)

R <- cor(post_draws)

corrplot(R,
         method      = "color",
         type        = "upper",
         addCoef.col = "black",
         number.cex  = 0.55,
         tl.col      = "black",
         tl.cex      = 0.75,
         tl.srt      = 45,
         mar         = c(0,0,2,0),
         title       = "Posterior correlations: β, Vα, α, ω, h")


############################### ESS

#───────────────────────────────────────────────────────────────
#  ESS (Effective Sample Size) de cada cadena resumida
#───────────────────────────────────────────────────────────────
if (!requireNamespace("coda", quietly = TRUE))
  install.packages("coda")
library(coda)

ess_vec <- sapply(as.data.frame(post_draws), effectiveSize)

print(round(ess_vec, 1))          # muestra el ESS de cada parámetro

## (opcional) guarda para comparar corridas
# saveRDS(ess_vec, file = "Output/ESS_vector_unaugmented.rds")

