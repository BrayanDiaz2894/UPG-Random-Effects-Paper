#-------------------------------------------------------------
# Gibbs Sampler with Random Effects for Logistic Posterior Estimation
#-------------------------------------------------------------

# Set the Preamble.

#1. Packages.
library(BayesLogit)   # For Polya-Gamma sampling (rpg)
library(MASS)         # For multivariate normal (mvrnorm)
library(MCMCpack)     # For inverse-Wishart sampling (riwish)
library(truncnorm)    # For truncated normal sampling
library(coda)         # For MCMC diagnostics (if needed)
library(ggplot2)      # For plotting (if needed)
library(profvis)

#2. Auxiliary Functions. 
source("Auxiliary Code/Gamma_FunctionV2.R") #Always use gamma V2. It is actualy used in the sampler with only gamma augmentation.
#source("Auxiliary Code/Delta_Function.R") #This one is outdated. I am not using delta anymore as a separate function. 



#II. Truncated Sampling.
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

#3. Invert a symmetric positive-definite matrix via Cholesky
chol_inverse <- function(M) {
  # Cholesky factor R:  M = Rᵀ R   (R is upper-triangular)
  R <- chol(M)
  
  # Fast inversion via chol2inv:  (Rᵀ R)^{-1} = R^{-1} R^{-ᵀ}
  M_inv <- chol2inv(R)
  
  return(M_inv)
}

#4. Fast multivariate-normal sampler via Cholesky
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

# 1. Simulation. 

# I. Set simulation parameters
n      <- 10000 # Number of individuals
t_obs  <- 10    # Observations per individual
p      <- 3     # Number of fixed-effect predictors (including intercept)
q      <- 2     # Number of random-effect covariates (random effects dimension)
unbalance <- 6
iterations <- 5000


# II. Set seed for reproducibility
set.seed(123)


# III. Simulate data

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
beta_true <- c(-unbalance, 2, -1.5)  # length p = 3 (intercept, X2, X3)

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
mean(Y)

#IV. Set up priors and initial values for the Gibbs sampler
# 1. Priors for fixed effects beta ~ N(mu0, Sigma0)
mu0     <- rep(0, p)
Sigma0  <- diag(p)          # relatively diffuse prior (variance 10 on diagonal)
Sigma0_inv <- solve(Sigma0)

# 2. Prior for and DELTA (global location parameter) ~ N(0, G0)
G0    <- 10                      # variance of gamma's prior (could be tuned as needed)
d0 <- 2            # prior IG(d0,D0)  – ajusta a gusto
D0 <- 1

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
burnin     <- floor(iterations * 0.1)    # 10% burn-in
Beta_save    <- matrix(NA, nrow = iterations, ncol = p)    # store beta draws
Alpha_save   <- array(NA, dim = c(iterations, n, q))       # store alpha draws
V_alpha_save <- vector("list", length = iterations)        # store V_alpha draws (as matrices)
gamma_save   <- numeric(iterations)                       # store gamma draws

#2. Run the gibbs sampler for the different algorithms.

#Save the elements in the environment in this point.
object_names <- ls()
object_names <- ls()

  # #iv. profile.
  # source("Additional Versions/UPG Gamma and Delta Integrated Alpha V2 Profiler.R")   # carga la función con una ruta concreta
  # profvis(run_gibbs_sampler())
  # rm(run_gibbs_sampler)


#I. Basic PG

#i. Fix stuff outside the loop.
ZZ_precomp <- array(NA, dim = c(n, q, q))
for (i in 1:n) {
  ZZ_precomp[i,,] <- tcrossprod(Z[i, ])  # outer product: Z_i Z_i^T
}
X_flat <- matrix(X, nrow = n * t_obs, ncol = p) #Used to vectorize beta.
X_flatt_cross <- t(X_flat)%*%X_flat
row_id    <- rep(1:n, times = t_obs)
V_alpha_inv <- chol_inverse(V_alpha)
omega_store <- array(1, dim = c(iterations, n, t_obs))  # Store omega samples

#ii. write store values.
Beta_save    <- matrix(NA, nrow = iterations, ncol = p)    # store beta draws
Alpha_save   <- array(NA, dim = c(iterations, n, q))       # store alpha draws
V_alpha_save <- array(1, dim = c(iterations, q, q))

#iii. run the loop.
source("Main Algorithms/Basic PG.R")

#iv. clean
object_names1 <- ls() #Store all the objects after we ran the first loop.
dropset1 <- setdiff(object_names1,object_names) #Objects that are now but were not creating in the simulation
rm(dropset1)



#II. Hybrid PG and Latent Variable

#i. Fix stuff outside the loop.
X_mat      <- matrix(X, nrow = n * t_obs, ncol = p)
ZZ_precomp <- array(NA, dim = c(n, q, q))
for(i in 1:n) {
  ZZ_precomp[i,,] <- tcrossprod(Z[i,])
}
omega <- matrix(1, nrow = n, ncol = t_obs)      # Polya-Gamma latent vars (n x t_obs)

## ii.fix initial values.
Beta_save    <- matrix(NA, nrow = iterations, ncol = p)    # store beta draws
Alpha_save   <- array(NA, dim = c(iterations, n, q))       # store alpha draws
V_alpha_save <- vector("list", length = iterations)        # store V_alpha draws (as matrices)

#iii. run the loop.
source("Main Algorithms/Hybrid PG and Latent Variable.R")

#iv. clean
object_names2 <- ls() #Store all the objects after we ran the first loop.
dropset2 <- setdiff(object_names2,object_names) #Objects that are now but were not creating in the simulation
rm(dropset2)

III. UPG Gamma and Delta Integrated Alpha

#i. Fix stuff outside the loop.

## ii.fix initial values.
Beta_save    <- matrix(NA, nrow = iterations, ncol = p)    # store beta draws
Alpha_save   <- array(NA, dim = c(iterations, n, q))       # store alpha draws
V_alpha_save <- vector("list", length = iterations)        # store V_alpha draws (as matrices)
gamma_save   <- numeric(iterations)                       # store gamma draws
delta_save <- numeric(iterations)   # reserva espacio

#iii. run the loop.
source("Main Algorithms/UPG Gamma and Delta Integrated Alpha V2.R")

#iv. clean
object_names3 <- ls() #Store all the objects after we ran the first loop.
dropset3 <- setdiff(object_names3,object_names) #Objects that are now but were not creating in the simulation
rm(dropset3)

#IV. UPG Gamma Integrated Alpha

#i. Fix stuff outside the loop.

## ii.fix initial values.
Beta_save    <- matrix(NA, nrow = iterations, ncol = p)    # store beta draws
Alpha_save   <- array(NA, dim = c(iterations, n, q))       # store alpha draws
V_alpha_save <- vector("list", length = iterations)        # store V_alpha draws (as matrices)
gamma_save   <- numeric(iterations)                       # store gamma draws

#iii. run the loop.
source("Main Algorithms/UPG Gamma Integrated Alpha V2.R")

#iv. clean
object_names4 <- ls() #Store all the objects after we ran the first loop.
dropset4 <- setdiff(object_names4,object_names) #Objects that are now but were not creating in the simulation
rm(dropset4)



