#-------------------------------------------------------------
# Gibbs Sampler with Random Effects for Logistic Posterior Estimation
#-------------------------------------------------------------

# Set the Preamble.
run_gibbs_sampler <- function(num_iter = iterations) {
#1. Packages.
library(BayesLogit)   # For Polya-Gamma sampling (rpg)
library(MASS)         # For multivariate normal (mvrnorm)
library(MCMCpack)     # For inverse-Wishart sampling (riwish)
library(truncnorm)    # For truncated normal sampling
library(coda)         # For MCMC diagnostics (if needed)
library(ggplot2)      # For plotting (if needed)
# library(doParallel)
# library(foreach)

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
rmvnorm_chol <- function(mu, Sigma, n_draw = 1) {
  d <- length(mu)
  R <- chol(Sigma)                 # Sigma = R' R (R upper)
  
  if (n_draw == 1L) {
    z <- rnorm(d)
    return(mu + drop(z %*% R))     # evita crear matrices grandes
  } else {
    Z <- matrix(rnorm(n_draw * d), n_draw, d)
    return(Z %*% R + matrix(mu, n_draw, d, byrow = TRUE))
  }
}


# 1. Simulation. 

# I. Set simulation parameters
n      <- 100 # Number of individuals
t_obs  <- 30    # Observations per individual
p      <- 3     # Number of fixed-effect predictors (including intercept)
q      <- 2     # Number of random-effect covariates (random effects dimension)
unbalance <- 6
iterations <- 1000


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
Sigma0_inv <- chol_inverse(Sigma0)

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


# Gibbs Sampler loop
X_flat <- matrix(X, nrow = n * t_obs, ncol = p)

matriz <- matrix(data = 1:(3 * iterations), nrow = iterations, ncol = 3)

start_time <- Sys.time()
#iter <- 2
for (iter in 1:iterations) {
  
  if (iter %% 100 == 0) cat(iter, "/", iterations, "\n")
  
  print("Delta")
  ## ---------------------------------------------------------------
  ## STEP 1 : draw latent utilities h_{ij}
  ## ---------------------------------------------------------------
  # 1) Parte fija: X_all %*% Beta  (N_total x 1)
  xb_vec <- drop(X_all %*% Beta)
  
  # 2) Reorganizar a (n x t_obs). OJO: como X_all se llenó por individuo,
  #    el primer bloque de 't_obs' filas es el individuo 1, etc.
  eta_X <- matrix(xb_vec, n, t_obs, byrow = TRUE)
  
  # 3) Escalar por individuo: Z[i,] %*% Alpha[i,]  (n)
  zalpha <- rowSums(Z * Alpha)
  
  # 4) Sumar ese escalar a cada fila (t_obs veces)
  eta_current <- sweep(eta_X, 1, zalpha, "+")   # o: eta_current <- eta_X + zalpha
  
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
  eps_mat <- log( V_mat / (1 - (V_mat)) )          # F^{-1}_ε (V)
  h_mat   <- eta_current + eps_mat               # h_ij = η_ij + ε_ij
  # ------------------------------------------------------------------------
  
  ## (Optional sanity check)
  all( h_mat[Y == 1]  > 0 - 1e-12 )   # should be TRUE
  all( h_mat[Y == 0] <= 0 + 1e-12 )   # should be TRUE
  
  
  
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
  # ---- Draw γ | ·  (truncated normal) ----
  V_alpha_inv <- chol_inverse(V_alpha)
  omega_vector <- as.vector(omega_mat)
  h_vec     <- as.vector(h_mat)
  
  ## I · sumatorios globales Σ ω   y   Σ ω h
  sum_omega   <- sum(omega_mat)
  sum_omega_h <- sum(omega_mat * h_mat)
  
  ## II and III ·Bloque optimizado (usa X_all ya existente; mismo output que tu loop)
  
  ## ----------  sumas por individuo ------------------------------------------
  sw  <- rowSums(omega_mat)              # Σ_j ω_ij        (n)
  swh <- rowSums(omega_mat * h_mat)      # Σ_j ω_ij h_ij   (n)
  
  ## ----------  b‑vectores ----------------------------------------------------
  b1_mat <- sweep(Z, 1, swh, `*`)        # n × q
  b2_mat <- sweep(Z, 1, sw , `*`)        # n × q
  
  m_alpha01 <- matrix(as.numeric(t(b1_mat)), n * q, 1)
  m_alpha02 <- matrix(as.numeric(t(b2_mat)), n * q, 1)
  
  ## ----------  swX  usando X_all --------------------------------------------
  omega_vector <- c(t(omega_mat))                     # (n·t_obs) — ROW‑wise
  row_id    <- rep(seq_len(n), each = t_obs)
  
  swX_mat <- t( rowsum(X_all * omega_vector, row_id) )   # p × n
  
  ## ----------  listas --------------------------------------------------------
  
  eta_list        <- vector("list", n)
  psi_list        <- lapply(seq_len(n), \(i) matrix(b1_mat[i, ], ncol = 1))
  psi_gamma_list  <- lapply(seq_len(n), \(i) matrix(b2_mat[i, ], ncol = 1))
  
  Q_alpha_inv <- matrix(0, n * q, n * q)   # (n q) × (n q)
  P_alpha_inv <- Q_alpha_inv               # idéntico
  
  
  for (i in seq_len(n)) {
    Zi  <- matrix(Z[i, ], ncol = 1)                    # q × 1
    Qi  <- V_alpha_inv + sw[i] * Z_outer[,, i]        # q × q
    
    Q_i_inv <- chol_inverse(Qi)
    
    ## índice de filas/columnas para el bloque i
    idx <- ((i - 1) * q + 1):(i * q)
    
    ## escribir el bloque inverso directamente en la diagonal
    Q_alpha_inv[idx, idx] <- Q_i_inv
    eta_list[[i]]  <- Zi %*% t(swX_mat[, i, drop = FALSE])  # q × p
  }
  
  ## ----------  bloques diagonales finales ------------------------------------
  P_alpha_inv <- Q_alpha_inv             # idéntico
  
  ## IV · η , ψ , ψγ  (apilados)
  
  eta       <- do.call(rbind, eta_list)        # (nq)×p
  psi       <- do.call(rbind, psi_list)        # (nq)×1
  psi_gamma <- do.call(rbind, psi_gamma_list)  # (nq)×1
  
  
  ## V · μ*  y  Σ*  de β   (incluyen η  y  ψ/ψγ correctos)
  
  eta_P_alpha <- t(eta) %*% P_alpha_inv
  
  Sigma_beta_star <- t(X_all) %*% (omega_vector * X_all) +
    Sigma0_inv -
    eta_P_alpha %*% eta
  Sigma_beta_star_inv <- chol_inverse(Sigma_beta_star)
  
  mu_beta1_star <- colSums(omega_vector * h_vec * X_all) +
    as.numeric(Sigma0_inv %*% mu0) -
    as.numeric(eta_P_alpha %*% psi)
  
  mu_beta2_star <- colSums(omega_vector * X_all) -
    as.numeric(t(psi_gamma) %*% P_alpha_inv %*% eta)
  
  ## VI · escalares de numerador y denominador  (con signos correctos)
  
  m_alpha_02_Q <- t(m_alpha02) %*% Q_alpha_inv
  term_m2_Qinv_m1 <- as.numeric(m_alpha_02_Q %*% m_alpha01)
  term_m2_Qinv_m2 <- as.numeric(m_alpha_02_Q %*% m_alpha02)
  
  
  mu_beta2_star_Sigma_beta <- t(mu_beta2_star) %*% Sigma_beta_star_inv 
  numer  <- sum_omega_h -
    term_m2_Qinv_m1 -
    (mu_beta2_star_Sigma_beta %*% mu_beta1_star)
  
  denom  <- sum_omega -
    term_m2_Qinv_m2 -
    (mu_beta2_star_Sigma_beta %*% mu_beta2_star) +
    1 / G0
  
  sigma2_gamma <- as.numeric(1 / denom)
  mu_gamma     <- as.numeric(numer * sigma2_gamma)
  
  ## VII ·  sorteo de γ  truncado
  
  gamma_new <- rtruncnorm(1,
                          a = L_bound, b = U_bound,
                          mean = mu_gamma, sd = sqrt(sigma2_gamma))
  
  # 6. **Shift utilities:** define new h = h_tilde - gamma_new (to center latent utilities back around 0 threshold)
  h_mat <- h_tilde - gamma_new
  #    Update gamma to the newly sampled value (store for output)
  gamma <- gamma_new
  
  # 7 : prior-δ  ,  posterior-δ  y re-escalado de  h
  ## ---------------------------------------------------------------
  
  # I ·  Draw δ_prior  ~  IG(d0 , D0)
  delta_prior <- 1 / rgamma(1, shape = d0, rate = D0)
  
  ## ----- DRAW δ | todo ----- ##
  
  ## Conveniencias locales -----------------------------------------
  # h_mat ya fue "recentred" (h_tilde - gamma_new).
  h_vec <- c(t(h_mat))               # N_total × 1 (row-wise flatten)
  
  ## (i). Σ ω h^2 -----------------------------------------------------
  sum_h2 <- sum(omega_vector * (h_vec^2))
  
  ## (ii). XtΩX y XtΩh SIN loop (reuso de estructuras gamma) ----------
  XtOmegaX <- t(X_all) %*% (omega_vector * X_all)        # p × p
  XtOmegah <- t(X_all) %*% (omega_vector * h_vec)        # p × 1
  
  ## (iii). Bloques colapsados por individuo (reuso total) -----------
  nu          <- matrix(as.numeric(t(b1_mat)), n * q, 1)  # (nq) × 1
  P_alpha_inv <- Q_alpha_inv                              # alias
  
  ## (iv). Componentes para β colapsada -------------------------------
  # eta_P_eta = t(eta) %*% P_alpha_inv %*% eta
  eta_P_eta       <- t(eta) %*% (P_alpha_inv %*% eta)
  Sigma_beta_prec <- Sigma0_inv + XtOmegaX - eta_P_eta    # precisión (p×p)
  Sigma_beta_cov  <- chol_inverse(Sigma_beta_prec)               # covarianza posterior
  
  # b_beta = XtΩh - t(eta) %*% P_alpha_inv %*% nu + Sigma0_inv %*% mu0
  b_beta       <- XtOmegah - t(eta) %*% (P_alpha_inv %*% nu) + Sigma0_inv %*% mu0
  mu_beta_star <- Sigma_beta_cov %*% b_beta                # media posterior colapsada
  
  ## (v). Términos cuadráticos ---------------------------------------
  # quad_beta   = b_beta' Σ_beta b_beta
  quad_beta <- as.numeric(t(b_beta) %*% (Sigma_beta_cov %*% b_beta))
  
  # quad_alpha  = ν' P_alpha_inv ν
  quad_alpha <- as.numeric(t(nu) %*% (P_alpha_inv %*% nu))
  
  ## (vi). Parámetros IG para δ --------------------------------------
  a_post <- d0 + N_total / 2
  
  b_post <- D0 + (delta_prior / 2) * (sum_h2 - quad_beta - quad_alpha)
  #                                  ^ signos:  suma global menos correcciones colapsadas
  if (b_post <= 0 || !is.finite(b_post))
    stop("draw_delta(): b_post ≤ 0 o no finito – revisa insumos (sum_h2, quad_beta, quad_alpha).")
  
  # δ | ·  ~ IG(a_post, b_post)   (parametrización: dens ∝ x^{-(a_post+1)} exp(-b_post/x))
  # Usando rgamma: si X ~ IG(a,b) ⇒ 1/X ~ Gamma(a, rate=b).
  delta_post <- 1 / rgamma(1, shape = a_post, rate = b_post)
  
  
  
  ## III ·  Proponer nuevo h  →  sqrt(δ_prior / δ_post) · h
  scale_h <- sqrt(delta_prior / delta_post)
  #h_actual   <- scale_h * h_mat          # re-escala TODA la matriz de utilities
  
  ## IV·  Guardar las dos δ si quieres monitorizarlas
  #delta_save [iter] <- delta_post
  
  
  # 8. **Sample Beta | h, ω, Alpha** from its full conditional (Gaussian).
  ## ---------------------------------------------------------------
  ## STEP 8 : draw  β | δ , h , ω , α      (versión h⋆  re-escalada)
  ## ---------------------------------------------------------------
  
  ## 1 · zα_i  (suma de random effects por individuo)
  zalpha_eff <- rowSums(Z * Alpha)             # longitud-n
  
  ## 2 · residuos  r⋆_ij = h⋆_ij − zα_i
  r_star_mat <- sweep(h_mat*scale_h, 1, zalpha_eff, "-")   # n × T
  r_star_vec <- as.vector(r_star_mat)                   # (N_tot)
  
  ## 3 · precisión  Q_beta  = Σ ω X Xᵀ  +  Σ0⁻¹
  weighted_X <- X_flat * omega_vec            # (N_tot × p)
  Q_beta <- Sigma0_inv + t(X_flat) %*% weighted_X   # (p × p)
  
  ## 4 · vector  s = Σ ω X r⋆  +  Σ0⁻¹ μ0
  s_vec <- as.numeric(
    Sigma0_inv %*% mu0 +
      t(X_flat) %*% (omega_vec * r_star_vec)
  )
  
  ## 5 · momentos posteriores
  Rbeta   <- chol(Q_beta)                          
  ybeta   <- forwardsolve(t(Rbeta), s_vec, upper.tri = FALSE)       
  m_ibeta <- backsolve(Rbeta,  ybeta, upper.tri = TRUE)      
  
  zbeta <- rnorm(p)                                 # z ~ N(0, I_q)
  Beta <- drop(m_ibeta + backsolve(Rbeta, zbeta, upper.tri = TRUE))
  
  

  ## ---------------------------------------------------------------
  ## STEP 9 : draw α_i | δ , β , h , ω         (usa h⋆ re-escalado)
  ## ---------------------------------------------------------------
  
  for (i in 1:n) {
    ## --- datos del individuo i -----------------------------------
    omega_i <- omega_mat[i, ]                  # (T)
    h_star_i <- h_mat[i, ]               # (T)
    X_i <- X[i, , ]                            # T × p
    Z_i <- matrix(Z[i, ], nrow = q, ncol = 1) # q × 1  (constante en j)
    
    ## --- precisión  Σ_{α_i}^{-1} ---------------------------------
    sum_omega_i <- sum(omega_i)
    Sigma_alpha_inv_i <- V_alpha_inv + sum_omega_i * Z_outer[,, i]
    
    
    ## --- media  m_{α_i} ------------------------------------------
    resid_star_i <- h_star_i*scale_h - (X_i %*% Beta)             # (T)
    S_i <- Z_i * sum(omega_i * resid_star_i)              # q × 1
    ## --- media  m_{α_i} ------------------------------------------
    
    #Let m_i to be the mean and Q_i to be the variance. Let S_i to be the linear term from the kernel.
    #Then the mean m_i = Q_i^-1 * S_i   or Q_i*m_i=S_i 
    #We dont want to invert Q_i. 
    #We decompose Q_i = R'R, then Q_i*m_i=S_i which implies that R'R*m_i=S_i
    #First we solve R' y = S_i (for y) so we have: y=R'^-1*S_i
    #Then we solve R m_i = y so we have  m_i = R^-1*y=(RR)'^-1*S_i
    #Backsolve is used for R triangular inferior, forward solve for R triangular sueperior#Backsolve is used for R triangular inferior, forward solve for R triangular sueperior
    
  
    R   <- chol(Sigma_alpha_inv_i)                          # Q_i = R'R
    y   <- forwardsolve(t(R), S_i, upper.tri = FALSE)       # resuelve R' y = S_i
    m_i <- backsolve(R,  y, upper.tri = TRUE)               # luego R m_i = y
    
    
    z <- rnorm(q)                                 # z ~ N(0, I_q)
    Alpha[i, ] <- drop(m_i + backsolve(R, z, upper.tri = TRUE))
    
  }
  
  
  
  # 10. **Sample V_alpha | Alpha** from Inverse-Wishart(ν0 + n, Λ0 + ∑_{i=1}^n α_i α_i^T).
  S_post <- Lambda0 + t(Alpha) %*% Alpha   # q x q scale matrix (sum of outer products of Alpha_i plus Lambda0)
  V_alpha <- riwish(nu0 + n, S_post)
  
  # -- Save current draws
  Beta_save[iter, ]    <- Beta
  Alpha_save[iter, , ] <- Alpha
  V_alpha_save[[iter]] <- V_alpha
  gamma_save[iter]     <- as.numeric(gamma)
}
end_time <- Sys.time()
timediff <- end_time - start_time
print(timediff)
}

