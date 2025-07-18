               
# Gibbs Sampler loop
X_flat <- matrix(X, nrow = n * t_obs, ncol = p)

matriz <- matrix(data = 1:(3 * iterations), nrow = iterations, ncol = 3)

start_time <- Sys.time()
#iter <- 2
for (iter in 1:iterations) {
  
  print(iter*100/iterations)
  print("Delta")
  ## ---------------------------------------------------------------
  ## STEP 1 : draw latent utilities h_{ij}
  ## ---------------------------------------------------------------
  # 1. linear predictor η_ij
  ## ----------  STEP (Z)  : sample latent utilities  z_ij  -----------------
  ## Linear predictor η_ij  (same as before)
  eta_current <- matrix(0, n, t_obs)
  for (i in 1:n)
    eta_current[i, ] <- X[i, , ] %*% Beta + as.numeric(Z[i, ] %*% Alpha[i, ])
  
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
  V_alpha_inv <- solve(V_alpha)
  omega_vec <- as.vector(omega_mat)
  h_vec     <- as.vector(h_mat)
  
  gamma_new <- draw_gamma(
    omega_mat   = omega_mat,
    h_mat       = h_tilde,
    X           = X,
    Z           = Z,
    Sigma0_inv  = Sigma0_inv,
    mu0         = mu0,
    V_alpha_inv = V_alpha_inv,
    G0          = G0,
    L_bound     = L_bound,
    U_bound     = U_bound
  )
  
  # 6. **Shift utilities:** define new h = h_tilde - gamma_new (to center latent utilities back around 0 threshold)
  h_mat <- h_tilde - gamma_new
  #    Update gamma to the newly sampled value (store for output)
  gamma <- gamma_new
  
  ## ---------------------------------------------------------------
  ## STEP 7 : prior-δ  ,  posterior-δ  y re-escalado de  h
  ## ---------------------------------------------------------------
  
  ## 7-A ·  Draw δ_prior  ~  IG(d0 , D0)
  delta_prior <- 1 / rgamma(1, shape = d0, rate = D0)
  
  ## 7-B ·  Draw δ_post   ~  IG(a_post , b_post)   (función autónoma)
  delta_post <- draw_delta(
    h_mat       = h_mat,
    omega_mat   = omega_mat,
    X_arr       = X,
    Z_mat       = Z,
    mu0         = mu0,
    Sigma0_inv  = Sigma0_inv,
    tilde_delta = delta_prior,        # o tu valor
    V_alpha     = V_alpha,
    d0          = d0,
    D0          = D0
  )
  
  ## 7-C ·  Proponer nuevo h  →  sqrt(δ_prior / δ_post) · h
  scale_h <- sqrt(delta_prior / delta_post)
  #h_actual   <- scale_h * h_mat          # re-escala TODA la matriz de utilities
  
  ## 7-D ·  Guardar las dos δ si quieres monitorizarlas
  delta_save [iter] <- delta_post
  
  
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
  Sigma_beta <- chol_inverse(Q_beta)          # Q_beta⁻¹
  mu_beta    <- Sigma_beta %*% s_vec          # m_beta
  
  ## 6 · draw β  ~  N( mu_beta , Sigma_beta )
  Beta <- as.numeric(mvrnorm(1, mu = mu_beta, Sigma = Sigma_beta))
  
  
  
  
  ## ---------------------------------------------------------------
  ## STEP 9 : draw α_i | δ , β , h , ω         (usa h⋆ re-escalado)
  ## ---------------------------------------------------------------
  
           # n × T   (h⋆)
  
  for (i in 1:n) {
    ## --- datos del individuo i -----------------------------------
    omega_i <- omega_mat[i, ]                  # (T)
    h_star_i <- h_mat[i, ]               # (T)
    X_i <- X[i, , ]                            # T × p
    Z_i <- matrix(Z[i, ], nrow = q, ncol = 1) # q × 1  (constante en j)
    
    ## --- precisión  Σ_{α_i}^{-1} ---------------------------------
    sum_omega_i <- sum(omega_i)
    Sigma_alpha_inv_i <- solve(V_alpha) +
      sum_omega_i * (Z_i %*% t(Z_i))   # q × q
    
    ## --- media  m_{α_i} ------------------------------------------
    resid_star_i <- h_star_i*scale_h - (X_i %*% Beta)             # (T)
    S_i <- Z_i * sum(omega_i * resid_star_i)              # q × 1
    m_alpha_i <- solve(Sigma_alpha_inv_i, S_i)            # q × 1
    
    ## --- draw α_i  ~  N(m_{α_i}, Σ_{α_i}) ------------------------
    cov_alpha_i <- solve(Sigma_alpha_inv_i)               # Σ_{α_i}
    Alpha[i, ] <- as.numeric(
      mvrnorm(1, mu = m_alpha_i  ,
              Sigma = cov_alpha_i))
  }
  

  
  # 9. **Sample V_alpha | Alpha** from Inverse-Wishart(ν0 + n, Λ0 + ∑_{i=1}^n α_i α_i^T).
  S_post <- Lambda0 + t(Alpha) %*% Alpha   # q x q scale matrix (sum of outer products of Alpha_i plus Lambda0)
  V_alpha <- riwish(nu0 + n, S_post)
  
  # -- Save current draws
  Beta_save[iter, ]    <- Beta
  Alpha_save[iter, , ] <- Alpha
  V_alpha_save[[iter]] <- V_alpha
  gamma_save[iter]     <- as.numeric(gamma)
}
end_time <- Sys.time()
print(end_time - start_time)
timesamples_upgd <- end_time - start_time

# VIII. Discard burn-in samples
Beta_samples    <- Beta_save[(burnin+1):iterations, , drop=FALSE]
Alpha_samples   <- Alpha_save[(burnin+1):iterations, , , drop=FALSE]
V_alpha_samples <- V_alpha_save[(burnin+1):iterations]
gamma_samples   <- gamma_save[(burnin+1):iterations]

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


## ------------------------------------------------------------------
## 1.  Extraer columnas de Beta_samples  →  vectores individuales
## ------------------------------------------------------------------
beta_1_samples <- Beta_samples[ , 1]                  # (iter × 1)
beta_2_samples <- Beta_samples[ , 2]
beta_3_samples <- Beta_samples[ , 3]

## Etiquetas “aug” para ser coherentes
beta_1_samplesupgd <- beta_1_samples
beta_2_samplesupgd <- beta_2_samples
beta_3_samplesupgd <- beta_3_samples

## ------------------------------------------------------------------
## 2.  Convertir lista de V_alpha_samples en array 3-D y sacar elementos
##      (q = 2 ⇒ cada matriz es 2×2, así que simplificamos a 3-D)
## ------------------------------------------------------------------
V_arr <- simplify2array(V_alpha_samples)   # dims: 2 × 2 × draws

V1_upgd  <- V_arr[1, 1, ]   # trayectoria de V[1,1]
V2_upgd  <- V_arr[2, 2, ]   # trayectoria de V[2,2]
V21_upgd <- V_arr[2, 1, ]   # trayectoria de V[2,1]  (equivalente a V[1,2])

## ------------------------------------------------------------------
## 3.  Sacar ejemplos de alpha:
##     a1_ag = componente 1 de α para el INDIVIDUO 1  a lo largo de draws
##     a2_ag = componente 2 de α para el INDIVIDUO 1  a lo largo de draws
##     (Ajusta el índice “1” si te interesa otro individuo)
## ------------------------------------------------------------------
a1_upgd <- Alpha_samples[ , 1, ]   # (iter × 1)
a2_upgd <- Alpha_samples[ , 2, ]

## Mantén el array completo por si lo necesitas
alpha_samples_upgd <- Alpha_samples



## ------------------------------------------------------------------
## 5.  Guardar TODO en un único .RData, idéntico a tu otro script
## ------------------------------------------------------------------
if (!dir.exists("Output")) dir.create("Output")

unbalance <- abs(beta_true[1])

filename <- paste0(
  "Output/V_samplesupgd_",
  n, "_", t_obs, "_", iterations, "_", unbalance, ".RData"
)


save(a1_upgd, a2_upgd,
     V21_upgd, V2_upgd, V1_upgd,
     beta_1_samplesupgd, beta_2_samplesupgd, beta_3_samplesupgd,
     alpha_samples_upgd, timesamples_upgd,
     file = filename)


