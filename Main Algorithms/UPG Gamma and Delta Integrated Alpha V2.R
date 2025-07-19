               
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
        Qi  <- V_alpha_inv + sw[i] * tcrossprod(Zi)        # q × q
        
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
      Sigma_beta_star_inv <- solve(Sigma_beta_star)
      
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
      Sigma_beta_cov  <- solve(Sigma_beta_prec)               # covarianza posterior
      
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


