###############################################################################
## draw_gamma()  – versión corregida (ψ ≠ ψγ  y factores en numer / denom)
###############################################################################
draw_gamma <- function(omega_mat, h_mat, X, Z,
                       Sigma0_inv, mu0,
                       V_alpha_inv, G0,
                       L_bound, U_bound)
{
  
  
  #############################################################################
  ## 1 · sumatorios globales Σ ω   y   Σ ω h
  #############################################################################
  sum_omega   <- sum(omega_mat)
  sum_omega_h <- sum(omega_mat * h_mat)
  
  #############################################################################
  ## 2‑3 ·Bloque optimizado (usa X_all ya existente; mismo output que tu loop)
  #############################################################################
  
  ## ----------  sumas por individuo ------------------------------------------
  sw  <- rowSums(omega_mat)              # Σ_j ω_ij        (n)
  swh <- rowSums(omega_mat * h_mat)      # Σ_j ω_ij h_ij   (n)
  
  ## ----------  b‑vectores ----------------------------------------------------
  b1_mat <- sweep(Z, 1, swh, `*`)        # n × q
  b2_mat <- sweep(Z, 1, sw , `*`)        # n × q
  
  m_alpha01 <- matrix(as.numeric(t(b1_mat)), n * q, 1)
  m_alpha02 <- matrix(as.numeric(t(b2_mat)), n * q, 1)
  
  ## ----------  swX  usando X_all --------------------------------------------
  omega_vec <- c(t(omega_mat))                     # (n·t_obs) — ROW‑wise
  row_id    <- rep(seq_len(n), each = t_obs)
  
  swX_mat <- t( rowsum(X_all * omega_vec, row_id) )   # p × n
  
  ## ----------  listas --------------------------------------------------------
  
  eta_list        <- vector("list", n)
  psi_list        <- lapply(seq_len(n), \(i) matrix(b1_mat[i, ], ncol = 1))
  psi_gamma_list  <- lapply(seq_len(n), \(i) matrix(b2_mat[i, ], ncol = 1))
  
  Q_alpha_inv <- matrix(0, n * q, n * q)   # (n q) × (n q)
  P_alpha_inv <- Q_alpha_inv               # idéntico
  
  
  for (i in seq_len(n)) {
    Zi  <- matrix(Z[i, ], ncol = 1)                    # q × 1
    Qi  <- V_alpha_inv + sw[i] * tcrossprod(Zi)        # q × q
    
    Q_i_inv <- tryCatch(
      chol_inverse(Qi),
      error = function(e) chol_inverse(Qi + 1e-8 * diag(q))  # ridge si Qi no es PD
    )
    
    ## índice de filas/columnas para el bloque i
    idx <- ((i - 1) * q + 1):(i * q)
    
    ## escribir el bloque inverso directamente en la diagonal
    Q_alpha_inv[idx, idx] <- Q_i_inv
    eta_list[[i]]  <- Zi %*% t(swX_mat[, i, drop = FALSE])  # q × p
  }
  
  ## ----------  bloques diagonales finales ------------------------------------
  P_alpha_inv <- Q_alpha_inv             # idéntico
  
  
  
  
  #############################################################################
  ## 4 · η , ψ , ψγ  (apilados)
  #############################################################################
  eta       <- do.call(rbind, eta_list)        # (nq)×p
  psi       <- do.call(rbind, psi_list)        # (nq)×1
  psi_gamma <- do.call(rbind, psi_gamma_list)  # (nq)×1
  
  #############################################################################
  ## 5 · μ*  y  Σ*  de β   (incluyen η  y  ψ/ψγ correctos)
  #############################################################################
  
  eta_P_alpha <- t(eta) %*% P_alpha_inv
  
  Sigma_beta_star <- t(X_all) %*% (omega_vec * X_all) +
    Sigma0_inv -
    eta_P_alpha %*% eta
  Sigma_beta_star_inv <- solve(Sigma_beta_star)
  
  mu_beta1_star <- colSums(omega_vec * h_vec * X_all) +
    as.numeric(Sigma0_inv %*% mu0) -
    as.numeric(eta_P_alpha %*% psi)
  
  mu_beta2_star <- colSums(omega_vec * X_all) -
    as.numeric(t(psi_gamma) %*% P_alpha_inv %*% eta)
  
  #############################################################################
  ## 6 · escalares de numerador y denominador  (con signos correctos)
  #############################################################################
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
  
  #############################################################################
  ## 7 ·  sorteo de γ  truncado
  #############################################################################
  gamma_result <- rtruncnorm(1,
                             a = L_bound, b = U_bound,
                             mean = mu_gamma, sd = sqrt(sigma2_gamma))
  return(gamma_result)
  
}
