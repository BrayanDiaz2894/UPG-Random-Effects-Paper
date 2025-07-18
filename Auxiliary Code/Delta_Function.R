###############################################################################
## draw_delta() – sortea δ | h , ω   (posterior Inverse-Gamma)
###############################################################################
# Entradas
#   h_mat        : n × T  matriz con h_{ij}
#   omega_mat    : n × T  matriz con ω_{ij}
#   X_arr        : n × T × p array de predictores fijos
#   Z_mat        : n × q matriz (covariables de efectos aleatorios)
#   mu0          : p-vector (media a priori de β)
#   Sigma0_inv   : p × p matriz  (Σ₀⁻¹)
#   tilde_delta  : escalar  δ̃  del modelo
#   V_alpha      : q × q matriz  (covarianza base para αᵢ)
#   d0 , D0      : hiperparámetros de la prior IG(d0 , D0)
#
# Devuelve
#   escalar  δ_new  ~ IG(a_post , b_post)
###############################################################################
draw_delta <- function(h_mat, omega_mat, X_arr, Z_mat,
                       mu0, Sigma0_inv,
                       tilde_delta, V_alpha,
                       d0, D0)
{
  require(Matrix)                       # para bdiag()
  
  ## --- 0 · Dimensiones ------------------------------------------------------
  n     <- nrow(h_mat)
  Tobs  <- ncol(h_mat)
  p     <- dim(X_arr)[3]
  q     <- ncol(Z_mat)
  Ntot  <- n * Tobs
  
  ## --- 1 · Suma  Σ ω h² -----------------------------------------------------
  sum_h2 <- sum(omega_mat * (h_mat ^ 2))
  
  ## --- 2 · Inicializaciones -------------------------------------------------
  V_alpha_inv <- solve(V_alpha)
  
  XtOmegaX <- matrix(0, p, p)
  XtOmegah <- matrix(0, p, 1)
  
  eta_list <- vector("list", n)
  nu_list  <- vector("list", n)
  P_blocks <- vector("list", n)
  
  quad_alpha <- 0                       # ∑ bᵢᵀ Qᵢ⁻¹ bᵢ
  
  ## --- 3 · Loop por individuo ----------------------------------------------
  for (i in seq_len(n)) {
    w_i <- omega_mat[i, ]                          # longitud Tobs
    Z_i <- matrix(Z_mat[i, ], ncol = 1)            # q × 1
    X_i <- matrix(X_arr[i, , ], nrow = Tobs, ncol = p)
    h_i <- h_mat[i, ]                              # longitud Tobs
    
    ## Sumas de verosimilitud
    XtOmegaX <- XtOmegaX + t(X_i) %*% (w_i * X_i)
    XtOmegah <- XtOmegah + t(X_i) %*% (w_i * h_i)
    
    ## Bloques para colapso
    eta_i <- Z_i %*% t(colSums(w_i * X_i))         # q × p
    nu_i  <- Z_i * sum(w_i * h_i)                  # q × 1
    
    ## Q_i  y su inversa (covarianza)
    Q_i    <- V_alpha_inv + sum(w_i) * (Z_i %*% t(Z_i))
    Q_i_inv <- chol2inv(chol(Q_i))
    
    ## bᵢ  y término cuadrático
    b_i <- nu_i
    quad_alpha <- quad_alpha + as.numeric(t(b_i) %*% Q_i_inv %*% b_i)
    
    ## Acumular bloques
    eta_list[[i]] <- eta_i
    nu_list [[i]] <- nu_i
    P_blocks[[i]] <- Q_i_inv
  }
  
  ## --- 4 · Construcciones globales -----------------------------------------
  eta         <- do.call(rbind, eta_list)                  # (nq) × p
  nu          <- do.call(rbind, nu_list)                   # (nq) × 1
  P_alpha_inv <- Matrix::bdiag(P_blocks)                   # (nq) × (nq)
  
  ## --- 5 · Posterior colapsado de β -------------
  eta_P_eta       <- t(eta) %*% (P_alpha_inv %*% eta)     # p × p
  Sigma_beta_prec <- Sigma0_inv + XtOmegaX - eta_P_eta    # precisión
  Sigma_beta_cov  <- solve(Sigma_beta_prec)               # cov
  
  b_beta          <- XtOmegah - t(eta) %*% (P_alpha_inv %*% nu) + Sigma0_inv %*% mu0
  mu_beta_star    <- Sigma_beta_cov %*% b_beta            # p × 1
  
  # Usar la forma bᵀ Σ b  (≡ mᵀ Σ⁻¹ m)
  quad_beta <- as.numeric(t(b_beta) %*% (Sigma_beta_cov %*% b_beta))
  #              ^^^^^^^                          ^^^^^
  
  ## --- 6 · Parámetros IG  ------------------------
  a_post <- d0 + Ntot / 2
  b_post <- D0 + (tilde_delta / 2) * (sum_h2 - quad_beta - quad_alpha)
  #                                         ^  signos cambiados  ^
  if (b_post <= 0 || !is.finite(b_post))
    stop("draw_delta(): b_post ≤ 0 – revisa insumos")
  
  delta_new <- 1 / rgamma(1, shape = a_post, rate = b_post) 
  return(delta_new)
}
