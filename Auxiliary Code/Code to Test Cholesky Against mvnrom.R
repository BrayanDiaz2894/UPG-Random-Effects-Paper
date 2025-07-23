# # ---------------------------------------------------------------
# # 0.  Paquetes y funciones auxiliares
# # ---------------------------------------------------------------
# library(MASS)            # para mvrnorm()
# # (pega aquí tus definiciones de chol_inverse() y rmvnorm_chol())
# 
# # ---------------------------------------------------------------
# # 1.  Parámetros de la distribución
# # ---------------------------------------------------------------
# mu    <- c(1, -2, 0.5, 0, 1)     # media (d = 5)
# Sigma <- matrix(c( 2,  .3,  .4,  .1,  .2,
#                    .3,  1,  .2,  .0,  .1,
#                    .4,  .2, 1.5,  .2,  .3,
#                    .1,  .0,  .2,  1,   .1,
#                    .2,  .1,  .3,  .1, 1.2), 5, 5)
# 
# iterations <- 100000
# 
# # ---------------------------------------------------------------
# # 2.  Bucle con rmvnorm_chol()
# # ---------------------------------------------------------------
# time_chol <- system.time({
#   for (iter in seq_len(iterations)) {
#     ignore <- rmvnorm_chol(mu, Sigma, n_draw = 1)
#   }
# })
# 
# # ---------------------------------------------------------------
# # 3.  Bucle con mvrnorm()
# # ---------------------------------------------------------------
# time_mvr <- system.time({
#   for (iter in seq_len(iterations)) {
#     ignore <- mvrnorm(1, mu = mu, Sigma = Sigma)
#   }
# })
# 
# # ---------------------------------------------------------------
# # 4.  Resultado comparativo
# # ---------------------------------------------------------------
# print(time_chol)   # tiempo total con Cholesky rápido
# print(time_mvr)    # tiempo total con mvrnorm clásico

# ---------------------------------------------------------------
# 0.  Helper: fast Cholesky inverse
# ---------------------------------------------------------------
chol_inverse <- function(M) {
  R     <- chol(M)          # M = Rᵀ R
  chol2inv(R)               # (Rᵀ R)^{-1}
}

# ---------------------------------------------------------------
# 1.  Build a d × d SPD matrix once
# ---------------------------------------------------------------
set.seed(20250719)
d <- 100                                     # dimension (feel free to tweak)
A <- matrix(rnorm(d * d), d, d)
M <- crossprod(A) + diag(d) * 1e-6           # force SPD & well‑conditioned

iterations <- 10000                           # number of inversions

# ---------------------------------------------------------------
# 2.  Benchmark chol_inverse()
# ---------------------------------------------------------------
time_chol <- system.time({
  for (iter in seq_len(iterations)) {
    inv_M <- chol_inverse(M)
  }
})

# ---------------------------------------------------------------
# 3.  Benchmark solve()
# ---------------------------------------------------------------
time_solve <- system.time({
  for (iter in seq_len(iterations)) {
    inv_M <- solve(M)
  }
})

# ---------------------------------------------------------------
# 4.  Report
# ---------------------------------------------------------------
cat("Cholesky-based inverse:\n"); print(time_chol)
cat("solve() inverse:\n");        print(time_solve)

