############################################################
# UPG LOGIT (efectos fijos) + STAN (HMC) + COMPARACIÓN ESS
############################################################

## ==== 0. LIBRERÍAS ====
req_pkgs <- c("pgdraw","truncnorm","cmdstanr","posterior",
              "ggplot2","ggrepel","coda")
to_install <- req_pkgs[!sapply(req_pkgs, require, character.only = TRUE)]
if(length(to_install)) install.packages(to_install)
invisible(lapply(req_pkgs, require, character.only = TRUE))
library(posterior)

## ==== helper: chequear cmdstan ====
.have_cmdstan <- function() {
  ok <- TRUE
  tryCatch({ cmdstanr::cmdstan_path() }, error = function(e) ok <<- FALSE)
  ok
}

set.seed(1)
N <- 10000
P <- 3  # Puedes cambiar este valor libremente

X <- cbind(1, matrix(rnorm(N * (P - 1)), N, P - 1))
beta_true <- rnorm(P, mean = 0, sd = sqrt(10))  # Vector de dimensión P con media 0 y varianza 10
eta <- X %*% beta_true
p  <- 1 / (1 + exp(-eta))
y  <- rbinom(N, 1, p)
iterations <- 2000
constant <- 0 #if 1 we consider ESS/s if zero is just ESS. 

time_to_seconds <- function(x) {
  # Vectoriza
  x <- unlist(x, use.names = TRUE)
  
  convert_one <- function(z) {
    # NA o vacío
    if (is.null(z) || is.na(z) || !nzchar(z)) return(NA_real_)
    
    # Si ya es difftime
    if (inherits(z, "difftime")) return(as.numeric(z, units = "secs"))
    
    # Si es numérico puro: asumimos segundos
    if (is.numeric(z)) return(as.numeric(z))
    
    # Limpia espacios extra
    z <- trimws(z)
    
    # Formato hh:mm:ss o mm:ss
    if (grepl("^\\d{1,2}:\\d{2}(:\\d{2})?$", z)) {
      parts <- as.numeric(strsplit(z, ":")[[1]])
      if (length(parts) == 2) {           # mm:ss
        return(parts[1] * 60 + parts[2])
      } else if (length(parts) == 3) {    # hh:mm:ss
        return(parts[1] * 3600 + parts[2] * 60 + parts[3])
      }
    }
    
    # Extraer pares número-unidad (ej: "1.2 sec", "3 min", "4h")
    pattern <- "([0-9]*\\.?[0-9]+)\\s*([a-zA-Z]+)"
    m <- gregexpr(pattern, z, perl = TRUE)
    matches <- regmatches(z, m)[[1]]
    if (length(matches) == 0) return(NA_real_)
    
    # Tabla de conversión a segundos
    unit_map <- c(
      "ms" = 1e-3, "msec" = 1e-3, "millisecond" = 1e-3, "milliseconds" = 1e-3,
      "s" = 1, "sec" = 1, "second" = 1, "seconds" = 1,
      "m" = 60, "min" = 60, "minute" = 60, "minutes" = 60,
      "h" = 3600, "hr" = 3600, "hour" = 3600, "hours" = 3600,
      "d" = 86400, "day" = 86400, "days" = 86400
    )
    
    total <- 0
    for (piece in matches) {
      num  <- as.numeric(sub(pattern, "\\1", piece, perl = TRUE))
      unit <- tolower(sub(pattern, "\\2", piece, perl = TRUE))
      if (!unit %in% names(unit_map)) next
      total <- total + num * unit_map[unit]
    }
    total
  }
  
  out <- vapply(x, convert_one, numeric(1))
  names(out) <- names(x)
  out
}


## ==== 1. GIBBS UPG (SOLO EFECTOS FIJOS) ====
upg_logit_fe <- function(y, X,
                         nsave       = 10000,
                         nburn       = 2000,
                         d0          = 2.5,
                         D0          = 1.5,
                         G0          = 100,
                         B0          = 4,
                         A0          = 4,
                         gamma_boost = TRUE,
                         delta_boost = TRUE,
                         beta_start  = NULL,
                         verbose     = TRUE,
                         seed        = 2025){
  
  set.seed(seed)
  
  N  <- nrow(X); P <- ncol(X)
  ntot <- nburn + nsave
  
  A0_inv <- diag(1 / c(A0, rep(B0, P - 1)), P)
  
  beta_store  <- matrix(NA, nsave, P)
  gamma_store <- numeric(nsave)
  delta_store <- numeric(nsave)
  
  beta_draw <- if (is.null(beta_start)) rep(0, P) else beta_start
  
  Fe_ <- function(p) log(p) - log(1 - p)
  Fe  <- function(x) exp(x) / (1 + exp(x))
  
  if (verbose) pb <- txtProgressBar(min = 0, max = ntot, style = 3)
  
  Sy0 <- sum(y) == 0
  SyN <- sum(y) == N
  y0  <- y == 0
  y1  <- y == 1
  
  start_time <- Sys.time()
  for (irep in 1:ntot) {
    
    ## 1. Latent utilities
    U_i    <- runif(N)
    lambda <- as.vector(X %*% beta_draw)
    z      <- lambda + Fe_(y + U_i * (1 - y - (1 - Fe(-lambda))))
    
    ## 2. Latent scales (PG)
    omega  <- pgdraw::pgdraw(2, abs(z - lambda))
    
    ## 3. SHIFT (gamma)
    gamma_star <- if (gamma_boost) rnorm(1, 0, sqrt(G0)) else 0
    z_tilde    <- z + gamma_star
    
    XW    <- X * omega
    tXW   <- t(XW)
    Bn    <- chol2inv(chol(A0_inv + tXW %*% X))
    mb    <- colSums(XW)
    mg    <- sum(z_tilde * omega)
    mn    <- tXW %*% z_tilde
    tmbBn <- t(mb) %*% Bn
    Gn    <- 1 / ((1 / G0) + sum(omega) - tmbBn %*% mb)
    gn    <- Gn * (mg - tmbBn %*% mn)
    
    Ucut <- if (Sy0)  Inf else min(z_tilde[y1])
    Lcut <- if (SyN) -Inf else max(z_tilde[y0])
    
    gamma_draw <- if (gamma_boost)
      truncnorm::rtruncnorm(1, a = Lcut, b = Ucut, mean = gn, sd = sqrt(Gn)) else 0
    
    z_shift <- z_tilde - gamma_draw
    
    ## 4. SCALE (delta)
    delta_star <- if (delta_boost) 1 / rgamma(1, d0, D0) else 1
    
    mn2 <- tXW %*% z_shift
    bn  <- Bn %*% mn2
    
    quad <- as.numeric(t(bn) %*% A0_inv %*% bn +
                         sum(omega * (z_shift - as.vector(X %*% bn))^2))
    
    delta_draw <- if (delta_boost)
      1 / rgamma(1, d0 + 0.5 * N, D0 + 0.5 * delta_star * quad) else 1
    
    ## 5. Beta
    sqrtBn    <- t(chol(Bn))
    beta_draw <- as.vector(sqrt(delta_star / delta_draw) * bn + sqrtBn %*% rnorm(P))
    
    ## Store
    if (irep > nburn) {
      idx <- irep - nburn
      beta_store[idx, ]  <- beta_draw
      gamma_store[idx]   <- gamma_star - gamma_draw
      delta_store[idx]   <- delta_star / delta_draw
    }
    
    if (verbose) setTxtProgressBar(pb, irep)
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  timesamples_upg <- end_time - start_time
  
  list(y = y, X = X,
       beta  = beta_store,
       gamma = gamma_store,
       delta = delta_store,
       time = timesamples_upg)
}

## ==== 2. MODELO STAN (β, γ, δ; HMC) ====
stan_code_upg_fe <- "
data {
  int<lower=1> N;
  int<lower=1> P;
  array[N] int<lower=0,upper=1> y;
  matrix[N,P] X;

  real<lower=0> A0;
  real<lower=0> B0;
  real<lower=0> G0;
  real<lower=0> d0;
  real<lower=0> D0;
}
parameters {
  vector[P] beta_raw;
  real gamma_raw;
  real<lower=0> delta_raw;
}
transformed parameters {
  vector[P] beta;
  real gamma;
  real delta;

  gamma = gamma_raw * sqrt(G0);
  delta = delta_raw;

  beta[1] = beta_raw[1] * sqrt(A0);
  for (p in 2:P)
    beta[p] = beta_raw[p] * sqrt(B0);
}
model {
  beta_raw  ~ normal(0,1);
  gamma_raw ~ normal(0,1);
  delta_raw ~ inv_gamma(d0, D0);

  vector[N] eta = (X * beta + gamma) / sqrt(delta);
  y ~ bernoulli_logit(eta);
}
generated quantities {
  vector[N] lp = (X * beta + gamma) / sqrt(delta);
}
"

run_stan_upg_fe <- function(y, X,
                            iter   = 12000,
                            warmup = 2000,
                            chains = 4,
                            seed   = 2025,
                            A0 = 4, B0 = 4, G0 = 100, d0 = 2.5, D0 = 1.5){
  
  if (!.have_cmdstan()) {
    message("Instalando CmdStan...")
    cmdstanr::install_cmdstan()
  }
  
  stan_file <- cmdstanr::write_stan_file(stan_code_upg_fe)
  mod <- cmdstanr::cmdstan_model(stan_file)
  
  data_list <- list(N = nrow(X), P = ncol(X),
                    y = as.integer(y),
                    X = X,
                    A0 = A0, B0 = B0, G0 = G0, d0 = d0, D0 = D0)
  
  start_time <- Sys.time()
  fit <- mod$sample(data = data_list,
                    iter_warmup     = warmup,
                    iter_sampling   = iter - warmup,
                    chains          = chains,
                    parallel_chains = chains,
                    seed            = seed,
                    refresh         = 500)
  end_time <- Sys.time()
  print(end_time - start_time)
  timesamples_hmc <- end_time - start_time
  fit <- list(fit = fit, time = timesamples_hmc)
}



# ==== 4. EJEMPLO RÁPIDO ====

fit_r    <- upg_logit_fe(y, X, nsave = iterations - 0.1*iterations, nburn = 0.1*iterations)
fit_stan <- run_stan_upg_fe(y, X, iter = iterations, warmup = 0.1*iterations, chains = 1)

timesamples_upg <- time_to_seconds(fit_r$time)^constant
timesamples_hmc <- time_to_seconds(fit_stan$time)^constant
fit_stan_after <- fit_stan$fit 




## ---------- UPG ----------
upg_mat <- as.matrix(fit_r$beta)                        # n_iter x P
colnames(upg_mat) <- paste0("beta[", 1:ncol(upg_mat), "]")

ess_upg <- apply(upg_mat, 2, function(v) ess_bulk(v))/timesamples_upg   # vector ESS por parámetro


## ---------- Stan ----------
# Saca solo beta y pásalo a matriz normal
stan_draws <- fit_stan_after$draws(variables = "beta")        # draws_array
stan_mat   <- as.matrix(as_draws_matrix(stan_draws))    # n_iter*chains x P
beta_cols  <- grep("^beta\\[", colnames(stan_mat), value = TRUE)
stan_mat   <- stan_mat[, beta_cols, drop = FALSE]

ess_stan <- apply(stan_mat, 2, function(v) ess_bulk(v))/timesamples_hmc # vector ESS por parámetro


## ---------- Mostrar ----------
cat("ESS UPG:\n");  print(ess_upg)
timesamples_upg
cat("\nESS Stan:\n"); print(ess_stan)
timesamples_hmc


