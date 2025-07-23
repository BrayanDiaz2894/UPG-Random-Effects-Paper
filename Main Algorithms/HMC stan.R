

############################################################
# Stan model: Logistic mixed-effects (q-dim random effects)
############################################################
stan_code_logit_re <- '
data {
  int<lower=1> N;                 // total obs
  int<lower=1> n_ind;             // individuals
  int<lower=1> p;                 // fixed effects dim
  int<lower=1> q;                 // RE dim
  array[N] int<lower=0,upper=1> y;
  matrix[N, p] X;
  matrix[n_ind, q] Z;
  array[N] int<lower=1, upper=n_ind> id;

  // IW hyperparameters
  real<lower=q-1> nu0;            // e.g. q + 3
  cov_matrix[q] Lambda0;          // e.g. diag_matrix(rep_vector(2,q))
}

parameters {
  vector[p] beta;
  matrix[n_ind, q] alpha;
  cov_matrix[q] V_alpha;
}

model {
  // Priors
  beta    ~ normal(0, 1);                     // matches Sigma0 = I_p
  V_alpha ~ inv_wishart(nu0, Lambda0);        // IW prior

  // Random effects
  for (i in 1:n_ind)
    alpha[i] ~ multi_normal(rep_vector(0, q), V_alpha);

  // Likelihood
  for (n_i in 1:N) {
    real eta_n = dot_product(row(X, n_i), beta)
               + dot_product(row(Z, id[n_i]), alpha[id[n_i]]);
    y[n_i] ~ bernoulli_logit(eta_n);
  }
}

'

stan_file_logit_re <- write_stan_file(stan_code_logit_re)

############################################################
# 3. COMPILE AND SAMPLE
############################################################

mod_logit_re <- cmdstan_model(stan_file_logit_re)

# IW hyperparameters
nu0     <- q + 3
Lambda0 <- diag(2, q)   # same as your Lambda0

start_time <- Sys.time()
fit <- mod_logit_re$sample(
  data = list(
    N       = N_total,
    n_ind   = n,
    p       = p,
    q       = q,
    y       = y_all,
    X       = X_all,
    Z       = Z,
    id      = id,
    nu0     = nu0,
    Lambda0 = Lambda0
  ),
  chains = 1,
  parallel_chains = 1,
  iter_warmup = iterations * 0.1,
  iter_sampling = (iterations)*0.9,
  seed = 123,
  refresh = 1
)

end_time <- Sys.time()
print(end_time - start_time)
timesamples_hmc <- time_to_seconds(end_time - start_time)
############################################################
# 3. EXTRACT DRAWS AS PURE MATRICES
############################################################
burnin <- floor(iterations * 0.1)

# 100% matrix base:
draws <- fit$draws(format = "matrix")
draws <- unclass(draws)                 # quita clase draws_matrix
storage.mode(draws) <- "double"         # asegura numeric

iter_total <- nrow(draws)               # usa lo que realmente salió
post       <- (burnin + 1):iter_total   # índices post-burnin

## 2a. Beta (iter × p)
Beta_save <- draws[, paste0("beta[", 1:p, "]"), drop = FALSE]

## 2b. Alpha  (iter × n × q)
Alpha_save <- array(NA_real_, dim = c(iter_total, n, q))
for (i in 1:n) {
  for (j in 1:q) {
    Alpha_save[, i, j] <- draws[, paste0("alpha[", i, ",", j, "]")]
  }
}

## 2c. V_alpha  (iter × q × q)
V_alpha_save <- array(NA_real_, dim = c(iter_total, q, q))
for (r in 1:q) {
  for (c in 1:q) {
    V_alpha_save[, r, c] <- draws[, paste0("V_alpha[", r, ",", c, "]")]
  }
}

############################################################
# 4. SLICE POST-BURN-IN AND NAME OBJECTS
############################################################
# Alpha components (q = 2)
a1_hmc <- as.matrix(Alpha_save[post, , 1])    # matrix: iter_post × n
a2_hmc <- as.matrix(Alpha_save[post, , 2])

# Beta
beta_1_sampleshmc <- as.matrix(Beta_save[post, 1])
beta_2_sampleshmc <- as.matrix(Beta_save[post, 2])
beta_3_sampleshmc <- as.matrix(Beta_save[post, 3])

# V_alpha components
V1_hmc  <- as.matrix(V_alpha_save[post, 1, 1])
V2_hmc  <- as.matrix(V_alpha_save[post, 2, 2])
V21_hmc <- as.matrix(V_alpha_save[post, 2, 1])

# Full Alpha block
alpha_samples_hmc <- Alpha_save[post, , ]

############################################################
# 5. SAVE
############################################################
unbalance <- abs(beta_true[1])
filename  <- paste0("Output/V_sampleshmc_",
                    n, "_", t_obs, "_", iter_total, "_", unbalance, ".RData")

save(a1_hmc, a2_hmc, V21_hmc, V2_hmc, V1_hmc,
     beta_1_sampleshmc, beta_2_sampleshmc, beta_3_sampleshmc,
     alpha_samples_hmc, timesamples_hmc,
     file = filename)


