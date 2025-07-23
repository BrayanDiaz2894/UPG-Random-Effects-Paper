#------------------------------------------------------------------
# Gibbs Sampler with Random Effects for Logistic Posterior Estimation
# Only alpha is sampled; beta and V_alpha are fixed to the true population values.
#Version with two augmentations latent and pg.
#------------------------------------------------------------------

rm(list = ls())  # Remove all variables
gc()             # Free memory

start_time <- Sys.time()

# I. Load required libraries
library(BayesLogit)  # For Polya-Gamma sampling (rpg)
library(MASS)        # For mvrnorm
library(ggplot2)
library(MCMCpack)    # For riwish() (not used now)
library(coda)
library(truncnorm)   # For truncated normal sampling
library(coda)
library(ggplot2)
library(tidyr)   
library(tidyverse)

constant = 0 ##if set to one, we use the normal time, if set to zero we only consider mixing 

#Pick a random individual
pick_random_real <- function(n = 1, min = 0, max = 100) {
  runif(n, min = min, max = max)
}
set.seed(123)
rnumb <- ceiling(pick_random_real())  # un número real


# II. Set seed for reproducibility
set.seed(123)

# III. Set parameters
n <- 100       # Number of individuals
t_obs <- 30    # Observations per individual
p <- 3         # Number of predictors (fixed effects)
q <- 2         # Number of individual-level covariates (random effects)
reg <- 0       # Regularization parameter (if needed)
num_iters <- 2000
unbalance <- 6

# IV. Simulate data
# i. Basic data:
X <- array(rnorm(n * t_obs * p), dim = c(n, t_obs, p))  # Predictor array (n x t_obs x p)
X[, , 1] <- 1   # Force intercept = 1
Z <- matrix(rnorm(n * q), nrow = n, ncol = q)          # Individual-level covariates (n x q)

# ii. Additional elements:
V_alpha_true <- matrix(c(1, 0.5,
                         0.5, 3), nrow = q, ncol = q)    # True covariance for random effects
alpha_true <- mvrnorm(n, mu = rep(0, q), Sigma = V_alpha_true)  # True random effects (n x q)
beta_true <- c(-unbalance, 2, -1.5)  

filename1 <- paste0(
  "Output/V_samplesnaug_",
  n, "_", t_obs, "_", num_iters, "_", unbalance, ".RData"
)

filename2 <- paste0(
  "Output/V_samplesaug_",
  n, "_", t_obs, "_", num_iters, "_", unbalance, ".RData"
)

num_hmc <- num_iters - num_iters*0.1
filename3 <- paste0(
  "Output/V_sampleshmc_",
  n, "_", t_obs, "_", num_hmc, "_", unbalance, ".RData"
)

filename4 <- paste0(
  "Output/V_samplesupgd_",
  n, "_", t_obs, "_", num_iters, "_", unbalance, ".RData"
)

load(filename1)
load(filename2)
load(filename3)
load(filename4)
#load("V_sampleslatentvarrealvlpha.RData") 
#load("V_sampleslatentvarrealalpha.RData") 
#load("V_samplesnlatentvarrealalpha.RData") 
#load("V_samplesauglatentvarrealvlpha.RData") 

#load("V_samplesnlatentvarrealalphaeqsize.RData") 
#load("V_sampleslatentvarrealalphaeqsize.RData") 


# Load required package

# Define the density_overlay function with a maximum of 4 input vectors
density_overlay <- function(..., 
                            vertical_line = NULL, 
                            vline_color = "black", 
                            vline_linetype = "dashed", 
                            plot_title = "Overlay Density Plot", 
                            xlab = "Value", 
                            ylab = "Density",
                            labels = NULL) {
  # Capture the list of vectors provided via ellipsis
  vec_list <- list(...)
  
  # Check that at least two vectors were provided
  if (length(vec_list) < 2) {
    stop("Please provide at least two vectors for comparison.")
  }
  
  # Limit the number of vectors to 4.
  if (length(vec_list) > 4) {
    warning("More than 4 vectors provided. Only the first 4 will be used.")
    vec_list <- vec_list[1:4]
  }
  
  # Determine names/labels for each vector
  if (!is.null(labels)) {
    if (length(labels) != length(vec_list)) {
      warning("Length of labels is not equal to the number of vectors provided. Using default names.")
      vec_names <- names(vec_list)
      if (is.null(vec_names)) {
        vec_names <- paste("Vector", seq_along(vec_list))
      }
    } else {
      vec_names <- labels
    }
  } else {
    vec_names <- names(vec_list)
    if (is.null(vec_names)) {
      vec_names <- paste("Vector", seq_along(vec_list))
    }
  }
  
  # Combine the vectors into one data frame.
  df <- data.frame(value = numeric(0), vector = character(0), stringsAsFactors = FALSE)
  for (i in seq_along(vec_list)) {
    tmp <- data.frame(value = vec_list[[i]], 
                      vector = rep(vec_names[i], length(vec_list[[i]])),
                      stringsAsFactors = FALSE)
    df <- rbind(df, tmp)
  }
  
  # Create the density plot with overlaid densities
  p <- ggplot(df, aes(x = value, color = vector)) +
    geom_density(size = 1) +
    labs(title = plot_title, x = xlab, y = ylab, color = "Vector") +
    theme_minimal()
  
  # Add vertical line(s) if the 'vertical_line' parameter is provided
  if (!is.null(vertical_line)) {
    for (v in vertical_line) {
      p <- p + geom_vline(xintercept = v, color = vline_color, 
                          linetype = vline_linetype, size = 1)
    }
  }
  
  return(p)
}


#### Some basic Statistics. 

#SE and posterior mean 

# Función para resumir una cadena MCMC
# Requiere el paquete coda
# Función para resumir una cadena MCMC e incluir prueba de Geweke
# Requiere el paquete coda
summary_mcmc <- function(draws) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Please install the 'coda' package to use this function.")
  }
  
  # Convertir a mcmc y extraer vector numérico (si fuera multivariante, toma la primera columna)
  draws_mcmc <- coda::as.mcmc(draws)
  draws_vec  <- as.numeric(draws_mcmc)
  
  # Estadísticos básicos
  m   <- mean(draws_vec)
  mn  <- min(draws_vec)
  mx  <- max(draws_vec)
  
  # Effective sample size
  n_eff <- coda::effectiveSize(draws_mcmc)
  # Desviación estándar de la cadena
  sd_chain <- sd(draws_vec)
  # Monte Carlo standard error
  mcse <- sd_chain / sqrt(n_eff)
  
  # Geweke diagnostic (z-score)
  gw <- coda::geweke.diag(draws_mcmc)
  z_geweke <- as.numeric(gw$z)[1]
  
  # Devolver lista con resultados
  list(
    mean        = m,
    minimum     = mn,
    maximum     = mx,
    mcse        = mcse,
    geweke_z    = z_geweke
  )
}


#### Plots for beta. 

#1. Compare Behaviors. 

density_plot <- density_overlay(vec1 = beta_1_samplesnaug, 
                                vec2 = beta_1_samplesaug,
                                vec3 = beta_1_sampleshmc,
                                vec4 = beta_1_samplesupgd,
                                vertical_line = beta_true[1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent", "HMC Stan" ,  "UPG Full Delta"))

print(density_plot)

density_plot <- density_overlay(vec1 = beta_2_samplesnaug, 
                                vec2 = beta_2_samplesaug,
                                vec3 = beta_2_sampleshmc,
                                vec4 = beta_2_samplesupgd,
                                vertical_line = beta_true[2],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent", "HMC Stan" , "UPG Gamma Delta"))

print(density_plot)

density_plot <- density_overlay(vec1 = beta_3_samplesnaug, 
                                vec2 = beta_3_samplesaug,
                                vec3 = beta_3_sampleshmc,
                                vec4 = beta_3_samplesupgd,
                                vertical_line = beta_true[3],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent", "HMC Stan" , "UPG Gamma Delta"))

print(density_plot)

#2. Get the Computational Efficiency. 

#Set computational time to 1, just in case: 

time_naug <- as.numeric(timesamples_naug, units = "secs")^constant 
time_aug <- as.numeric(timesamples_aug, units = "secs")^constant  
#time_upg <- as.numeric(timesamples_hmc, units = "secs")^constant  
time_hmc <- as.numeric(timesamples_hmc, units = "secs")^constant 
time_upgd <-as.numeric(timesamples_upgd, units = "secs")^constant  

#ESS. 
beta_1_samplesnaug <- as.mcmc(beta_1_samplesnaug)
essb1naug <- effectiveSize(beta_1_samplesnaug)/time_naug
beta_1_samplesaug <- as.mcmc(beta_1_samplesaug)
essb1aug <- effectiveSize(beta_1_samplesaug)/time_aug
beta_1_sampleshmc <- as.mcmc(beta_1_sampleshmc)
essb1hmc <- effectiveSize(beta_1_sampleshmc)/time_hmc
beta_1_samplesupgd <- as.mcmc(beta_1_samplesupgd)
essb1upgd <- effectiveSize(beta_1_samplesupgd)/time_upgd



beta_2_samplesnaug <- as.mcmc(beta_2_samplesnaug)
essb2naug <- effectiveSize(beta_2_samplesnaug)/time_naug
beta_2_samplesaug <- as.mcmc(beta_2_samplesaug)
essb2aug <- effectiveSize(beta_2_samplesaug)/time_aug
beta_2_sampleshmc <- as.mcmc(beta_2_sampleshmc)
essb2hmc <- effectiveSize(beta_2_sampleshmc)/time_hmc
beta_2_samplesupgd <- as.mcmc(beta_2_samplesupgd)
essb2upgd <- effectiveSize(beta_2_samplesupgd)/time_upgd

beta_3_samplesnaug <- as.mcmc(beta_3_samplesnaug)
essb3naug <- effectiveSize(beta_3_samplesnaug)/time_naug
beta_3_samplesaug <- as.mcmc(beta_3_samplesaug)
essb3aug <- effectiveSize(beta_3_samplesaug)/time_aug
beta_3_sampleshmc <- as.mcmc(beta_3_sampleshmc)
essb3hmc <- effectiveSize(beta_3_sampleshmc)/time_hmc
beta_3_samplesupgd <- as.mcmc(beta_3_samplesupgd)
essb3upgd <- effectiveSize(beta_3_samplesupgd)/time_upgd

# build a data frame with correct factor levels
ess_df <- data.frame(
  beta = rep(paste0("Beta", 1:3), each = 4),   # 3 betas × 4 métodos = 12
  method = factor(
    rep(c("Polya Gamma Augmentation",
          "Polya Gamma + Latent Variable Augmentation",
          "HMC Stan",
          "UPG Gamma Delta"), times = 3),
    levels = c("Polya Gamma Augmentation",
               "Polya Gamma + Latent Variable Augmentation",
               "HMC Stan",
               "UPG Gamma Delta")                 # orden en la leyenda
  ),
  ess = c(
    essb1naug, essb1aug, essb1hmc, essb1upgd,
    essb2naug, essb2aug, essb2hmc, essb2upgd,
    essb3naug, essb3aug, essb3hmc, essb3upgd
  )
)


# grouped bar chart with legend control
ggplot(ess_df, aes(x = beta, y = ess, fill = method)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7)) +
  scale_fill_manual(
    values = c(
      "Polya Gamma Augmentation"                       = "#1f77b4",  # azul
      "Polya Gamma + Latent Variable Augmentation"     = "#d62728",  # rojo
      "HMC Stan"                                       = "#2ca02c",  # verde
      "UPG Gamma Delta"                                    = "#ff7f0e"   # naranja
    ),
    breaks = c("Polya Gamma Augmentation",
               "Polya Gamma + Latent Variable Augmentation",
               "HMC Stan",
               "UPG Gamma Delta")                          # orden en la leyenda
  ) +
  scale_x_discrete(labels = c(
    Beta1 = expression(beta[1]),
    Beta2 = expression(beta[2]),
    Beta3 = expression(beta[3])
  )) +
  labs(
    x    = NULL,
    y    = "ESS por segundo",
    fill = "Sampler"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = "Output/Results/essbeta.pdf", width = 8, height = 4, dpi = 300)

# #3. Tracepots.
# 
# #Plot for PG and Hybrid.
# plot_beta_trace <- function(idx, k = NULL) {
#   # build the object names
#   name_naug <- paste0("beta_", idx, "_samplesnaug")
#   name_aug  <- paste0("beta_", idx, "_samplesaug")
#   
#   # grab them from the environment
#   b_naug <- get(name_naug)
#   b_aug  <- get(name_aug)
#   
#   # optionally restrict to first k draws
#   if (!is.null(k)) {
#     b_naug <- head(b_naug, k)
#     b_aug  <- head(b_aug,  k)
#   }
#   
#   # build a data.frame
#   df <- data.frame(
#     iteration = c(seq_along(b_naug), seq_along(b_aug)),
#     value     = c(b_naug, b_aug),
#     sampler   = factor(
#       c(
#         rep("Polya Gamma Augmentation", length(b_naug)),
#         rep("Polya Gamma + Latent Variable", length(b_aug))
#       ),
#       levels = c("Polya Gamma Augmentation", "Polya Gamma + Latent Variable")
#     )
#   )
#   
#   # plot with horizontal line for true value
#   ggplot(df, aes(x = iteration, y = value, color = sampler)) +
#     geom_line(alpha = 0.7) +
#     geom_hline(yintercept = beta_true[idx], 
#                linetype = "dashed", 
#                color = "black", 
#                size = 0.9) +
#     scale_color_manual(
#       values = c(
#         "Polya Gamma Augmentation"         = "blue", 
#         "Polya Gamma + Latent Variable"    = "red"
#       )
#     ) +
#     labs(
#       x     = "Iteration",
#       y     = bquote(beta[.(idx)]),
#       color = "Sampler"
#     ) +
#     theme_minimal() +
#     theme(legend.position = "bottom")
# }
# 
# #Plot for UPG. 
# plot_beta_traceupg <- function(idx, k = NULL) {
#   # build the object names
#   name_upg <- paste0("beta_", idx, "_samplesupg")
#   name_upgd  <- paste0("beta_", idx, "_samplesupgd")
#   
#   # grab them from the environment
#   b_upg <- get(name_upg)
#   b_upgd  <- get(name_upgd)
#   
#   # optionally restrict to first k draws
#   if (!is.null(k)) {
#     b_upg <- head(b_hmc, k)
#     b_upgd  <- head(b_upgd,  k)
#   }
#   
#   # build a data.frame
#   df <- data.frame(
#     iteration = c(seq_along(b_upg), seq_along(b_upgd)),
#     value     = c(b_hmc, b_upgd),
#     sampler   = factor(
#       c(
#         rep("Ultimate Polya Gamma Full", length(b_upg)),
#         rep("Ultimate Polya Gamma Partial", length(b_upgd))
#       ),
#       levels = c("Ultimate Polya Gamma Full", "Ultimate Polya Gamma Partial")
#     )
#   )
#   
#   # plot with horizontal line for true value
#   ggplot(df, aes(x = iteration, y = value, color = sampler)) +
#     geom_line(alpha = 0.7) +
#     geom_hline(yintercept = beta_true[idx], 
#                linetype = "dashed", 
#                color = "black", 
#                size = 0.9) +
#     scale_color_manual(
#       values = c(
#         "Ultimate Polya Gamma Full"         = "blue", 
#         "Ultimate Polya Gamma Partial"    = "red"
#       )
#     ) +
#     labs(
#       x     = "Iteration",
#       y     = bquote(beta[.(idx)]),
#       color = "Sampler"
#     ) +
#     theme_minimal() +
#     theme(legend.position = "bottom")
# }
# 
# 
# k <- num_iters
# plot_beta_trace(1, k) 
# ggsave(filename = "Output/Results/tracebeta1.pdf", width = 8, height = 4, dpi = 300)
# plot_beta_traceupg(1, k) 
# ggsave(filename = "Output/Results/tracebeta1upg.pdf", width = 8, height = 4, dpi = 300)
# 
# plot_beta_trace(2, k) 
# ggsave(filename = "Output/Results/tracebeta2.pdf", width = 8, height = 4, dpi = 300)
# plot_beta_traceupg(2, k) 
# ggsave(filename = "Output/Results/tracebeta2upg.pdf", width = 8, height = 4, dpi = 300)
# 
# plot_beta_trace(3, k) 
# ggsave(filename = "Output/Results/tracebeta3.pdf", width = 8, height = 4, dpi = 300)
# plot_beta_traceupg(3, k) 
# ggsave(filename = "Output/Results/tracebeta3upg.pdf", width = 8, height = 4, dpi = 300)

#### Plots for V_alpha

#1. Compare Behaviors.
density_plot <- density_overlay(vec1 = V1_naug,
                                vec2 = V1_aug,
                                vec3 = V1_hmc,
                                vec4 = V1_upgd,
                                vertical_line = V_alpha_true[1,1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent" , "HMC Stan" , "UPG Gamma Delta"))

print(density_plot)

density_plot <- density_overlay(vec1 = V2_naug,
                                vec2 = V2_aug,
                                vec3 = V2_hmc,
                                vec4 = V2_upgd,
                                vertical_line = V_alpha_true[1,1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent" , "HMC Stan" , "UPG Gamma Delta"))

print(density_plot)

density_plot <- density_overlay(vec1 = V21_naug,
                                vec2 = V21_aug,
                                vec3 = V21_hmc,
                                vec4 = V21_upgd,
                                vertical_line = V_alpha_true[2,1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent" , "HMC Stan" , "UPG Gamma Delta"))

print(density_plot)


#2. Get the Computational Efficiency.

#ESS.
V1_naug <- as.mcmc(V1_naug)
essV1_naug <- effectiveSize(V1_naug)/time_naug
V1_aug <- as.mcmc(V1_aug)
essV1_aug <- effectiveSize(V1_aug)/time_aug
V1_hmc <- as.mcmc(V1_hmc)
essV1_hmc <- effectiveSize(V1_hmc)/time_hmc
V1_upgd <- as.mcmc(V1_upgd)
essV1_upgd <- effectiveSize(V1_upgd)/time_upgd

V2_naug <- as.mcmc(V2_naug)
essV2_naug <- effectiveSize(V2_naug)/time_naug
V2_aug <- as.mcmc(V2_aug)
essV2_aug <- effectiveSize(V2_aug)/time_aug
V2_hmc <- as.mcmc(V2_hmc)
essV2_hmc <- effectiveSize(V2_hmc)/time_hmc
V2_upgd <- as.mcmc(V2_upgd)
essV2_upgd <- effectiveSize(V2_upgd)/time_upgd

V21_naug <- as.mcmc(V21_naug)
essV21_naug <- effectiveSize(V21_naug)/time_naug
V21_aug <- as.mcmc(V21_aug)
essV21_aug <- effectiveSize(V21_aug)/time_aug
V21_hmc <- as.mcmc(V21_hmc)
essV21_hmc <- effectiveSize(V21_hmc)/time_hmc
V21_upgd <- as.mcmc(V21_upgd)
essV21_upgd <- effectiveSize(V21_upgd)/time_upgd

## 1) Construir data-frame ----------------------------------------------
essV_df <- data.frame(
  param  = rep(c("V11", "V22", "V21"), each = 4),   # 3 parámetros × 4 métodos = 12 filas
  method = factor(
    rep(c("Polya Gamma Augmentation",
          "Polya Gamma + Latent Variable Augmentation",
          "HMC Stan",
          "UPG Gamma Delta"),
        times = 3),
    levels = c("Polya Gamma Augmentation",
               "Polya Gamma + Latent Variable Augmentation",
               "HMC Stan",
               "UPG Gamma Delta")                       # orden fijo en la leyenda
  ),
  ess = c(
    essV1_naug,  essV1_aug,  essV1_hmc,  essV1_upgd,
    essV2_naug,  essV2_aug,  essV2_hmc,  essV2_upgd,
    essV21_naug, essV21_aug, essV21_hmc, essV21_upgd
  )
)

## 2) Gráfico -------------------------------------------------------------
library(ggplot2)

ggplot(essV_df, aes(x = param, y = ess, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  scale_fill_manual(
    values = c(
      "Polya Gamma Augmentation"                       = "#1f77b4", # azul
      "Polya Gamma + Latent Variable Augmentation"     = "#d62728", # rojo
      "HMC Stan"                                       = "#2ca02c", # verde
      "UPG Gamma Delta"                                    = "#ff7f0e"  # naranja
    )
  ) +
  scale_x_discrete(labels = c(
    V11 = expression(V[alpha][11]),
    V22 = expression(V[alpha][22]),
    V21 = expression(V[alpha][21])
  )) +
  labs(
    x    = NULL,
    y    = "ESS por segundo",
    fill = "Sampler"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("Output/Results/essvalpha.pdf", width = 8, height = 4, dpi = 300)


# #3) traceplots
# 
# plot_V_trace <- function(idx, k = NULL) {
#   # map idx→object names
#   naug_names <- c("V1_naug", "V2_naug", "V21_naug")
#   aug_names  <- c("V1_aug",  "V2_aug",  "V21_aug")
#   
#   # grab the correct vectors
#   v_naug <- as.numeric(get(naug_names[idx]))
#   v_aug  <- as.numeric(get(aug_names[idx]))
#   
#   # optional head()
#   if (!is.null(k)) {
#     v_naug <- head(v_naug, k)
#     v_aug  <- head(v_aug,  k)
#   }
#   
#   # assemble data
#   df <- data.frame(
#     iteration = c(seq_along(v_naug), seq_along(v_aug)),
#     value     = c(v_naug, v_aug),
#     sampler   = factor(
#       c(rep("Polya Gamma Augmentation", length(v_naug)),
#         rep("Polya Gamma + Latent Variable", length(v_aug))),
#       levels = c("Polya Gamma Augmentation","Polya Gamma + Latent Variable")
#     )
#   )
#   
#   # nice y-labels for each entry of V_alpha
#   ylabels <- list(
#     expression(V[alpha][1]),       # idx=1
#     expression(V[alpha][2]),       # idx=2
#     expression(V[alpha][2][1])     # idx=3 (off-diagonal)
#   )
#   
#   # true population value from V_alpha_true
#   pop_val <- switch(idx,
#                     V_alpha_true[1, 1],  # idx = 1
#                     V_alpha_true[2, 2],  # idx = 2
#                     V_alpha_true[2, 1]   # idx = 3
#   )
#   
#   ggplot(df, aes(x = iteration, y = value, color = sampler)) +
#     geom_line(alpha = 0.7) +
#     geom_hline(yintercept = pop_val,
#                linetype = "dashed",
#                color = "black",
#                size = 0.9) +
#     scale_color_manual(
#       values = c(
#         "Polya Gamma Augmentation"      = "blue",
#         "Polya Gamma + Latent Variable" = "red"
#       )
#     ) +
#     labs(
#       x     = "Iteration",
#       y     = ylabels[[idx]],
#       color = "Sampler"
#     ) +
#     theme_minimal() +
#     theme(legend.position = "bottom")
# }
# 
# 
# # Examples:
# plot_V_trace(1, k = 5000)         # V[alpha][1,1]
# ggsave(filename = "Output/Results/tracevalpha1.pdf", width = 8, height = 4, dpi = 300)
# plot_V_trace(2, k = 5000) # first 200 draws for V[alpha][2,2]
# ggsave(filename = "Output/Results/tracevalpha2.pdf", width = 8, height = 4, dpi = 300)
# plot_V_trace(3, k = 5000)      # the 'third' V_alpha entry (correlation)
# ggsave(filename = "Output/Results/tracevalpha21.pdf", width = 8, height = 4, dpi = 300)



#### Plots for alpha

#1. Compare Behaviors.


#Componente 1 individuo 1
density_plot <- density_overlay(vec1 = a1_nag[,1],
                                vec2 = a1_ag[,1],
                                vec3 = a1_hmc[,1],
                                vec4 = a1_upgd[,1],
                                vertical_line = alpha_true[1,1],
                                labels = c("Augmented by Gamma", "Augmented by Gamma and latent", "HMC Stan", "UPG Gamma Delta"))

print(density_plot)


density_plot <- density_overlay(
  vec1 = a2_nag[,2],   # alpha[1,2] from PG
  vec2 = a2_ag[,2],    # alpha[1,2] from PG+Latent
  vec3 = a2_hmc[,2],       # alpha[1,2] from UPG Full
  vec4 = a2_upgd[,2],      # alpha[1,2] from UPG Partial
  vertical_line = alpha_true[1,2],  # <- corregido
  labels = c("Augmented by Gamma",
             "Augmented by Gamma and latent",
             "HMC Stan",
             "UPG Gamma Delta")
)

print(density_plot)


#2. Get the Computational Efficiency.

#ESS.
test <- as.matrix(alpha_samples_aug[,1,1])
dim(test)


## ────────────────────────────────────────────────────────────────
## 0. Paquetes
library(coda)            # effectiveSize()
library(dplyr)           # manipulación
library(tidyr)           # pivot_longer()
library(ggplot2)         # gráficos

## ────────────────────────────────────────────────────────────────
## 1. Inicializar df_ess con las nuevas columnas
df_ess <- data.frame(
  ind             = integer(n),
  ess_alpha1naug  = numeric(n),
  ess_alpha1aug   = numeric(n),
  ess_alpha1hmc   = numeric(n),
  ess_alpha1upgd  = numeric(n),
  ess_alpha2naug  = numeric(n),
  ess_alpha2aug   = numeric(n),
  ess_alpha2hmc   = numeric(n),
  ess_alpha2upgd  = numeric(n),
  meanalpha1_max  = numeric(n),
  meanalpha2_max  = numeric(n)
)

## ────────────────────────────────────────────────────────────────
## 2. Rellenar fila j  (bucle sobre individuos)
for (j in seq_len(n)) {
  
  ## ─ cadenas del individuo j ────────────────────────────────────
  chain1aug   <-  alpha_samples_aug[ , j, 1]
  chain2aug   <-  alpha_samples_aug[ , j, 2]
  chain1naug  <- alpha_samples_naug[ , j, 1]
  chain2naug  <- alpha_samples_naug[ , j, 2]
  chain1hmc   <-  alpha_samples_hmc[ , j, 1]
  chain2hmc   <-  alpha_samples_hmc[ , j, 2]
  chain1upgd  <- alpha_samples_upgd[ , j, 1]
  chain2upgd  <- alpha_samples_upgd[ , j, 2]
  
  ## ─ ESS / segundo ──────────────────────────────────────────────
  ess1aug  <- effectiveSize(as.mcmc(chain1aug))  /
    time_aug
  ess2aug  <- effectiveSize(as.mcmc(chain2aug))  /
    time_aug
  ess1nau  <- effectiveSize(as.mcmc(chain1naug)) /
    time_naug
  ess2nau  <- effectiveSize(as.mcmc(chain2naug)) /
    time_naug
  ess1hmc  <- effectiveSize(as.mcmc(chain1hmc))  /
    time_hmc
  ess2hmc  <- effectiveSize(as.mcmc(chain2hmc))  /
    time_hmc
  ess1upgd <- effectiveSize(as.mcmc(chain1upgd)) /
    time_upgd
  ess2upgd <- effectiveSize(as.mcmc(chain2upgd)) /
    time_upgd
  
  ## ─ medias ─────────────────────────────────────────────────────
  m1 <- c(mean(chain1naug), mean(chain1aug),
          mean(chain1hmc),  mean(chain1upgd))
  m2 <- c(mean(chain2naug), mean(chain2aug),
          mean(chain2hmc),  mean(chain2upgd))
  
  ## ─ rellenar fila j ────────────────────────────────────────────
  df_ess[j, ] <- list(
    ind             = j,
    ess_alpha1naug  = ess1nau,
    ess_alpha1aug   = ess1aug,
    ess_alpha1hmc   = ess1hmc,
    ess_alpha1upgd  = ess1upgd,
    ess_alpha2naug  = ess2nau,
    ess_alpha2aug   = ess2aug,
    ess_alpha2hmc   = ess2hmc,
    ess_alpha2upgd  = ess2upgd,
    meanalpha1_max  = max(m1),
    meanalpha2_max  = max(m2)
  )
}

## ────────────────────────────────────────────────────────────────
## 3. Función de plot con cuatro barras
plot_ess <- function(df, alpha = 1, n = nrow(df)) {
  stopifnot(alpha %in% 1:2,
            n >= 1, n <= nrow(df))
  
  ## columnas de ESS por algoritmo --------------------------------
  col_pat <- sprintf("^ess_alpha%d", alpha)
  col_mean <- sprintf("meanalpha%d_max", alpha)
  
  ## ordenar y recortar -------------------------------------------
  df_sub <- df[order(df[[col_mean]]), ][1:n, ]
  df_sub$ind <- seq_len(n)               # re-indexar eje x
  
  ## pasar a formato largo ----------------------------------------
  df_long <- df_sub |>
    select(ind, matches(col_pat)) |>
    pivot_longer(-ind,
                 names_to   = "Sampler",
                 values_to  = "ESSps") |>
    mutate(Sampler = recode(Sampler,
                            !!sprintf("ess_alpha%dnaug", alpha)  := "Polya Gamma",
                            !!sprintf("ess_alpha%daug",  alpha)  := "PG + Utility",
                            !!sprintf("ess_alpha%dhmc",  alpha)  := "HMC Stan",
                            !!sprintf("ess_alpha%dupgd", alpha)  := "UPG Gamma Delta"
    ))
  
  ## gráfico -------------------------------------------------------
  ggplot(df_long,
         aes(x = ind, y = ESSps, fill = Sampler)) +
    geom_col(position = position_dodge()) +
    scale_fill_manual(
      name   = "Sampler",
      values = c("Polya Gamma"     = "blue",
                 "PG + Utility"    = "red",
                 "HMC Stan"        = "forestgreen",
                 "UPG Gamma Delta"     = "orange")
    ) +
    scale_x_continuous(breaks = seq(0, n, by = max(1, floor(n/10)))) +
    labs(x = "Índice ordenado por media posterior",
         y = "ESS / s") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

## ────────────────────────────────────────────────────────────────
## Uso de ejemplo
## plot_ess(df_ess, alpha = 1)     # componente α₁
## plot_ess(df_ess, alpha = 2, n = 50)



plot_ess(df_ess, alpha = 1, n = 50)
ggsave("Output/Results/ess_alpha1.pdf", width = 8, height = 4, dpi = 300)

plot_ess(df_ess, alpha = 2, n = 100)
ggsave("Output/Results/ess_alpha2.pdf", width = 8, height = 4, dpi = 300)

plot_ess_diff(df_ess, alpha = 1, n = 100)
ggsave("Output/Results/ess_alpha1diff.pdf", width = 8, height = 4, dpi = 300)
plot_ess_diff(df_ess, alpha = 2, n = 100)
ggsave("Output/Results/ess_alpha2diff.pdf", width = 8, height = 4, dpi = 300)


#3) Traceplots.


plot_alpha_trace <- function(rnumb, chain_num = 1, n_iter = NULL) {
  # Validación del argumento
  if (!(chain_num %in% c(1, 2))) {
    stop("El argumento 'chain_num' debe ser 1 o 2.")
  }
  
  # Seleccionar las cadenas correctas
  chain_aug  <- alpha_samples_aug[, rnumb, chain_num]
  chain_naug <- alpha_samples_naug[, rnumb, chain_num]
  
  # Limitar número de iteraciones si se especifica
  if (!is.null(n_iter)) {
    chain_aug  <- head(chain_aug, n_iter)
    chain_naug <- head(chain_naug, n_iter)
  }
  
  # Valor verdadero de α
  alpha_true_val <- alpha_true[rnumb, chain_num]
  
  # Crear dataframe
  df_trace <- tibble(
    Iteration     = 1:length(chain_aug),
    `Polya Gamma Augmentation`         = chain_naug,
    `Polya Gamma + Latent Variable`    = chain_aug
  )
  
  # Formato largo
  df_long <- df_trace %>%
    pivot_longer(cols = -Iteration, names_to = "Sampler", values_to = "Value") %>%
    mutate(Sampler = factor(Sampler, levels = c(
      "Polya Gamma Augmentation", "Polya Gamma + Latent Variable"
    )))
  
  # Graficar
  ggplot(df_long, aes(x = Iteration, y = Value, color = Sampler)) +
    geom_line(alpha = 0.8) +
    geom_hline(yintercept = alpha_true_val, 
               linetype = "dashed", color = "black", size = 0.9) +
    scale_color_manual(
      values = c(
        "Polya Gamma Augmentation"      = "blue",
        "Polya Gamma + Latent Variable" = "red"
      )
    ) +
    labs(
      title = paste(""),
      y     = bquote(alpha[.(chain_num)]),
      x     = "Iteration",
      color = "Sampler"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_alpha_trace(rnumb,chain_num = 1, n_iter = 1000)
ggsave(filename = "Output/Results/alpharandomtrace.pdf", width = 8, height = 4, dpi = 300)

