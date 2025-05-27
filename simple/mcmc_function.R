# simple/mcmc_function.R 
# Author: Gina Cuomo-Dannenburg
# Date: 2025-04-25
# Purpose: A functional form of mcmc.R to enable cluster integration

# libraries are loaded and lookup.RDS is read in as dependencies
# function of i where i is a row of lookup.RDS

# set parameters of the MCMC process -- increase these once we have working {glam}
lookup <- readRDS("simple/lookup.RDS")

burnin_iterations <- 1e2
sampling_iterations <- 5e2
num_chains <- 1

# set priors for all of the mcmc fittings
# define priors
lambda_prior <- function(x) dlnorm_reparam(x, mean = 0.5, sd = 0.2, return_log = TRUE)
theta_prior <- function(x) dlnorm_reparam(x, mean = 3, sd = 3, return_log = TRUE)
decay_rate_prior <- function(x) dlnorm_reparam(x, mean = 0.5, sd = 0.2, return_log = TRUE)
sens_prior <- function(x) dbeta(x, shape1 = 99, shape2 = 1, log = TRUE)

# fix haplo freqs 
haplo_freqs <- rep(0.05, 20) # equal frequency of all haplotypes
max_infections <- 20

mcmc_fitting <- function(i) { 
  set.seed(i)
  id <- lookup$sim_id[i]
  sim <- readRDS(paste0("simple/data/sim", id, ".RDS"))
  
  # create a new MCMC object and load data
  g <- GLAM::glam_mcmc$new(df_data = sim$df_data)
  
  # initialize. If parameters are set to NULL, they are estimated. If they are
  # given a value, they take this fixed value. Useful for giving the inference
  # method the correct true value of some parameters, to diagnose how well it
  # estimates others.
  g$init(lambda_prior = lambda_prior,
         theta_prior = theta_prior,
         decay_rate_prior = decay_rate_prior,
         sens_prior = sens_prior,
         haplo_freqs = haplo_freqs,
         chains = num_chains,
         max_infections = max_infections)
  
  # run burn-in and sampling
  g$burn(iterations = burnin_iterations) |>
    system.time()
  
  g$sample(iterations = sampling_iterations)  |>
    system.time()
  
  # save output
  saveRDS(g, paste0("simple/mcmc_outputs/out", id,".RDS"))
}

# create a function to fit a block of simulations with defined size
mcmc_fitting_block <- function(block_num, block_size) {
  iterations <- seq(((block_num-1) * block_size)+1, block_num*block_size)
  
  lapply(iterations, mcmc_fitting)
}
