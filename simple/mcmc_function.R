# simple/mcmc_function.R 
# Author: Gina Cuomo-Dannenburg
# Date: 2025-04-25
# Purpose: A functional form of mcmc.R to enable cluster integration

# libraries are loaded and lookup.RDS is read in as dependencies
# function of i where i is a row of lookup.RDS

# set parameters of the MCMC process -- increase these once we have working {glam}
burnin_iterations <- 1e2
sampling_iterations <- 5e2
num_chains <- 1
num_rungs <- 1

mcmc_fitting <- function(i) { 
  set.seed(i)
  id <- lookup$sim_id[i]
  sim <- readRDS(paste0("simple/data/sim", id, ".RDS"))
  
  # create a new MCMC object and load data
  g <- glam_mcmc$new(df_data = sim$df_data)
  
  # initialize. If parameters are set to NULL, they are estimated. If they are
  # given a value, they take this fixed value. Useful for giving the inference
  # method the correct true value of some parameters, to diagnose how well it
  # estimates others.
  g$init(start_time = 0,
         end_time = 10,
         haplo_freqs = haplo_freqs,
         theta = NULL,
         decay_rate = NULL,
         lambda = NULL,
         sens = NULL,
         n_infections = NULL,
         infection_times = NULL,
         max_infections = 20,
         chains = num_chains, 
         rungs = num_rungs)
  
  # run burn-in and sampling
  g$burn(iterations = burnin_iterations) |>
    system.time()
  
  g$sample(iterations = sampling_iterations)  |>
    system.time()
  
  # save output
  saveRDS(g, paste0("simple/mcmc_outputs/out", id,".RDS"))
}
