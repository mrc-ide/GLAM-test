# simple/mcmc.R 
# Author: Gina Cuomo-Dannenburg
# Date: 23/04/2025
# Purpose: simulate a simple bank of datasets for testing

# installation for testing of specified version of {GLAM} - SHA1 (c9cdc36f)
# devtools::install_github("mrc-ide/GLAM@test/gina_working_example")
# this specific version does not currently work at actually fitting the data
library(GLAM)

# load other libraries
library(dplyr)
library(ggplot2)

# read source R script of additional functions
source("functions.R")

# read in lookup table
lookup <- readRDS("simple/lookup.RDS")

# set parameters of the MCMC process -- increase these once we have working {glam}
burnin_iterations <- 1e2
sampling_iterations <- 5e2
num_chains <- 1
num_rungs <- 1

# choose parameter values that we want to select
# chosen parameters in alignment with the deploy script for the moment
samp_time <- seq(0, 10, 1) # keep the same for all in this case
haplo_freqs <- rep(0.05, 20) # equal frequency of all haplotypes
theta_val <- 2
decay_val <- 0.1
lambda_val <- 1/10 # set to n_inf/max(samp_time)
sens <- 1
n_inf <- 1
cohort_size <- 10

# need to update when we have multiple samp_time to ensure that we filter by this
sim_ids <- lookup |>
  dplyr::filter(theta == theta_val,
                decay == decay_val,
                lambda == lambda_val,
                sens == sens,
                n == n_inf,
                cohort_size == cohort_size) |>
  dplyr::pull(sim_id)

for (i in 1:nrow(sim_ids)) {
  id <- sim_ids[i]
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
         theta = theta_val,
         decay_rate = decay_val,
         sens = sens,
         n_infections = mapply(function(x) length(x$t_inf), sim$raw_list),
         infection_times = NULL,
         max_infections = 20,
         chains = num_chains, rungs = num_rungs)
  
  # run burn-in and sampling
  g$burn(iterations = burnin_iterations) |>
    system.time()
  
  g$sample(iterations = sampling_iterations)  |>
    system.time()
  
  # save output
  saveRDS(g, paste0("simple/mcmc_outputs/out", id,".RDS"))
  
  # summary outputs
  
  
}


