# simple/simulate_function.R 
# Author: Gina Cuomo-Dannenburg
# Date: 2025-04-25
# Purpose: Generate the lookup table - included here for ease of duplicating for 
# banks of simulations and to allow it to also be run on the cluster

# add this to hipercow environment
source("functions.R")
require(dplyr)
require(tidyr)

# define the vector of options that we want to cycle over for the simulations
samp_time <- seq(0, 10, 1) # keep the same for all in this case
haplo_freqs <- rep(0.05, 20) # equal frequency of all haplotypes
theta_vec <- c(0.5, 1, 2, 3, 4)
decay_vec <- c(0.01,0.1,0.3)
sens <- 1
n_inf <- NULL
cohort_size <- 10
repetitions <- 10 # set repetitions to be 100 so we can test how often MCMC 
# returns true parameters

# generate a lookup table with all the parameter values
lookup <- 
  # start with cycling over the vectors
  tidyr::expand_grid(theta = theta_vec, 
                     decay = decay_vec, 
                     repetition = 1:repetitions) |>
  dplyr::mutate(samp_time = list(samp_time),
                haplo_freqs = list(haplo_freqs),
                sens = sens,
                cohort_size = cohort_size) |>
  dplyr::select(cohort_size, samp_time, haplo_freqs, 
                theta, decay, sens, repetition) |>
  dplyr::mutate(sim_id = zero_pad_fixed(row_number(), 5)) # generate a simulation id number
saveRDS(lookup, "null/lookup_null.RDS")
