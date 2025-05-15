# simple/simulate_function.R 
# Author: Gina Cuomo-Dannenburg
# Date: 2025-04-25
# Purpose: A functional form of simulate.R to enable cluster integration

# libraries are loaded and lookup_null.RDS is read in as dependencies
# function of i where i is a row of lookup_null.RDS

lookup_null <- readRDS("null/lookup_null.RDS")

simulate_data_null <- function(i) {
  set.seed(i)
  samp_time <- unlist(lookup_null$samp_time[i])
  id <- lookup_null$sim_id[i]
  
  # simulate the cohort corresponding to the parameter values
  sim <- GLAM::sim_cohort(
    n = lookup_null$cohort_size[i],
    samp_time = samp_time,
    haplo_freqs = unlist(lookup_null$haplo_freqs[i]),
    theta = lookup_null$theta[i],
    lambda = 0.1,
    decay_rate = lookup_null$decay[i],
    sens = lookup_null$sens[i],
    n_inf = rep(NULL, lookup_null$cohort_size[i]),
    return_full = TRUE
  )
  saveRDS(sim, paste0("null/data/sim", id, ".RDS"))
  
}

# create a function to run a block of simulations with defined size
simulate_block_null <- function(block_num, block_size) {
  iterations <- seq(((block_num-1) * block_size)+1, block_num*block_size)
  
  lapply(iterations, simulate_data_null)
}
