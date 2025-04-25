# simple/simulate_function.R 
# Author: Gina Cuomo-Dannenburg
# Date: 2025-04-25
# Purpose: A functional form of simulate.R to enable cluster integration

# libraries are loaded and lookup.RDS is read in as dependencies

# function of i where i is a row of lookup.RDS

simulate_data <- function(i) {
  samp_time <- unlist(lookup$samp_time[i])
  n_inf <- lookup$n[i]
  id <- lookup$sim_id[i]
  
  # simulate the cohort corresponding to the parameter values
  sim <- GLAM::sim_cohort(
    n = lookup$cohort_size[i],
    samp_time = samp_time,
    haplo_freqs = unlist(lookup$haplo_freqs[i]),
    theta = lookup$theta[i],
    lambda = lookup$lambda[i],
    decay_rate = lookup$decay[i],
    sens = lookup$sens[i],
    n_inf = rep(n_inf, lookup$cohort_size[i]),
    return_full = TRUE
  )
  saveRDS(sim, paste0("simple/data/sim", id, ".RDS"))
  
}
