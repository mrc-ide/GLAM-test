# simulate_data.R
# 
# Author: Gina Cuomo-Dannenburg
# Date: 03/04/2025
# 
# Purpose: simulate datasets for testing of the {mrc-ide/GLAM} package
# Script is intended to produce core examples of datasets, from basic to more
# complex to test the full bredth of functionality of {GLAM}, to validate methods etc.
# 
# Initial thoughts:
# - Cohorts with the same # total infections
#     * Does the MCMC accurately fit the timing of infections?
#     * Are there issues with saturation at higher #s of infections?
# - Low and high transmission 
#     * Is there and issue identifying specific parameters at high and low lambda?
# - Fixing individual level parameters and ability to fit global parameters
# - Inconsistent sampling time
# - Scalability with # of haplotypes
# - Frequency of specific haplotypes
# - Collinearity between specific parameters where they are likely to create a 
#   gradient given the likelihood and other factors? (~ shape and scale in some 
#   distributions)
# - How to deal with the fact that you may have multiple infections between sampling?
#   Is part of the issue fixing the # infections in this circumstance?

# installation for testing of specified version of {GLAM}
# devtools::install_github("mrc-ide/GLAM@test/gina_working_example")

library(GLAM)
set.seed(2)

# -----------------------
# SIMULATE DATA

# set simulation parameters
samp_time <- seq(0, 10, 1)
haplo_freqs <- rep(0.05, 20)
lambda <- 0.2
theta <- 2
decay_rate <- 0.1
sens <- 0.9
max_infections <- 20

# cohorts with n_inf equal - use the relationship between lambda and median n_inf
# from build_intuition.R: median(n_inf) ~ max(samp_time) * lambda

# simple example but scalable
cohort_size <- 10
max_time <- 10
n_inf <- 1:20

for(i in n_inf) {
  print(paste0("Simulate cohort with ", i, " infections"))
  samp_time <- seq(0, max_time, 1)
  sim <- sim_cohort(n = cohort_size,
                    samp_time = samp_time,
                    haplo_freqs = haplo_freqs,
                    lambda = i/max_time,
                    theta = theta,
                    decay_rate = decay_rate,
                    sens = sens,
                    n_inf = rep(i, cohort_size),
                    return_full = TRUE)
  
  assign(paste0("n_",i), sim)
  saveRDS(sim, paste0("simulated_data/cohort_n_infections/cohort_",i,"_infections.RDS"))
}

# plot_ind(filter_cohort(n_1, 5))

n_inf <- 1:5
sens_vec <- seq(0.8, 1, by = 0.05)
for(i in n_inf) {
  for(j in sens_vec) {  
    print(paste0("Simulate cohort with ", i, " infections and ", j, " sensitivity"))
    samp_time <- seq(0, max_time, 1)
    sim <- sim_cohort(n = cohort_size,
                      samp_time = samp_time,
                      haplo_freqs = haplo_freqs,
                      lambda = i/max_time,
                      theta = theta,
                      decay_rate = decay_rate,
                      sens = j,
                      n_inf = rep(i, cohort_size),
                      return_full = TRUE)
    
    # save so easily accessible in the local environment
    assign(paste0("n_",i,"_",j), sim)
    saveRDS(sim, paste0("simulated_data/sensitivity/cohort_",i,"_inf_",j,"_sens.RDS"))
  }
}

theta_vec <- c(seq(0.1, 1, 0.1), seq(1.5, 5, 0.5))
for(i in n_inf) {
  for(j in theta_vec) {  
    print(paste0("Simulate cohort with ", i, " infections and ", j, " COI parameter"))
    samp_time <- seq(0, max_time, 1)
    sim <- sim_cohort(n = cohort_size,
                      samp_time = samp_time,
                      haplo_freqs = haplo_freqs,
                      lambda = i/max_time,
                      theta = j,
                      decay_rate = decay_rate,
                      sens = sens,
                      n_inf = rep(i, cohort_size),
                      return_full = TRUE)
    
    # save so easily accessible in the local environment
    assign(paste0("n_",i,"_",j), sim)
    saveRDS(sim, paste0("simulated_data/coi/cohort_",i,"_inf_",j,"_coi.RDS"))
  }
}
