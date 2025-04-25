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
sens <- 1
n_inf <- 2
lambda_val <- n_inf/max(samp_time) # set to n_inf/max(samp_time)
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

if(length(sim_ids) == 0) {
  stop("Parameter value set does not exist in the lookup table")
}

# create a stop message for if the parameter value set exists but the simulation doesn't

# I think this still needs a lot of work. Check with Bob if this is even the approach we want
# I was imagining having this as a way of checking all reps with the same parameter set 
# And extracting some information about the sensitiviy of the MCMC in detecting the true values
# Outputs are organised in subfolders with the same id ref so it's easy to pull in everything together

for (i in 1:length(sim_ids)) {
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
         lambda = lambda_val,
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
  # have queried Bob for what outputs would be helpful 
  # current thoughts:
  # credible intervals of the posterior of each parameter
  # credible intervals of the posterior of each infection time
  # dataframe showing whether the true parameter value is within the CrI for all params matching the defined set
  # dataframe showing whether the true infection time is within the CrI 
  
  # somewhat redundant at the moment b/c we fix the parameter values
  param_cri <- g$get_output_global() |>
    dplyr::filter(phase == "sampling") |> 
    tidyr::pivot_longer(lambda:sens, names_to = "parameter") |>
    dplyr::select(parameter, value) |>
    dplyr::group_by(parameter) |>
    dplyr::reframe(lower_cri = quantile(value, 0.025), 
                   upper_cri = quantile(value, 0.975),
                   mean = mean(value),
                   median = median(value)) |>
    dplyr::mutate(repetition = lookup$repetition[which(lookup$sim_id == id)])
  saveRDS(param_cri, paste0("simple/summary/metrics/param_cri/cri", id, ".RDS"))
  
  # infection time is static for some reason
  t_inf_cri <- g$get_output_infection_times() |> 
    dplyr::filter(phase == "sampling") |> 
    dplyr::group_by(individual, infection) |>
    dplyr::reframe(lower_cri = quantile(value, 0.025), 
                   upper_cri = quantile(value, 0.975))
  saveRDS(t_inf_cri, paste0("simple/summary/metrics/t_inf_cri/cri", id, ".RDS"))
  
  # real times of infection
  df_inf_time_true <- mapply(function(i) {
    t <- sim$raw_list[[i]]$t_inf
    if (length(t) == 0) {
      return(NULL)
    }
    data.frame(individual = i,
               infection = seq_along(t),
               time = t)
  }, seq_along(sim$raw_list), SIMPLIFY = FALSE) |>
    bind_rows()
  
  
}


