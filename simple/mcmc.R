# simple/mcmc.R 
# Author: Gina Cuomo-Dannenburg
# Date: 23/04/2025
# Purpose: simulate a simple bank of datasets for testing

# Notes: this script should run an MCMC on every simulation based on the lookup table
# and save the results of the MCMC to another folder (namely mcmc_outputs) as .RDS


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
sampling_iterations <- 1e3
num_chains <- 1
num_rungs <- 1

lookup <- readRDS("simple/lookup.RDS")

# I think this still needs a lot of work. Check with Bob if this is even the approach we want
# I was imagining having this as a way of checking all reps with the same parameter set 
# And extracting some information about the sensitiviy of the MCMC in detecting the true values
# Outputs are organised in subfolders with the same id ref so it's easy to pull in everything together

for (i in 1:nrow(lookup)) {
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
  
  
}


