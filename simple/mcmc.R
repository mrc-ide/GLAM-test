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

# Loops over every row of the lookup table and runs the MCMC and saves the output
# start <- Sys.time()
for (i in 1:nrow(lookup)) {
  set.seed(i)
  # for now, print out the progress so that I can see how far through we are
  if(i %% 10 == 0) {
    print(paste("MCMC on simulation set", i))
  }
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
# Sys.time() - start

