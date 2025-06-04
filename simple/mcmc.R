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

# set priors for all of the mcmc fittings
# define priors
lambda_prior <- function(x) dlnorm_reparam(x, mean = 0.5, sd = 0.2, return_log = TRUE)
theta_prior <- function(x) dlnorm_reparam(x, mean = 3, sd = 3, return_log = TRUE)
decay_rate_prior <- function(x) dlnorm_reparam(x, mean = 0.5, sd = 0.2, return_log = TRUE)
sens_prior <- function(x) dbeta(x, shape1 = 99, shape2 = 1, log = TRUE)

# fixed parameters for all models
haplo_freqs <- rep(0.05, 20) # equal frequency of all haplotypes
max_infections <- 20

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
  g <- GLAM::glam_mcmc$new(df_data = sim$df_data)
  
  # initialize. If parameters are set to NULL, they are estimated. If they are
  # given a value, they take this fixed value. Useful for giving the inference
  # method the correct true value of some parameters, to diagnose how well it
  # estimates others.
  # copied from GLAM/R_scripts/deploy.R
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
# Sys.time() - start

