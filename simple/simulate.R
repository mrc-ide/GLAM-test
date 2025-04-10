# simple/simulate.R 
# Author: Gina Cuomo-Dannenburg
# Date: 10/04/2025
# Purpose: simulate a simple bank of datasets for testing

# set seed to ensure reproducibility
set.seed(2)

# source functions
source("functions.R")

# installation for testing of specified version of {GLAM} - SHA1 (c9cdc36f)
# devtools::install_github("mrc-ide/GLAM@test/gina_working_example")
library(GLAM)

# load libraries
library(ggplot2)
library(tidyverse)

# parameters involved: 
# samp_time; haplo_freqs; lambda; theta; decay; sens; n_inf; t_inf
# also need the repetition number 
# indexing becomes simulation number

# define the vector of options that we want to cycle over for the simulations
samp_time <- seq(0, 10, 1) # keep the same for all in this case
haplo_freqs <- rep(0.05, 20) # equal frequency of all haplotypes
theta_vec <- c(seq(0.1, 0.5, 0.1), seq(1, 5, 0.5))
decay_vec <- c(seq(0.1, 1, 0.1))
sens <- 1
n_inf <- 1:3
cohort_size <- 10
repetitions <- 10

# estimating necessary zero padding -- use 9 digits 
# 1*1*1*14*10*5*10*2*100

# generate a lookup table with all the parameter values
lookup <- 
  # start with cycling over the vectors
  tidyr::expand_grid(theta = theta_vec, 
                             decay = decay_vec, 
                             n = n_inf,
                             repetition = 1:repetitions)|>
  dplyr::mutate(lambda = n/max(samp_time),
                samp_time = list(samp_time),
                haplo_freqs = list(haplo_freqs),
                sens = sens,
                cohort_size = cohort_size) |>
  dplyr::select(cohort_size, samp_time, haplo_freqs, lambda,
                  theta, decay, sens, n, repetition) |>
  dplyr::mutate(sim_id = zero_pad_fixed(row_number(), 9)) # generate a simulation id number
saveRDS(lookup, "simple/lookup.RDS")

# generate the simulated datasets
for (i in 1:nrow(lookup)) {
  id <- lookup$sim_id[i]
  samp_time <- unlist(lookup$samp_time[i])
  n_inf <- lookup$n[i]
  
  # must be a better way of imputing t_inf with multiple options 
  # this is the case where they are equally spaced
  # t_inf <- seq(from = max(samp_time)/(n_inf+1),
  #              by = max(samp_time)/(n_inf+1),
  #              length.out = n_inf) - 0.3 # stagger to be sure that it isn't at same time as sampling
  # 
  # simulate the cohort corresponding to the parameter values
  sim <- GLAM::sim_cohort(
    n = lookup$cohort_size[i],
    samp_time = samp_time,
    haplo_freqs = unlist(lookup$haplo_freqs[i]),
    theta = lookup$theta[i],
    lambda = lookup$lambda[i],
    decay_rate = lookup$decay[i],
    sens = lookup$sens[i],
    n_inf = rep(n_inf, lookup$cohort_size[i])
    )
  saveRDS(sim, paste0("simple/data/sim", id, ".RDS"))
}



