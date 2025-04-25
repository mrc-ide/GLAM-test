# functions.R
# 
# Author: Gina Cuomo-Dannenburg
# Date: 08/04/2025
# 
# Purpose: Functions to help with the testing and understanding of outputs of {glam}

# load libraries
library(tidyverse)
library(GLAM)
library(ggplot2)
library(ggpubr)

# Function to zero-pad integer ensuring fixed length
zero_pad_fixed <- function(num, places) {
  sprintf(paste0("%0", places, "d"), as.integer(num))
}


# # test functions
# # set simulation parameters
# samp_time <- seq(0, 10, 1)
# haplo_freqs <- rep(0.05, 20)
# lambda <- 0.2
# theta <- 2
# decay_rate <- 0.1
# sens <- 0.9
# max_infections <- 20
# cohort_size <- 10
# 
# # simulate cohort
# sim <- sim_cohort(n = cohort_size,
#                    samp_time = samp_time,
#                    haplo_freqs = haplo_freqs,
#                    lambda = lambda,
#                    theta = theta,
#                    decay_rate = decay_rate,
#                    sens = sens,
#                    n_inf = NULL,
#                    return_full = TRUE)

# ind1 <- filter_cohort(sim, 1)
# plot_ind(ind1)
# 
# ind2 <- filter_cohort(sim, 2)
# plot_ind(ind2)

# plot MCMC posterior and individual data in one panel
