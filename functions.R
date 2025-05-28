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

# adding functions from GLAM that aren't exported but are needed to run the code

#------------------------------------------------
# log-normal distribution parametrized in terms of mean and SD
#' @importFrom stats dlnorm
#' @noRd
dlnorm_reparam <- function(x, mean, sd, return_log = FALSE) {
  sigma_sq <- log(sd^2 / mean^2 + 1)
  dlnorm(x, meanlog = log(mean) - sigma_sq / 2, sdlog = sqrt(sigma_sq), log = return_log)
}

#------------------------------------------------
# gamma distribution parametrized in terms of mean and SD
#' @importFrom stats dgamma
#' @noRd
dgamma_reparam <- function(x, mean, sd, return_log = FALSE) {
  dgamma(x, shape = mean^2 / sd^2, rate = mean / sd^2, log = return_log)
}

#------------------------------------------------
# zero-truncated Poisson with rate parameter lambda (expectation
# lambda/(1-exp(-lambda)))
#' @importFrom stats qpois
#' @noRd
rztpois <- function(n, lambda) {
  if (lambda <= 0) {
    stop("lambda must be positive")
  }
  qpois(runif(n, exp(-lambda), 1), lambda)
}

#------------------------------------------------
# create nested list of phases in chains
#' @noRd
create_chain_phase_list <- function(chains, base) {
  lapply(1:chains, function(x) {
    list(
      tune = base,
      burn = base,
      sample = base
    )
  })
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
