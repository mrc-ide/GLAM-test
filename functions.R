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

# filter a cohort to look at one individual and enable plotting
filter_cohort <- function(cohort_sim, ind) {
  return(cohort_sim$raw_list[[ind]])
}

# plot an individual without information on when a bite was and without the 
# information on what haplotypes were introduced 
# seeing the individuala as the MCMC would

mcmc_sees <- function (ind) {
  start <- haplo <- end <- time <- NULL
  # this throws an error despite it being in plot_ind and that rendering?!
  # assert_in(c("state_obs", "t_inf", "w_init", "clear_init", 
  #             "w_inf", "clear_inf", "state_true", "samp_time"), names(ind))
  state_obs <- ind$state_obs
  t_inf <- ind$t_inf
  w_init <- ind$w_init
  clear_init <- ind$clear_init
  w_inf <- ind$w_inf
  clear_inf <- ind$clear_inf
  state_true <- ind$state_true
  samp_time <- ind$samp_time
  n_haplo <- length(w_init)
  n_samp <- length(samp_time)
  start_time <- samp_time[1]
  end_time <- samp_time[n_samp]
  n_inf <- length(t_inf)
  df_plot <- data.frame(start = rep(start_time, sum(w_init)), 
                        end = clear_init[w_init == 1], haplo = which(w_init == 
                                                                       1))
  if (n_inf > 0) {
    df_plot <- bind_rows(df_plot, data.frame(start = (w_inf * 
                                                        t_inf)[w_inf == 1], end = clear_inf[w_inf == 1], 
                                             haplo = which(w_inf == 1, arr.ind = TRUE)[, 2]))
  }

  df_plot <- mutate(data.frame(time = samp_time, 
                               haplo = rep(1:n_haplo, each = n_samp), 
                               state_true = c("Absent", "Present")[state_true + 1], 
                               state_obs = c("Unobserved", "Observed")[state_obs + 1]), 
                    state_true = factor(state_true, levels = c("Absent",  "Present")), 
                    state_obs = factor(state_obs, levels = c("Unobserved", "Observed"))) |>
    dplyr::mutate(positive = if_else(state_true == "Present" & state_obs == "Observed", "Positive", "Negative"),
                  positive = factor(positive, levels = c("Negative","Positive")))
  ggplot(df_plot) + theme_bw() + 
    geom_point(aes(x = time, y = haplo, alpha = positive, col = positive), size = 3, data = df_plot) + 
    scale_color_discrete(name = "Test result") + 
    scale_alpha_manual(values = c(0.1, 1), name = "Test result") +
    xlab("Time") + 
    ylab("Haplotype") + 
    xlim(c(start_time, end_time)) + 
    ylim(c(1, n_haplo)) + 
    theme(panel.grid = element_blank())
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
