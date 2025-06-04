# simple/summary.R 
# Author: Gina Cuomo-Dannenburg
# Date: 2025-05-28
# Purpose: Create summary statistics of the MCMC outputs
# Current thought: separate files for each of the MCMC outputs identified by the id

# This script will need a function that creates summary values for every simulated
# dataset and the associated mcmc
# start with a script and we will adapt to a function for cluster running

# ideas for inclusion here:
# - posterior 95% CrI for each parameter
# - posterior 95% CrI for each infection time
# - % of the time that the MCMC estimate each of these correctly
#   i.e. true val w/i 95% CrI
# - Sophie's suggestions - posterior predictive check and percentile-based residuals

library(GLAM)
library(tidyverse)
library(ggplot2)

lookup <- readRDS("simple/lookup.RDS")

# write a function that takes the row number i and outputs the posteriors

posterior_summary <- function(i) {
  id <- lookup$sim_id[i]
  n_ind <- lookup$cohort_size[i]
  
  sim <- readRDS(paste0("simple/data/sim", id, ".RDS"))
  mcmc <- readRDS(paste0("simple/mcmc_outputs/out", id, ".RDS"))
  
  # parameter names
  param_true <- lookup[i,] |>
    dplyr::select(lambda:n) |>
    dplyr::rename(n_inf = n) |>
    tidyr::pivot_longer(lambda:n_inf, names_to = "parameter_name", values_to = "true_val")
  
  # get the t_inf for all individuals
  t_inf <- mapply(function(x, i) {
    if (length(x$t_inf) == 0) {
      return(NULL)
    }
    data.frame(ind = i,
               t_inf = x$t_inf)
  }, sim$raw_list, 1:n_ind, SIMPLIFY = FALSE) |>
    bind_rows() |> 
    # add infection number
    group_by(ind) |> 
    dplyr::mutate(infection = row_number()) |>
    dplyr::rename(individual = ind,
                  true_val = t_inf) |>
    dplyr::select(individual, infection, true_val) |>
    dplyr::mutate(infection = factor(infection, levels = seq(1, 20, 1)))
    
  
  # find the 95% CrI for all of the paraemters
  params_pred <- mcmc$get_output_global() |>
    dplyr::filter(phase == "sampling") |> 
    dplyr::rename(decay = decay_rate) |>
    tidyr::pivot_longer(lambda:sens, names_to = "parameter_name") |>
    dplyr::group_by(parameter_name) |>
    dplyr::reframe(median = median(value),
                   lower_cri = quantile(value, 0.025),
                   upper_cri = quantile(value, 0.975))
  
  # this includes even those with such low probabilities that they are not super meaningful
  # ask BV how to proceed with these..
  t_inf_pred <- mcmc$get_output_infection_times() |> 
    dplyr::filter(phase == "sampling") |> 
    dplyr::group_by(individual, infection) |> 
    dplyr::reframe(median = median(value),
                   lower_cri = quantile(value, 0.025),
                   upper_cri = quantile(value, 0.975))
  
  n_inf_pred <- mcmc$get_output_n_infections() |> 
    dplyr::filter(phase == "sampling") |> 
    dplyr::group_by(ind) |> 
    dplyr::reframe(median = median(value),
                   lower_cri = quantile(value, 0.025),
                   upper_cri = quantile(value, 0.975))
  
  # return if the values are within the real values
  param_correct <- param_true |>
    dplyr::filter(parameter_name != "n_inf") |>
    dplyr::left_join(params_pred) |>
    dplyr::mutate(correct = if_else(true_val >= lower_cri & true_val <= upper_cri,
                                    TRUE, FALSE))
  n_inf <- param_true |>
    dplyr::filter(parameter_name == "n_inf") |>
    dplyr::pull(true_val)
  
  n_inf_correct <- n_inf_pred |> 
    dplyr::mutate(true_val = n_inf) |> 
    dplyr::mutate(correct = if_else(true_val >= lower_cri & true_val <= upper_cri,
                                    TRUE, FALSE))
  
  # t_inf is hard because sometimes there is a disconnect with the # of true inf and predicted
  
  
  
}


