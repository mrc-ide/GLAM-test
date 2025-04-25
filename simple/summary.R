# simple/summary.R 
# Author: Gina Cuomo-Dannenburg
# Date: 2025-04-25
# Purpose: Create summary statistics of the MCMC outputs
# Current thought: separate files for each of the MCMC outputs identified by the id

# Code currently included is the code that was copied from mcmc.R but now moved over

# # summary outputs
# # have queried Bob for what outputs would be helpful 
# # current thoughts:
# # credible intervals of the posterior of each parameter
# # credible intervals of the posterior of each infection time
# # dataframe showing whether the true parameter value is within the CrI for all params matching the defined set
# # dataframe showing whether the true infection time is within the CrI 
# 
# # somewhat redundant at the moment b/c we fix the parameter values
# param_cri <- g$get_output_global() |>
#   dplyr::filter(phase == "sampling") |> 
#   tidyr::pivot_longer(lambda:sens, names_to = "parameter") |>
#   dplyr::select(parameter, value) |>
#   dplyr::group_by(parameter) |>
#   dplyr::reframe(lower_cri = quantile(value, 0.025), 
#                  upper_cri = quantile(value, 0.975),
#                  mean = mean(value),
#                  median = median(value)) |>
#   dplyr::mutate(repetition = lookup$repetition[which(lookup$sim_id == id)])
# saveRDS(param_cri, paste0("simple/summary/metrics/param_cri/cri", id, ".RDS"))
# 
# # infection time is static for some reason
# t_inf_cri <- g$get_output_infection_times() |> 
#   dplyr::filter(phase == "sampling") |> 
#   dplyr::group_by(individual, infection) |>
#   dplyr::reframe(lower_cri = quantile(value, 0.025), 
#                  upper_cri = quantile(value, 0.975))
# saveRDS(t_inf_cri, paste0("simple/summary/metrics/t_inf_cri/cri", id, ".RDS"))
# 
# # real times of infection
# df_inf_time_true <- mapply(function(i) {
#   t <- sim$raw_list[[i]]$t_inf
#   if (length(t) == 0) {
#     return(NULL)
#   }
#   data.frame(individual = i,
#              infection = seq_along(t),
#              time = t)
# }, seq_along(sim$raw_list), SIMPLIFY = FALSE) |>
#   bind_rows()
# 
