# test_mcmc.R
# 
# Author: Gina Cuomo-Dannenburg
# Date: 07/04/2025
# 
# Purpose: test glam MCMC methods on simulated datasets

# test how well the MCMC works on cohorts with the same n_inf
source("functions.R")

n_inf <- seq(1:20)

for(i in n_inf) {
  print(paste("Performing MCMC on cohort with ", i, " infections"))
  sim <- readRDS(paste0("simulated_data/cohort_n_infections/cohort_",i,"_infections.RDS"))
  
  # so it's easy to access cohorts afterwards
  assign(paste0("n_",i), sim)
  
  # create a new MCMC object and load data
  g <- glam_mcmc$new(df_data = sim$df_data)
  
  # initialize. If parameters are set to NULL, they are estimated. If they are
  # given a value, they take this fixed value. Useful for giving the inference
  # method the correct true value of some parameters, to diagnose how well it
  # estimates others.
  g$init(start_time = 0,
         end_time = 10,
         haplo_freqs = haplo_freqs,
         lambda = lambda,
         theta = theta,
         decay_rate = decay_rate,
         sens = sens,
         n_infections = mapply(function(x) length(x$t_inf), sim$raw_list),
         infection_times = NULL,
         max_infections = 20,
         chains = 1, rungs = 1)
  
  # run burn-in and sampling
  g$burn(iterations = 1e3) |>
    system.time()
  
  g$sample(iterations = 5e3) |>
    system.time()
  
  # -----------------------
  # EXTRACT OUTPUTS AND PLOT RESULTS
  
  # extract outputs into data.frames
  df_global <- g$get_output_global()
  df_n_infections <- g$get_output_n_infections()
  df_infection_times <- g$get_output_infection_times()
  
  # dataframe of the times of infections in the simulated cohort
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
  
  # trace plot of a single infection time
  df_infection_times |>
    filter(infection == 1) |>
    ggplot() + theme_bw() +
    geom_point(aes(x = iteration, y = value, col = chain)) +
    ylim(c(0, max_infections)) + facet_wrap(~individual)
  ggsave(paste0("figures/cohort_mcmc/infection_chain_", i,".png"), dpi = 300)
  
  # density plot of all infection times
  df_infection_times |>
    filter(phase == "sampling") |>
    group_by(individual, infection) |>
    do({
      d <- density(.$value)
      data.frame(x = d$x, y = d$y)
    }) |>
    ungroup() |>
    ggplot(aes(x = x, y = y, color = infection)) +
    geom_vline(aes(xintercept = time), col = grey(0.5), 
               data = df_inf_time_true) +
    geom_line() +
    facet_wrap(~ individual) +
    labs(x = "Value", y = "Density", color = "Infection") +
    theme_bw() +
    theme(panel.grid.major = element_blank())
  
  ggsave(paste0("figures/cohort_mcmc/infection_time_", i,".png"), dpi = 300)
}

# # look at specific individuals to try and see why the MCMC is fitting the incorrect time of infection
# plot_ind(filter_cohort(sim, 7))

# impact of sensitivity
n_inf <- 1:5
sens_vec <- seq(0.8, 1, by = 0.05)
for(i in n_inf) {
  for(j in sens_vec) {
    print(paste("Performing MCMC on cohort with ", i, " infections and a sensitivity of ", j))
    sim <- readRDS(paste0("simulated_data/sensitivity/cohort_",i,"_inf_",j,"_sens.RDS"))
    
    # so it's easy to access cohorts afterwards
    assign(paste0("n_",i,"_",j), sim)
    
    # create a new MCMC object and load data
    g <- glam_mcmc$new(df_data = sim$df_data)
    
    # initialize. If parameters are set to NULL, they are estimated. If they are
    # given a value, they take this fixed value. Useful for giving the inference
    # method the correct true value of some parameters, to diagnose how well it
    # estimates others.
    g$init(start_time = 0,
           end_time = 10,
           haplo_freqs = haplo_freqs,
           lambda = lambda,
           theta = theta,
           decay_rate = decay_rate,
           sens = sens,
           n_infections = mapply(function(x) length(x$t_inf), sim$raw_list),
           infection_times = NULL,
           max_infections = 20,
           chains = 1, rungs = 1)
    
    # run burn-in and sampling
    g$burn(iterations = 1e3) |>
      system.time()
    
    g$sample(iterations = 5e3) |>
      system.time()
    
    # -----------------------
    # EXTRACT OUTPUTS AND PLOT RESULTS
    
    # extract outputs into data.frames
    df_global <- g$get_output_global()
    df_n_infections <- g$get_output_n_infections()
    df_infection_times <- g$get_output_infection_times()
    
    # dataframe of the times of infections in the simulated cohort
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
    
    # trace plot of a single infection time
    df_infection_times |>
      filter(infection == 1) |>
      ggplot() + theme_bw() +
      geom_point(aes(x = iteration, y = value, col = chain)) +
      ylim(c(0, max_infections)) + facet_wrap(~individual)
    ggsave(paste0("figures/sensitivity/infection_chain_", i,"_",j,".png"), dpi = 300)
    
    # density plot of all infection times
    df_infection_times |>
      filter(phase == "sampling") |>
      group_by(individual, infection) |>
      do({
        d <- density(.$value)
        data.frame(x = d$x, y = d$y)
      }) |>
      ungroup() |>
      ggplot(aes(x = x, y = y, color = infection)) +
      geom_vline(aes(xintercept = time), col = grey(0.5), 
                 data = df_inf_time_true) +
      geom_line() +
      facet_wrap(~ individual) +
      labs(x = "Value", y = "Density", color = "Infection") +
      theme_bw() +
      theme(panel.grid.major = element_blank())
    
    ggsave(paste0("figures/sensitivity/infection_time_", i,"_",j,".png"), dpi = 300)
  }
}


