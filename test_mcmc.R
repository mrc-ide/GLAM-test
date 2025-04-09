# test_mcmc.R
# 
# Author: Gina Cuomo-Dannenburg
# Date: 07/04/2025
# 
# Purpose: test glam MCMC methods on simulated datasets

# test how well the MCMC works on cohorts with the same n_inf
source("functions.R")

# set parameters of the MCMC process
burnin_iterations <- 1e3
sampling_iterations <- 5e3
num_chains <- 1
num_rungs <- 1

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
         chains = num_chains, rungs = num_rungs)
  
  # run burn-in and sampling
  g$burn(iterations = burnin_iterations) |>
    system.time()
  
  g$sample(iterations = sampling_iterations)  |>
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
           chains = num_chains, rungs = num_rungs)
    
    # run burn-in and sampling
    g$burn(iterations = burnin_iterations) |>
      system.time()
    
    g$sample(iterations = sampling_iterations)  |>
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

# impact of sensitivity
n_inf <- 1:5
theta_vec <- c(seq(0.1, 1, 0.1), seq(1.5, 5, 0.5))
for(i in n_inf) {
  for(j in theta_vec) {
    print(paste("Performing MCMC on cohort with ", i, " infections and a COI parameter of ", j))
    sim <- readRDS(paste0("simulated_data/coi/cohort_",i,"_inf_",j,"_coi.RDS"))
    
    # so it's easy to access cohorts afterwards 
    # come up with more robust naming convention
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
           chains = num_chains, rungs = num_rungs)
    
    # run burn-in and sampling
    g$burn(iterations = burnin_iterations) |>
      system.time()
    
    g$sample(iterations = sampling_iterations)  |>
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
    ggsave(paste0("figures/coi/infection_chain_", i,"_",j,".png"), dpi = 300)
    
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
    
    ggsave(paste0("figures/coi/infection_time_", i,"_",j,".png"), dpi = 300)
  }
}


### ----  EXPLORE WHEN MCMC SUCCEEDS AND FAILS ---- OMIT WHEN GLAM WORKS ---- ###

# look at specific output
# systematically look at individuals where the MCMC succeeds and fails to see if there is somethings connecting them
# successes
sim_1_0.1 <- readRDS("~/code/GLAM-test/simulated_data/coi/cohort_1_inf_0.1_coi.RDS")
mcmc_sees(filter_cohort(sim_1_0.1, 7)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_1_0.1, 7)) + theme(legend.position = "bottom")
mcmc_sees(filter_cohort(sim_1_0.1, 1)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_1_0.1, 1)) + theme(legend.position = "bottom")

# fails
# infected between 6 and 7 -- posteriof implies very early on in sampling
mcmc_sees(filter_cohort(sim_1_0.1, 4)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_1_0.1, 4)) + theme(legend.position = "bottom")

# identifies an infection well after the haplotype infected has decayed
mcmc_sees(filter_cohort(sim_1_0.1, 8)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_1_0.1, 8)) + theme(legend.position = "bottom")

# here, the posterior distribution is > 0 above true_inf and == 0 below 
mcmc_sees(filter_cohort(sim_1_0.1, 10)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_1_0.1, 10)) + theme(legend.position = "bottom")

# now with more infections
sim_2_0.1 <- readRDS("~/code/GLAM-test/simulated_data/coi/cohort_2_inf_0.1_coi.RDS")

# individual 7 succeeds for infection 1 and fails for 2 
# honestly, totally sensible here!
mcmc_sees(filter_cohort(sim_2_0.1, 7)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_0.1, 7)) + theme(legend.position = "bottom")

# again, also looks sensible here
mcmc_sees(filter_cohort(sim_2_0.1, 1)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_0.1, 1)) + theme(legend.position = "bottom")

# mistake #2 as #1 
# eyeballing I think it's fairly obvious
# is there something wrong with the likelihood that the absense of a positive for two sampling and then a positive
# is more likely than being infected later? (posterior: inf 1 ~ 0.2, inf 2 ~ 1.5 (actually when # 1 happened))
mcmc_sees(filter_cohort(sim_2_0.1, 6)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_0.1, 6)) + theme(legend.position = "bottom")

# flattish posterior for infection #1
# again, pretty obvious when the first infection is 
mcmc_sees(filter_cohort(sim_2_0.1, 2)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_0.1, 2)) + theme(legend.position = "bottom")

# when both succeed in fitting
sim_2_0.3 <- readRDS("~/code/GLAM-test/simulated_data/coi/cohort_2_inf_0.3_coi.RDS")
mcmc_sees(filter_cohort(sim_2_0.3, 1)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_0.3, 1)) + theme(legend.position = "bottom")

sim_2_0.5 <- readRDS("~/code/GLAM-test/simulated_data/coi/cohort_2_inf_0.5_coi.RDS")
mcmc_sees(filter_cohort(sim_2_0.5, 1)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_0.5, 1)) + theme(legend.position = "bottom")

sim_2_3.5 <- readRDS("~/code/GLAM-test/simulated_data/coi/cohort_2_inf_3.5_coi.RDS")
mcmc_sees(filter_cohort(sim_2_3.5, 8)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_3.5, 8)) + theme(legend.position = "bottom")

##  reflections: how much of being able to "eyeball" when the timing should be is due to 
#   knowing the number of infections --> this is also inherently (in these simulations) related to 
#   the value of lambda i.e. transmission intensity. 
#   i.e. in this last example, you could estimate potentially 5/6 distinct infections based on what 
#   the MCMC 'sees'
# this kind of notion of what is more likely -- a decay and new infection or absence due to sensitivity

# posterior of both infections does not cover the true timing
# challenges of multiple haplotypes getting false negatives at the same sampling
# I guess this is more likely than only 1? due to issues with amplification perhaps?
mcmc_sees(filter_cohort(sim_2_3.5, 5)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_3.5, 5)) + theme(legend.position = "bottom")

# what about when sensitivity is perfect
sim_2_perfect <- readRDS("~/code/GLAM-test/simulated_data/sensitivity/cohort_2_inf_1_sens.RDS")
mcmc_sees(filter_cohort(sim_2_perfect, 8)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_perfect, 8)) + theme(legend.position = "bottom")

# I don't know what caused it to struggle here 
mcmc_sees(filter_cohort(sim_2_perfect, 9)) + theme(legend.position = "bottom")
plot_ind(filter_cohort(sim_2_perfect, 9)) + theme(legend.position = "bottom")
