# build_intuition.R
# 
# Author: Gina Cuomo-Dannenburg
# Date: 03/04/2025
# 
# Purpose: script to try and build intuition for how different parameters influence
# the simulated data and think about potential challenges and limitations of methodologies

# load libraries
library(GLAM)
library(tidyverse)
library(ggplot2)

# default parameters
# set simulation parameters
samp_time <- seq(0, 10, 1)
haplo_freqs <- rep(0.05, 20)
lambda <- 0.2
theta <- 2
decay_rate <- 0.1
sens <- 0.9
max_infections <- 20

set.seed(1234)

# explore different sampling times
# larger spacing
samp_time <- seq(0, 10, 2)

sim1 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 5)
plot_ind(sim1)

samp_time <- seq(0, 10, 3)

sim2 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 5)
plot_ind(sim2)

samp_time <- c(0, 1, 3, 5, 7, 10)


sim3 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 5)
plot_ind(sim3)

# return to default
samp_time <- seq(0, 10, 1)

# explore asymmetric haplotype frequencies
denom <- sum(1:20)
haplo_freqs <- c(1:20/denom)

sim1 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 5)
plot_ind(sim1)

haplo_freqs <- c(rep(0.01, 3), rep((1-0.03)/17,17))


sim2 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 5)
plot_ind(sim2)

haplo_freqs <- c(rep(0.3, 2), rep((1-0.6)/18,18))
sim3 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 5)
plot_ind(sim3)

sim4 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 2)
plot_ind(sim4)

sim5 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 12)
plot_ind(sim5)

# higher frequency mutations are just reintroduced many times when the number of 
# n_inf is higher -- hard to distinguish what is high frequency and what is repeated
# reinfections of the same haplotypes
# e.g. in sim5, the two haplotypes with 0.3 frequency are introduced at nearly all infections
# and therefore detected ant most time points.
# haplotype 13 was only introduced once, but remained for all samplings (although not 
# always observed)

# revert to default
haplo_freqs <- rep(0.05, 20)

# explore lambda ~ FOI -- set n_inf to NULL to reflect lambda
lambda <- 0.01

sim1 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = NULL)
plot_ind(sim1)

lambda <- 0.05

sim2 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = NULL)
plot_ind(sim2)

lambda <- 0.1

sim3 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = NULL)
plot_ind(sim3)

lambda <- 0.2

sim4 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = NULL)
plot_ind(sim4)

lambda <- 0.5

sim5 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = NULL)
plot_ind(sim5)

# examine the relationship between n_inf and lambda 
lambda_vec <- c(0.001, 0.005, seq(0.01, 0.99, by = 0.01))

lambda_inf <- data.frame(lambda = rep(lambda_vec, each = 100),
                         individual = rep(1:100, times = length(lambda_vec)),
                         n_inf = NA)

# increase number of sampling times to smooth the output
samp_time <- seq(0, 100, 1)
for(i in 1:nrow(lambda_inf)) {
  sim <- sim_ind(samp_time = samp_time,
                 haplo_freqs = haplo_freqs,
                 lambda = lambda_inf$lambda[i],
                 theta = theta, 
                 decay_rate = decay_rate,
                 sens = sens,
                 return_full = TRUE,
                 n_inf = NULL)
  lambda_inf$n_inf[i] <- length(sim$t_inf)
}

lambda_summary <- lambda_inf |>
  dplyr::group_by(lambda) |>
  dplyr::reframe(median = median(n_inf),
                 lower = quantile(n_inf, 0.025),
                 upper = quantile(n_inf, 0.975))

ggplot(lambda_summary, aes(x = lambda, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) + theme_bw()

max(lambda_inf$n_inf)
# maximum # infections > max time at high FOI

# fit a linear model
lm_lambda <- lm(median ~ lambda, data = lambda_summary)
summary(lm_lambda)
lm_upper <- lm(upper ~ lambda, data = lambda_summary)
lm_lower <- lm(lower ~ lambda, data = lambda_summary)

# interesting -- could we use relationship between n_inf ~ lambda to help
# with identifiability issues?
ggplot(lambda_summary, aes(x = lambda, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) + theme_bw() +
  geom_abline(slope = lm_lambda$coefficients[2],
              intercept = lm_lambda$coefficients[1],
              lwd = 1, lty = 4, col = "blue") + 
  geom_abline(slope = lm_upper$coefficients[2],
              intercept = lm_upper$coefficients[1],
              lwd = 1, lty = 4, col = "blue") + 
  geom_abline(slope = lm_lower$coefficients[2],
              intercept = lm_lower$coefficients[1],
              lwd = 1, lty = 4, col = "blue")

# return to default values
lambda <- 0.2
samp_time <- seq(0, 10, 1)
theta <- 2
# theta = COI intensity parameter 

sim1 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 3)
plot_ind(sim1)

theta <- 0.5
sim2 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 3)
plot_ind(sim2)


theta <- 20
sim3 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = 3)
plot_ind(sim3)

# higher theta, higher likelihood of multiple haplotypes being introduced simultaneously
# a little confused because this doesn't actually change COI because we have the same number
# of simultaneous infections

theta <- 2

# explore decay rate 
decay_rate <- 0.01
sim1 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = NULL)
plot_ind(sim1)

decay_rate <- 0.1
sim2 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = NULL)
plot_ind(sim2)

decay_rate <- 0.5
sim3 <- sim_ind(samp_time = samp_time,
                haplo_freqs = haplo_freqs,
                lambda = lambda,
                theta = theta, 
                decay_rate = decay_rate,
                sens = sens,
                return_full = TRUE,
                n_inf = NULL)
plot_ind(sim3)

# a little confused by the relationship between theta and decay rate
# very low decay ~ high theta --> I guess that makes sense but discuss with Bob
