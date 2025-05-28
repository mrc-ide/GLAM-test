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
