# simulate_data.R
# 
# Author: Gina Cuomo-Dannenburg
# Date: 03/04/2025
# 
# Purpose: simulate datasets for testing of the {mrc-ide/GLAM} package
# Script is intended to produce core examples of datasets, from basic to more
# complex to test the full bredth of functionality of {GLAM}, to validate methods etc.
# 
# Initial thoughts:
# - Cohorts with the same # total infections
#     * Does the MCMC accurately fit the timing of infections?
#     * Are there issues with saturation at higher #s of infections?
# - Low and high transmission 
#     * Is there and issue identifying specific parameters at high and low lambda?
# - Fixing individual level parameters and ability to fit global parameters
# - Inconsistent sampling time
# - Scalability with # of haplotypes
# - Frequency of specific haplotypes
# - Collinearity between specific parameters where they are likely to create a 
#   gradient given the likelihood and other factors? (~ shape and scale in some 
#   distributions)

