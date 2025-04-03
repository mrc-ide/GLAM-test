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
#     * Are there issues with saturations with higher #s of infections?
# - Low and high transmission 
#     * Is there identifiability issues at high and low lambda?
# - Fixing individual level parameters and ability to fit global parameters
