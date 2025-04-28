library(hipercow)
hipercow_init()
hipercow_configure(driver = "dide-windows")
hipercow_provision(method = "pkgdepends")
hipercow_environment_create(sources = c("simple/simulate_function.R", 
                                        "simple/mcmc_function.R"))
hipercow_resources(cores = 1)

lookup <- readRDS("simple/lookup.RDS")
# 
# # try and debug with running one task outside of a batch call
# simulate_task <- task_create_expr(simulate_data(2))
# task_log_show(simulate_task)
# 
# mcmc_task <- task_create_expr(mcmc_fitting(2))
# task_log_show(mcmc_task)

rows <- 1:100 #nrow(lookup)
bundle <- task_create_bulk_call(fn = simulate_data, 
                                data = rows, 
                                bundle_name = "simulate-data")
bundle <- task_create_bulk_call(fn = mcmc_fitting, 
                                data = rows, 
                                bundle_name = "fit-simulated-data")