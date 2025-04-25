library(hipercow)
hipercow_init()
hipercow_configure(driver = "windows")
hipercow_provision(method = "pkgdepends")
hipercow_environment_create(sources = c("simple/simulate_function.R", 
                                        "mcmc_function.R"))
hipercow_resources(cores = 20)

lookup <- readRDS("simple/lookup.RDS")
rows <- 1:nrow(lookup)

bundle <- task_create_bulk_call(simulate_data, rows, 
                                bundle_name = "simulate-data")
bundle <- task_create_bulk_call(mcmc_fitting, rows, 
                                bundle_name = "fit-simulated-data")
