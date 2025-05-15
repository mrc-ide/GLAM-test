library(hipercow)
# library(rrq) # enable workers

# hipercow set up
hipercow_init()
hipercow_configure(driver = "dide-windows")
hipercow_provision(method = "pkgdepends")
hipercow_environment_create(sources = c("simple/simulate_function.R", 
                                        "simple/mcmc_function.R",
                                        "null/simulate_function.R", 
                                        "null/mcmc_function.R"))
hipercow_resources(cores = 15)

# running sims and MCMC with same # of n_inf
lookup <- readRDS("simple/lookup.RDS")

rows <- 1:nrow(lookup)

# try and send them in batches based on the functions in the scripts
block_size <- 50
num_blocks <- nrow(lookup)/block_size
blocks <- 1:num_blocks
bundle <- task_create_bulk_call(fn = simulate_block, 
                                data = blocks,
                                args = list(block_size = 50),
                                bundle_name = "simulate-block")
bundle <- task_create_bulk_call(fn = mcmc_fitting_block, 
                                data = blocks, 
                                args = list(block_size = 50),
                                bundle_name = "fit-block")

# running sims and MCMC with same # of n_inf
lookup_null <- readRDS("null/lookup_null.RDS")

rows <- 1:nrow(lookup_null)

# try and send them in batches based on the functions in the scripts
block_size <- 30
num_blocks <- nrow(lookup_null)/block_size
blocks <- 1:num_blocks
bundle <- task_create_bulk_call(fn = simulate_block_null, 
                                data = blocks,
                                args = list(block_size = block_size),
                                bundle_name = "simulate-null")
bundle <- task_create_bulk_call(fn = mcmc_fitting_block_null, 
                                data = blocks, 
                                args = list(block_size = block_size),
                                bundle_name = "fit-block-null")
