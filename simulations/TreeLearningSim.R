source("Tree_functions.R")
library(graphicalExtremes)
library(igraph)
library(matrixcalc)
library(mvtnorm)
library(tidyverse)
library(rngtools)
library(doFuture)
library(here)
library(jsonlite)

# change from here ----
sim_setting <- "sim_study_1"
strategy <-  "sequential"  # "sequential" or "parallel"
n_workers <- 96
sim_function <- "sim_study"
# to here ----

# import json
json <- fromJSON(here("/simulations/config", paste0(sim_setting, ".json")))

# process parameters
param <- json$params %>% 
  expand.grid(stringsAsFactors = FALSE) %>% 
  as_tibble()

if (!has_name(param, "p")){
  param <- param %>% 
    mutate(k = floor(n ** .8),
           p = 1 - k / n)
}

# join with simulation repetitions
sims_args <- rep_tibble(param, json$other_params$nsim) %>% 
  rename(nsim = rep_id) %>% 
  assign_random_seed(grouping_vars = json$other_params$random_seed_for_each,
                     seed = json$other_params$seed) %>% 
  ungroup() %>% 
  rownames_to_column(var = "rowid")



# split arguments between fun_args and other_sargs
other_args_nms <- names(sims_args)[
  !(names(sims_args) %in% names(formals(eval(as.name(sim_function)))))]
fun_args <- sims_args %>% select(!any_of(other_args_nms))
other_args <- sims_args %>% select(any_of(other_args_nms))
m <- nrow(sims_args)


# set up file names
dttime <- gsub(pattern = " |:", x = Sys.time(), replacement = "_")
file_log <- here("simulations/output", "progress.txt")
file_rds <- here("simulations/output", paste0(sim_setting,  "-", dttime, ".rds"))


# set up cluster
registerDoFuture()
if(strategy == "parallel"){
  
  plan(multisession, workers = n_workers)
  
} else if (strategy == "sequential"){
  
  plan(sequential)
  
} else if (strategy == "mixed"){
  
  strategy_1 <- tweak(sequential)
  strategy_2 <- tweak(multisession, workers = n_workers)
  plan(list(strategy_1, strategy_2))
  
} else {
  
  stop ("wrong parallelization strategy")
  
}


# iterate over simulations
ptm <- proc.time()
cat("**** Simulation ---", sim_setting , "**** \n", file = file_log)
ll <- foreach(i = 1:m, .combine = bind_rows, 
              .options.future = list(scheduling = FALSE),
              .errorhandling = "remove") %dopar% {
                
                # Run simulation
                cat("Simulation", i, "out of", m, "\n", file = file_log, append = TRUE)
                wrapper_sim(i, other_args$rowid[i], eval(as.name(sim_function)), fun_args)
                
              }
sink(file = file_log, append = TRUE)
cat("\n Time \n")
print(proc.time() - ptm)
sink()


# collect results
ll <- ll %>% 
  nest(perf = c(type, value)) %>% 
  left_join(sims_args, by = "rowid")


# Error log
rowid_errors <- which(!(sims_args$rowid %in% ll$rowid))
if (length(rowid_errors) > 0){
  sink(file = file_log, append = TRUE)
  cat("\nError handling \n")
  cat(paste0("Error occured in iteration with rowid: ", rowid_errors, "\n"))
  sink()
}

# send email

# save results
saveRDS(ll, file = file_rds)

