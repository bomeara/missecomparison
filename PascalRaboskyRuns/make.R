source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan


config <- drake_config(plan)
vis_drake_graph(config)

#future::plan(future::multiprocess)
future::plan(cluster, workers = c("omearaclustera.nomad.utk.edu", "localhost"))



make(
  plan, # defined in R/plan.R
  verbose = 2,
  parallelism = "future", jobs = 2
)
