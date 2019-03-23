source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan


config <- drake_config(plan)
vis_drake_graph(config)

#future::plan(future::multiprocess)

good_clusters = c("a", "b", "c", "e", "i", "h", "l")
future::plan(cluster, workers = c(paste0("omearacluster", good_clusters , ".nomad.utk.edu")))



make(
  plan, # defined in R/plan.R
  verbose = 2,
  parallelism = "future", jobs = length(good_clusters)
)
