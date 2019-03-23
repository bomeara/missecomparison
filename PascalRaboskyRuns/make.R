source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan


#config <- drake_config(plan)
#vis_drake_graph(config)

#future::plan(future::multiprocess)

good_cluster_nodes = c("a", "b", "c", "e", "i", "h", "l")
all_nodes <- c()
for (i in seq_along(good_cluster_nodes)) {
  if(good_cluster_nodes[i] %in% c("a", "b")) {
    all_nodes <- append(all_nodes, rep(paste0("omearacluster", good_cluster_nodes[i] , ".nomad.utk.edu"), 16))
  } else {
    all_nodes <- append(all_nodes, rep(paste0("omearacluster", good_cluster_nodes[i] , ".nomad.utk.edu"), 24))

  }
}
future::plan(cluster, workers = all_nodes)



make(
  plan, # defined in R/plan.R
  verbose = 2,
  parallelism = "future", jobs = length(all_nodess)
)
