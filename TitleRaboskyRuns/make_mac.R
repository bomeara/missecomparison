source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan


#config <- drake_config(plan)
#vis_drake_graph(config)

#future::plan(future::multicore)

#good_cluster_nodes = c(paste0("omearacluster",c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), ".nomad.utk.edu"), paste0("omearalab",c(7,8,9,11,18,22), ".nomad.utk.edu"), "omearashiny1.desktop.utk.edu")


#good_cluster_nodes = c(paste0("omearacluster",c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"), ".nomad.utk.edu"), paste0("omearalab",c(9,11,22), ".nomad.utk.edu"), "omearashiny1.desktop.utk.edu") #removing 7,8, 18, since they often went unreachable

#good_cluster_nodes = c("omearaclusterc.nomad.utk.edu", "omearaclusterf.nomad.utk.edu", "omearalab22.nomad.utk.edu") #not running condor at the moment

good_cluster_nodes = c( "omearaclustera.local", "omearaclusterb.local", "omearatc1.local", "omearatc2.local", "omearaclustere.local", "omearaclusterl.local", "omearaclusterg.local") #not running condor at the moment


all_nodes <- good_cluster_nodes
# all_nodes <- c()
# free_nodes <- 10
# for (i in seq_along(good_cluster_nodes)) {
#   if(good_cluster_nodes[i] %in% c("a", "b")) {
#     all_nodes <- append(all_nodes, rep(paste0("omearacluster", good_cluster_nodes[i] , ".nomad.utk.edu"), 1))
#   } else {
#     all_nodes <- append(all_nodes, rep(paste0("omearacluster", good_cluster_nodes[i] , ".nomad.utk.edu"), 1))
#   }
# }


#cl <- parallel::makeCluster(all_nodes, rscript="/usr/bin/Rscript")
#future::plan(cluster, workers = cl)
future::plan(future::multicore)


make(
  plan_mac, # defined in R/plan.R
  verbose = 2,
  #parallelism = "future", jobs = length(all_nodes)
  parallelism = "future", jobs = 6

)
