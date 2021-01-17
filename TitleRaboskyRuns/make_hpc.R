source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan



workers <- c(rep(c("omearaclustera.local", "omearaclusterb.local", "omearaclusterl.local", "omearaclusterg.local"), 23), rep(c("omearatc1.local", "omearatc2.local"),11))
workers <- sample(workers, length(workers), replace=FALSE)
cl <- parallel::makeCluster(workers)

future::plan(cluster, workers=cl)


make(plan_hpc, parallelism = "future", jobs = 114, verbose=4)

