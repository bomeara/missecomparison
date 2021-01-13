source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan



workers <- c(rep(c("omearaclustera.local", "omearaclusterb.local", "omearaclusterl.local", "omearaclusterg.local"), 24), rep(c("omearatc1.local", "omearatc2.local"),12))
cl <- makeClusterPSOCK(workers)

future::plan(cluster, workers=cl)


make(plan_hpc, parallelism = "future", jobs = 120)

