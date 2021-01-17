source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan

drake::drake_cache("/Users/bomeara/Documents/MyDocuments/GitClones/missecomparison/TitleRaboskyRuns/.drake")$unlock()

workers <- c(rep(c("omearaclustera.local", "omearaclusterb.local", "omearaclusterl.local", "omearaclusterg.local"), 12), rep(c("omearatc1.local", "omearatc2.local"),6))

#workers <- c(rep(c("omearaclustera.local", "omearaclusterb.local", "omearaclusterl.local", "omearaclusterg.local"), 1), rep(c("omearatc1.local", "omearatc2.local"),1))
workers <- sample(workers, length(workers), replace=FALSE)
cl <- parallel::makeCluster(workers,  rscript="/usr/bin/Rscript", setup_strategy = "sequential")

future::plan(cluster, workers=cl)


make(plan_hpc, parallelism = "future", jobs = 60, verbose=4)

