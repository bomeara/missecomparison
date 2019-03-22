source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan
make(
  plan, # defined in R/plan.R
  verbose = 2
)
