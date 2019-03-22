source("R/packages.R")  # loads packages
source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan


tree_info <- read.csv(file=file_in("data/pascal_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
tree_names <- unique(tree_info$treeName)[1:3] #limit to 1:3 for trial


evaluate_plan(
  plan,
  rules = list(
    dir__ = tree_names
  )
)
