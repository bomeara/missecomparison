tree_info <- read.csv(file=file_in("data/pascal_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
tree_names <- unique(tree_info$treeName)[1:3] #limit to 1:3 for trial
tree_names <- paste0("evolvingRates_", sequence(4))
trees <- list()
for (i in seq_along(tree_names)) {
  trees[[i]] <- ape::read.tree(paste0("data/pascal_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
}
names(trees) <- tree_names


plan <- drake_plan(
  #hisse_out = DoSingleRun("Rabosky2014_DD_k1_1")
   hisse_out = target(
     DoSingleRun(phy=phy[[1]], neps=neps),
     transform = cross(phy=trees, neps=c(2))
   )
)
