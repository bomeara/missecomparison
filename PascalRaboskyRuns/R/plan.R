tree_info <- read.csv(file=file_in("data/pascal_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
tree_names <- unique(tree_info$treeName)[1:3] #limit to 1:3 for trial
tree_names <- paste0("evolvingRates_", sequence(4))
trees <- list()
for (i in seq_along(tree_names)) {
  trees[[i]] <- ape::read.tree(paste0("data/pascal_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
}
names(trees) <- tree_names

# # Works but single variable
# plan <- drake_plan(
#   #hisse_out = DoSingleRun("Rabosky2014_DD_k1_1")
#    hisse_out = target(DoSingleRun(dir), transform = map(dir = !!tree_names))
# )

plan <- drake_plan(
   hisse_out = target(DoSingleRun(dir, nturnover=turnover_states, neps_same=neps_same_states, root_type=root_type_states), transform = cross(dir = !!tree_names, turnover_states=!!sequence(10), neps_same_states=c(TRUE,FALSE), root_type_states=!!c('madfitz')))
)
