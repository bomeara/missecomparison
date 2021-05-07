tree_info <- read.csv(file=file_in("data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
tree_names <- unique(tree_info$treeName)
tree_names <- tree_names[grepl("1$", tree_names)] # for speed, only take a tenth of the trees: those ending in a 1.
trees <- list()
for (i in seq_along(tree_names)) {
  trees[[i]] <- ape::read.tree(paste0("data/title_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
}
names(trees) <- tree_names

print(paste0("There are ", length(trees), " to analyze"))

# # Works but single variable
# plan <- drake_plan(
#   #hisse_out = DoSingleRun("Rabosky2014_DD_k1_1")
#    hisse_out = target(DoSingleRun(dir), transform = map(dir = !!tree_names))
# )
#
# plan <- drake_plan(
#    hisse_out = target(DoSingleRun(dir, nturnover=turnover_states, neps_same=neps_same_states, root_type=root_type_states), transform = cross(dir = !!tree_names, turnover_	states=!!sequence(10), neps_same_states=c(TRUE), root_type_states=!!c('madfitz')))
# )


#Working

# plan <- drake_plan(
#    hisse_out = target(DoSingleRun(dir=tree_names[tree_index], phy=trees[[tree_index]], nturnover=turnover_states, neps_same=neps_same_states, root_type=root_type_states), transform = cross(tree_index=!!sequence(length(trees)), turnover_states=!!sequence(12), neps_same_states=c(TRUE), root_type_states=!!c('madfitz')))
# )


plan <- drake_plan(
   possible_combos = hisse::generateMiSSEGreedyCombinations(max.param=50, vary.both=TRUE),
   hisse_out = target(DoSingleRun(dir=tree_names[tree_index], phy=trees[[tree_index]], root_type=root_type_states, possibilities=possible_combos, tree_index=tree_index), transform = cross(tree_index=!!sequence(length(tree_names)),  root_type_states=!!c('madfitz'))),
   combined_df = target(
      dplyr::bind_rows(hisse_out, .id = "id"),
      transform = combine(hisse_out)
   )
)


plan_mac <- drake_plan(
   possible_combos = hisse::generateMiSSEGreedyCombinations(max.param=10, vary.both=TRUE),
   hisse_out = target(DoSingleRun(dir=tree_names[tree_index], phy=trees[[tree_index]], root_type=root_type_states, possibilities=possible_combos, tree_index=tree_index), transform = cross(tree_index=!!sequence(length(tree_names)),  root_type_states=!!c('madfitz'))),
   combined_df = target(
      dplyr::bind_rows(hisse_out, .id = "id"),
      transform = combine(hisse_out)
   )
)


plan_hpc <- drake_plan(
   possible_combos = hisse::generateMiSSEGreedyCombinations(max.param=10, vary.both=TRUE),
   hisse_out = target(DoSingleRun(dir=tree_names[tree_index], phy=trees[[tree_index]], root_type=root_type_states, possibilities=possible_combos, tree_index=tree_index, n.cores=1), transform = cross(tree_index=!!sequence(length(tree_names)),  root_type_states=!!c('madfitz'))),
   combined_df = target(
      dplyr::bind_rows(hisse_out, .id = "id"),
      transform = combine(hisse_out)
   )
)