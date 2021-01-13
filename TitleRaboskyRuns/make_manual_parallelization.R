# ansible linux -a 'nohup Rscript /share/missecomparison/TitleRaboskyRuns/make_manual_parallelization.R &'

setwd("/share/missecomparison/TitleRaboskyRuns")


source("R/packages.R")
source("R/functions.R")
Sys.setenv('R_MAX_VSIZE'=32000000000)

set.seed(as.integer(round(runif(1, min=1, max=1e5)) + sqrt(as.numeric(gsub("\\.", "", as.character(ipify::get_ip()))))))

tree_info <- read.csv(file=file_in("data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
tree_names <- unique(tree_info$treeName)
tree_names <- tree_names[grepl("1$", tree_names)] # for speed, only take a tenth of the trees: those ending in a 1.
trees <- list()
for (i in seq_along(tree_names)) {
  trees[[i]] <- ape::read.tree(paste0("data/title_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
}
names(trees) <- tree_names


possible_combos = hisse::generateMiSSEGreedyCombinations(max.param=50, vary.both=TRUE)

tree_indices <- sample(sequence(length(tree_names)), replace=FALSE) # randomize order

results <- list()
for (i in seq_along(tree_indices)) {
	results[[i]] <- DoSingleRun(dir=tree_names[tree_index], phy=trees[[tree_index]], root_type=root_type_states, possibilities=possible_combos, tree_index=tree_index)
	save(results, tree_indices, file=paste0("manual", Sys.info()['nodename'], ".rda"))
}