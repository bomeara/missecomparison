# ansible linux -a 'nohup Rscript /share/missecomparison/TitleRaboskyRuns/make_manual_parallelization.R &' -f 10

setwd("/share/missecomparison/TitleRaboskyRuns")


source("R/packages.R")
source("R/functions.R")

library(foreach)
library(doParallel)
registerDoParallel(floor(parallel::detectCores()/10)) # so on machines with >20 cores it runs multiple

Sys.setenv('R_MAX_VSIZE'=32000000000)
ipifyseed <- 100
try(ipifyseed <- sqrt(as.numeric(gsub("\\.", "", as.character(ipify::get_ip())))))
set.seed(as.integer(round(runif(1, min=1, max=1e5)) + ipifyseed))

trees_for_me_to_run <- read.csv("trees_for_Brian.csv")

tree_info <- read.csv(file="data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv",stringsAsFactors=FALSE)
tree_names <- unique(tree_info$treeName)
#tree_names <- tree_names[grepl("1$", tree_names)] # for speed, only take a tenth of the trees: those ending in a 1.
tree_names <- tree_names[tree_names%in%trees_for_me_to_run[,2]]
trees <- list()
for (i in seq_along(tree_names)) {
  trees[[i]] <- ape::read.tree(paste0("data/title_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
}
names(trees) <- tree_names



tree_indices <- sample(sequence(length(tree_names)), replace=FALSE) # randomize order

results <- list()
#for (i in seq_along(tree_indices)) {
foreach (i=seq_along(tree_indices)) %dopar% { 
	started_runs <- list.files(path="results", pattern="starting.*.rda")
	tree_index <- tree_indices[i]

	if(!any(grepl(paste0("starting_", tree_index, "_.rda"), started_runs))) { # so we skip ones already started
		starting_session <- sessionInfo()
		node <- unname(Sys.info()["nodename"])
		save(tree_index, starting_session, node, file=paste0("results/starting_", tree_index, "_.rda"))
		local_result <- NULL
		possible_combos = hisse::generateMiSSEGreedyCombinations(max.param=max(4,round(ape::Ntip(trees[[tree_index]])/10)), vary.both=TRUE)
		try(local_result <- DoSingleRun(dir=tree_names[tree_index], phy=trees[[tree_index]], root_type="madfitz", possibilities=possible_combos, tree_index=tree_index, n.cores=parallel::detectCores(), chunk.size=10))
		#save(local_result, file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, "_local_result_newrun.rda"))
		if(!is.null(local_result)) {
			results[[i]] <- local_result
			#save(results, tree_indices, file=paste0("manual", Sys.info()['nodename'], "_newrun.rda"))
		}
	}
}