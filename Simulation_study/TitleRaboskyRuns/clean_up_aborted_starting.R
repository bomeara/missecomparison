setwd("/share/missecomparison/TitleRaboskyRuns")

source("R/packages.R")
source("R/functions.R")

library(foreach)
library(doParallel)
registerDoParallel(parallel::detectCores())

Sys.setenv('R_MAX_VSIZE'=32000000000)
ipifyseed <- 100
try(ipifyseed <- sqrt(as.numeric(gsub("\\.", "", as.character(ipify::get_ip())))))
set.seed(as.integer(round(runif(1, min=1, max=1e5)) + ipifyseed))

tree_info <- read.csv(file="data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv",stringsAsFactors=FALSE)
tree_names <- unique(tree_info$treeName)
tree_names <- tree_names[grepl("1$", tree_names)] # for speed, only take a tenth of the trees: those ending in a 1.
trees <- list()
for (i in seq_along(tree_names)) {
  trees[[i]] <- ape::read.tree(paste0("data/title_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
}
names(trees) <- tree_names

started_runs <- list.files(path="results", pattern="starting.*.rda")
finished_runs <- list.files(path="results", pattern="*done_recon.*.rda")

for (tree_index in seq_along(trees)) {
	matching <- sum(grepl(paste0("done_recon_", tree_index, "_newrun.rda"), finished_runs))
	if(matching==0) {
		system(paste0("rm results/starting_", tree_index, "_.rda"))
	}	
}