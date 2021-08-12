
# recovering trees that crashed
#setwd("~/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns") # TV local
# setwd("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns")
# rm(list=ls())

library(tidyverse)
library(magrittr)
library(stats)   
library(ggplot2)
library(progress)
library(drake)
library(hisse)
library(ape)
library(purrr)
library(dplyr)
library(future)
library(future.batchtools)
library(batchtools)
library(parallel)
library(data.table)
library(Metrics)
library(viridis)

################################################################
# Loading all trees
################################################################
tree_info <- read.csv(file=file_in("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns/data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
#tree_info <- read.csv(file=file_in("/home/tvasconcelos/missecomparison/TitleRaboskyRuns/data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
tree_names <- unique(tree_info$treeName)
tree_names <- tree_names[grepl("1$", tree_names)] # for speed, only take a tenth of the trees: those ending in a 1.
trees <- list()
for (i in seq_along(tree_names)) {
  trees[[i]] <- ape::read.tree(paste0(getwd(),"/data/title_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
}
names(trees) <- tree_names

################################################################
# Loading preliminary results
################################################################
# This is useful because Brian's code also loads and organizes
# the true values
################################################################
treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

all_results <- data.frame()
model_average_results <- data.frame()
best_results <- data.frame()
dones <- c(list.files("results", pattern="done_recon", full.names=TRUE), list.files("results_March2021", pattern="done_recon", full.names=TRUE))
for (i in seq_along(dones)) {
  print(paste0("Loading ", i, " of ", length(dones)))
  try(load(dones[i]))
  if(exists("summary_df")) {
    all_results <- plyr::rbind.fill(all_results, summary_df)
    local_weighted <- summary_df %>% group_by(treeName, taxon_id_in_phy, tipName) %>% summarise(
      turnover = weighted.mean(turnover, AICc_weight),
      extinction.fraction = weighted.mean(extinction.fraction, AICc_weight),
      net.div = weighted.mean(net.div, AICc_weight), 
      speciation = weighted.mean(speciation, AICc_weight),
      extinction = weighted.mean(extinction, AICc_weight)
    )
    model_average_results <- plyr::rbind.fill(model_average_results, local_weighted)
    best_results <- plyr::rbind.fill(best_results, subset(summary_df, deltaAICc==0))
    
    rm(summary_df)
    rm(delta_AICc)
    rm(local_weighted)
  }
  
}
# Loading all results as they may be useful
rates.ours <- treeTipMerge(all_results)
rates.ours.model.average <- treeTipMerge(model_average_results)
rates.theirs <- treeTipMerge(read.csv("data/title_rabosky_dryad/tipRates_dryad/dataFiles/estimatedTipRates.csv", stringsAsFactors=FALSE))
rates.true <- treeTipMerge(read.csv("data/title_rabosky_dryad/tipRates_dryad/dataFiles/trueTipRates.csv", stringsAsFactors=FALSE))
rates.ours$unique_string <- paste(rates.ours$treeName, "nturnover", rates.ours$nturnover, "neps", rates.ours$neps, "root", rates.ours$root_type, sep="_")
rates.ours.best <- best_results
rates.theirs <- rates.theirs[rates.theirs$treeName %in% rates.ours.best$treeName,]
rates.true <- rates.true[rates.true$treeName %in% rates.ours.best$treeName,]

colnames(rates.true)[which(colnames(rates.true)=="tipLambda")] <- "lambdaTRUE"
colnames(rates.true)[which(colnames(rates.true)=="tipMu")] <- "muTRUE"
colnames(rates.true)[which(colnames(rates.true)=="tipNetDiv")] <- "netDivTRUE"
rates.true$turnoverTRUE <- rates.true$lambdaTRUE + rates.true$muTRUE
rates.true$extinctionFractionTRUE <- rates.true$muTRUE / rates.true$lambdaTRUE

colnames(rates.ours.best)[which(colnames(rates.ours.best)=="turnover")] <- "turnoverMiSSEbest"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="extinction.fraction")] <- "extinctionFractionMiSSEbest"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="net.div")] <- "netDivMiSSEbest"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="speciation")] <- "lambdaMiSSEbest"
colnames(rates.ours.best)[which(colnames(rates.ours.best)=="extinction")] <- "muMiSSEbest"

colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="turnover")] <- "turnoverMiSSEavg"
colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="extinction.fraction")] <- "extinctionFractionMiSSEavg"
colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="net.div")] <- "netDivMiSSEavg"
colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="speciation")] <- "lambdaMiSSEavg"
colnames(rates.ours.model.average)[which(colnames(rates.ours.model.average)=="extinction")] <- "muMiSSEavg"

rates.combined <- merge(rates.ours.best, rates.ours.model.average)
rates.combined <- merge(rates.combined, rates.theirs)
rates.combined <- merge(rates.combined, rates.true)

rates.combined$muBAMM <- rates.combined$lambdaBAMM - rates.combined$netDivBAMM
rates.combined$turnoverBAMM <- rates.combined$muBAMM + rates.combined$lambdaBAMM
rates.combined$extinctionFractionBAMM <- rates.combined$muBAMM / rates.combined$lambdaBAMM

# get rid of NAs and weird values for now?
rates.cleaned <- rates.combined
rates.cleaned <- rates.cleaned[!is.na(rates.cleaned$netDivBAMM),]
rates.cleaned <- rates.cleaned[(rates.cleaned$lambdaTRUE>0),] #yeah, not sure why there'd be no speciation in reality for a tree with >2 taxa

################################################################
# Do single run to test specific things
################################################################

all_files <- list.files(paste0(getwd(),"/results"))
include <- c("omeara")
exclude <- c("recon","donefitting",".log")
raw_files <- all_files[!grepl(paste(exclude, collapse="|"), all_files)]
raw_files <- unique(grep(paste(include, collapse="|"), raw_files, value=TRUE))
all_lik_results <- list()
for(result_index in seq_along(raw_files)) {
  try(load(paste0(getwd(),"/results/",raw_files[result_index])))   
  if(exists("possible.combos")) {
    all_lik_results[[result_index]] <- possible.combos
    #names(all_lik_results)[result_index] <- dir # adding tree name
  }
}
# load(file="results_fit.rda")
# starting.vals

sub <- subset(rates.cleaned, rates.cleaned$treeName %in% names(random_tree))

View(sub)

misse.list
trees_done <- trees[names(trees) %in% rates.cleaned$treeName]


random_tree <- trees_done[sample(1:length(trees_done), 1)]
names(all_lik_results) == names(random_tree)


# get.true.values 



  phy <- random_tree[[1]]
  possibilities = hisse::generateMiSSEGreedyCombinations(max.param=round(ape::Ntip(phy)/10), vary.both=TRUE, fixed.eps.tries = NA)[1:5,]
  

    hisse_result_all <- hisse::MiSSEGreedy(phy, f=1, possible.combos=possibilities, chunk.size=4, n.cores=4, save.file=paste0(getwd(),"/", "test_lik.rda"), stop.deltaAICc=10, sann=TRUE)
    
    
    
    hisse_result_nonredundant <- PruneRedundantModels(hisse_result_all)
    AIC_weights <- hisse::GetAICWeights(hisse_result_nonredundant, criterion="AIC")
    delta_AIC <- sapply(hisse_result_nonredundant, "[[", "AIC") - min(sapply(hisse_result_nonredundant, "[[", "AIC"))
    AICc_weights <- hisse::GetAICWeights(hisse_result_nonredundant, criterion="AICc")
    delta_AICc <- sapply(hisse_result_nonredundant, "[[", "AICc") - min(sapply(hisse_result_nonredundant, "[[", "AICc"))
 


