# recovering trees that crashed
# setwd("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns")
# rm(list=ls())

library(Metrics)
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

# Loading all trees
load.all.trees <- function() {
  tree_info <- read.csv(file=file_in("data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
  tree_names <- unique(tree_info$treeName)
  tree_names <- tree_names[grepl("1$", tree_names)] # for speed, only take a tenth of the trees: those ending in a 1.
  trees <- list()
  for (i in seq_along(tree_names)) {
    trees[[i]] <- ape::read.tree(paste0("data/title_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
  }
  names(trees) <- tree_names
return(trees)
}
  
treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

# Checking which ones didn't crashed
get.all.results <- function() {
  all_results <- data.frame()
  model_average_results <- data.frame()
  best_results <- data.frame() 
  dones <- c(list.files("results", pattern="done_recon", full.names=TRUE)) #list.files("results_March2021", pattern="done_recon", full.names=TRUE)
  file_to_tree_number <- c()
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
  rates.ours <- treeTipMerge(all_results)
  rates.ours.model.average <- treeTipMerge(model_average_results)
  done_trees_number <- unique(model_average_results$treeName)
  return(done_trees_number)
}

# 
# Selecting those that did crash
get.crashed.trees <- function(trees, done_trees_number) {
  tree_names <- names(trees)
  crashed_trees_number <- tree_names[!tree_names %in% done_trees_number]
  crashed_trees <- trees[which(names(trees) %in% crashed_trees_number)]
  return(crashed_trees)
}


# do single run for crashed trees


DoSingleRun_crashed <- function(dir, phy, root_type="madfitz", n.cores=NULL) {
  phy <- phy[which(names(phy)==dir)]
  possibilities = hisse::generateMiSSEGreedyCombinations(max.param=round(ape::Ntip(phy)/10), vary.both=TRUE)
  
  start_time <- Sys.time()
  if(!is.null(phy)) {
    summary_df <- data.frame()
    hisse_result_all <- hisse::MiSSEGreedy(phy, f=1, root.type=root_type, possible.combos=possibilities, chunk.size=1, n.cores=1, save.file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, ".rda"), stop.deltaAICc=1000, sann=TRUE)
    hisse_result_nonredundant <- PruneRedundantModels(hisse_result_all)
    AIC_weights <- hisse::GetAICWeights(hisse_result_nonredundant, criterion="AIC")
    delta_AIC <- sapply(hisse_result_nonredundant, "[[", "AIC") - min(sapply(hisse_result_nonredundant, "[[", "AIC"))
    AICc_weights <- hisse::GetAICWeights(hisse_result_nonredundant, criterion="AICc")
    delta_AICc <- sapply(hisse_result_nonredundant, "[[", "AICc") - min(sapply(hisse_result_nonredundant, "[[", "AICc"))
    save(list=ls(), file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, "_donefitting.rda"))
    model_fit_time <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
    for(model_index in sequence(length(hisse_result_nonredundant))) {
      if(delta_AICc[model_index]<20) {
        cat(paste0("Doing recon on model ", model_index, " at ", Sys.time(), "\n"), file=paste0("results/",unname(Sys.info()["nodename"]), "_",tree_index, ".log"), append=TRUE)
        
        start_time <- Sys.time()
        nturnover <- length(unique(hisse_result_nonredundant[[model_index]]$turnover))
        neps <- length(unique(hisse_result_nonredundant[[model_index]]$eps)) - ifelse(is.null(hisse_result_nonredundant[[model_index]]$fixed.eps), 0,1)
        
        
        hisse_recon <- hisse::MarginReconMiSSE(phy=hisse_result_nonredundant[[model_index]]$phy, f=1, hidden.states=nturnover, fixed.eps=hisse_result_nonredundant[[model_index]]$fixed.eps, pars=hisse_result_nonredundant[[model_index]]$solution, AIC=hisse_result_nonredundant[[model_index]]$AIC, root.type=root_type, get.tips.only=TRUE, n.cores=1)
        
        save(hisse_recon, hisse_result_nonredundant, hisse_result_all, file=paste0("results/", unname(Sys.info()["nodename"]), "_pre_summarizing_recon_",tree_index, "_model_", model_index, "_raw_.rda"))
        
        
        tip_mat_transformed <- hisse_recon$tip.mat[,-1]
        if(max(tip_mat_transformed) == 0) {
          tip_mat_transformed[,1] <- 1 #deal with misse bug of no weight if no hidden
        }
        summary_df_local <- t(hisse_recon$rates.mat %*% t(tip_mat_transformed))
        summary_df_local <- data.frame(summary_df_local)
        summary_df_local$treeName <- dir
        summary_df_local$tree_index <- tree_index
        summary_df_local$model_index <- model_index
        summary_df_local$length_all_models <- length(hisse_result_all)
        summary_df_local$length_nonredundant_models <- length(hisse_result_nonredundant)
        summary_df_local$length_nonredundant_delta_below_20 <- length(which(delta_AICc<20))
        
        summary_df_local$taxon_id_in_phy <- sequence(nrow(summary_df_local))
        summary_df_local$tipName <- phy$tip.label
        summary_df_local$nturnover <- nturnover
        summary_df_local$neps <- neps
        summary_df_local$nparam <- neps + nturnover
        summary_df_local$fixed_eps <- hisse_result_nonredundant[[model_index]]$fixed.eps
        summary_df_local$loglik <- hisse_result_nonredundant[[model_index]]$loglik
        
        summary_df_local$AIC <- hisse_result_nonredundant[[model_index]]$AIC
        summary_df_local$AICc <- hisse_result_nonredundant[[model_index]]$AICc
        summary_df_local$AIC_weight <- AIC_weights[model_index]
        summary_df_local$AICc_weight <- AICc_weights[model_index]
        summary_df_local$deltaAIC <- delta_AIC[model_index]
        summary_df_local$deltaAICc <- delta_AICc[model_index]
        summary_df_local$root_type <- hisse_result_nonredundant[[model_index]]$root.type
        summary_df_local$elapsed_mins_recon <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
        summary_df_local$elapsed_mins_params_all_models_together <- model_fit_time
        summary_df_local$ntip <- ape::Ntip(hisse_result_nonredundant[[model_index]]$phy)
        summary_df_local$treeHeight <- max(branching.times(hisse_result_nonredundant[[model_index]]$phy))
        summary_df_local$min_brlen <- min(hisse_result_nonredundant[[model_index]]$phy$edge.length)
        if(nrow(summary_df)==0) {
          try(summary_df <- summary_df_local)
        } else {
          try(summary_df <- plyr::rbind.fill(summary_df, summary_df_local))
        }
        save(summary_df, hisse_result_nonredundant, AICc_weights, delta_AICc, file=paste0("results/", unname(Sys.info()["nodename"]), "_post_summarizing_recon_",tree_index,"_model_", model_index, "_newrun.rda"))
      }
    }
    save(summary_df, hisse_result_nonredundant, hisse_result_all, AICc_weights, delta_AICc, model_fit_time, file=paste0("results/", unname(Sys.info()["nodename"]), "_done_recon_",tree_index, "_newrun.rda"))
    return(summary_df)
  } else {
    return(list(failure=dir))
  }
}


# Rerunning trees
Recover_crashed_trees <- drake_plan(
    all_trees = load.all.trees(),
    done_tree_numbers = get.all.results(),
    crashed_trees = get.crashed.trees(all_trees, done_trees_number),
    crashed_trees_names = names(crashed_trees),
    target(DoSingleRun_crashed(dir=crashed_trees_names, phy=crashed_trees, root_type="madfitz", n.cores=NULL), 
           dynamic=map(crashed_trees_names))
  )
    

