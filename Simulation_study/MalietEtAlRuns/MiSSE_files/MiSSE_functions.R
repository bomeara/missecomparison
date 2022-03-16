# setwd("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison/")
# rm(list=ls())

library(tidyverse)
library(magrittr)
library(stats)
library(ggplot2)
library(progress)
library(drake)
#library(hisse)
devtools::install_github("thej022214/hisse")
library(hisse)
library(ape)
library(purrr)
library(dplyr)
library(future)
library(future.batchtools)
library(batchtools)
library(parallel)
library(data.table)


# Loading all trees
load.all.trees <- function(base.dir, where=c("local","labcomputer"), ref=c("title_rabosky","maliet")) {
  if(ref=="title_rabosky") {
    if(where=="labcomputer") {
    tree_info <- read.csv(file=file_in("/home/tvasconcelos/missecomparison/TitleRaboskyRuns/data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
    }
    if(where=="local") {
    tree_info <- read.csv(file=file_in("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns/data/title_rabosky_dryad/tipRates_dryad/dataFiles/treeSummary.csv"),stringsAsFactors=FALSE)
    }
    tree_names <- unique(tree_info$treeName)
    tree_names <- tree_names[grepl("1$", tree_names)] # for speed, only take a tenth of the trees: those ending in a 1.
    trees <- list()
    for (i in seq_along(tree_names)) {
      trees[[i]] <- ape::read.tree(paste0(base.dir,"/data/title_rabosky_dryad/trees/", tree_names[i], "/", tree_names[i], ".tre"))
    }
    names(trees) <- tree_names
  return(trees)
  }
  if(ref=="maliet"){
    trees <- list()
    if(where=="local"){
      path <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/Simulation_study/MalietEtAlRuns/MalietEtAl_ClaDS2_trees/"
    }
    if(where=="labcomputer") {
      path <- "/home/tvasconcelos/missecomparison/Simulation_study/MalietEtAlRuns/MalietEtAl_ClaDS2_trees/"
    }
    tree_files <- list.files(path)
    tree_files <- tree_files[grep(".Rdata", tree_files)]
    for(i in 1:length(tree_files)){
      load(paste0(path, tree_files[i]))
      if(exists("tree")){
        trees[[i]] <- tree
        rm("tree")
      }
    }
    tree_names <- gsub(".Rdata","", tree_files)
    names(trees) <- tree_names
    return(trees)
  }
}
  
treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

# Checking which ones didn't crashed
get.all.results <- function(base.dir, x=NULL) {
  all_results <- data.frame()
  model_average_results <- data.frame()
  best_results <- data.frame() 
  dones <- c(list.files(paste0(base.dir,"/results"), pattern="done_recon", full.names=TRUE)) #list.files("results_March2021", pattern="done_recon", full.names=TRUE)
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

get.crashed.trees2 <- function(base.dir, where=c("local","labcomputer")) {
    all_results <- data.frame()
    model_average_results <- data.frame()
    best_results <- data.frame()
    dones <- c(list.files("results", pattern="done_recon", full.names=TRUE), list.files("new_results", pattern="done_recon", full.names=TRUE))
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

    rates.cleaned <- rates.combined

    dones <- unique(rates.cleaned$treeName)
    
    all_trees = load.all.trees(base.dir=base.dir, where=where)
    to_rerun <- names(all_trees)[!names(all_trees) %in% dones]
    to_rerun <- subset(to_rerun, !grepl("MeyerWiens2017",to_rerun))
  return(to_rerun)
}


# 
# Selecting those that did crash
get.crashed.trees <- function(trees, done_trees_number) {
  tree_names <- names(trees)
  crashed_trees_number <- tree_names[!tree_names %in% done_trees_number]
  crashed_trees <- trees[which(names(trees) %in% crashed_trees_number)]
  return(crashed_trees)
}


# Do single run 
DoSingleRun_new <- function(dir, phy, root_type="madfitz", n.cores=NULL) {
  dir <- unname(dir)
  tree_index <- dir
  phy <- phy[which(names(phy)%in%dir)][[1]]
  possibilities = generateMiSSEGreedyCombinations(max.param=max(4, round(ape::Ntip(phy)/10)), vary.both=TRUE, fixed.eps.tries=NA, shuffle.start=TRUE)
  
  start_time <- Sys.time()
  if(!is.null(phy)) {
    summary_df <- data.frame()
    hisse_result_all <- hisse::MiSSEGreedy(phy, f=1, root.type=root_type, possible.combos=possibilities, chunk.size=10, n.cores=1, save.file=paste0("/home/tvasconcelos/missecomparison/Simulation_study/TitleRaboskyRuns/new_results/",unname(Sys.info()["nodename"]), "_",tree_index, ".rda"), stop.deltaAICc=10, sann=TRUE)
    hisse_result_nonredundant <- PruneRedundantModels(hisse_result_all)
    AIC_weights <- hisse::GetAICWeights(hisse_result_nonredundant, criterion="AIC")
    delta_AIC <- sapply(hisse_result_nonredundant, "[[", "AIC") - min(sapply(hisse_result_nonredundant, "[[", "AIC"))
    AICc_weights <- hisse::GetAICWeights(hisse_result_nonredundant, criterion="AICc")
    delta_AICc <- sapply(hisse_result_nonredundant, "[[", "AICc") - min(sapply(hisse_result_nonredundant, "[[", "AICc"))
    save(list=ls(), file=paste0("/home/tvasconcelos/missecomparison/Simulation_study/TitleRaboskyRuns/new_results/",unname(Sys.info()["nodename"]), "_",tree_index, "_donefitting.rda"))
    model_fit_time <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
    for(model_index in sequence(length(hisse_result_nonredundant))) {
      if(delta_AICc[model_index]<20) {
        cat(paste0("Doing recon on model ", model_index, " at ", Sys.time(), "\n"), file=paste0("/home/tvasconcelos/missecomparison/Simulation_study/TitleRaboskyRuns/new_results/",unname(Sys.info()["nodename"]), "_",tree_index, ".log"), append=TRUE)
        
        start_time <- Sys.time()
        nturnover <- length(unique(hisse_result_nonredundant[[model_index]]$turnover))
        neps <- length(unique(hisse_result_nonredundant[[model_index]]$eps)) - ifelse(is.null(hisse_result_nonredundant[[model_index]]$fixed.eps), 0,1)
        hisse_recon <- hisse::MarginReconMiSSE(phy=hisse_result_nonredundant[[model_index]]$phy, f=1, hidden.states=nturnover, fixed.eps=hisse_result_nonredundant[[model_index]]$fixed.eps, pars=hisse_result_nonredundant[[model_index]]$solution, AIC=hisse_result_nonredundant[[model_index]]$AIC, root.type=root_type, get.tips.only=TRUE, n.cores=1)
        save(hisse_recon, hisse_result_nonredundant, hisse_result_all, file=paste0("/home/tvasconcelos/missecomparison/Simulation_study/TitleRaboskyRuns/new_results/", unname(Sys.info()["nodename"]), "_pre_summarizing_recon_",tree_index, "_model_", model_index, "_raw_.rda"))
        
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
        save(summary_df, hisse_result_nonredundant, AICc_weights, delta_AICc, file=paste0("/home/tvasconcelos/missecomparison/Simulation_study/TitleRaboskyRuns/new_results/", unname(Sys.info()["nodename"]), "_post_summarizing_recon_",tree_index,"_model_", model_index, "_newrun.rda"))
      }
    }
    save(summary_df, hisse_result_nonredundant, hisse_result_all, AICc_weights, delta_AICc, model_fit_time, file=paste0("/home/tvasconcelos/missecomparison/Simulation_study/TitleRaboskyRuns/new_results/", unname(Sys.info()["nodename"]), "_done_recon_",tree_index, "_newrun.rda"))
    return(summary_df)
  } else {
    return(list(failure=dir))
  }
}


# Rerunning trees
# Dividing into two runs to use less cores 
# First run
# Recover_crashed_trees1 <- drake_plan(
#    base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
#    #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
#    all_trees = load.all.trees(base.dir),
#    done_tree_numbers = get.all.results(base.dir, x=all_trees),
#    crashed_trees = get.crashed.trees(all_trees, done_tree_numbers),
#    crashed_trees_names = names(crashed_trees),
#    subset_crashed = crashed_trees_names[1:40],
#    target(DoSingleRun_new(dir=subset_crashed, 
#                               phy=crashed_trees, 
#                               root_type="madfitz", n.cores=NULL), dynamic=map(subset_crashed))
#  )

# Second run
# Recover_crashed_trees2 <- drake_plan(
#  base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
#  #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
#  all_trees = load.all.trees(base.dir),
#  done_tree_numbers = get.all.results(base.dir),
#  crashed_trees = get.crashed.trees(all_trees, done_tree_numbers),
#  crashed_trees_names = names(crashed_trees),
#  target(DoSingleRun_new(dir=crashed_trees_names, 
#                             phy=crashed_trees, 
#                             root_type="madfitz", n.cores=NULL), dynamic=map(dir))
# )

#############################################################
# New runs
New_sim_runs_beaulieu3 <- drake_plan(
  base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
  #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
  all_trees = load.all.trees(base.dir, where="labcomputer"),
  tree_names = names(all_trees),
  subset_trees1 = tree_names[1:96],
  target(DoSingleRun_new(dir=subset_trees1, 
                             phy=all_trees, 
                             root_type="madfitz", n.cores=NULL), dynamic=map(subset_trees1))
)

New_sim_runs_beaulieu4 <- drake_plan(
  base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
  #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
  all_trees = load.all.trees(base.dir, where="labcomputer"),
  tree_names = names(all_trees),
  subset_trees2 = tree_names[97:192],
  target(DoSingleRun_new(dir=subset_trees2, 
                         phy=all_trees, 
                         root_type="madfitz", n.cores=NULL), dynamic=map(subset_trees2))
)

#write.csv(tree_names[193:312], file="trees_for_Brian.csv") # These are going to Brian's 

Second_sim_runs_beaulieu4 <- drake_plan(
  base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
  #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
  all_trees = load.all.trees(base.dir, where="labcomputer"),
  tree_names = names(all_trees),
  subset_trees3 = tree_names[313:372],
  target(DoSingleRun_new(dir=subset_trees3, 
                         phy=all_trees, 
                         root_type="madfitz", n.cores=NULL), dynamic=map(subset_trees3))
)


Second_sim_runs_beaulieu3 <- drake_plan( # all done
  base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
  #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
  all_trees = load.all.trees(base.dir, where="labcomputer"),
  tree_names = names(all_trees),
  subset_trees4 = tree_names[373:412],
  target(DoSingleRun_new(dir=subset_trees4, 
                         phy=all_trees, 
                         root_type="madfitz", n.cores=NULL), dynamic=map(subset_trees4))
)


#write.csv(tree_names[413:482], file="more_trees_for_Brian.csv")

Third_sim_runs_beaulieu4 <- drake_plan(
  base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
  #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
  all_trees = load.all.trees(base.dir, where="local"),
  tree_names = names(all_trees),
  subset_trees5 = tree_names[483:521],
  target(DoSingleRun_new(dir=subset_trees5, 
                         phy=all_trees, 
                         root_type="madfitz", n.cores=NULL), dynamic=map(subset_trees5))
)

crashed_runs_beaulieu4 <- drake_plan(
    dones = c(list.files("results", pattern="done_recon", full.names=TRUE), list.files("new_results", pattern="done_recon", full.names=TRUE)),
    base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
    all_trees = load.all.trees(base.dir, where="labcomputer"),
#---- possible crashes in labcomputer4
    dones_beaulieu4 = dones[grep("beaulieulab4", dones)],
    names_done_b4 = gsub(paste0(c("_newrun.rda","new_results/beaulieu3_done_recon_","new_results/beaulieulab4_done_recon_"), collapse="|"),"", dones_beaulieu4),
    subset_b4 = names(all_trees)[c(97:192, 313:372, 483:521)], # trees that ran at beaulieu4
    possible_crashes_b4 = subset_b4[which(!subset_b4 %in% names_done_b4)],
    target(DoSingleRun_new(dir=possible_crashes_b4, 
                       phy=all_trees, 
                       root_type="madfitz", n.cores=NULL), dynamic=map(possible_crashes_b4))
)

##########################################
# 2021-09-13
##########################################

# Had to run this out of drake to save names of crashed trees
#dones = c(list.files("results", pattern="done_recon", full.names=TRUE), list.files("new_results", pattern="done_recon", full.names=TRUE))
#base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns"
#all_trees = load.all.trees(base.dir, where="labcomputer")
#crashed_trees = get.crashed.trees2(base.dir, where="labcomputer")
#write.csv(crashed_trees, file="crashed.trees.csv")


rerun_crashed_runs_beaulieu4_final <- drake_plan(
  dones = c(list.files("results", pattern="done_recon", full.names=TRUE), list.files("new_results", pattern="done_recon", full.names=TRUE)),
  #base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
  base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
  all_trees = load.all.trees(base.dir, where="local"),
  crashed_trees = get.crashed.trees2(base.dir, where="local"),
  target(DoSingleRun_new(dir=crashed_trees, 
                         phy=all_trees, 
                         root_type="madfitz", n.cores=NULL), dynamic=map(crashed_trees))
)


##########################################
# 2022-03-02
##########################################
# Running MiSSE on ClaDS simulated data

run_misse_on_clads <- drake_plan(
  base.dir = "/home/tvasconcelos/missecomparison/TitleRaboskyRuns",
  #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
  all_trees = load.all.trees(base.dir, where="labcomputer", ref="maliet"),
  tree_names = names(all_trees),
  subset_trees1 = tree_names[1:40],
  target(DoSingleRun_new(dir=subset_trees1, 
                         phy=all_trees, 
                         root_type="madfitz", n.cores=NULL), dynamic=map(subset_trees1))
)


run_misse_on_clads_crashed <- drake_plan(
  base.dir = "/home/tvasconcelos/missecomparison/Simulation_study/TitleRaboskyRuns",
  #base.dir = "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns",
  all_trees = load.all.trees(base.dir, where="labcomputer", ref="maliet"),
  #all_trees = load.all.trees(base.dir, where="local", ref="maliet"),
  #tree_names = names(all_trees),
  subset_trees1 = c("ClaDS2tree_100_3_2","ClaDS2tree_100_4_1","ClaDS2tree_100_6_2","ClaDS2tree_100_7_4"),
  target(DoSingleRun_new(dir=subset_trees1, 
                         phy=all_trees, 
                         root_type="madfitz", n.cores=NULL), dynamic=map(subset_trees1))
)


