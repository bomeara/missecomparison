# rm(list=ls())
#  setwd("~/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns")
library(Metrics)
library(tidyverse)
library(magrittr)
library(stats)
library(ggplot2)
library(progress)
library(viridis)
library(hisse)
library(PMCMR)

treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}

#base.dir <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison"
base.dir <- "/home/tvasconcelos/missecomparison/ClaDScomparison"
############################################################
# Load previous results
############################################################
all_results <- data.frame()
model_average_results <- data.frame()
best_results <- data.frame()
dones <- c(list.files("results", pattern="done_recon", full.names=TRUE), list.files("new_results", pattern="done_recon", full.names=TRUE))


# Getting trees for ClaDS
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
    
    # saving trees for ClaDS
    
    #write.tree(hisse_result_all[[1]]$phy, file=paste0(base.dir, "/trees/", unique(best_results$treeName), ".tre"))
    
    rm(summary_df)
    rm(delta_AICc)
    rm(local_weighted)
  }
}

tree_names <- unique(all_results$treeName)
original_tree_path <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/TitleRaboskyRuns/data/title_rabosky_dryad/trees"

for(i in 1:length(tree_names)) {
  one_tree_label <- list.files(paste0(original_tree_path, "/", tree_names[i]))
  one_tree_label <- one_tree_label[grep(".tre", one_tree_label)]
  one_tree <- read.tree(paste0(original_tree_path, "/", tree_names[i],"/",one_tree_label))
  write.tree(one_tree, paste0(base.dir, "/trees/", tree_names[i], ".tre"))
}


# Making ClaDS script for manual parallelization in Julia


new_tree_path <- paste0(base.dir, "/trees")

groups <- seq(from=1, to=length(tree_names), by=50)
for(k in 1:length(groups)) {
  if(k==length(groups)){
    one_group <- tree_names[groups[k]:length(tree_names)]
  } else { one_group <- tree_names[groups[k]:(groups[k]+50)] }
  sink(paste0("group",k,"_script.jl"))
  cat("# mode: julia","\n")
  cat("\t","using PANDA","\n")
  cat("# mode: julia","\n")
  cat("\t","using JLD2","\n")
  for(i in 1:length(one_group)) {
    cat("# mode: julia","\n")
    cat("\t",paste0('my_tree = load_tree("', new_tree_path, '/', one_group[i], '.tre")'),"\n")
    cat("# mode: julia","\n")
    cat("\t","output = infer_ClaDS(my_tree, print_state = 100)","\n")
    cat("# mode: julia","\n")
    #cat("@save 'output_clads' output", "\n")
    cat("\t",paste0('save_ClaDS_in_R(output, "',base.dir, '/results/', one_group[i], '_clads_results.Rdata")'),"\n")
  }
    sink()
}


# individual scripts
for(i in 1:length(tree_names)) {
  sink(paste0(base.dir, "/manual_parallel/", tree_names[i],"_script.jl"))
  cat("# mode: julia","\n")
  cat("\t","using PANDA","\n")
  cat("# mode: julia","\n")
  cat("\t",paste0('my_tree = load_tree("', new_tree_path, '/', tree_names[i], '.tre")'),"\n")
  cat("# mode: julia","\n")
  cat("\t","output = infer_ClaDS(my_tree, print_state = 100)","\n")
  cat("# mode: julia","\n")
  cat("\t","using JLD2","\n")
  cat("# mode: julia","\n")
  #cat("@save 'output_clads' output", "\n")
  cat("\t",paste0('save_ClaDS_in_R(output, "',base.dir, '/results/', tree_names[i], '_clads_results.Rdata")'),"\n")
  sink()
}

# making more script for slow trees on Feb 20

#base.dir <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison"
new_tree_path <- paste0(base.dir, "/trees")
labels_all <- gsub(".tre","",list.files(new_tree_path))

clads.dir <- "/Volumes/NO NAME/results"
clads.results <- list.files(clads.dir)
labels_partial_results <- gsub("_clads_results.Rdata","", clads.results)

undones <- setdiff(labels_all, labels_partial_results)

tree_names = undones 
groups <- seq(from=1, to=length(tree_names), by=3)

#base.dir <- "/home/tvasconcelos/missecomparison/ClaDScomparison"
new_tree_path <- paste0(base.dir, "/trees")

for(k in 1:length(groups)) {
  if(k==length(groups)){
    one_group <- tree_names[groups[k]:length(tree_names)]
  } else { one_group <- tree_names[groups[k]:(groups[k]+50)] }
  sink(paste0("group",k,"_script.jl"))
  cat("# mode: julia","\n")
  cat("\t","using PANDA","\n")
  cat("# mode: julia","\n")
  cat("\t","using JLD2","\n")
  for(i in 1:length(one_group)) {
    cat("# mode: julia","\n")
    cat("\t",paste0('my_tree = load_tree("', new_tree_path, '/', one_group[i], '.tre")'),"\n")
    cat("# mode: julia","\n")
    cat("\t","output = infer_ClaDS(my_tree, print_state = 100)","\n")
    cat("# mode: julia","\n")
    #cat("@save 'output_clads' output", "\n")
    cat("\t",paste0('save_ClaDS_in_R(output, "',base.dir, '/results/', one_group[i], '_clads_results.Rdata")'),"\n")
  }
  sink()
}


###################################################
# Making ClaDS script for manual parallelization in Julia - 2nd run
###################################################
new_tree_path <- paste0(base.dir, "/trees")
labels_all <- gsub(".tre","",list.files(new_tree_path))

base.dir <- "/home/tvasconcelos/missecomparison/ClaDScomparison"
new_tree_path <- paste0(base.dir, "/trees")

groups <- seq(from=1, to=length(labels_all), by=10)
for(k in 1:length(groups)) {
  if(k==length(groups)){
    one_group <- labels_all[groups[k]:length(labels_all)]
  } else { one_group <- labels_all[groups[k]:(groups[k]+50)] }
  sink(paste0("group",k,"_script.jl"))
  cat("# mode: julia","\n")
  cat("\t","using PANDA","\n")
  cat("# mode: julia","\n")
  cat("\t","using JLD2","\n")
  for(i in 1:length(one_group)) {
    cat("# mode: julia","\n")
    cat("\t",paste0('my_tree = load_tree("', new_tree_path, '/', one_group[i], '.tre")'),"\n")
    cat("# mode: julia","\n")
    cat("\t","output = infer_ClaDS(my_tree, print_state = 100)","\n")
    cat("# mode: julia","\n")
    #cat("@save 'output_clads' output", "\n")
    cat("\t",paste0('save_ClaDS_in_R(output, "',base.dir, '/results_2nd_run/', one_group[i], '_clads_results.Rdata")'),"\n")
  }
  sink()
}


###################################################
# Making ClaDS script for manual parallelization in Julia - 2nd run
###################################################
base.dir <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison"
new_tree_path <- paste0(base.dir, "/original_trees/trees/ClaDS2")
tree_files <- list.files(new_tree_path)[grep(".tre$", list.files(new_tree_path))]
labels_all <- gsub(".tre$","",tree_files)

base.dir <- "/home/tvasconcelos/missecomparison/ClaDScomparison"
new_tree_path <- paste0(base.dir, "/original_trees/trees/ClaDS2")

groups <- seq(from=1, to=length(labels_all), by=5)
for(k in 1:length(groups)) {
  if(k==length(groups)){
    one_group <- labels_all[groups[k]:length(labels_all)]
  } else { one_group <- labels_all[groups[k]:(groups[k]+50)] }
  sink(paste0("group",k,"_script.jl"))
  cat("# mode: julia","\n")
  cat("\t","using PANDA","\n")
  cat("# mode: julia","\n")
  cat("\t","using JLD2","\n")
  for(i in 1:length(one_group)) {
    cat("# mode: julia","\n")
    cat("\t",paste0('my_tree = load_tree("', new_tree_path, '/', one_group[i], '.tre")'),"\n")
    cat("# mode: julia","\n")
    cat("\t","output = infer_ClaDS(my_tree, print_state = 100)","\n")
    cat("# mode: julia","\n")
    #cat("@save 'output_clads' output", "\n")
    cat("\t",paste0('save_ClaDS_in_R(output, "',base.dir, '/MiSSE_ClaDS_comparison/clads/', one_group[i], '_clads_results.Rdata")'),"\n")
  }
  sink()
}
