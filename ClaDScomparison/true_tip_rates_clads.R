
# rm(list=ls())

# setwd("/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison")
# Extracting trees from simulated data (ClaDS2) to run MiSSE
library(ape)
tree_files <- list.files("original_trees/trees/ClaDS2/")
tree_files <- tree_files[grep("Rdata",tree_files)]
for(i in 1:length(tree_files)) {
  load(paste0(getwd(),"/original_trees/trees/ClaDS2/", tree_files[i]))
  if(exists("tree")){
   label <- gsub(".Rdata","",tree_files[i])
    write.tree(tree, file=paste0(getwd(),"/original_trees/trees/ClaDS2/",label,".tre"))
    rm("tree")
  }
}

# Get the true rates from the ClaDS2 model
tree_files <- list.files("original_trees/trees/ClaDS2/")
tree_files <- tree_files[grep("Rdata",tree_files)]
true_clads_tiprates <- list()
for(i in 1:length(tree_files)) {
  load(paste0(getwd(),"/original_trees/trees/ClaDS2/", tree_files[i]))
  if(exists("tree")){
    label <- gsub(".Rdata","",tree_files[i])
    n_tip <- Ntip(tree)
    tmp_tiprates <- as.data.frame(matrix(nrow=n_tip, ncol=7))
    tip_speciation_rates <- c()
    tip_extinction_rates <- c()
    for(tip_index in 1:n_tip){
      tip_speciation_rates[tip_index] <- speciation_rates[which(tree$edge[,2]==tip_index)]
      tip_extinction_rates[tip_index] <- extinction_rates[which(tree$edge[,2]==tip_index)]
    }
    tmp_tiprates[,1] <- label
    tmp_tiprates[,2] <- tree$tip.label
    #library(RPANDA)
    #plot_ClaDS_phylo(tree, speciation_rates) # have to make sure the order is right
    tmp_tiprates[,3] <- tip_speciation_rates
    tmp_tiprates[,4] <- tip_extinction_rates
    tmp_tiprates[,5] <- tip_speciation_rates + tip_extinction_rates 
    tmp_tiprates[,6] <- tip_speciation_rates - tip_extinction_rates 
    tmp_tiprates[,7] <- tip_extinction_rates / tip_speciation_rates 
    true_clads_tiprates[[i]] <- tmp_tiprates
    rm("tree")
  }
}
true_clads_tiprates <- do.call(rbind, true_clads_tiprates)
colnames(true_clads_tiprates) <- c("treeName","tipName","lambdaTRUE","muTRUE",
                                   "netDivTRUE","turnoverTRUE","extinctionFractionTRUE")


# Add MiSSE rates:
library(Metrics)
library(tidyverse)
library(magrittr)
library(stats)
library(ggplot2)
library(progress)
library(viridis)
library(PMCMR)

treeTipMerge <- function(x) {
  x$treeTipString <- paste0(x$treeName, "_", x$tipName)
  return(x)
}
############################################################
# Load results
############################################################
clads.dir <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison"
misse_files <- list.files(paste0(clads.dir, "/MiSSE_ClaDS_comparison/misse"))

clads.dir <- "/Volumes/NO NAME/new_results/"
clads.results <- list.files(clads.dir)
clads.results <- clads.results[grep("_done_recon", clads.results)]
clads.results <- clads.results[grep("ClaDS2tree", clads.results)]

all_results <- data.frame()
model_average_results <- data.frame()
best_results <- data.frame()
dones <- misse_files

for (i in seq_along(dones)) {
  print(paste0("Loading ", i, " of ", length(dones)))
  try(load(paste0(clads.dir, "/MiSSE_ClaDS_comparison/misse/", dones[i])))
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
rates.ours$unique_string <- paste(rates.ours$treeName, "nturnover", rates.ours$nturnover, "neps", rates.ours$neps, "root", rates.ours$root_type, sep="_")
rates.ours.best <- best_results

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
rates.combined <- merge(rates.combined, true_clads_tiprates)

# get rid of NAs and weird values for now
#rates.cleaned <- rates.combined
#rates.cleaned <- rates.cleaned[(rates.cleaned$lambdaTRUE>0),] #yeah, not sure why there'd be no speciation in reality for a tree with >2 taxa


# Retriving ClaDS results:
clads.dir <- "/Users/thaisvasconcelos/Desktop/misse_mme_paper/missecomparison/ClaDScomparison"
clads.results <- list.files(paste0(clads.dir, "/MiSSE_ClaDS_comparison/clads"))

clads_tip_rates <- list()
for(result_index in 1:length(clads.results)) {
  load(paste0(clads.dir, "/MiSSE_ClaDS_comparison/clads/", clads.results[result_index]))
  tmp_tiprates <- data.frame(CladsOutput$tree$tip.label, CladsOutput$lambdatip_map, CladsOutput$eps_map)
  colnames(tmp_tiprates) <- c("tipName","lambdaClaDS", "extinctionFractionClaDS")
  tmp_tiprates$muClaDS <- tmp_tiprates$lambdaClaDS * tmp_tiprates$extinctionFractionClaDS  
  tmp_tiprates$turnoverClaDS <- tmp_tiprates$lambdaClaDS + tmp_tiprates$muClaDS
  tmp_tiprates$netDivClaDS <- tmp_tiprates$lambdaClaDS - tmp_tiprates$muClaDS
  tmp_tiprates$treeName <- gsub("_results.Rdata","", clads.results[result_index])
  clads_tip_rates[[result_index]] <- tmp_tiprates
  print(result_index)
}
all_clads_tip_rates <- do.call(rbind, clads_tip_rates)

#subset_results <- subset(rates.cleaned, rates.cleaned$treeName %in% tree_labels)
subset_results <- merge(rates.combined, all_clads_tip_rates, by=c("treeName","tipName"))



############################################################
# Organize results by type of tree, metric and parameter
############################################################
approaches <- c("ClaDS", "MiSSEbest", "MiSSEavg")
parameters <- c("mu", "lambda", "netDiv", "turnover", "extinctionFraction")
error.measurements <- c("RMSE","absoluteError.mean","absoluteError.median")
all_trees <- unique(subset_results$treeName)

results_all_trees <- matrix(nrow=0, ncol=4)
for (i in seq_along(all_trees)) {
  one_tree <- all_trees[i]
  tips_of_one_tree <- subset_results[subset_results$treeName==one_tree,]
  results_one_approach <- matrix(nrow=0, ncol=4)
  for (approach.index in seq_along(approaches)) {
    for (parameter.index in seq_along(parameters)) {
      truth <- tips_of_one_tree[,paste0(parameters[parameter.index],"TRUE")]
      estimate_name <- paste0(parameters[parameter.index], approaches[approach.index])
      if(estimate_name %in% colnames(tips_of_one_tree)) {
        estimate <- tips_of_one_tree[,estimate_name]
        RMSE.results <- Metrics::rmse(estimate, truth)
        absoluteError.mean.results <- mean(abs(estimate-truth))
        absoluteError.median.results <- median(abs(estimate-truth))
        tmp_errors <- c(RMSE.results, absoluteError.mean.results, absoluteError.median.results)
        results_one_approach <- rbind(results_one_approach, cbind(rep(one_tree, 3), c("RMSE", "absoluteError.mean", "absoluteError.median"), unname(estimate_name), round(unname(tmp_errors),5)))
        #names(tmp_errors) <- c(paste0(c("RMSE", "absoluteError.mean", "absoluteError.median"), "_",estimate_name))  
        #results_one_approach <- c(results_one_approach, tmp_errors)
      }
    }
  }
  results_all_trees <- rbind(results_all_trees, results_one_approach)
  print(i)
}

#rownames(results_all_trees) <- all_trees
results_all_trees <- as.data.frame(results_all_trees)
colnames(results_all_trees) <- c("tree_name","error","parameter","value")

results_all_trees$value <- as.numeric(results_all_trees$value)

types_of_trees <- types_of_trees <- unique(unlist(lapply(strsplit(rates.combined$treeTipString, "_"),"[[",1)))
length_tree_type <- c()
for(i in types_of_trees){
  #n <- length(grep(i, rownames(results_all_trees)))
  n <- length(grep(i, all_trees))
  length_tree_type <- c(length_tree_type, n)  
}
names(length_tree_type) <- types_of_trees



############################################################
# Organize results by type of tree, metric and parameter
############################################################
load("ClaDS-master/Simulations/BAMM_Cl2.Rdata")
BAMM_mse <- data_BAMM[,c("name","seed","MSE_ext","MSE","MSE_div")]



############################################################
# Plot results
############################################################

require(ggplot2)
library(gridExtra)

organize_table <- function(results_all_trees, tree_type, parameter, error.measurement){
  results_one_tree_type <- subset(results_all_trees, grepl(tree_type, results_all_trees$tree_name))
  results_one_tree_type_one_parameter <- subset(results_one_tree_type, !is.na(results_one_tree_type))
  results_one_tree_type_one_parameter <- results_one_tree_type_one_parameter[grep(parameter, results_one_tree_type_one_parameter$parameter),]
  results_one_tree_type_one_parameter <- results_one_tree_type_one_parameter[grep(error.measurement, results_one_tree_type_one_parameter$error),]
  #results_one_tree_type_one_parameter$error <- as.factor(results_one_tree_type_one_parameter$error)
  results_one_tree_type_one_parameter$parameter <- gsub(parameter,"", results_one_tree_type_one_parameter$parameter)
  
  results_one_tree_type_one_parameter$parameter <- factor(results_one_tree_type_one_parameter$parameter, levels=c("TB", "ND", "DR", "BAMM", "ClaDS","MiSSEbest", "MiSSEavg"))
  return(results_one_tree_type_one_parameter)
}

one_tree_type <- c("ClaDS2tree")
approaches <- c("ClaDS", "MiSSEbest", "MiSSEavg")
parameters <- c("mu", "lambda", "netDiv", "turnover", "extinctionFraction")
error.measurements <- c("RMSE","absoluteError.mean","absoluteError.median")


############################################################
# Define measure of error to plot
############################################################
one.error.measurement = "absoluteError.mean" # any of c("RMSE","absoluteError.mean","absoluteError.median")

{
  pal.name <- "Viridis"
  pal <- rev(hcl.colors(7, palette = pal.name, alpha = 0.75))
  
  pdf(paste0("../Supplementary_Material/",one.error.measurement, "_results_comparison_10.pdf"), width=15, height=20)
  #colnames(results_one_tree_type_one_parameter)
  lambda1 <- organize_table(results_all_trees, one_tree_type[1], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda1 <- ggplot(lambda1, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    #geom_violin(trim=FALSE)+
    theme_bw() +
    theme(legend.position = "none") + 
    #xlab("") +
    #ggtitle(paste0(one_tree_type[1], " (lambda)", " (N=",length_tree_type[1],")")) +
    ggtitle(paste0(one_tree_type[1], " (N=",length_tree_type[1],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    #geom_boxplot(width=0.1) +
    coord_flip(ylim = c(0.001,10)) #+
    #theme(axis.title.y=element_blank(),
    #      axis.text.y=element_blank(),
    #      axis.ticks.y=element_blank())
  
  
  mu1 <- organize_table(results_all_trees, one_tree_type[1], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu1 <- ggplot(mu1, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[1], " (mu)", " (N=",length_tree_type[1],")")) +
    ggtitle("") +
    scale_x_discrete("parameter",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    #coord_cartesian(ylim = c(-0.25,0.75)) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv1 <- organize_table(results_all_trees, one_tree_type[1], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv1 <- ggplot(netDiv1, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[1], " (netDiv)", " (N=",length_tree_type[1],")")) +
    ggtitle("") +
    scale_x_discrete("parameter",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ##
  turnover1 <- organize_table(results_all_trees, one_tree_type[1], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover1 <- ggplot(turnover1, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[1], " (turnover)", " (N=",length_tree_type[1],")")) +
    ggtitle("") +
    scale_x_discrete("parameter",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction1 <- organize_table(results_all_trees, one_tree_type[1], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction1 <- ggplot(extinctionFraction1, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[1], " (extinctionFraction)", " (N=",length_tree_type[1],")")) +
    ggtitle("") +
    scale_x_discrete("parameter",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  grid.arrange(plot_lambda1, plot_mu1, plot_netDiv1, plot_turnover1, plot_extinctionFraction1,
               plot_lambda2, plot_mu2, plot_netDiv2, plot_turnover2, plot_extinctionFraction2,
               plot_lambda3, plot_mu3, plot_netDiv3, plot_turnover3, plot_extinctionFraction3,
               plot_lambda4, plot_mu4, plot_netDiv4, plot_turnover4, plot_extinctionFraction4,
               plot_lambda5, plot_mu5, plot_netDiv5, plot_turnover5, plot_extinctionFraction5,
               plot_lambda6, plot_mu6, plot_netDiv6, plot_turnover6, plot_extinctionFraction6,
               plot_lambda7, plot_mu7, plot_netDiv7, plot_turnover7, plot_extinctionFraction7,
               plot_lambda8, plot_mu8, plot_netDiv8, plot_turnover8, plot_extinctionFraction8,
               ncol=5, nrow = 8)
  
  
  dev.off()
}



