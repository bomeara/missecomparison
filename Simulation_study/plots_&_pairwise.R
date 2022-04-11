############################################################
# (1) Boxplots of error measurements
# (2) Accuracy comparison using Kruskal-Wallis 
############################################################
# rm(list=ls())
library(data.table)
library(PMCMR)

# (1) Loading tables


subset_results_titlerabosky <- as.data.frame(fread("2022-03-29_preliminary_results_titlerabosky.csv"))
subset_results_malietetal <- as.data.frame(fread("2022-03-28_preliminary_results_malietetal.csv"))

subset_results_titlerabosky <- subset_results_titlerabosky[,-which(!colnames(subset_results_titlerabosky) %in% colnames(subset_results_malietetal))]
subset_results_titlerabosky <- subset_results_titlerabosky[order(match(colnames(subset_results_titlerabosky), colnames(subset_results_malietetal)))]

all_results <- rbind(subset_results_titlerabosky, subset_results_malietetal)

############################################################
# Organize results by type of tree, metric and parameter
############################################################
approaches <- c("TB", "ND", "DR", "BAMM", "ClaDS", "MiSSEbest", "MiSSEavg")
parameters <- c("mu", "lambda", "netDiv", "turnover", "extinctionFraction")
error.measurements <- c("RMSE","absoluteError.mean","absoluteError.median")
all_trees <- unique(all_results$treeName)

results_all_trees <- matrix(nrow=0, ncol=4)
for (i in seq_along(all_trees)) {
  one_tree <- all_trees[i]
  tips_of_one_tree <- all_results[all_results$treeName==one_tree,]
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
  cat(i, "\r")
}

#rownames(results_all_trees) <- all_trees
results_all_trees <- as.data.frame(results_all_trees)
colnames(results_all_trees) <- c("tree_name","error","parameter","value")
results_all_trees$value <- as.numeric(results_all_trees$value)

types_of_trees <- unique(unlist(lapply(strsplit(results_all_trees$tree_name, "_"),"[[",1)))
length_tree_type <- c()
for(i in types_of_trees){
  #n <- length(grep(i, rownames(results_all_trees)))
  n <- length(grep(i, all_trees))
  length_tree_type <- c(length_tree_type, n)  
}
names(length_tree_type) <- types_of_trees

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

#one_tree_type <- c("evolvingRates","fossilBAMM","lambdaConstantVariance","MitchellRabosky2016","MooreEtAl2016","netDivConstantVariance","Rabosky2014","RaboskyEtAl2017")
one_tree_type <- types_of_trees

approaches <- c("TB", "ND", "DR", "BAMM", "ClaDS", "MiSSEbest", "MiSSEavg")
parameters <- c("mu", "lambda", "netDiv", "turnover", "extinctionFraction")
error.measurements <- c("RMSE","absoluteError.mean","absoluteError.median")

############################################################
# Define measure of error to plot
############################################################
one.error.measurement = "absoluteError.mean" # any of c("RMSE","absoluteError.mean","absoluteError.median")

############################################################
# Make median+CI plots
############################################################


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
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
  
  
  lambda2 <- organize_table(results_all_trees, one_tree_type[2], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda2 <- ggplot(lambda2, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[2], " (lambda)", " (N=",length_tree_type[2],")")) +
    ggtitle(paste0(one_tree_type[2], " (N=",length_tree_type[2],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    coord_flip(ylim = c(0.001,10))  +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  lambda3 <- organize_table(results_all_trees, one_tree_type[3], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda3 <- ggplot(lambda3, aes(x=parameter, y=value, fill=parameter)) +
    #geom_boxplot() + 
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[3], " (lambda)", " (N=",length_tree_type[3],")")) +
    ggtitle(paste0(one_tree_type[3], " (N=",length_tree_type[3],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    coord_flip(ylim = c(0.001,10))  +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  lambda4 <- organize_table(results_all_trees, one_tree_type[4], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda4 <- ggplot(lambda4, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[4], " (lambda)", " (N=",length_tree_type[4],")")) +
    ggtitle(paste0(one_tree_type[4], " (N=",length_tree_type[4],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    coord_flip(ylim = c(0.001,10))  +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  lambda5 <- organize_table(results_all_trees, one_tree_type[5], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda5 <- ggplot(lambda5, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[5], " (lambda)", " (N=",length_tree_type[5],")")) +
    ggtitle(paste0(one_tree_type[5], " (N=",length_tree_type[5],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    coord_flip(ylim = c(0.001,10))  +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  lambda6 <- organize_table(results_all_trees, one_tree_type[6], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda6 <- ggplot(lambda6, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[6], " (lambda)", " (N=",length_tree_type[6],")")) +
    ggtitle(paste0(one_tree_type[6], " (N=",length_tree_type[6],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    coord_flip(ylim = c(0.001,10))  +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  lambda7 <- organize_table(results_all_trees, one_tree_type[7], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda7 <- ggplot(lambda7, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=7, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[7], " (lambda)", " (N=",length_tree_type[7],")")) +
    ggtitle(paste0(one_tree_type[7], " (N=",length_tree_type[7],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    coord_flip(ylim = c(0.001,10))  +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  lambda8 <- organize_table(results_all_trees, one_tree_type[8], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda8 <- ggplot(lambda8, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=8, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (lambda)", " (N=",length_tree_type[8],")")) +
    ggtitle(paste0(one_tree_type[8], " (N=",length_tree_type[8],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    coord_flip(ylim = c(0.001,10))  +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  lambda9 <- organize_table(results_all_trees, one_tree_type[9], parameters[parameters=="lambda"], error.measurements[error.measurements==one.error.measurement])
  plot_lambda9 <- ggplot(lambda9, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=8, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (lambda)", " (N=",length_tree_type[8],")")) +
    ggtitle(paste0(one_tree_type[9], " (N=",length_tree_type[9],")")) +
    theme(plot.title = element_text(size = 10)) +
    scale_x_discrete("",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal) +
    coord_flip(ylim = c(0.001,10))  +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
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
  
  ##
  mu2 <- organize_table(results_all_trees, one_tree_type[2], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu2 <- ggplot(mu2, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[2], " (mu)", " (N=",length_tree_type[2],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv2 <- organize_table(results_all_trees, one_tree_type[2], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv2 <- ggplot(netDiv2, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[2], " (netDiv)", " (N=",length_tree_type[2],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  turnover2 <- organize_table(results_all_trees, one_tree_type[2], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover2 <- ggplot(turnover2, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[2], " (turnover)", " (N=",length_tree_type[2],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction2 <- organize_table(results_all_trees, one_tree_type[2], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction2 <- ggplot(extinctionFraction2, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[2], " (extinctionFraction)", " (N=",length_tree_type[2],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ##
  
  mu3 <- organize_table(results_all_trees, one_tree_type[3], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu3 <- ggplot(mu3, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[3], " (mu)", " (N=",length_tree_type[3],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv3 <- organize_table(results_all_trees, one_tree_type[3], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv3 <- ggplot(netDiv3, aes(x=parameter, y=value, fill=parameter)) +
    #geom_boxplot() + 
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    ggtitle("") +
    xlab("") +
    #ggtitle(paste0(one_tree_type[3], " (netDiv)", " (N=",length_tree_type[3],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  turnover3 <- organize_table(results_all_trees, one_tree_type[3], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover3 <- ggplot(turnover3, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[3], " (turnover)", " (N=",length_tree_type[3],")")) +
    ggtitle("") +
    scale_x_discrete("parameter",  drop=FALSE) +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction3 <- organize_table(results_all_trees, one_tree_type[3], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction3 <- ggplot(extinctionFraction3, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[3], " (extinctionFraction)", " (N=",length_tree_type[3],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ##
  
  mu4 <- organize_table(results_all_trees, one_tree_type[4], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu4 <- ggplot(mu4, aes(x=parameter, y=value, fill=parameter)) +
    geom_boxplot(outlier.size = 0.8)  +
    scale_y_continuous(trans='log10') +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[4], " (mu)", " (N=",length_tree_type[4],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv4 <- organize_table(results_all_trees, one_tree_type[4], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv4 <- ggplot(netDiv4, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[4], " (netDiv)", " (N=",length_tree_type[4],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  turnover4 <- organize_table(results_all_trees, one_tree_type[4], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover4 <- ggplot(turnover4, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[4], " (turnover)", " (N=",length_tree_type[4],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction4 <- organize_table(results_all_trees, one_tree_type[4], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction4 <- ggplot(extinctionFraction4, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[4], " (extinctionFraction)", " (N=",length_tree_type[4],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ##
  
  mu5 <- organize_table(results_all_trees, one_tree_type[5], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu5 <- ggplot(mu5, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[5], " (mu)", " (N=",length_tree_type[5],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv5 <- organize_table(results_all_trees, one_tree_type[5], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv5 <- ggplot(netDiv5, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[5], " (netDiv)", " (N=",length_tree_type[5],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  turnover5 <- organize_table(results_all_trees, one_tree_type[5], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover5 <- ggplot(turnover5, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[5], " (turnover)", " (N=",length_tree_type[5],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction5 <- organize_table(results_all_trees, one_tree_type[5], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction5 <- ggplot(extinctionFraction5, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[5], " (extinctionFraction)", " (N=",length_tree_type[5],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ##
  
  mu6 <- organize_table(results_all_trees, one_tree_type[6], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu6 <- ggplot(mu6, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[6], " (mu)", " (N=",length_tree_type[6],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv6 <- organize_table(results_all_trees, one_tree_type[6], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv6 <- ggplot(netDiv6, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[6], " (netDiv)", " (N=",length_tree_type[6],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  turnover6 <- organize_table(results_all_trees, one_tree_type[6], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover6 <- ggplot(turnover6, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[6], " (turnover)", " (N=",length_tree_type[6],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction6 <- organize_table(results_all_trees, one_tree_type[6], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction6 <- ggplot(extinctionFraction6, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[6], " (extinctionFraction)", " (N=",length_tree_type[6],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ##
  
  mu7 <- organize_table(results_all_trees, one_tree_type[7], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu7 <- ggplot(mu7, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[7], " (mu)", " (N=",length_tree_type[7],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv7 <- organize_table(results_all_trees, one_tree_type[7], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv7 <- ggplot(netDiv7, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[7], " (netDiv)", " (N=",length_tree_type[7],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  turnover7 <- organize_table(results_all_trees, one_tree_type[7], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover7 <- ggplot(turnover7, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[7], " (turnover)", " (N=",length_tree_type[7],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction7 <- organize_table(results_all_trees, one_tree_type[7], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction7 <- ggplot(extinctionFraction7, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[7], " (extinctionFraction)", " (N=",length_tree_type[7],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ##
  
  mu8 <- organize_table(results_all_trees, one_tree_type[8], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu8 <- ggplot(mu8, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (mu)", " (N=",length_tree_type[8],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv8 <- organize_table(results_all_trees, one_tree_type[8], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv8 <- ggplot(netDiv8, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (netDiv)", " (N=",length_tree_type[8],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  turnover8 <- organize_table(results_all_trees, one_tree_type[8], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover8 <- ggplot(turnover8, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (turnover)", " (N=",length_tree_type[8],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction8 <- organize_table(results_all_trees, one_tree_type[8], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction8 <- ggplot(extinctionFraction8, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (extinctionFraction)", " (N=",length_tree_type[8],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  ##
  
  mu9 <- organize_table(results_all_trees, one_tree_type[9], parameters[parameters=="mu"], error.measurements[error.measurements==one.error.measurement])
  plot_mu9 <- ggplot(mu9, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (mu)", " (N=",length_tree_type[8],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  netDiv9 <- organize_table(results_all_trees, one_tree_type[9], parameters[parameters=="netDiv"], error.measurements[error.measurements==one.error.measurement])
  plot_netDiv9 <- ggplot(netDiv9, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (netDiv)", " (N=",length_tree_type[8],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  turnover9 <- organize_table(results_all_trees, one_tree_type[9], parameters[parameters=="turnover"], error.measurements[error.measurements==one.error.measurement])
  plot_turnover9 <- ggplot(turnover9, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (turnover)", " (N=",length_tree_type[8],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
    #scale_fill_brewer(palette=pal) +
    scale_fill_manual(values=pal[4:7]) +
    coord_flip(ylim = c(0.001,10)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  extinctionFraction9 <- organize_table(results_all_trees, one_tree_type[9], parameters[parameters=="extinctionFraction"], error.measurements[error.measurements==one.error.measurement])
  plot_extinctionFraction9 <- ggplot(extinctionFraction9, aes(x=parameter, y=value, fill=parameter)) +
    scale_y_continuous(trans='log10') +
    geom_boxplot(outlier.size = 0.8)  +
    #geom_boxplot() + 
    theme_bw() +
    theme(legend.position = "none") + 
    #annotate(geom="text", x=6, y=7, size=5, hjust=0, label=one.error.measurement) + 
    xlab("") +
    #ggtitle(paste0(one_tree_type[8], " (extinctionFraction)", " (N=",length_tree_type[8],")")) +
    scale_x_discrete("parameter",  drop=FALSE) +
    ggtitle("") +
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
               plot_lambda9, plot_mu9, plot_netDiv9, plot_turnover9, plot_extinctionFraction9,
               ncol=5, nrow = 9)
  dev.off()
}


######################################################
# (2) Kruskal-wallis
######################################################
all_tables <- list(lambda1, mu1, netDiv1, turnover1, extinctionFraction1,
                   lambda2, mu2, netDiv2, turnover2, extinctionFraction2,
                   lambda3, mu3, netDiv3, turnover3, extinctionFraction3,
                   lambda4, mu4, netDiv4, turnover4, extinctionFraction4,
                   lambda5, mu5, netDiv5, turnover5, extinctionFraction5,
                   lambda6, mu6, netDiv6, turnover6, extinctionFraction6,
                   lambda7, mu7, netDiv7, turnover7, extinctionFraction7,
                   lambda8, mu8, netDiv8, turnover8, extinctionFraction8,
                   lambda9, mu9, netDiv9, turnover9, extinctionFraction9)

rate_names <- rep(c("lambda","mu","netDiv","turnover","extinctionFraction"),9)

sink(paste0("../Supplementary_Material/",one.error.measurement,"_kruskal_wallis_results.txt"))
for(i in 1:length(all_tables)) {
  one_table <- all_tables[[i]]
  tree_name <- gsub("\\_.*","", one_table$tree_name[1])
  result <- posthoc.kruskal.conover.test(one_table$value, one_table$parameter)
  print(tree_name)
  print(rate_names[i])
  print(result)
  print("------------------------")
}
sink()
