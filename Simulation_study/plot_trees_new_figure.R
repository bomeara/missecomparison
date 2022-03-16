
# Some plots for new figure:

# Random plots
# rm(list=ls())
# Selecting random tree
library(ape)
# setwd("~/Desktop/misse_mme_paper/missecomparison/Simulation_study/")

tree_files_tr <- list.files("TitleRaboskyRuns/trees",full.names = T)
tree_files_m <- list.files("MalietEtAlRuns/MalietEtAl_ClaDS2_trees", full.names = T)[grep(".tre$", list.files("MalietEtAlRuns/MalietEtAl_ClaDS2_trees/"))]

all_trees <- lapply(c(tree_files_m, tree_files_tr), read.tree)
names(all_trees) <- gsub(".tre$","",unlist(lapply(strsplit(c(tree_files_m, tree_files_tr),"/"), "[[", 3)))

#lapply(all_trees[grep("ClaDS", names(all_trees))], Ntip)


library(data.table)
subset_results_titlerabosky <- as.data.frame(fread("Simulation_study2022-03-07_preliminary_results_titlerabosky.csv"))
subset_results_malietetal <- as.data.frame(fread("Simulation_study2022-03-08_preliminary_results_malietetal.csv"))


one_set <- data.frame()
while(nrow(one_set)==0) {
  one_tree <- all_trees[sample(names(all_trees), 1)]
  one_set <- subset(subset_results_malietetal, subset_results_malietetal$treeName == names(one_tree))
}

#all_trees <- all_trees[c("ClaDS2tree_100_5_2","ClaDS2tree_100_5_4","fossilBAMM_101","fossilBAMM_31","fossilBAMM_71",
#                         "lambdaConstantVariance_201","lambdaConstantVariance_21","lambdaConstantVariance_281",
#                         "lambdaConstantVariance_311","MooreEtAl2016_1","MooreEtAl2016_91","Rabosky2014_DD_k1_1",
#                         "Rabosky2014_DD_k1_211","Rabosky2014_DD_k1_221","Rabosky2014_DD_k1_251",
#                         "Rabosky2014_DD_k1_291","Rabosky2014_DD_k2_231", "RaboskyEtAl2017_71",
#                         "evolvingRates_1131", "evolvingRates_111","evolvingRates_1041","evolvingRates_101")]
#all_trees <- all_trees[c("ClaDS2tree_100_2_2","evolvingRates_921","Rabosky2014_DD_k1_261","lambdaConstantVariance_281")]

#clads_trees <- names(all_trees)[grep("ClaDS2tree", names(all_trees))]

all_trees <- all_trees[c("ClaDS2tree_100_5_4","evolvingRates_711","lambdaConstantVariance_281","Rabosky2014_DD_k1_1")]


pdf("lambda_plot.pdf", width=4, height=5)
for(i in 1:length(all_trees)) {
  one_tree <- all_trees[[i]]
  one_label <- names(all_trees)[i]
  if(any(subset_results_malietetal$treeName %in% one_label)){
    one_set <- subset(subset_results_malietetal, subset_results_malietetal$treeName == one_label)
  } else 
    if(any(subset_results_titlerabosky$treeName %in% one_label)) {
      one_set <- subset(subset_results_titlerabosky, subset_results_titlerabosky$treeName == one_label)
    } else {next}
  
  if(nrow(one_set)>0) {
    
    tree_pruned <- one_tree
    # Small corrections so that tips appear in the right order in the plot
    is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
    ordered_tips <- tree_pruned$edge[is_tip, 2]
    right_order <- as.character(tree_pruned$tip.label[ordered_tips])
    
    # Organizing so tip rates are in the same order as tips of the tree
    cleaned_table <- one_set[match(as.character(right_order), as.character(one_set$tipName)),]
    
    if(length(which(is.na(cleaned_table$lambdaTRUE)))==0) {
      
      # Plotting
      color_breaks = 4 
      layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
      layout(mat = layout.matrix,widths = c(2,1,1,1))
      par(mar=c(5,0.5,3,0.5))
      # Tree
      plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
      title(main=one_label)
      ape::axisPhylo()
      
      palette_name = "Inferno"
      
      # Load rates
      lambdaTRUE <- as.numeric(cleaned_table$lambdaTRUE)
      lambdaBAMM <- as.numeric(cleaned_table$lambdaBAMM)
      lambdaClaDS <- as.numeric(cleaned_table$lambdaClaDS)
      lambdaMiSSEavg <- as.numeric(cleaned_table$lambdaMiSSEavg)
      
      set_true <- 1:length(lambdaTRUE)
      set_rate_one <- (tail(set_true,1)+1):(length(set_true)*2)
      set_rate_two <- (tail(set_rate_one,1)+1):(length(set_true)*3)
      set_rate_three <- (tail(set_rate_two,1)+1):(length(set_true)*4)
      
      rounded_rates <- round(c(lambdaTRUE, lambdaBAMM, lambdaClaDS, lambdaMiSSEavg), color_breaks)
      pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name, alpha = 0.75)
      pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
      
      
      #segments(lambdaTRUE, 1:length(lambdaTRUE), lambdaTRUE[1:length(lambdaTRUE)], 1:length(lambdaTRUE), col= pal[set_true],lwd = 1)
      
      # Lambda BAMM
      x <- 1:length(lambdaTRUE)
      plot(lambdaTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(lambdaTRUE), max(lambdaTRUE)))
      segments(lambdaBAMM, 1:length(lambdaBAMM), lambdaTRUE[1:length(lambdaTRUE)], 1:length(lambdaTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(lambdaBAMM)
      plot(lambdaBAMM, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="BAMM lambda", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_one]) #xlim=range(c(min(lambdaBAMM), max(lambdaBAMM)))
      
      # Lambda ClaDS
      x <- 1:length(lambdaTRUE)
      plot(lambdaTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(lambdaTRUE), max(lambdaTRUE)))
      segments(lambdaClaDS, 1:length(lambdaClaDS), lambdaTRUE[1:length(lambdaTRUE)], 1:length(lambdaTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(lambdaClaDS)
      plot(lambdaClaDS, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="ClaDS lambda", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_two]) #xlim=range(c(min(lambdaBAMM), max(lambdaBAMM)))
      
      # Lambda MiSSE
      x <- 1:length(lambdaTRUE)
      plot(lambdaTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(lambdaTRUE), max(lambdaTRUE)))
      segments(lambdaMiSSEavg, 1:length(lambdaMiSSEavg), lambdaTRUE[1:length(lambdaTRUE)], 1:length(lambdaTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(lambdaMiSSEavg)
      plot(lambdaMiSSEavg, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="MiSSE lambda", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_three]) #xlim=range(c(min(lambdaBAMM), max(lambdaBAMM)))
      

    }
  }
}
dev.off()

##################
##################

pdf("mu_plot.pdf", width=4, height=5)
for(i in 1:length(all_trees)) {
  one_tree <- all_trees[[i]]
  one_label <- names(all_trees)[i]
  if(any(subset_results_malietetal$treeName %in% one_label)){
    one_set <- subset(subset_results_malietetal, subset_results_malietetal$treeName == one_label)
  } else 
    if(any(subset_results_titlerabosky$treeName %in% one_label)) {
      one_set <- subset(subset_results_titlerabosky, subset_results_titlerabosky$treeName == one_label)
    } else {next}
  
  if(nrow(one_set)>0) {
    
    tree_pruned <- one_tree
    # Small corrections so that tips appear in the right order in the plot
    is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
    ordered_tips <- tree_pruned$edge[is_tip, 2]
    right_order <- as.character(tree_pruned$tip.label[ordered_tips])
    
    # Organizing so tip rates are in the same order as tips of the tree
    cleaned_table <- one_set[match(as.character(right_order), as.character(one_set$tipName)),]
    
    if(length(which(is.na(cleaned_table$muTRUE)))==0) {
      
      # Plotting
      color_breaks = 4 
      layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
      layout(mat = layout.matrix,widths = c(2,1,1,1))
      par(mar=c(5,0.5,3,0.5))
      # Tree
      plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
      title(main=one_label)
      ape::axisPhylo()
      
      palette_name = "Inferno"
      
      # Load rates
      muTRUE <- as.numeric(cleaned_table$muTRUE)
      muBAMM <- as.numeric(cleaned_table$muBAMM)
      muClaDS <- as.numeric(cleaned_table$muClaDS)
      muMiSSEavg <- as.numeric(cleaned_table$muMiSSEavg)
      
      set_true <- 1:length(muTRUE)
      set_rate_one <- (tail(set_true,1)+1):(length(set_true)*2)
      set_rate_two <- (tail(set_rate_one,1)+1):(length(set_true)*3)
      set_rate_three <- (tail(set_rate_two,1)+1):(length(set_true)*4)
      
      rounded_rates <- round(c(muTRUE, muBAMM, muClaDS, muMiSSEavg), color_breaks)
      pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name, alpha = 0.75)
      pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
      
      
      #segments(muTRUE, 1:length(muTRUE), muTRUE[1:length(muTRUE)], 1:length(muTRUE), col= pal[set_true],lwd = 1)
      
      # mu BAMM
      x <- 1:length(muTRUE)
      plot(muTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(muTRUE), max(muTRUE)))
      segments(muBAMM, 1:length(muBAMM), muTRUE[1:length(muTRUE)], 1:length(muTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(muBAMM)
      plot(muBAMM, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="BAMM mu", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_one]) #xlim=range(c(min(muBAMM), max(muBAMM)))
      
      # mu ClaDS
      x <- 1:length(muTRUE)
      plot(muTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(muTRUE), max(muTRUE)))
      segments(muClaDS, 1:length(muClaDS), muTRUE[1:length(muTRUE)], 1:length(muTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(muClaDS)
      plot(muClaDS, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="ClaDS mu", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_two]) #xlim=range(c(min(muBAMM), max(muBAMM)))
      
      # mu MiSSE
      x <- 1:length(muTRUE)
      plot(muTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(muTRUE), max(muTRUE)))
      segments(muMiSSEavg, 1:length(muMiSSEavg), muTRUE[1:length(muTRUE)], 1:length(muTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(muMiSSEavg)
      plot(muMiSSEavg, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="MiSSE mu", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_three]) #xlim=range(c(min(muBAMM), max(muBAMM)))
      
      
    }
  }
}
dev.off()

#######
#######
pdf("turnover_plot.pdf", width=4, height=5)
for(i in 1:length(all_trees)) {
  one_tree <- all_trees[[i]]
  one_label <- names(all_trees)[i]
  if(any(subset_results_malietetal$treeName %in% one_label)){
    one_set <- subset(subset_results_malietetal, subset_results_malietetal$treeName == one_label)
  } else 
    if(any(subset_results_titlerabosky$treeName %in% one_label)) {
      one_set <- subset(subset_results_titlerabosky, subset_results_titlerabosky$treeName == one_label)
    } else {next}
  
  if(nrow(one_set)>0) {
    
    tree_pruned <- one_tree
    # Small corrections so that tips appear in the right order in the plot
    is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
    ordered_tips <- tree_pruned$edge[is_tip, 2]
    right_order <- as.character(tree_pruned$tip.label[ordered_tips])
    
    # Organizing so tip rates are in the same order as tips of the tree
    cleaned_table <- one_set[match(as.character(right_order), as.character(one_set$tipName)),]
    
    if(length(which(is.na(cleaned_table$turnoverTRUE)))==0) {
      
      # Plotting
      color_breaks = 4 
      layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
      layout(mat = layout.matrix,widths = c(2,1,1,1))
      par(mar=c(5,0.5,3,0.5))
      # Tree
      plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
      title(main=one_label)
      ape::axisPhylo()
      
      palette_name = "Inferno"
      
      # Load rates
      turnoverTRUE <- as.numeric(cleaned_table$turnoverTRUE)
      turnoverBAMM <- as.numeric(cleaned_table$turnoverBAMM)
      turnoverClaDS <- as.numeric(cleaned_table$turnoverClaDS)
      turnoverMiSSEavg <- as.numeric(cleaned_table$turnoverMiSSEavg)
      
      set_true <- 1:length(turnoverTRUE)
      set_rate_one <- (tail(set_true,1)+1):(length(set_true)*2)
      set_rate_two <- (tail(set_rate_one,1)+1):(length(set_true)*3)
      set_rate_three <- (tail(set_rate_two,1)+1):(length(set_true)*4)
      
      rounded_rates <- round(c(turnoverTRUE, turnoverBAMM, turnoverClaDS, turnoverMiSSEavg), color_breaks)
      pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name, alpha = 0.75)
      pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
      
      
      #segments(turnoverTRUE, 1:length(turnoverTRUE), turnoverTRUE[1:length(turnoverTRUE)], 1:length(turnoverTRUE), col= pal[set_true],lwd = 1)
      
      # turnover BAMM
      x <- 1:length(turnoverTRUE)
      plot(turnoverTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(turnoverTRUE), max(turnoverTRUE)))
      segments(turnoverBAMM, 1:length(turnoverBAMM), turnoverTRUE[1:length(turnoverTRUE)], 1:length(turnoverTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(turnoverBAMM)
      plot(turnoverBAMM, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="BAMM turnover", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_one]) #xlim=range(c(min(turnoverBAMM), max(turnoverBAMM)))
      
      # turnover ClaDS
      x <- 1:length(turnoverTRUE)
      plot(turnoverTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(turnoverTRUE), max(turnoverTRUE)))
      segments(turnoverClaDS, 1:length(turnoverClaDS), turnoverTRUE[1:length(turnoverTRUE)], 1:length(turnoverTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(turnoverClaDS)
      plot(turnoverClaDS, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="ClaDS turnover", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_two]) #xlim=range(c(min(turnoverBAMM), max(turnoverBAMM)))
      
      # turnover MiSSE
      x <- 1:length(turnoverTRUE)
      plot(turnoverTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(turnoverTRUE), max(turnoverTRUE)))
      segments(turnoverMiSSEavg, 1:length(turnoverMiSSEavg), turnoverTRUE[1:length(turnoverTRUE)], 1:length(turnoverTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(turnoverMiSSEavg)
      plot(turnoverMiSSEavg, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="MiSSE turnover", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_three]) #xlim=range(c(min(turnoverBAMM), max(turnoverBAMM)))
      
      
    }
  }
}
dev.off()


#######
#######
pdf("netDiv_plot.pdf", width=4, height=5)
for(i in 1:length(all_trees)) {
  one_tree <- all_trees[[i]]
  one_label <- names(all_trees)[i]
  if(any(subset_results_malietetal$treeName %in% one_label)){
    one_set <- subset(subset_results_malietetal, subset_results_malietetal$treeName == one_label)
  } else 
    if(any(subset_results_titlerabosky$treeName %in% one_label)) {
      one_set <- subset(subset_results_titlerabosky, subset_results_titlerabosky$treeName == one_label)
    } else {next}
  
  if(nrow(one_set)>0) {
    
    tree_pruned <- one_tree
    # Small corrections so that tips appear in the right order in the plot
    is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
    ordered_tips <- tree_pruned$edge[is_tip, 2]
    right_order <- as.character(tree_pruned$tip.label[ordered_tips])
    
    # Organizing so tip rates are in the same order as tips of the tree
    cleaned_table <- one_set[match(as.character(right_order), as.character(one_set$tipName)),]
    
    if(length(which(is.na(cleaned_table$netDivTRUE)))==0) {
      
      # Plotting
      color_breaks = 4 
      layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
      layout(mat = layout.matrix,widths = c(2,1,1,1))
      par(mar=c(5,0.5,3,0.5))
      # Tree
      plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
      title(main=one_label)
      ape::axisPhylo()
      
      palette_name = "Inferno"
      
      # Load rates
      netDivTRUE <- as.numeric(cleaned_table$netDivTRUE)
      netDivBAMM <- as.numeric(cleaned_table$netDivBAMM)
      netDivClaDS <- as.numeric(cleaned_table$netDivClaDS)
      netDivMiSSEavg <- as.numeric(cleaned_table$netDivMiSSEavg)
      
      set_true <- 1:length(netDivTRUE)
      set_rate_one <- (tail(set_true,1)+1):(length(set_true)*2)
      set_rate_two <- (tail(set_rate_one,1)+1):(length(set_true)*3)
      set_rate_three <- (tail(set_rate_two,1)+1):(length(set_true)*4)
      
      rounded_rates <- round(c(netDivTRUE, netDivBAMM, netDivClaDS, netDivMiSSEavg), color_breaks)
      pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name, alpha = 0.75)
      pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
      
      
      #segments(netDivTRUE, 1:length(netDivTRUE), netDivTRUE[1:length(netDivTRUE)], 1:length(netDivTRUE), col= pal[set_true],lwd = 1)
      
      # netDiv BAMM
      x <- 1:length(netDivTRUE)
      plot(netDivTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(netDivTRUE), max(netDivTRUE)))
      segments(netDivBAMM, 1:length(netDivBAMM), netDivTRUE[1:length(netDivTRUE)], 1:length(netDivTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(netDivBAMM)
      plot(netDivBAMM, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="BAMM netDiv", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_one]) #xlim=range(c(min(netDivBAMM), max(netDivBAMM)))
      
      # netDiv ClaDS
      x <- 1:length(netDivTRUE)
      plot(netDivTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(netDivTRUE), max(netDivTRUE)))
      segments(netDivClaDS, 1:length(netDivClaDS), netDivTRUE[1:length(netDivTRUE)], 1:length(netDivTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(netDivClaDS)
      plot(netDivClaDS, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="ClaDS netDiv", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_two]) #xlim=range(c(min(netDivBAMM), max(netDivBAMM)))
      
      # netDiv MiSSE
      x <- 1:length(netDivTRUE)
      plot(netDivTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(netDivTRUE), max(netDivTRUE)))
      segments(netDivMiSSEavg, 1:length(netDivMiSSEavg), netDivTRUE[1:length(netDivTRUE)], 1:length(netDivTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(netDivMiSSEavg)
      plot(netDivMiSSEavg, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="MiSSE netDiv", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_three]) #xlim=range(c(min(netDivBAMM), max(netDivBAMM)))
      
      
    }
  }
}
dev.off()

#######
pdf("extinctionFraction_plot.pdf", width=4, height=5)
for(i in 1:length(all_trees)) {
  one_tree <- all_trees[[i]]
  one_label <- names(all_trees)[i]
  if(any(subset_results_malietetal$treeName %in% one_label)){
    one_set <- subset(subset_results_malietetal, subset_results_malietetal$treeName == one_label)
  } else 
    if(any(subset_results_titlerabosky$treeName %in% one_label)) {
      one_set <- subset(subset_results_titlerabosky, subset_results_titlerabosky$treeName == one_label)
    } else {next}
  
  if(nrow(one_set)>0) {
    
    tree_pruned <- one_tree
    # Small corrections so that tips appear in the right order in the plot
    is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
    ordered_tips <- tree_pruned$edge[is_tip, 2]
    right_order <- as.character(tree_pruned$tip.label[ordered_tips])
    
    # Organizing so tip rates are in the same order as tips of the tree
    cleaned_table <- one_set[match(as.character(right_order), as.character(one_set$tipName)),]
    
    if(length(which(is.na(cleaned_table$extinctionFractionTRUE)))==0) {
      
      # Plotting
      color_breaks = 4 
      layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
      layout(mat = layout.matrix,widths = c(2,1,1,1))
      par(mar=c(5,0.5,3,0.5))
      # Tree
      plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
      title(main=one_label)
      ape::axisPhylo()
      
      palette_name = "Inferno"
      
      # Load rates
      extinctionFractionTRUE <- as.numeric(cleaned_table$extinctionFractionTRUE)
      extinctionFractionBAMM <- as.numeric(cleaned_table$extinctionFractionBAMM)
      extinctionFractionClaDS <- as.numeric(cleaned_table$extinctionFractionClaDS)
      extinctionFractionMiSSEavg <- as.numeric(cleaned_table$extinctionFractionMiSSEavg)
      
      set_true <- 1:length(extinctionFractionTRUE)
      set_rate_one <- (tail(set_true,1)+1):(length(set_true)*2)
      set_rate_two <- (tail(set_rate_one,1)+1):(length(set_true)*3)
      set_rate_three <- (tail(set_rate_two,1)+1):(length(set_true)*4)
      
      rounded_rates <- round(c(extinctionFractionTRUE, extinctionFractionBAMM, extinctionFractionClaDS, extinctionFractionMiSSEavg), color_breaks)
      pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name, alpha = 0.75)
      pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
      
      
      #segments(extinctionFractionTRUE, 1:length(extinctionFractionTRUE), extinctionFractionTRUE[1:length(extinctionFractionTRUE)], 1:length(extinctionFractionTRUE), col= pal[set_true],lwd = 1)
      
      # extinctionFraction BAMM
      x <- 1:length(extinctionFractionTRUE)
      plot(extinctionFractionTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(extinctionFractionTRUE), max(extinctionFractionTRUE)))
      segments(extinctionFractionBAMM, 1:length(extinctionFractionBAMM), extinctionFractionTRUE[1:length(extinctionFractionTRUE)], 1:length(extinctionFractionTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(extinctionFractionBAMM)
      plot(extinctionFractionBAMM, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="BAMM extinctionFraction", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_one]) #xlim=range(c(min(extinctionFractionBAMM), max(extinctionFractionBAMM)))
      
      # extinctionFraction ClaDS
      x <- 1:length(extinctionFractionTRUE)
      plot(extinctionFractionTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(extinctionFractionTRUE), max(extinctionFractionTRUE)))
      segments(extinctionFractionClaDS, 1:length(extinctionFractionClaDS), extinctionFractionTRUE[1:length(extinctionFractionTRUE)], 1:length(extinctionFractionTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(extinctionFractionClaDS)
      plot(extinctionFractionClaDS, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="ClaDS extinctionFraction", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_two]) #xlim=range(c(min(extinctionFractionBAMM), max(extinctionFractionBAMM)))
      
      # extinctionFraction MiSSE
      x <- 1:length(extinctionFractionTRUE)
      plot(extinctionFractionTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="", ylab="", frame.plot=T, cex=0.75, col="lightgrey") #xlim=range(c(min(extinctionFractionTRUE), max(extinctionFractionTRUE)))
      segments(extinctionFractionMiSSEavg, 1:length(extinctionFractionMiSSEavg), extinctionFractionTRUE[1:length(extinctionFractionTRUE)], 1:length(extinctionFractionTRUE), col="lightgrey",lwd = 1)
      par(new=TRUE)
      x <- 1:length(extinctionFractionMiSSEavg)
      plot(extinctionFractionMiSSEavg, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="MiSSE extinctionFraction", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_three]) #xlim=range(c(min(extinctionFractionBAMM), max(extinctionFractionBAMM)))
      
      
    }
  }
}
dev.off()