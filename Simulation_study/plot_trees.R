# Random plots
# rm(list=ls())
# Selecting random tree
library(ape)
# setwd("~/Desktop/misse_mme_paper/missecomparison/Simulation_study/")

tree_files_tr <- list.files("TitleRaboskyRuns/trees",full.names = T)
tree_files_m <- list.files("MalietEtAlRuns/MalietEtAl_ClaDS2_trees", full.names = T)[grep(".tre$", list.files("MalietEtAlRuns/MalietEtAl_ClaDS2_trees/"))]

all_trees <- lapply(c(tree_files_m, tree_files_tr), read.tree)
names(all_trees) <- gsub(".tre$","",unlist(lapply(strsplit(c(tree_files_m, tree_files_tr),"/"), "[[", 3)))


library(data.table)
subset_results_titlerabosky <- as.data.frame(fread("Simulation_study2022-03-07_preliminary_results_titlerabosky.csv"))
subset_results_malietetal <- as.data.frame(fread("Simulation_study2022-03-07_preliminary_results_malietetal.csv"))


one_set <- data.frame()
while(nrow(one_set)==0) {
  one_tree <- all_trees[sample(names(all_trees), 1)]
  one_set <- subset(subset_results_malietetal, subset_results_malietetal$treeName == names(one_tree))
}


#clads_trees <- names(all_trees)[grep("ClaDS2tree", names(all_trees))]
pdf("lambda_plots.pdf", width=6, height=5)
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
  layout.matrix <- matrix(c(1,2,3,4,5), nrow = 1, ncol = 5)
  layout(mat = layout.matrix,widths = c(2,1,1,1,1))
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
  
  
  x <- 1:length(lambdaTRUE)
  plot(lambdaTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
       pch=19, yaxt = "n", xlab="true lambda", ylab="", frame.plot=T, cex=0.75, col=pal[set_true]) #xlim=range(c(min(lambdaTRUE), max(lambdaTRUE)))
  segments(min(lambdaTRUE), 1:length(lambdaTRUE), lambdaTRUE[1:length(lambdaTRUE)], 1:length(lambdaTRUE), col= pal[set_true],lwd = 1)
  title(main=Ntip(tree_pruned))
  
  # Lambda BAMM
  x <- 1:length(lambdaBAMM)
  plot(lambdaBAMM, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
       pch=19, yaxt = "n", xlab="BAMM lambda", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_one]) #xlim=range(c(min(lambdaBAMM), max(lambdaBAMM)))
  segments(min(lambdaBAMM), 1:length(lambdaBAMM), lambdaBAMM[1:length(lambdaBAMM)], 1:length(lambdaBAMM), col= pal[set_rate_one],lwd = 1)
  
  # Lambda ClaDS
  x <- 1:length(lambdaClaDS)
  plot(lambdaClaDS, x,  lwd = 1,xlim=c(min(rounded_rates), max(rounded_rates)), 
       pch=19, yaxt = "n", xlab="ClaDS lambda", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_two]) #xlim=range(c(min(lambdaClaDS), max(lambdaClaDS)))
  segments(min(lambdaClaDS), 1:length(lambdaClaDS), lambdaClaDS[1:length(lambdaClaDS)], 1:length(lambdaClaDS), col= pal[set_rate_two],lwd = 1)
  
  # Lambda MiSSE
  x <- 1:length(lambdaMiSSEavg)
  plot(lambdaMiSSEavg, x,  lwd = 1,xlim=c(min(rounded_rates), max(rounded_rates)), 
       pch=19, yaxt = "n", xlab="MiSSE lambda", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_three]) #xlim=range(c(min(lambdaTRUE), max(lambdaTRUE)))
  segments(min(lambdaMiSSEavg), 1:length(lambdaMiSSEavg), lambdaMiSSEavg[1:length(lambdaMiSSEavg)], 1:length(lambdaMiSSEavg), col= pal[set_rate_three],lwd = 1)
  }
}
}
dev.off()

##################
##################
pdf("mu_plots.pdf", width=6, height=5)
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
      layout.matrix <- matrix(c(1,2,3,4,5), nrow = 1, ncol = 5)
      layout(mat = layout.matrix,widths = c(2,1,1,1,1))
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
      
      
      x <- 1:length(muTRUE)
      plot(muTRUE, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="true mu", ylab="", frame.plot=T, cex=0.75, col=pal[set_true]) #xlim=range(c(min(muTRUE), max(muTRUE)))
      segments(min(muTRUE), 1:length(muTRUE), muTRUE[1:length(muTRUE)], 1:length(muTRUE), col= pal[set_true],lwd = 1)
      
      # mu BAMM
      x <- 1:length(muBAMM)
      plot(muBAMM, x,  lwd = 1, xlim=c(min(rounded_rates), max(rounded_rates)),
           pch=19, yaxt = "n", xlab="BAMM mu", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_one]) #xlim=range(c(min(muBAMM), max(muBAMM)))
      segments(min(muBAMM), 1:length(muBAMM), muBAMM[1:length(muBAMM)], 1:length(muBAMM), col= pal[set_rate_one],lwd = 1)
      
      # mu ClaDS
      x <- 1:length(muClaDS)
      plot(muClaDS, x,  lwd = 1,xlim=c(min(rounded_rates), max(rounded_rates)), 
           pch=19, yaxt = "n", xlab="ClaDS mu", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_two]) #xlim=range(c(min(muClaDS), max(muClaDS)))
      segments(min(muClaDS), 1:length(muClaDS), muClaDS[1:length(muClaDS)], 1:length(muClaDS), col= pal[set_rate_two],lwd = 1)
      
      # mu MiSSE
      x <- 1:length(muMiSSEavg)
      plot(muMiSSEavg, x,  lwd = 1,xlim=c(min(rounded_rates), max(rounded_rates)), 
           pch=19, yaxt = "n", xlab="MiSSE mu", ylab="", frame.plot=T, cex=0.75, col=pal[set_rate_three]) #xlim=range(c(min(muTRUE), max(muTRUE)))
      segments(min(muMiSSEavg), 1:length(muMiSSEavg), muMiSSEavg[1:length(muMiSSEavg)], 1:length(muMiSSEavg), col= pal[set_rate_three],lwd = 1)
    }
  }
}
dev.off()
