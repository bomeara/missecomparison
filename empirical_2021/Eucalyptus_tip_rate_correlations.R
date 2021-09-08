#----------------------------
# Post-run 
# setwd("~/Desktop/misse_mme_paper/missecomparison/empirical_2021")
# rm(list=ls())
library(ape)
library(hisse)

################################################
# BAMM
################################################
library(BAMMtools)
library(coda)
eucalypts <- read.tree("attachments/ML1_modified.tre")

# checking mcmc convergence
mcmcout <- read.csv("attachments/mcmc_out.txt", header=T)
#plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ] # 10% burnin
effectiveSize(postburn$N_shifts) # over 200
effectiveSize(postburn$logLik) # over 200

table(postburn$N_shifts ) / nrow(postburn)
# plotPrior(postburn)

ed <- getEventData(eucalypts, "attachments/event_data.txt", nsamples = 45000)

summary(ed)

pal.name = "Viridis" 
pal <- hcl.colors(30, palette = pal.name, alpha = 1)

plot <- plot.bammdata(ed, pal=pal)
addBAMMlegend(plot)
addBAMMshifts(plot)



#z <- plot.bammdata(ed, tau = 0.002, lwd=2)
#addBAMMlegend(z)

bamm_tiprates <- getTipRates(ed)

lambda <- bamm_tiprates$lambda.avg
mu <- bamm_tiprates$mu.avg
bamm_tiprates_final <- data.frame(species=names(lambda), lambda=unname(lambda), mu=unname(mu))
write.csv(bamm_tiprates_final, file="bamm_tiprates_final.csv")

################################################
# MiSSE 
################################################
load("Eucalypts_example_ml1_tree.Rsave") # load results


MiSSENet(model.set_pruned_ml1)

PlotMisseSpace(model.set_pruned_ml1, size=10, size2=5, width=2, arrow.size=0.2, arrow.width=0.2)



?plot.igraph

#
load("Eucalypts_example_bayes_tree.Rsave") # load results
model.recons_ml1 <- model.recons_bayes
tip.rates_ml1 <- tip.rates_bayes


tree <- model.recons_ml1[[1]]$phy

#View(possible.combos_ml1)
     #tip.rates_best <- hisse::GetModelAveRates(model.recons_ml1[[which.min(unlist(lapply(model.recons_ml1, "[[","AIC")))]], type = "tips")
#tip.rates_avg <- hisse::GetModelAveRates(model.recons_ml1, type = "tips")
#tip.rates_ml1


###############################################
# Some exploratory things
###############################################

# Getting tip age for plots:
tip_age <- c()
for(i in 1:length(tree$tip.label)) {
  tip_age[i] <- tree$edge.length[which(tree$edge[,2] == i)]
}
names(tip_age) <- tree$tip.label

#tip.rates_avg$taxon[which.max(tip.rates_avg$turnover)]
#tip_age[names(tip_age)==tip.rates_avg$taxon[which.max(tip.rates_avg$turnover)]]

painted.tree <- hisse::plot.misse.states(x = model.recons_ml1, 
                                         rate.param = "speciation", type = "phylo", show.tip.label = F) 

###############################################
# Tip-correlation with plant height
###############################################


heights <- read.csv("Eucalypt_heights.csv")[,c(2,3)]
heights <- heights[heights$max_height_m!="no_info_yet",]

#plot(tree, show.tip.label = T, cex=0.2) 


# Combining heights and rates

full_table_ml1 <- merge(tip.rates_ml1, heights, by.x="taxon",by.y="species")
rownames(full_table_ml1) <- full_table_ml1$taxon

#full_table_best <- merge(tip.rates_best, heights, by.x="taxon",by.y="species")
#rownames(full_table_best) <- full_table_best$taxon

#full_table_avg <- merge(tip.rates_avg, heights, by.x="taxon",by.y="species")
#rownames(full_table_avg) <- full_table_avg$taxon


# plot(full_table$turnover, full_table$max_height_m)

#----------------------------
# Some plots
#----------------------------
# Eucalypts Height 
full_table_ml1$max_height_m <- as.numeric(full_table_ml1$max_height_m)
# removing outliers?
full_table_ml1 <- subset(full_table_ml1, full_table_ml1$turnover < 100)


#hist(exp(full_table_ml1$max_height_m), breaks = 50)

tree_pruned <- ape::drop.tip(tree, setdiff(tree$tip.label, full_table_ml1$taxon))

# Small corrections so that tips appear in the right order
is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
ordered_tips <- tree_pruned$edge[is_tip, 2]
right_order <- as.character(tree_pruned$tip.label[ordered_tips])

# Organizing so tip rates are in the same order as tips of the tree
cleaned_table <- full_table_ml1[match(as.character(right_order), as.character(full_table_ml1$taxon)),]



#cleaned_table <- full_table_avg[match(as.character(right_order), as.character(full_table_avg$taxon)),]

# cleaned_table_heights <- tip.rates_1[match(as.character(right_order), as.character(tip.rates_1$species)),]
height_max <- as.numeric(cleaned_table$max_height_m)

# Plotting
color_breaks = 4 # an interger to be defined as argument in the main function
layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,widths = c(3,1,1,1))
par(mar=c(5,0.5,3,0.5))

# Tree
plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
title(main="Eucalypts")
ape::axisPhylo()

x <- 1:length(height_max)
rounded_rates <- round(height_max, color_breaks)
pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = "Viridis", alpha = 0.75)
pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
plot(height_max, x,  lwd = 0.2, xlim=range(c(min(height_max), max(height_max))),
     pch=19, yaxt = "n", xlab="max_height_m", ylab="", frame.plot=T, cex=0.75, col=pal)
segments(min(height_max), 1:length(height_max), height_max[1:length(height_max)], 1:length(height_max), col= pal,lwd = 0.2)


# Some plots:

#plot(x=1:length(run_1), y=unlist(lapply(run_1[1:length(run_1)], "[[", "loglik")), xlab="comb", ylab="loglik")
#plot(x=1:length(run_1), y=unlist(lapply(run_1[1:length(run_1)], "[[", "AICc")), xlab="comb", ylab="AICc")

#boxplot(tip.rates_1$turnover, tip.rates_1$net.div, 
#        tip.rates_1$speciation,tip.rates_1$extinct.frac, 
#        tip.rates_1$extinction, main="tip rates run 1", 
#        names=c("turnover", "netdiv", "speciation", "eps", "extinction"))

#painted_tree1 <- hisse::plot.misse.states(x = model.recons1, 
#                                          rate.param = "turnover", type = "phylogram",  
#                                          show.tip.label = F, fsize=0.9, ) 


turnover.mean <- as.numeric(cleaned_table$turnover)

#net.div.mean <- as.numeric(cleaned_table$net.div)


x <- 1:length(turnover.mean)
rounded_rates <- round(turnover.mean, color_breaks)
pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = "Viridis", alpha = 0.75)
pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
plot(turnover.mean, x,  lwd = 0.2, xlim=range(c(min(turnover.mean), max(turnover.mean))),
     pch=19, yaxt = "n", xlab="turnover", ylab="", frame.plot=T, cex=0.75, col=pal)
segments(min(turnover.mean), 1:length(turnover.mean), turnover.mean[1:length(turnover.mean)], 1:length(turnover.mean), col= pal,lwd = 0.2)

#x <- 1:length(net.div.mean)
#rounded_rates <- round(net.div.mean, color_breaks)
#pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = "Viridis", alpha = 0.75)
#pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
#plot(net.div.mean, x, lwd = 0.2, xlim=range(c(min(net.div.mean), max(net.div.mean))),
#     pch=19, yaxt = "n", xlab="net_div", ylab="", frame.plot=T, cex=0.75, col=pal)
#segments(min(net.div.mean), 1:length(net.div.mean), net.div.mean[1:length(net.div.mean)], 1:length(net.div.mean), col= pal,lwd = 0.2)


###
dev.off()

library(phytools)

matrix_netdiv <- cleaned_table[,c("net.div", "max_height_m")]
phylomorphospace(tree_pruned, matrix_netdiv, label=c("off"), col.edge=pal)

matrix_turnover <- cleaned_table[,c("turnover", "max_height_m")]
phylomorphospace(tree_pruned, matrix_turnover, label=c("off"), col.edge=pal)


#----- PGLS

full_table_avg$max_height_m <- log(as.numeric(full_table_avg$max_height_m))
full_table_avg$turnover <- log(as.numeric(full_table_avg$turnover))
full_table_avg$net.div <- log(as.numeric(full_table_avg$net.div))

tree_pruned <- ape::keep.tip(tree, full_table_avg$taxon) 

#result_best <- lm(full_table_avg$net.div~full_table_avg$max_height_m)
#result_best <- lm(full_table_avg$turnover~full_table_avg$max_height_m)

  result1 <- phylolm::phylolm(turnover~max_height_m, data=full_table_avg, phy=tree_pruned, model="BM")
summary(result1)


summary(result_best)

result_best <- phylolm::phylolm(turnover~max_height_m, data=full_table_avg, phy=tree_pruned, model="BM")


# Combining heights and rates for PGLS
full_table_avg <- merge(tip.rates_avg, heights, by.x="taxon",by.y="species")
rownames(full_table_avg) <- full_table_avg$taxon
tree_pruned <- ape::keep.tip(tree, full_table_avg$taxon) 
result_avg <- phylolm::phylolm(turnover~max_height_m, data=full_table_avg, phy=tree_pruned, model="BM")
summary(result_avg)

full_table_avg <- subset(full_table_avg, full_table_avg$turnover < 100)

plot(full_table_avg$turnover, full_table_avg$max_height_m)

  



pgls_table <- full_table_ml1


pgls_table$max_height_m <- log(as.numeric(pgls_table$max_height_m))
pgls_table$turnover <- log(as.numeric(pgls_table$turnover))
tree_pruned <- ape::keep.tip(tree, pgls_table$taxon) 


result_ml1 <- phylolm::phylolm(turnover~max_height_m, data=pgls_table, phy=tree_pruned, model="BM")
summary(result_ml1)

plot(pgls_table$turnover, pgls_table$max_height_m)

###############################################
# Some plots
###############################################

