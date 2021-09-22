################################################
# Post-run 
################################################
# rm(list=ls())
library(ape)
library(hisse)
library(BAMMtools)
library(gridExtra)
library(viridis)
library(ggplot2)

################################################
# Load MiSSE results
################################################
load("MiSSE_files/Eucalypts_example_ml1_tree.Rsave") # load results
tree_misse <- model.recons_ml1[[1]]$phy

################################################
# Load BAMM results
################################################
tree_bamm <- read.tree("BAMM_files/ML1_modified-Thornhill_et_al_2019.tre")
event_data <- getEventData(tree_bamm, "BAMM_files/event_data.txt", nsamples = 45000) 
bamm_tiprates <- getTipRates(event_data)

# Getting BAMM net-turnover rates
lambda <- bamm_tiprates$lambda.avg
mu <- bamm_tiprates$mu.avg
bamm_tiprates_final <- data.frame(species=names(lambda), lambda=unname(lambda), mu=unname(mu))
colnames(bamm_tiprates_final) <- c("species","speciation.bamm","extinction.bamm")
bamm_tiprates_final$turnover.bamm <- bamm_tiprates_final$speciation.bamm + bamm_tiprates_final$extinction.bamm
#write.csv(bamm_tiprates_final, file="bamm_tiprates_final.csv")

###############################################
# Tip-correlation with plant height
###############################################
heights <- read.csv("Eucalypt_height-EUCLID.csv")[,c(2,3)] # Loading data from EUCLID
heights <- heights[heights$max_height_m!="no_info_yet",]

# Combining heights and rates
full_table_ml1 <- merge(tip.rates_ml1, heights, by.x="taxon",by.y="species")
full_table_ml1 <- merge(full_table_ml1, bamm_tiprates_final, by.x="taxon",by.y="species")
full_table_ml1$max_height_m <- as.numeric(full_table_ml1$max_height_m)
rownames(full_table_ml1) <- full_table_ml1$taxon

###############################################
# PICs
###############################################
tree_pruned <- ape::keep.tip(tree_misse, full_table_ml1$taxon) 

# Small corrections so that tips appear in the right order in the plot
is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
ordered_tips <- tree_pruned$edge[is_tip, 2]
right_order <- as.character(tree_pruned$tip.label[ordered_tips])

# Organizing so tip rates are in the same order as tips of the tree
cleaned_table <- full_table_ml1[match(as.character(right_order), as.character(full_table_ml1$taxon)),]
height_max <- as.numeric(cleaned_table$max_height_m)

# PICs
t0 <- log(cleaned_table$max_height_m)
r0 <- log(cleaned_table$turnover)

t0Pic <- pic(t0, tree_pruned, scaled =T)
r0Pic <- pic(r0, tree_pruned, scaled =T)

names(t0Pic) <- (ape::Ntip(tree_pruned) + 1) : (ape::Ntip(tree_pruned) + (ape::Ntip(tree_pruned) - 1))
names(r0Pic) <- (ape::Ntip(tree_pruned) + 1) : (ape::Ntip(tree_pruned) + (ape::Ntip(tree_pruned) - 1))

# Identifying cherries
internals <- ape::Ntip(tree_pruned) + sequence(ape::Nnode(tree_pruned))
ndescendants <- rep(NA, length(internals))
for (i in seq_along(internals)) {
        ndescendants[i] <- length(phytools::getDescendants(tree_pruned, node=internals[i]))
}
cherry_internals <- internals[ndescendants==2]
t0Pic_nocherry <- t0Pic[!(as.numeric(names(t0Pic)) %in% cherry_internals)]
r0Pic_nocherry <- r0Pic[!(as.numeric(names(r0Pic)) %in% cherry_internals)]

# Positivizing PICs and Regression through the origin
r0Pic[which(t0Pic < 0)] <- -1 * r0Pic[which(t0Pic < 0)]
t0Pic <- abs(t0Pic)
picModel <- lm(r0Pic ~ 0 + t0Pic)
regression.origin_w_cherry <- lmorigin(picModel, picModel$model, origin=TRUE, nperm=999)

r0Pic_nocherry[which(t0Pic_nocherry < 0)] <- -1 * r0Pic[which(t0Pic_nocherry < 0)]
t0Pic_nocherry <- abs(t0Pic_nocherry)
picModel <- lm(r0Pic_nocherry ~ 0 + t0Pic_nocherry)
regression.origin_no_cherry <- lmorigin(picModel, picModel$model, origin=TRUE, nperm=999)

sink("PIC_results.txt")
cat("Results including cherries")
cat("------------------------------")
print(regression.origin_w_cherry)
cat("Results removing cherries")
cat("------------------------------")
print(regression.origin_no_cherry)
sink()

###############################################
###############################################
# Plots for Figure 5
###############################################
# Insert (d) (PICs)
all_pic_table_no_cherries <- as.data.frame(cbind(t0Pic_nocherry, r0Pic_nocherry))
all_pic_table_w_cherries <- as.data.frame(cbind(t0Pic, r0Pic))

pdf("PICs.pdf", width=9, height=4)
no_cherries <- ggplot(all_pic_table_no_cherries, aes(t0Pic_nocherry, r0Pic_nocherry, color = t0Pic_nocherry)) +
        geom_point(shape = 19, size = 3, show.legend = FALSE, alpha=0.75) +
        theme_minimal() +
        scale_color_viridis(option = "B") +
        geom_abline(intercept=0, slope=regression.origin_no_cherry$reg$coefficients, color="red",  linetype="dashed")  +
        geom_smooth(se=FALSE, formula=y~x-1, col="red", lwd=0.5) +
        ylim(-2, 2) +
        xlim(0,8.5) +
        xlab("PIC log(plant height) (m)") 

w_cherries <- ggplot(all_pic_table_w_cherries, aes(t0Pic, r0Pic, color = t0Pic)) +
        geom_point(shape = 19, size = 3, show.legend = FALSE, alpha=0.75) +
        theme_minimal() +
        scale_color_viridis(option = "B") +
        geom_abline(intercept=0, slope=regression.origin_w_cherry$reg$coefficients, color="red",  linetype="dashed")  +
        geom_smooth(se=FALSE, formula=y~x-1, col="red", lwd=0.5) +
        ylim(-2, 2) +
        xlim(0,8.5) +
        xlab("PIC log(plant height) (m)") 
grid.arrange(w_cherries, no_cherries, nrow=1, ncol=2)

dev.off()

        
###############################################
# Main figure with Eucalypts phylogeny, tree height and rates
###############################################
# Small corrections so that tips appear in the right order in the plot
is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
ordered_tips <- tree_pruned$edge[is_tip, 2]
right_order <- as.character(tree_pruned$tip.label[ordered_tips])

# Organizing so tip rates are in the same order as tips of the tree
cleaned_table <- full_table_ml1[match(as.character(right_order), as.character(full_table_ml1$taxon)),]
height_max <- as.numeric(cleaned_table$max_height_m)

# Plotting
{
pdf("Eucalypts_comparison.pdf", width=3.5, height=8)
color_breaks = 4 
layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
layout(mat = layout.matrix,widths = c(2,1,1,1))
par(mar=c(5,0.5,3,0.5))
# Tree
plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
title(main="Eucalypts")
ape::axisPhylo()

palette_name = "Inferno"

# Plant height
rounded_rates <- round(height_max, color_breaks)
pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name, alpha = 0.75)
pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
x <- 1:length(height_max)
plot(height_max, x,  lwd = 0.2, xlim=range(c(min(height_max), max(height_max))),
     pch=19, yaxt = "n", xlab="max_height_m", ylab="", frame.plot=T, cex=0.75, col=pal)
segments(min(height_max), 1:length(height_max), height_max[1:length(height_max)], 1:length(height_max), col= pal,lwd = 0.2)

# Turnover MiSSE
turnover.mean <- as.numeric(cleaned_table$turnover)
rounded_rates <- round(turnover.mean, color_breaks)
pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name, alpha = 0.75)
pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
x <- 1:length(turnover.mean)
plot(turnover.mean, x,  lwd = 0.2, xlim=c(0,6), 
     pch=19, yaxt = "n", xlab="turnover.misse", ylab="", frame.plot=T, cex=0.75, col=pal) #xlim=range(c(min(turnover.mean), max(turnover.mean)))
segments(min(turnover.mean), 1:length(turnover.mean), turnover.mean[1:length(turnover.mean)], 1:length(turnover.mean), col= pal,lwd = 0.2)

# Turnover BAMM
turnover.bamm <- as.numeric(cleaned_table$turnover.bamm)
rounded_rates <- round(turnover.bamm, color_breaks)
pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = palette_name, alpha = 0.75)
pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
x <- 1:length(turnover.bamm)
plot(turnover.bamm, x,  lwd = 0.2, xlim=c(0,6),
     pch=19, yaxt = "n", xlab="turnover.bamm", ylab="", frame.plot=T, cex=0.75, col=pal)
segments(min(turnover.bamm), 1:length(turnover.bamm), turnover.bamm[1:length(turnover.bamm)], 1:length(turnover.bamm), col= pal,lwd = 0.2)

dev.off()
}

