#tree <- ape::read.tree(file.choose())
#tree_pruned <- tree
# tip.rates_1 <- heights

load(file.choose()) # load results
run_1 <- model.set
model.recons1 <- recon_models
# tip.rates_1 <- hisse::GetModelAveRates(model.recons1, type = c("tips"))

tip.rates_1 <- tip.rates
tree_pruned <- model.recons1[[1]]$phy

plot(x=1:length(run_1), y=unlist(lapply(run_1[1:length(run_1)], "[[", "loglik")), xlab="comb", ylab="loglik")
plot(x=1:length(run_1), y=unlist(lapply(run_1[1:length(run_1)], "[[", "AICc")), xlab="comb", ylab="AICc")

# comb 1 - both turnover and eps vary, 32 models - averaging models
boxplot(tip.rates_1$turnover, tip.rates_1$net.div, 
        tip.rates_1$speciation,tip.rates_1$extinct.frac, 
        tip.rates_1$extinction, main="tip rates run 1", 
        names=c("turnover", "netdiv", "speciation", "eps", "extinction"))

painted_tree1 <- hisse::plot.misse.states(x = model.recons1, 
                                          rate.param = "turnover", type = "phylogram",  
                                          show.tip.label = F, fsize=0.9) 
   

tree_pruned <- ape::drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, tip.rates_1$taxon))

# Small corrections so that tips appear in the right order
is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
ordered_tips <- tree_pruned$edge[is_tip, 2]
right_order <- as.character(tree_pruned$tip.label[ordered_tips])
    


# Organizing so tip rates are in the same order as tips of the tree
cleaned_table <- tip.rates_1[match(as.character(right_order), as.character(tip.rates_1$taxon)),]
# cleaned_table <- tip.rates_1[match(as.character(right_order), as.character(tip.rates_1$species)),]
#turnover.mean <- as.numeric(cleaned_table$max_height_m)

turnover.mean <- as.numeric(cleaned_table$turnover)
# turnover.mean <- 
extinction.frac.mean <- as.numeric(cleaned_table$extinct.frac)
net.div.mean <- as.numeric(cleaned_table$net.div)
    
    # Plotting
    color_breaks = 4 # an interger to be defined as argument in the main function
    layout.matrix <- matrix(c(1,2,3,4), nrow = 1, ncol = 4)
    layout(mat = layout.matrix,widths = c(4,1,1,1))
    par(mar=c(5,0.5,3,0.5))
    
    # Tree
    plot(tree_pruned, show.tip.label=F, edge.width=0.2, adj=1, cex=0.08)
    title(main="Eucalypts")
    ape::axisPhylo()
    
    x <- 1:length(turnover.mean)
    rounded_rates <- round(turnover.mean, color_breaks)
    pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = "Viridis", alpha = 0.75)
    pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
    plot(turnover.mean, x,  lwd = 0.2, xlim=range(c(min(turnover.mean), max(turnover.mean))),
         pch=19, yaxt = "n", xlab="turnover", ylab="", frame.plot=T, cex=0.75, col=pal)
    segments(min(turnover.mean), 1:length(turnover.mean), turnover.mean[1:length(turnover.mean)], 1:length(turnover.mean), col= pal,lwd = 0.2)
    
    x <- 1:length(net.div.mean)
    rounded_rates <- round(net.div.mean, color_breaks)
    pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = "Viridis", alpha = 0.75)
    pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
    plot(net.div.mean, x, lwd = 0.2, xlim=range(c(min(net.div.mean), max(net.div.mean))),
         pch=19, yaxt = "n", xlab="net.div", ylab="", frame.plot=T, cex=0.75, col=pal)
    segments(min(net.div.mean), 1:length(net.div.mean), net.div.mean[1:length(net.div.mean)], 1:length(net.div.mean), col= pal,lwd = 0.2)
    
    x <- 1:length(extinction.frac.mean)
    rounded_rates <- round(extinction.frac.mean, color_breaks)
    pal <- hcl.colors(length(levels(as.factor(rounded_rates))), palette = "Viridis", alpha = 0.75)
    pal <- pal[match(rounded_rates, as.numeric(levels(as.factor(rounded_rates))))] 
    plot(extinction.frac.mean, x, lwd = 0.2, xlim=range(c(min(extinction.frac.mean), max(extinction.frac.mean))),
         pch=19, yaxt = "n", xlab="extinction.frac", ylab="", frame.plot=T, cex=0.75, col=pal)
    segments(min(extinction.frac.mean), 1:length(extinction.frac.mean), extinction.frac.mean[1:length(extinction.frac.mean)], 1:length(extinction.frac.mean), col= pal,lwd = 0.2)

    
#----------------------------
#----------------------------
#----------------------------
    
# Height plots

# Plots
heights <- read.csv("Eucalypt_heights.csv")[,c(2,3)]
heights <- heights[heights$max_height_m!="no_info_yet",]
heights$max_height_m <- as.numeric(heights$max_height_m)
    
hist(heights$max_height_m, breaks = 50)

tree <- ape::read.tree(file.choose())
tree_pruned <- tree
tip.rates_1 <- heights
 
tree_pruned <- ape::drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, tip.rates_1$species))
    
# Small corrections so that tips appear in the right order
is_tip <- tree_pruned$edge[,2] <= length(tree_pruned$tip.label)
ordered_tips <- tree_pruned$edge[is_tip, 2]
right_order <- as.character(tree_pruned$tip.label[ordered_tips])
    
# Organizing so tip rates are in the same order as tips of the tree
cleaned_table <- tip.rates_1[match(as.character(right_order), as.character(tip.rates_1$species)),]

# cleaned_table <- tip.rates_1[match(as.character(right_order), as.character(tip.rates_1$species)),]
height_max <- as.numeric(cleaned_table$max_height_m)
    
# Plotting
color_breaks = 4 # an interger to be defined as argument in the main function
layout.matrix <- matrix(c(1,2), nrow = 1, ncol = 2)
layout(mat = layout.matrix,widths = c(3,1))
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



    
