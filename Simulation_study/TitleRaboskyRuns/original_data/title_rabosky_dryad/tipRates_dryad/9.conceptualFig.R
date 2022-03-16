# conceptual figure showing true vs BAMM vs DR tip rates, on a simple rate shift tree, and on an evolving rates tree


require(BAMMtools)
require(phytools)

source('sourceScripts/tipRateFunctions.R')
source('sourceScripts/sourceFxns.R')



# ------------------------------------

setwd('conceptualFigTrees/toytree')

toytree <- read.tree('tree.tre')
toytree <- ladderize(toytree)

# create bamm object based on true parameters
true <- as.data.frame(matrix(nrow = 2, ncol = 8)) 
colnames(true) <- c('generation', 'leftchild', 'rightchild', 'abstime', 'lambdainit', 'lambdashift', 'muinit', 'mushift')
true[, 'generation'] <- 1
true[, 'abstime'] <- c(0, 75)
true[, 'lambdainit'] <- c(0.02995732, 0.1360479)
true[, 'lambdashift'] <- 0
true[, 'muinit'] <- 0
true[, 'mushift'] <- 0
true[, 'leftchild'] <- c('t38', 's1')
true[, 'rightchild'] <- c('t35', 's30')

trueEd <- getEventData(toytree, true)
trueTips <- getTipRates(trueEd)$lambda.avg

checkBAMMconvergence(burnin=0.25)

ed <- getEventData(toytree, 'event_data.txt', nsamples=1000, burnin=0.25)

bamm <- getTipRates(ed)$lambda.avg
dr <- DRstat(toytree)

tipOrder <- sapply(toytree$tip.label, function(x) which(toytree$edge[,2] == which(toytree$tip.label == x)))
tipOrder <- names(sort(tipOrder))

trueTips <- trueTips[tipOrder]
bamm <- bamm[tipOrder]
dr <- dr[tipOrder]

# shift-specific means
tipRegimes <- rep(0, length(toytree$tip.label))
names(tipRegimes) <- tipOrder
tipRegimes[names(tipRegimes) %in% tips(toytree, getMRCA(toytree, c('s1','s30')))] <- 1

regime0MeanDR <- mean(dr[names(tipRegimes)[which(tipRegimes == 0)]])
regime1MeanDR <- mean(dr[names(tipRegimes)[which(tipRegimes == 1)]])
regime0MeanBAMM <- mean(bamm[names(tipRegimes)[which(tipRegimes == 0)]])
regime1MeanBAMM <- mean(bamm[names(tipRegimes)[which(tipRegimes == 1)]])



pdf('fig6.pdf', width=6, height=7)
# other option involving error in rates
layout(matrix(1:3, nrow=1, ncol=3), widths=c(0.3, 0.35, 0.35))
par(mar = c(5,4,4,0))
plot.phylo(toytree, show.tip.label = FALSE, edge.width=1.5)
addBAMMshifts(trueEd, par.reset=F, bg='orange', cex=2.5)
par(mar = c(5,0,4,2))
plot.new()
plot.window(xlim = c(0, 0.3), ylim = c(0, 69))
axis(1)
points(trueTips, 1:length(toytree$tip.label), pch=20)
lines(trueTips, 1:length(toytree$tip.label), lwd=1)
points(bamm, 1:length(toytree$tip.label), col='blue', pch=20)
lines(bamm, 1:length(toytree$tip.label), lwd=0.5, col='blue')
points(dr, 1:length(toytree$tip.label), col='red', pch=20)
lines(dr, 1:length(toytree$tip.label), lwd=0.5, col='red')
mtext('speciation rate', side = 1, line = 2.5)

legend(0.15, 69, legend = c(expression(paste(lambda['TRUE'])), expression(paste(lambda['BAMM'])), expression(paste(lambda['DR']))), fill=c('black', 'blue','red'), bty='n', cex=1.5)

plot.new()
plot.window(xlim = c(0, 0.2), ylim = c(0, 69))
axis(1)
segments(abs(bamm - trueTips), 1:length(toytree$tip.label), abs(dr - trueTips), 1:length(toytree$tip.label), lwd=0.5)
points(abs(bamm - trueTips), 1:length(toytree$tip.label), pch=20, col='blue')
points(abs(dr - trueTips), 1:length(toytree$tip.label), pch=20, col='red')
abline(v=0, lty=3)

# add means
points(mean(abs(bamm - trueTips)), -1, col='blue', pch=8)
points(mean(abs(dr - trueTips)), -1, col='red', pch=8)

mtext('absolute error', side = 1, line = 2.5)
dev.off()

mean(abs(bamm - trueTips))
mean(abs(dr - trueTips))

var(abs(bamm - trueTips))
var(abs(dr - trueTips))
# ----------------------------------------
# do the same with a evolvingRates tree
setwd('../../conceptualFigTrees/evolvingRates_250')
# sigma = 0.12

alltipRates <- readRDS('../../dataFiles/allRatesList.rds')
trueTips <- alltipRates[['evolvingRates_250']]$tipLambda
names(trueTips) <- alltipRates[['evolvingRates_250']]$tipName

tr <- read.tree(grep('\\.tre$', list.files(), value=TRUE))
tr <- ladderize(tr)

ed <- getEventData(tr, grep('event_data', list.files(), value=TRUE), nsamples=1000, burnin=0.5)
bamm <- getTipRates(ed)$lambda.avg
dr <- DRstat(tr)

tipOrder <- sapply(tr$tip.label, function(x) which(tr$edge[,2] == which(tr$tip.label == x)))
tipOrder <- names(sort(tipOrder))

trueTips <- trueTips[tipOrder]
bamm <- bamm[tipOrder]
dr <- dr[tipOrder]

ptsize <- 0.5

pdf('fig7.pdf', width=6, height=7)
# other option involving error in rates
layout(matrix(1:3, nrow=1, ncol=3), widths=c(0.3, 0.35, 0.35))
par(mar = c(5,4,4,0))
plot.phylo(tr, show.tip.label = FALSE, edge.width=1.5)
par(mar = c(5,0,4,2))
plot.new()
plot.window(xlim = c(0, max(c(max(trueTips), max(bamm), max(dr)))), ylim = c(0, length(tr$tip.label)))
axis(1)
points(trueTips, 1:length(tr$tip.label), pch=20, cex= ptsize)
lines(trueTips, 1:length(tr$tip.label), lwd=1)
points(bamm, 1:length(tr$tip.label), col='blue', pch=20, cex= ptsize)
lines(bamm, 1:length(tr$tip.label), lwd=0.5, col='blue')
points(dr, 1:length(tr$tip.label), col='red', pch=20, cex= ptsize)
lines(dr, 1:length(tr$tip.label), lwd=0.5, col='red')
mtext('speciation rate', side = 1, line = 2.5)

legend(0.4, 207, legend = c(expression(paste(lambda['TRUE'])), expression(paste(lambda['BAMM'])), expression(paste(lambda['DR']))), fill=c('black', 'blue','red'), bty='n', cex=1.5)

plot.new()
plot.window(xlim = c(0, 0.7), ylim = c(0, length(tr$tip.label)))
axis(1)
segments(abs(bamm - trueTips), 1:length(tr$tip.label), abs(dr - trueTips), 1:length(tr$tip.label), lwd=0.5)
points(abs(bamm - trueTips), 1:length(tr$tip.label), pch=20, col='blue', cex= ptsize)
points(abs(dr - trueTips), 1:length(tr$tip.label), pch=20, col='red', cex= ptsize)
abline(v=0, lty=3)

# add means
points(mean(abs(bamm - trueTips)), -5, col='blue', pch=8)
points(mean(abs(dr - trueTips)), -5, col='red', pch=8)
mtext('absolute error', side = 1, line = 2.5)

dev.off()

mean(abs(bamm - trueTips))
mean(abs(dr - trueTips))

var(abs(bamm - trueTips))
var(abs(dr - trueTips))

