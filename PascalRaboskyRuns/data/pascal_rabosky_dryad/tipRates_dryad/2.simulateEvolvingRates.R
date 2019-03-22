# Simulation of evolving rates phylogenies, using simulation code from:
# Beaulieu JM, O'Meara BC (2015) Extinction can be estimated from moderately sized molecular phylogenies. Evolution, 69, 1036–1043.
# https://doi.org/10.1111/evo.12614 


# to simulate phylogenies that have an expected diversity of N=100
# pair lambda = 0.078, 0.103, 0.145, 0.249
# with eps = 0, 0.25, 0.5, 0.75
# across sigma (diffusion parameter) = 0.03, 0.06, 0.12
# and bound the number of taxa at 25 - 2500

# lambda and eps are selected so as to achieve an expected diversity of 100
# This is done by optimizing equation 5 of Magallon and Sanderson 2001:
expectedN <- function(a, lambda, eps, t) {
	mu <- eps * lambda
	r <- lambda - mu
	beta_t <- (exp(r * t) - 1) / (exp(r * t) - eps)
	alpha_t <- eps * beta_t
	(a * exp(r * t)) / (1 - alpha_t ^ a)	
}

lambdaVec <- c(0.078, 0.103, 0.145, 0.249)
epsVec <- c(0, 0.25, 0.5, 0.75)
sigmaVec <- c(0.03, 0.06, 0.12)

expectedN(2, lambdaVec[1], epsVec[1], t = 50)
expectedN(2, lambdaVec[2], epsVec[2], t = 50)
expectedN(2, lambdaVec[3], epsVec[3], t = 50)
expectedN(2, lambdaVec[4], epsVec[4], t = 50)

setname <- 'evolvingRates'

require(ape)
require(geiger)
require(diversitree)
require(BAMMtools)

source('sourceScripts/sourceFxns.R')
source('sourceScripts/evolveRateTree_JMB_sqT.R')

# randomly sample from the 4 lambda/eps options in order to get on average 100 taxa
# As the number of exant taxa can be very small, sometimes 1, we will only keep simulations that lead to 50 or more extant taxa.

treeList <- vector('list', 1200)
trueparams <- as.data.frame(matrix(nrow = 1200, ncol = 4))
colnames(trueparams) <- c('treeName', 'lambda', 'eps', 'sigma')

counter <- 1
for (i in 1:length(lambdaVec)) {
	
	for (j in 1:length(sigmaVec)) {
		
		cat('\t', i, '--', j, '\n')
		for (k in 1:100) {
		
			while (1) {
				tr <- evolveRateTree.eps(b=lambdaVec[i], eps = epsVec[i], stdev = sigmaVec[j], mintax = 25, maxtax = 2500)
				if (length(edgesFromTime(tr, time = 0)) >= 25) {
					break
				}
			}
			
			treeList[[counter]] <- tr
			trueparams[counter,] <- c(NA, lambdaVec[i], epsVec[i], sigmaVec[j])
			
			counter <- counter + 1
		}
			
	}
}

trueparams$treeName <- paste0(setname, '_', 1:1200)

# list element 'lambda' in phylo object is a value of lambda for each branch
## keep terminal branch rates and prune tree to extant only

getTerminalLambda <- function(phy) {
	rates <- phy$lambda
	names(rates) <- phy$edge[,2]   # name them after corresponding node
	phy <- drop.extinct(phy) 
	return(rates[as.character(phy$tip.label)])	
}

tiprateList <- vector('list', 1200)
for (i in 1:length(treeList)) {
	
	tiprateList[[i]] <- getTerminalLambda(treeList[[i]])
	
}		

names(tiprateList) <- trueparams$treeName

# drop extinct lineages	
treeList_alltips <- treeList
treeList <- lapply(treeList, drop.extinct)

class(treeList_alltips) <- 'multiPhylo'
class(treeList) <- 'multiPhylo'


