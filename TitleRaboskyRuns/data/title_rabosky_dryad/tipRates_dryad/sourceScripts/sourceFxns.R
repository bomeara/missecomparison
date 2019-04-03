require(ape)
require(geiger)
require(diversitree)
require(coda)

# externalize BAMMtools NU.branching.times() function for non-ultrametric trees
NU.branching.times <- BAMMtools:::NU.branching.times

# function to return the branch indices that exist at a certain time
# time = 0 at present
edgesFromTime <- function(tree, time, tol = NULL) {
	
	btimes <- NU.branching.times(tree)
	
	if (is.null(tol)) {
		tol <- min(tree$edge.length)/100
	}
	
	branchTimes <- t(apply(tree$edge, 1, function(x) btimes[as.character(x)]))
	
	# NA values are terminal branches
	## if tree were extant only, ultrametric, we would set all NA values to 0, as terminal branches terminate at time 0. 
	## However, we will make this functional with non-ultrametric trees as well.
	## So for terminal branches, we will pull the branch length and subtract from the age of the parent node.
	
	terminalEdge <- which(is.na(branchTimes[,2]))
	branchTimes[terminalEdge, 2] <- branchTimes[terminalEdge, 1] - tree$edge.length[terminalEdge]
	
	# some terminal branches should reach to 0, but don't quite make it
	extant <- which((branchTimes[,2] - 0) < tol)
	if (length(extant) > 0) {
		branchTimes[extant, 2] <- 0
	}

	if (anyNA(branchTimes)) {stop('NA in branchTimes.')}
	
	branchInd <- which(apply(branchTimes, 1, function(x) time <= x[1] & time >= x[2]) == TRUE)
	
	return(branchInd)
}


# function to take simulation parameters and create a BAMMdata object.
modelDat2edata <- function(modelDat, tree) {
	
	res <- as.data.frame(matrix(nrow = nrow(modelDat), ncol = 8)) 
	colnames(res) <- c('generation', 'leftchild', 'rightchild', 'abstime', 'lambdainit', 'lambdashift', 'muinit', 'mushift')
	
	res[, 'generation'] <- 1
	res[, 'abstime'] <- modelDat$regimeTime
	res[, 'lambdainit'] <- modelDat$lambda0
	res[, 'lambdashift'] <- 0
	res[, 'muinit'] <- modelDat$mu0
	res[, 'mushift'] <- 0
	
	for (i in 1:nrow(modelDat)) {
	
		res[i, c('leftchild', 'rightchild')] <- BAMMtools:::getSpanningTaxonPair(tree, geiger::tips(tree, modelDat[i, 'regimeMRCA']))
 	
	}
	
	return(getEventData(tree, res))
}



# Computes the Colless imbalance statistic 
#	across an entire tree.
# from Rabosky 2016, Evolution
colless <- function(phy){
	
	bb <- balance(phy);
	ss <- sum(abs(bb[,1] - bb[,2]));
	n <- length(phy$tip.label);
	return((2 / ((n-1)*(n-2))) * ss);	
}

# calculate tip rates from a diversity-dependent model
tipLambdaDD <- function(lambda0, nTips, K) {
	res <- lambda0 * (1 - nTips / K)
	ifelse(res < 0, 0, res)
}

# function to return the edge indices, in order, from root to tip
root2tipEdges <- function(phy, sp) {

	if (!sp %in% phy$tip.label) {
		stop('sp not a tip label.')
	}

	res <- c()
	childnode <- which(phy$tip.label == sp)
	branch <- which(phy$edge[,2] == childnode)
	parentnode <- phy$edge[branch, 1]
	res <- c(res, branch)
	childnode <- parentnode
	
	while (childnode != length(phy$tip.label) + 1) {
		branch <- which(phy$edge[,2] == childnode)
		parentnode <- phy$edge[branch, 1]
		res <- c(res, branch)
		childnode <- parentnode
	}
	return(rev(res))
}


# intended for evolving rates model
# extract rates for terminal branches
getTerminalLambda <- function(phy) {
	rates <- phy$lambda
	names(rates) <- phy$edge[,2]   # name them after corresponding node
	phy <- drop.extinct(phy) 
	return(rates[as.character(phy$tip.label)])	
}

getTerminalMu <- function(phy) {
	rates <- phy$mu
	names(rates) <- phy$edge[,2]   # name them after corresponding node
	phy <- drop.extinct(phy) 
	return(rates[as.character(phy$tip.label)])	
}

getTerminalEdgeIndex <- function(tree) {
	
	setNames(sapply(1:length(tree$tip.label), function(x) which(tree$edge[,2] == x)), tree$tip.label)
}

# return convergence statistics from BAMM mcmc file
# returns effective size and Geweke z-score for the log likelihood and the number of shifts
# effective size > 200 is generally adequate
# Geweke z-score between -2 and 2 is good.
checkBAMMconvergence <- function(mcmcFile = NULL, burnin = 0.5) {
	if (is.null(mcmcFile)) {
		files <- list.files(pattern='mcmc_out')
		if (length(files) == 1) {
			mcmcFile <- files
		} else {
			stop('Please specify an mcmc_out file. Several were detected.')
		}
	}
	
	if (class(mcmcFile) == 'character') {
		# mcmc <- read.csv(mcmcFile)
		mcmc <- data.table::fread(mcmcFile, data.table = FALSE)
	}
	mcmc2 <- mcmc[floor(burnin * nrow(mcmc)):nrow(mcmc), ]
	
	res <- numeric(4)
	names(res) <- c('loglik_effSize', 'nShifts_effSize', 'loglik_geweke', 'nShift_geweke')

	res[1] <- coda::effectiveSize(mcmc2[,4]) #Nshifts
	res[2] <- coda::effectiveSize(mcmc2[,2]) #loglik
	res[3] <- coda::geweke.diag(coda::as.mcmc(mcmc2[,4]))$z
	res[4] <- coda::geweke.diag(coda::as.mcmc(mcmc2[,2]))$z
	
	return(res)
}



