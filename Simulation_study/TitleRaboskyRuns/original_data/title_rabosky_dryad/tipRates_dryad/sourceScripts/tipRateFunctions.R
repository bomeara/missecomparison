require(ape)
require(diversitree)
require(phangorn)


# DR metric / inverse equal splits
DRstat <- function(tree) {
	
	spRate <- function(sp, tree) {
		#get branch lengths from root to tip
		edges <- vector()
		daughterNode <- match(sp, tree$tip.label)
		while (daughterNode != (length(tree$tip.label) + 1)) {
			parentNode <- tree$edge[which(tree$edge[,2] == daughterNode), 1]
			edges <- c(edges, tree$edge.length[which(tree$edge[,1] == parentNode & tree$edge[,2] == daughterNode)])
			daughterNode <- parentNode
		}
		
		res <- sum(sapply(1:length(edges), function(x) edges[x] * (1/(2 ^ (x-1)))))
		res <- res ^ (-1)
		
		return(res)
	}
	
	rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
	names(rates) <- tree$tip.label
	
	return(rates)
}




# ratio of number of speciation events / age of clade
nodeDensity <- function(tree) {
	maxBT <- max(branching.times(tree))
	nodeCounts <- sapply(1:length(tree$tip.label), function(x) {
		n <- 0
		childnode <- x
		parentNode <- tree$edge[which(tree$edge[,2] == childnode), 1]
		while(parentNode != (length(tree$tip.label) + 1)) {
			n <- n + 1
			childnode <- parentNode
			parentNode <- tree$edge[which(tree$edge[,2] == childnode), 1]
		}
		return(n)
	})
	
	# avoid rates of 0 by adding the root node
	nodeCounts <- nodeCounts + 1
	
	return(setNames(nodeCounts / maxBT, tree$tip.label))
}




# inverse of terminal branch lengths
inverseTerminalBranches <- function(tree) {
	tb <- sapply(1:length(tree$tip.label), function(x) {
		tree$edge.length[which(tree$edge[,2] == x)]	
	})
	return(setNames((1 / (2*tb)), tree$tip.label))	
}




# Fit constant-rate birth-death process function
#    note bounds on lambda optimization; may want to increase these or check for boundary fails
fitCRBD <- function(phy, nopt=5, lmin=0.001, lmax=5.0, MAXBAD = 200){
	
	if (length(phy$tip.label) < 3){
		pars <- c(0.0001,0)
		names(pars) <- c("lambda", "mu")
		return(pars)
	}
	
	fx <- make.bd(phy)
	
	for (i in 1:nopt){
	
		lam <- runif(1, 0, 0.5)	
	 	mu <- lam * runif(1, 0, 1)
	 
		badcount <- 0
 
		resx <- try(optim(c(lam, mu) ,fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T)
		while (class(resx) == 'try-error'){

			lam <- runif(1, 0, 0.5)	
	 		mu <- lam * runif(1, 0, 1)
			
			resx <- try(optim(c(lam, mu) , fx, method='L-BFGS-B', control=list(maxit=1000, fnscale=-1), lower=lmin, upper=lmax), silent=T);
			
			badcount <- badcount + 1;
			if (badcount > MAXBAD){
				stop("Too many fails in fitDiversitree\n")
			}
		}
		
		if (i == 1){
			best <- resx
		}else{
			if (best$value < resx$value){
				best <- resx
			}
		}
		
	}
	
	fres <- list(pars=best$par, loglik=best$value)
	fres$AIC <- -2*fres$loglik + 2*length(argnames(fx))
	fres$counts <- best$counts
	#fres$like_function <- fx
	fres$convergence <- best$convergence
	fres$message <- best$message
	
	pars <- fres$pars
	names(pars) <- c("lambda", "mu")
	
	return(pars)
}



