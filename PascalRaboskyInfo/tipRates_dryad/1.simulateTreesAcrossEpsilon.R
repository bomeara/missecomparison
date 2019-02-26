# Tree simulation for testing whether DR is tracking speciation rate or net div rate. 

# issue: If we simulate trees by sampling lambda and epsilon from uniform distributions, then net diversification will tend to decrease with epsilon, and skew the relationship between net div and epsilon. 
# So we will simulate a set of trees such that the variance in lambda is constant across values of epsilon (set 1), and another set of trees such that variance in net div is constant across epsilon (set 2).

# simulate two sets of trees:
# (1) 1000 100-tip trees where lambda and eps are drawn from independent uniform distributions
# (2) 1000 100-tip trees where net div and eps are drawn from independent uniform distributions

setname1 <- 'lambdaConstantVariance'
setname2 <- 'netDivConstantVariance'


require(TreeSim)
require(BAMMtools)

set1Trees <- vector('list', 1000)
set1params <- matrix(nrow=1000, ncol=3)
colnames(set1params) <- c('lambda','mu','eps')
for (i in 1:1000) {
	
	cat(i, '\n')
	Lambda <- runif(1, 0.05, 0.3)
	Eps <- runif(1, 0, 1)
	Mu <- Eps * Lambda
	set1params[i, ] <- c(Lambda, Mu, Eps)
	set1Trees[[i]] <- sim.bd.taxa(100, 1, lambda=Lambda, mu=Mu, complete=FALSE)[[1]]
}


set2Trees <- vector('list', 1000)
set2params <- matrix(nrow=1000, ncol=3)
colnames(set2params) <- c('lambda','mu','eps')
for (i in 1:1000) {	
	
	cat(i, '\n')
	NetDiv <- runif(1, 0.05, 0.3)
	Eps <- runif(1, 0, 1)
	Lambda <- NetDiv / (1 - Eps)
	Mu <- Eps * Lambda
	set2params[i, ] <- c(Lambda, Mu, Eps)
	set2Trees[[i]] <- sim.bd.taxa(100, 1, lambda=Lambda, mu=Mu, complete=FALSE)[[1]]
}


set1params <- as.data.frame(set1params)
set2params <- as.data.frame(set2params)
set1params$netdiv <- set1params$lambda - set1params$mu
set2params$netdiv <- set2params$lambda - set2params$mu


class(set1Trees) <- 'multiPhylo'
class(set2Trees) <- 'multiPhylo'


# set 1 should have lambda variance constant across epsilon, set 2 should have net div variance constant across epsilon




pdf('figS1.pdf', width=7, height=7)


m <- matrix(c(1,1,2,3,4,4,5,6), nrow=4, ncol=2, byrow=TRUE)
layout(m, heights=c(0.1, 0.4, 0.1, 0.4))

par(mar=c(4,4,0,1))

plot.new()
text(0.5, 0.2, 'simulations for evaluating speciation rate', pos=1, cex=1.5, font=2, xpd=NA)
plot.new()
plot.window(xlim = c(0,1), ylim = c(0, 0.3))
axis(1, at = c(-0.5, axTicks(1)))
axis(2, at = c(-0.5, axTicks(2)))
points(set1params$eps, set1params$lambda, cex=0.5)
mtext(expression(mu~'/'~lambda), side=1, line=2.5)
mtext(expression(lambda), side=2, line=2.5)

plot.new()
plot.window(xlim = c(0,1), ylim = c(0, 0.3))
axis(1, at = c(-0.5, axTicks(1)))
axis(2, at = c(-0.5, axTicks(2)))
points(set1params$eps, set1params$netdiv, cex=0.5)
mtext(expression(mu~'/'~lambda), side=1, line=2.5)
mtext(expression(lambda~'-'~mu), side=2, line=2.5)

plot.new()
text(0.5, 0.2, 'simulations for evaluating net div. rate', pos=1, cex=1.5, font=2, xpd=NA)
plot.new()
plot.window(xlim = c(0,1), ylim = c(0, 10))
axis(1, at = c(-0.5, axTicks(1)))
axis(2, at = c(-0.5, axTicks(2)))
points(set2params$eps, set2params$lambda, cex=0.5)
mtext(expression(mu~'/'~lambda), side=1, line=2.5)
mtext(expression(lambda), side=2, line=2.5)

plot.new()
plot.window(xlim = c(0,1), ylim = c(0, 0.3))
axis(1, at = c(-0.5, axTicks(1)))
axis(2, at = c(-0.5, axTicks(2)))
points(set2params$eps, set2params$netdiv, cex=0.5)
mtext(expression(mu~'/'~lambda), side=1, line=2.5)
mtext(expression(lambda~'-'~mu), side=2, line=2.5)

dev.off()


