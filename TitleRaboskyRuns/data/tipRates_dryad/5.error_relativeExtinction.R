# analyze error in rates as a function of relative extinction
# compare speciation rate and net diversification rate

require(TeachingDemos)

dat <- read.csv('dataFiles/rateErrorMetrics.csv', stringsAsFactors=FALSE)
treeDat <- read.csv('dataFiles/treeSummary.csv', stringsAsFactors=FALSE)
# ---------------------------------------------------------
# ---------------------------------------------------------
## error of lambda and lambda - mu against eps
# first for single regime trees


# keep only single regime trees
dat_lambda <- dat[which(dat$setname == 'lambdaConstantVariance'),]
dat_netdiv <- dat[which(dat$setname == 'netDivConstantVariance'),]

# get true relative extinction rate
epsVec_lambda <- treeDat[which(treeDat$treeName %in% dat_lambda$treeName),]
epsVec_lambda <- setNames(epsVec_lambda$mu0 / epsVec_lambda$lambda0, epsVec_lambda$treeName)
epsVec_lambda <- epsVec_lambda[dat_lambda$treeName]

epsVec_netdiv <- treeDat[which(treeDat$treeName %in% dat_netdiv$treeName),]
epsVec_netdiv <- setNames(epsVec_netdiv$mu0 / epsVec_netdiv$lambda0, epsVec_netdiv$treeName)
epsVec_netdiv <- epsVec_netdiv[dat_netdiv$treeName]


tipMetrics <- c('lambdaND', 'lambdaDR', 'lambdaBAMM')
tipMetricsNetDiv <- c('lambdaND', 'lambdaDR', 'netDivBAMM')
errorMetricNames <- c(PE1='Prop. Error 1', PE2='Prop. Error 2', PE3='Prop. Error 3', PE4='Prop. 4', r2='OLS R2', slope='OLS slope', absoluteError='mean absolute error')

tipMetricLabel <- c(expression(bold(paste(lambda['ND']))), expression(bold(paste(lambda['DR']))),expression(bold(paste(lambda['BAMM']))))


###################################################

# mean absolute error
errorMetric <- 'absoluteError'

yrange <- c(0,5)
insetYrange <- c(0, 0.5)

pdf('fig1.pdf', width=9, height=5)

layout(matrix(1:6, nrow=2,ncol=3, byrow=T))
par(mar = c(3,3,0,0), oma=c(1,1,1,0))


for (i in 1:3) {
	
	tipMetric <- tipMetrics[i]
	errorLabel <- expression(paste("mean absolute error in ", lambda))
	errorCol <- Reduce(intersect, list(grep('tipLambda', colnames(dat)), grep(tipMetric, colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(0,1), ylim = yrange)
	axis(1, at=c(-0.1, axTicks(1)), cex.axis = 0.7)
	axis(2, at=c(-5, axTicks(2)), cex.axis = 0.7)
	points(epsVec_lambda, dat_lambda[, errorCol], cex=0.3)
	if (i == 1) mtext(errorLabel, side = 2, cex=0.7, line = 2.5)
	title(main = tipMetricLabel[i], line = -1, cex.main=1.5)
	abline(h=0, lty=2, col=gray(0.5), lwd=1)
	
	insetMap <- function() {
		plot.new()
		plot.window(xlim = c(0,1), ylim = insetYrange)
		axis(1, at=c(-0.1, axTicks(1)), cex.axis = 0.7, tcl = -0.2, mgp = c(3,0.3,0), las = 1)
		axis(2, at=c(-5, axTicks(2)), cex.axis = 0.7, tcl = -0.2, mgp = c(3,0.3,0), las = 1)
		points(epsVec_lambda, dat_lambda[, errorCol], cex=0.3)
	}
	
	subplot(insetMap(), x=0.5, y=1.5, vadj=0, hadj=0, size=c(1,1), type='plt')
}

# second row: net div
for (i in 1:3) {
	
	tipMetric <- tipMetricsNetDiv[i]
	errorLabel <- expression(paste("mean absolute error in ", italic('r')))
	errorCol <- Reduce(intersect, list(grep('tipNetDiv', colnames(dat)), grep(tipMetric, colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(0,1), ylim = yrange)
	axis(1, at=c(-0.1, axTicks(1)), cex.axis=0.7)
	axis(2, at=c(-5, axTicks(2)), cex.axis=0.7)
	points(epsVec_netdiv, dat_netdiv[, errorCol], cex=0.3)
	mtext(expression(mu~'/'~lambda), side = 1, cex=0.8, line = 2.5)
	if (i == 1) mtext(errorLabel, side = 2, cex=0.7, line = 2.5)
	if (tipMetric == 'netDivBAMM') {
		title(main = expression(bold(paste(italic(r)['BAMM']))), line = -2, cex.main=1.5)
	}
	abline(h=0, lty=2, col=gray(0.5), lwd=1)		
}

dev.off()


###################################################

# RMSE: Root mean square error
errorMetric <- 'rmse'

yrange <- c(0,5)
insetYrange <- c(0, 0.5)

pdf('figS2.pdf', width=9, height=5)

layout(matrix(1:6, nrow=2,ncol=3, byrow=T))
par(mar = c(3,3,0,0), oma=c(1,1,1,0))


for (i in 1:3) {
	
	tipMetric <- tipMetrics[i]
	errorLabel <- expression(paste("RMSE in ", lambda))
	errorCol <- Reduce(intersect, list(grep('tipLambda', colnames(dat)), grep(tipMetric, colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(0,1), ylim = yrange)
	axis(1, at=c(-0.1, axTicks(1)), cex.axis = 0.7)
	axis(2, at=c(-5, axTicks(2)), cex.axis = 0.7)
	points(epsVec_lambda, dat_lambda[, errorCol], cex=0.3)
	if (i == 1) mtext(errorLabel, side = 2, cex=0.7, line = 2.5)
	title(main = tipMetricLabel[i], line = -1, cex.main=1.5)
	abline(h=0, lty=2, col=gray(0.5), lwd=1)
	
	insetMap <- function() {
		plot.new()
		plot.window(xlim = c(0,1), ylim = insetYrange)
		axis(1, at=c(-0.1, axTicks(1)), cex.axis = 0.7, tcl = -0.2, mgp = c(3,0.3,0), las = 1)
		axis(2, at=c(-5, axTicks(2)), cex.axis = 0.7, tcl = -0.2, mgp = c(3,0.3,0), las = 1)
		points(epsVec_lambda, dat_lambda[, errorCol], cex=0.3)
	}
	
	subplot(insetMap(), x=0.5, y=1.5, vadj=0, hadj=0, size=c(1,1), type='plt')
}

# second row: net div
for (i in 1:3) {
	
	tipMetric <- tipMetricsNetDiv[i]
	errorLabel <- expression(paste("RMSE in ", italic('r')))
	errorCol <- Reduce(intersect, list(grep('tipNetDiv', colnames(dat)), grep(tipMetric, colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(0,1), ylim = yrange)
	axis(1, at=c(-0.1, axTicks(1)), cex.axis=0.7)
	axis(2, at=c(-5, axTicks(2)), cex.axis=0.7)
	points(epsVec_netdiv, dat_netdiv[, errorCol], cex=0.3)
	mtext(expression(mu~'/'~lambda), side = 1, cex=0.8, line = 2.5)
	if (i == 1) mtext(errorLabel, side = 2, cex=0.7, line = 2.5)
	if (tipMetric == 'netDivBAMM') {
		title(main = expression(bold(paste(italic(r)['BAMM']))), line = -2, cex.main=1.5)
	}
	abline(h=0, lty=2, col=gray(0.5), lwd=1)		
}

dev.off()

#######################
# PE1 on log scale

errorMetric <- 'PE1'

yrangeTop <- c(-0.5,4)
yrangeBottom <- c(-1.5,4)
insetYrange <- c(-0.4,0.6)
constant <- 4
width <- 0.03
correction <- log(constant)

pdf('figS3.pdf', width=9, height=5)

layout(matrix(1:6, nrow=2,ncol=3, byrow=T))
par(mar = c(3,3,0,0), oma=c(1,1,1,0))


for (i in 1:3) {
	
	tipMetric <- tipMetrics[i]
	errorLabel <- errorLabel <- expression(paste("log prop. error in ", lambda))
	errorCol <- Reduce(intersect, list(grep('tipLambda', colnames(dat)), grep(tipMetric, colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(0,1), ylim = yrangeTop)
	axis(1, at=c(-0.1, axTicks(1)), cex.axis = 0.7)
	axis(2, at=c(-5, axTicks(2)), cex.axis = 0.7)
	points(epsVec_lambda, log(dat_lambda[, errorCol] + constant) - correction, cex=0.3)
	if (i == 1) mtext(errorLabel, side = 2, cex=0.7, line = 2.5)
	title(main = tipMetricLabel[i], line = -1, cex.main=1.5)
	abline(h=log(0 + constant) - correction, lty=2, col=gray(0.5), lwd=1)
	
	insetMap <- function() {
		plot.new()
		plot.window(xlim = c(0,1), ylim = insetYrange)
		axis(1, at=c(-0.1, axTicks(1)), cex.axis = 0.7, tcl = -0.2, mgp = c(3,0.3,0), las = 1)
		axis(2, at=c(-5, axTicks(2)), cex.axis = 0.7, tcl = -0.2, mgp = c(3,0.3,0), las = 1)
		points(epsVec_lambda, log(dat_lambda[, errorCol] + constant) - correction, cex=0.3)
		abline(h=log(0 + constant) - correction, lty=2, col=gray(0.5), lwd=1)
	}
	
	subplot(insetMap(), x=0.5, y=1, vadj=0, hadj=0, size=c(1,1), type='plt')
	
}

# second row: net div
for (i in 1:3) {
	
	tipMetric <- tipMetricsNetDiv[i]
	errorLabel <- expression(paste("log prop. error in ", italic('r')))
	errorCol <- Reduce(intersect, list(grep('tipNetDiv', colnames(dat)), grep(tipMetric, colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(0,1), ylim = yrangeBottom)
	axis(1, at=c(-0.1, axTicks(1)), cex.axis = 0.7)
	axis(2, at=c(-5, axTicks(2)), cex.axis = 0.7)
	points(epsVec_netdiv, log(dat_netdiv[, errorCol] + constant) - correction, cex=0.3)
	mtext(expression(mu~'/'~lambda), side = 1, cex=0.8, line = 2.5)
	if (i == 1) mtext(errorLabel, side = 2, cex=0.7, line = 2.5)
	if (tipMetric == 'netDivBAMM') {
		title(main = expression(bold(paste(italic('r')['BAMM']))), line = -2, cex.main=1.5)
	}
	abline(h=log(0 + constant) - correction, lty=2, col=gray(0.5), lwd=1)
}

dev.off()


#######################

# Separate figure for lambdaTB, compared to lambdaDR and lambdaBAMM
tipMetrics <- c('lambdaTB', 'lambdaDR', 'lambdaBAMM')
tipMetricsNetDiv <- c('lambdaTB', 'lambdaDR', 'netDivBAMM')
errorMetricNames <- c(PE1='Prop. Error 1', PE2='Prop. Error 2', PE3='Prop. Error 3', PE4='Prop. 4', r2='OLS R2', slope='OLS slope', absoluteError='mean absolute error')

tipMetricLabel <- c(expression(bold(paste(lambda['TB']))), expression(bold(paste(lambda['DR']))),expression(bold(paste(lambda['BAMM']))))

# mean absolute error
errorMetric <- 'absoluteError'

yrange <- c(0,50)
insetYrange <- c(0, 30)
width <- 0.03

pdf('figS4.pdf', width=9, height=5)
	
layout(matrix(1:6, nrow=2,ncol=3, byrow=T))
par(mar = c(3,3,0,0), oma=c(1,1,1,0))


for (i in 1:3) {
	
	tipMetric <- tipMetrics[i]
	errorLabel <- expression(paste("mean absolute error in ", lambda))
	errorCol <- Reduce(intersect, list(grep('tipLambda', colnames(dat)), grep(tipMetric, colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(0,1), ylim = yrange)
	axis(1, at=c(-0.1, axTicks(1)), cex.axis = 0.7)
	axis(2, at=c(-5, axTicks(2)), cex.axis = 0.7)
	points(epsVec_lambda, dat_lambda[, errorCol], cex=0.3)
	if (i == 1) mtext(errorLabel, side = 2, cex=0.7, line = 2.5)
	title(main = tipMetricLabel[i], line = -1, cex.main=1.5)
	abline(h=0, lty=2, col=gray(0.5), lwd=1)
	
}

# second row: net div
for (i in 1:3) {
	
	tipMetric <- tipMetricsNetDiv[i]
	errorLabel <- expression(paste("mean absolute error in ", italic('r')))
	errorCol <- Reduce(intersect, list(grep('tipNetDiv', colnames(dat)), grep(tipMetric, colnames(dat)), grep(errorMetric, colnames(dat))))
	plot.new()
	plot.window(xlim = c(0,1), ylim = yrange)
	axis(1, at=c(-0.1, axTicks(1)), cex.axis=0.7)
	axis(2, at=c(-5, axTicks(2)), cex.axis=0.7)
	points(epsVec_netdiv, dat_netdiv[, errorCol], cex=0.3)
	mtext(expression(mu~'/'~lambda), side = 1, cex=0.8, line = 2.5)
	if (i == 1) mtext(errorLabel, side = 2, cex=0.7, line = 2.5)
	if (tipMetric == 'netDivBAMM') {
		title(main = expression(bold(paste(italic(r)['BAMM']))), line = -2, cex.main=1.5)
	}
	abline(h=0, lty=2, col=gray(0.5), lwd=1)	
	
}

dev.off()











# for speciation rate, some statistics
mean(dat_lambda$tipLambda_lambdaTB_absoluteError)
mean(dat_lambda$tipLambda_lambdaND_absoluteError)
mean(dat_lambda$tipLambda_lambdaDR_absoluteError)
mean(dat_lambda$tipLambda_lambdaBAMM_absoluteError)

sd(dat_lambda$tipLambda_lambdaTB_absoluteError)
sd(dat_lambda$tipLambda_lambdaND_absoluteError)
sd(dat_lambda$tipLambda_lambdaDR_absoluteError)
sd(dat_lambda$tipLambda_lambdaBAMM_absoluteError)

quantile(dat_lambda$tipLambda_lambdaDR_absoluteError, 0.95)
quantile(dat_lambda$tipLambda_lambdaBAMM_absoluteError, 0.95)		

quantile(dat_netdiv$tipNetDiv_lambdaDR_absoluteError, 0.95)
quantile(dat_netdiv$tipNetDiv_netDivBAMM_absoluteError, 0.95)		

