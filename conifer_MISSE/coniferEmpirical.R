

library(hisse)
library(parallel)


GetModelPars <- function(shift.no) {
    turnover <- 1:shift.no
    eps <- 1:shift.no
    par.settings <- NULL
    par.settings$turnover <- turnover
    par.settings$eps <- eps
    return(par.settings)
}


RunMiSSE <- function(phy, f, n.cores=10){

    #Step 1: Get a set of 20 random starting points.
    start.vals <- hisse:::starting.point.generator(phy, samp.freq.tree=f, k=1)
    turn.tries <- c(sum(start.vals[1:2]), exp(rnorm(19, log(sum(start.vals)), 0.25)))
    eps.tries <- c(start.vals[2]/start.vals[1], runif(19, 0.01, 0.99))

    print(turn.tries)
    print(eps.tries)

    #This begins looping through the model space
    for(model.index in 2:12){
        print(model.index)
        RandomStartVal <- function(iteration, model.index, phy, f, turnover, eps, turn.tries, eps.tries, trans.tries){
            print(iteration)
            tmp <- MiSSE(phy=phy, f=f, turnover=turnover, eps=eps, starting.vals=c(turn.tries[iteration], eps.tries[iteration], trans.tries), root.type="herr_als")
            save(tmp, file=paste(model.index, iteration, "e.Rsave", sep="."))
        }
        par.set <- GetModelPars(model.index)
        startTries <- mclapply(1:20, RandomStartVal, model.index=model.index, phy=phy, f=f, turnover=par.set$turnover, eps=par.set$eps, turn.tries=turn.tries, eps.tries=eps.tries, trans.tries=0.001, mc.cores=n.cores)
    }

}

#phy <- read.tree("Leslieetal2018.tre")
#RunMiSSE(phy=phy, f=0.9, n.cores=20)


RunMiSSERecon <- function(with.eps=FALSE, num.models, n.cores){
    
    for(rate.class in 2:num.models){
        
        if(with.eps == TRUE){
            model.restarts <- system(paste("ls -1 ", rate.class, "*..e.Rsave", sep = ""), intern = TRUE)
        }else{
            model.restarts.weps <- system(paste("ls -1 ", rate.class, "*..e.Rsave", sep = ""), intern = TRUE)
            model.restarts.all <- system(paste("ls -1 ", rate.class, "*.Rsave", sep = ""), intern = TRUE)
            model.restarts <- model.restarts.all[which(!model.restarts.all %in% model.restarts.weps)]
        }
	print(model.restarts)        
        res <- c()
        for(file.index in 1:length(model.restarts)){
            load(model.restarts[file.index])
            res <- c(res, tmp$loglik)
        }
        
        load(model.restarts[which.max(res)])
        pp <- MarginReconMiSSE(phy=tmp$phy, f=tmp$f, pars=tmp$solution, hidden.states=tmp$hidden.states, condition.on.survival=TRUE, root.type="herr_als", aic=tmp$AIC, verbose=TRUE, n.cores=n.cores)
        save(pp, file=paste("recon.bestmodel.weps", rate.class, ".Rsave", sep=""))
    }
}
RunMiSSERecon(with.eps=TRUE, num.models=12, n.cores=20)




#logliks <- c()
#for(shift.index in 2:20){
#model.restarts <- system(paste("ls -1 ", "", paste(shift.index, ".*", sep = ""), sep=""), intern = TRUE)

#res <- c()
#for(file.index in 1:length(model.restarts)){
#  load(model.restarts[file.index])
#  res <- c(res, tmp$loglik)
#}
#logliks <-c(logliks, res[which.max(res)])
#}



