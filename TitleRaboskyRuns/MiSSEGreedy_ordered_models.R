MiSSEGreedy <- function(phy, f=1, possible.combos = generateMiSSEGreedyCombinations(), stop.deltaAICc=10, save.file=NULL, n.cores=NULL, chunk.size=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils=FALSE, k.samples=NULL, sann=TRUE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0) {
  
  misse.list <- list()
  chunk.size <- ifelse(is.null(chunk.size),ifelse(is.null(n.cores),1,n.cores), chunk.size)
  total.chunks <- ceiling(nrow(possible.combos)/chunk.size)
  possible.combos$lnL <- NA
  possible.combos$AIC <- NA
  possible.combos$AICc <- NA
  possible.combos$deltaAICc <- NA
  possible.combos$elapsedMinutes <- NA
  possible.combos$predictedMinutes <- NA
  
  for (batch_index in sequence(total.chunks)) { # So, if we can do parallel, we do it in chunks so all cores are busy
    starting.time <- Sys.time()
    local.combos <- possible.combos[(1+chunk.size*(batch_index-1)):min(nrow(possible.combos), chunk.size*batch_index) ,]
    
    #cat("\nNow starting run with", paste(range(local.combos$turnover), collapse="-"), "turnover categories and", paste(range(local.combos$eps), collapse="-"), "extinction fraction categories", "\n")
    cat("Starting at ", as.character(starting.time), "\n running on ", n.cores, " cores.", sep="")
    cat("\n")
    print(local.combos[,1:3])
    cat("\n")
    
    misse.list <- append(misse.list, parallel::mcmapply(
      MiSSE,
      eps=local.combos$eps,
      turnover=local.combos$turnover,
      fixed.eps=local.combos$fixed.eps,
      MoreArgs=list(
        phy=phy,
        f=f,
        condition.on.survival=condition.on.survival,
        root.type=root.type,
        root.p=root.p,
        includes.fossils=includes.fossils,
        k.samples=k.samples,
        sann=sann,
        sann.its=sann.its,
        bounded.search=bounded.search,
        max.tol=max.tol,
        starting.vals=starting.vals,
        turnover.upper=turnover.upper,
        eps.upper=eps.upper,
        trans.upper=trans.upper,
        restart.obj=restart.obj,
        ode.eps=ode.eps,
        expand.mode=TRUE
      ),
      mc.cores=ifelse(is.null(n.cores),1,n.cores),
      SIMPLIFY=FALSE
    ))
    
    AICc <- unlist(lapply(misse.list, "[[", "AICc"))
    deltaAICc <- AICc-min(AICc)
    min.deltaAICc.this.chunk <- min(deltaAICc[(1+chunk.size*(batch_index-1)):min(nrow(possible.combos), chunk.size*batch_index)])
    
    possible.combos$lnL[1:min(nrow(possible.combos), chunk.size*batch_index)] <- unlist(lapply(misse.list, "[[", "loglik"))
    possible.combos$AIC[1:min(nrow(possible.combos), chunk.size*batch_index)] <- unlist(lapply(misse.list, "[[", "AIC"))
    possible.combos$AICc[1:min(nrow(possible.combos), chunk.size*batch_index)] <- AICc
    possible.combos$deltaAICc[1:min(nrow(possible.combos), chunk.size*batch_index)] <- deltaAICc
    possible.combos$elapsedMinutes[1:min(nrow(possible.combos), chunk.size*batch_index)] <- unlist(lapply(misse.list, "[[", "elapsed.minutes"))
    
    data.for.fit <- data.frame(nparam=(possible.combos$eps+possible.combos$turnover)[1:min(nrow(possible.combos), chunk.size*batch_index)], logmin=log(possible.combos$elapsedMinutes[1:min(nrow(possible.combos), chunk.size*batch_index)]))
    data.for.prediction <- data.frame(nparam=(possible.combos$eps+possible.combos$turnover))
    suppressWarnings(possible.combos$predictedMinutes <- exp(predict(lm(logmin ~ nparam, data=data.for.fit), newdata=data.for.prediction)))
    
    cat("\nResults so far\n")
    print(round(possible.combos,2))
    
    if(!is.null(save.file)) {
      save(misse.list, possible.combos, file=save.file)
    }
    
    
    if(batch_index<total.chunks) {
      if(stop.deltaAICc>min.deltaAICc.this.chunk) {
        print(paste0("Best AICc in this set of parallel runs was ", round(min.deltaAICc.this.chunk,2), " which is less than the cutoff to stop running (",stop.deltaAICc,"), so starting another set of parallel runs"))
      } else {
        print(paste0("Best AICc in this set of parallel runs was ", round(min.deltaAICc.this.chunk,2), " which is greater than the cutoff to stop running (",stop.deltaAICc,"), so stopping here"))
        break()
      }
    }
    
    # print("\n")
    # print(local.combos)
    # misse.list <- append(misse.list, MiSSE(
    #     eps=local.combos$eps[1],
    #     turnover=local.combos$turnover[1],
    #     fixed.eps=local.combos$fixed.eps[1],
    #
    #         phy=phy,
    #         f=f,
    #         condition.on.survival=condition.on.survival,
    #         root.type=root.type,
    #         root.p=root.p,
    #         sann=sann,
    #         sann.its=sann.its,
    #         bounded.search=bounded.search,
    #         max.tol=max.tol,
    #         starting.vals=starting.vals,
    #         turnover.upper=turnover.upper,
    #         eps.upper=eps.upper,
    #         trans.upper=trans.upper,
    #         restart.obj=restart.obj,
    #         ode.eps=ode.eps,
    #         expand.mode=TRUE
    # ))
  }
  return(misse.list)
}