source("recover_runs.R") 
setDTthreads(threads=1)
cl <- future::makeClusterPSOCK(workers=c(rep(c("10.4.9.45", "10.4.9.34"),48)), rscript="/usr/bin/Rscript")
#cl <- future::makeClusterPSOCK(workers=c(rep(c("10.4.8.174"),48)), rscript="/usr/bin/Rscript")

sessionInfo() 
future::plan(cluster, workers = cl)
cache_rerun_crashed <- drake::new_cache(path = "drake_cache_rerun_crashed", hash_algorithm = "md5")
make(Recover_crashed_trees, cache=cache_rerun_crashed, parallelism="future", jobs=150)
parallel::stopCluster(cl)

