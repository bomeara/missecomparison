source("new_runs.R") 
setDTthreads(threads=1)
print("passed loading packages")
#cl <- future::makeClusterPSOCK(workers=c(rep(c("10.4.9.34", "10.4.9.45"),48)), rscript="/usr/bin/Rscript")
cl <- future::makeClusterPSOCK(workers=c(rep(c("10.4.9.45"),96)), rscript="/usr/bin/Rscript")


#sessionInfo() 
#future::plan(cluster, workers = cl)
#cache_rerun_b3 <- drake::new_cache(path = "drake_cache_rerun_b3", hash_algorithm = "md5")
#make(New_sim_runs_beaulieu3, cache=cache_rerun_b3, parallelism="future", jobs=96)
#parallel::stopCluster(cl)

sessionInfo() 
future::plan(cluster, workers = cl)
#cache_rerun_b3 <- drake::new_cache(path = "drake_cache_rerun_b3", hash_algorithm = "md5")
make(Second_sim_runs_beaulieu3, cache=cache_rerun_b3, parallelism="future", jobs=40)
parallel::stopCluster(cl)

