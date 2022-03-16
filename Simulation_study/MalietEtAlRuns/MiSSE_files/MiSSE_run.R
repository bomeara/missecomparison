source("misse_clads_comparison.R") 
setDTthreads(threads=1)
print("passed loading packages")
#cl <- future::makeClusterPSOCK(workers=c(rep(c("10.4.9.34", "10.4.9.45"),48)), rscript="/usr/bin/Rscript")
cl <- future::makeClusterPSOCK(workers=c(rep(c("10.4.9.34"),48)), rscript="/usr/bin/Rscript")


#sessionInfo() 
#future::plan(cluster, workers = cl)
#cache_rerun_b4 <- drake::new_cache(path = "drake_cache_rerun_b4", hash_algorithm = "md5")
#make(New_sim_runs_beaulieu4, cache=cache_rerun_b4, parallelism="future", jobs=96)
#parallel::stopCluster(cl)

#sessionInfo() 
#future::plan(cluster, workers = cl)
#cache_rerun_b4 <- drake::drake_cache(path = "drake_cache_rerun_b4")
#make(Second_sim_runs_beaulieu4, cache=cache_rerun_b4, parallelism="future", jobs=60)
#parallel::stopCluster(cl)

#sessionInfo() 
#future::plan(cluster, workers = cl)
#cache_rerun_b4 <- drake::drake_cache(path = "drake_cache_rerun_b4")
#make(Third_sim_runs_beaulieu4, cache=cache_rerun_b4, parallelism="future", jobs=60)
#parallel::stopCluster(cl)

#sessionInfo() 
#future::plan(cluster, workers = cl)
#cache_rerun_b4 <- drake::drake_cache(path = "drake_cache_rerun_b4")
#make(crashed_runs_beaulieu4, cache=cache_rerun_b4, parallelism="future", jobs=60)
#parallel::stopCluster(cl)

#sessionInfo() 
#future::plan(cluster, workers = cl)
#cache_run_crashed <- drake::new_cache(path = "drake_cache_rerun_crashed", hash_algorithm = "md5")
#make(rerun_crashed_runs_beaulieu4_final, cache=cache_run_crashed, parallelism="future", jobs=96)
#parallel::stopCluster(cl)

sessionInfo() 
future::plan(cluster, workers = cl)
cache_run_crashed <- drake::new_cache(path = "drake_cache_misse_on_clads", hash_algorithm = "md5")
make(run_misse_on_clads, cache=cache_run_crashed, parallelism="future", jobs=96)
parallel::stopCluster(cl)

