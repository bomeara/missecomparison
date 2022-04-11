# install.packages("BAMMtools")
library(BAMMtools)
library(phytools)

one_tree <- read.tree("Lupinus.tre")
priors <- setBAMMpriors(one_tree, outfile=NULL)
generateControlFile(file = "Lupinus_divcontrol.txt", 
                      params = list(treefile = paste0(getwd(),"/Lupinus.tre"),
                                     runInfoFilename = paste0(getwd(), "/Lupinus_run_info.txt"),
                                     mcmcOutfile = paste0(getwd(), "/Lupinus_mcmc_out.txt"),
                                     eventDataOutfile = paste0(getwd(), "/Lupinus_event_data.txt"),
                                     globalSamplingFraction = '1',
                                     numberOfGenerations = '1000000',
                                     overwrite = '1',
                                     lambdaInitPrior = as.numeric(priors['lambdaInitPrior']), lambdaShiftPrior = as.numeric(priors['lambdaShiftPrior']), muInitPrior = as.numeric(priors['muInitPrior']), expectedNumberOfShifts = '1'),
)

library(coda)
mcmcout <- read.csv("Lupinus_mcmc_out.txt", header=T)
#plot(mcmcout$logLik ~ mcmcout$generation)  
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
tmp1 <- effectiveSize(postburn$N_shifts)
if(tmp1 < 200) {
  print(paste0("rerun tree ", label))
}
tmp2 <- effectiveSize(postburn$logLik)
if(tmp2 < 200) {
  print(paste0("rerun tree ", label))
}
  

