# Not require password
# git config --global credential.helper cache
# git config --global credential.helper 'cache --timeout=360000'

library(drake)

while (1<2) {
  drake::loadd(drake::cached())
  result.df <- dplyr::bind_rows(lapply(drake::cached(), get))
  write.csv(result.df, file="result.csv")
  system(paste0('git commit -m"updating with ', length(drake::cached()), ' runs done" result.csv'))
  try(system("git push"))
  Sys.sleep(60*60)
}
