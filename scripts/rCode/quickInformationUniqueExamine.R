min.vals <- seq(1, 76800, by = 768)
min.vals <- c(min.vals, 76800)

## now go through the first 10 and launch these to run in parallel
command.list <- list()
for(i in 1:9){
  ## prep the system call
  max.val <- min.vals[i] + 767
  command.val <- paste("for i in `seq ", min.vals[i], max.val, " | shuf` ; do Rscript ./scripts/rCode/informationUniqueExamine.R ${i} ; done")
  command.list[[i]] <- command.val
}

## Now launch these in parallel
parallel::mclapply(command.list, function(x) system(x), mc.cores = 3)