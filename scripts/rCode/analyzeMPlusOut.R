## First load packages
library(MplusAutomation)
library(foreach)
library(doParallel)
library(tidyverse)
source("./scripts/rCode/hurdleFunctions.r")
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
orig.dir <- getwd()
## Make all sim vals
sim.param.nitems <- c(5,7,10)
sim.param.ncates <- c(3,5,7)
sim.param.discri <- c(2)
sim.param.diffic <- c(-1,1)
sim.param.sample <- c(5000)
sim.param.faccor <- c(.3, .6)
sim.param.critco <- c(0,.3,.6)
sim.param.critco2 <- c(0, .3, .6)
sim.param.difgrm <- c(-1, 0)
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.diffic, sim.param.sample, sim.param.faccor, 
                            sim.param.critco, sim.param.critco2,sim.param.difgrm)


## Now load all of the parameters for an examination of error here
cl <- makeCluster(2) #not to overload your computer
registerDoParallel(cl)
timevals <- foreach(seedVal=1:10, .packages = c("tidyverse", "MplusAutomation")) %dopar% {
  out.directory <- paste(orig.dir, '/data/simdir_',seedVal ,sep='')
  repTime <- list()
  for(i in 1:nrow(all.sim.vals)){
    ## Check if files exist
    out.file3 <- paste(out.directory, '/fileVal_', i,'.RDS', sep='')
    if(file.exists(out.file3)){
      out.dat <- readRDS(out.file3)
      ## now grab the parameter estimates which will be used to estimate error
      all.params <- cbind(out.dat$a_z, out.dat$b_z, out.dat$a, out.dat$b)
      prepVec <- c("true_z_discrim", "est_z_discrim", "true_z_diff", "est_z_diff",
                   "true_grm_discrim", "est_grm_discrim")
      grmDiffVec <- (dim(all.params)[2] - length(prepVec)) / 2
      grmTruVec <- paste("true_grm_diff_", 1:grmDiffVec, sep='')
      grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
      prepVec <- c(prepVec, grmTruVec, grmEstVec)
      colnames(all.params) <- prepVec
      all.params <- data.frame(all.params)
      all.params$seedVal = seedVal
      all.params$repVal = i
      ## Now do the inter factor correlation
      est.fac.cor <- data.frame("TrueInter" = out.dat$varCovMat[2,1], "estInter" = out.dat$estRho)
      est.fac.cor$seedVal <- seedVal
      est.fac.cor$repVal <- i
      ## Now do the regressions
      tmp.dat <- out.dat$theta
      mod <- lm(X3 ~ eap2PL_Hurdle + eapGRM_Hurdle, data = tmp.dat)
      est.crit.cor <- data.frame(summary(mod)$coefficients)
      est.crit.cor$true <- c(0,out.dat$varCovMat[1,3],out.dat$varCovMat[2,3])
      est.crit.cor$seedVal <- seedVal
      est.crit.cor$repVal <- i
      all.out <- list(paramError = all.params, rhoError = est.fac.cor, critError = est.crit.cor)
      repTime[[i]] <- all.out
    }else{
      out.dat <- "ERROR"
    }
  }
  repTime
}
stopCluster(cl)

## Now grab all of the param error vals
step0 <- lapply(timevals, function(x) lapply(x, function(y) y$paramError))
step0 <- bind_rows(step0)
step1 <- lapply(timevals, function(x) lapply(x, function(y) y$rhoError))
step1 <- bind_rows(step1)
step2 <- lapply(timevals, function(x) lapply(x, function(y) y$critError))
step2 <- bind_rows(step2)
## rm intercept
step2 <- step2[grep(pattern = "Hurdle", x = rownames(step2)),]
## Now add rowname as a variable to step2
step2$parameterVal <- rownames(step2)
step2$parameterVal <- strSplitMatrixReturn(charactersToSplit = step2$parameterVal, splitCharacter = "_")[,1]
## Now create the factor corresponding to each of these
all.sim.vals <- data.frame(all.sim.vals)
colnames(all.sim.vals) <- c("nitems", "ncates", "discrim", "diffic2PL", "sample", "faccor", "critco", "critco2", "difgrm")
all.sim.vals$repVal <- 1:nrow(all.sim.vals)
for(i in c("nitems", "ncates", "discrim", "diffic2PL", "sample", "faccor", "critco", "critco2", "difgrm")){
  all.sim.vals[,i] <- factor(all.sim.vals[,i])
}
## Now combine these
step0 <- merge(step0, all.sim.vals, by=c("repVal"), all=TRUE)
step1 <- merge(step1, all.sim.vals, by=c("repVal"), all=TRUE)
step2 <- merge(step2, all.sim.vals, by=c("repVal"), all=TRUE)

## Now do an anova for step0 here
step0$absError <- abs(step0$true_z_discrim - step0$est_z_discrim)
mod <- lm(absError ~ (nitems+ncates+diffic2PL+faccor+difgrm)^2, data = step0)
step0$absError <- abs(step0$true_grm_discrim - step0$est_grm_discrim)
mod <- lm(absError ~ (nitems+ncates+diffic2PL+faccor+difgrm)^2, data = step0)

## Now do the anova for step1 here
step1$absError <- abs(step1$TrueInter - step1$estInter)
mod <- lm(absError ~ (nitems+ncates+diffic2PL+faccor+difgrm)^2, data = step1)


## Now do the anovas for step2 here
step2$absError <- abs(step2$Estimate - as.numeric(as.character(step2$critco)))
step2$absError[which(step2$parameterVal=="eapGRM")] <- abs(step2$Estimate[which(step2$parameterVal=="eapGRM")] - as.numeric(as.character(step2$critco2[which(step2$parameterVal=="eapGRM")])))
mod <- lm(absError ~ -1 + (nitems+ncates+diffic2PL+faccor+difgrm+critco+critco2)^4, data = step2)
