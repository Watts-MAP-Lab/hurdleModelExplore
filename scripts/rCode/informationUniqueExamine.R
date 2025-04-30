## This script will be used to run the winosrization hurdle model explore
## It is going to require simulation a zero-inflated test response pattern
## and changing response patterns by collapsing response categories into one another
## This will be performed using three separate strategies:
##  1. Winsorization: Taking the highest response category and collapsing it into the next highest
##  2. Rec Winsor: Collapsing the lowest response category into the 0 response
##  3. Total information preservation: Identify the difficulty parameter that by removing would preserve the total information; total information will be taken as the integral of the test informaiton function

##### --load-library-------------
source("./scripts/rCode/hurdleFunctions.r")
suppressPackageStartupMessages(library("tidyverse"))
library("numDeriv")
suppressPackageStartupMessages(library(MplusAutomation))
library("mirt")
library("psych")

##### --declare-sim-params-------
## Sim params will need to be modified at a later time point
sim.param.nitems <- c(6,16)
sim.param.ncates <- c(3,5)
sim.param.discri <- c(.3,1.5)
sim.param.2pl.spread <- c(1.2)
sim.param.sample <- c(15000)
sim.param.faccor <- c(.2,.8)
sim.param.difgrmF <- c(-3,-.5)
sim.param.difgrmC <- c(2)
sim.param.dif2pl <- c(-3,-1,1)
sim.param.discri2 <- sim.param.discri
sim.iter <- 1:50
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.2pl.spread,sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.param.discri2,sim.param.dif2pl, sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCat", "discrim2pl", "diffSpread", "n", "facCor", "difGrmF","difGrmC","grmDiscrim","dif2PL","iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)
seedVal <- as.integer(commandArgs(TRUE))
#seedVal <- 1
set.seed(seedVal)

### --check-run-status----------
out.file = paste("./data/hurdleCollapse/seedVal_", seedVal, "_take3.RDS", sep='')
if(file.exists(out.file)){
  print("Done")
  q()
}
## --run-test-sim ------------
## I will run through a single test sim here and see how I can best, and most efficiently run this
## I really need to focus on the memory leakage problem the R has, especially when running large loops
## First create the data -- this will start with the maximum full dataset, 9 total response categories, full range of difficulty parameters
## This will also showcase where I need to streamline code with custom functions
## Run a spread check for discrim values
add.val.2pl <- 1.2
add.val.grm <- 1.2
if(all.sim.vals$grmDiscrim[seedVal]==1.5){add.val.grm=2}
if(all.sim.vals$discrim2pl[seedVal]==1.5){add.val.2pl=2}
a = runif(n = all.sim.vals$nItem[seedVal], min = all.sim.vals$grmDiscrim[seedVal], all.sim.vals$grmDiscrim[seedVal] + add.val.grm)
b = genDiffGRM(num_items = all.sim.vals$nItem[seedVal], num_categories = all.sim.vals$nCat[seedVal], min = all.sim.vals$difGrmF[seedVal], max = all.sim.vals$difGrmF[seedVal] + 2.5, rnorm_var = .1)
a_z = runif(n = all.sim.vals$nItem[seedVal], min = all.sim.vals$discrim2pl[seedVal], all.sim.vals$discrim2pl[seedVal] + add.val.2pl)
## Need to generate 4 separate b_z levels
b_z1 = runif(all.sim.vals$nItem[seedVal], min = all.sim.vals$dif2PL[seedVal], max=all.sim.vals$dif2PL[seedVal]+all.sim.vals$diffSpread[seedVal])

muVals = c(0,0)
rho <- all.sim.vals$facCor[seedVal]
varCovMat = matrix(c(1,rho,rho,1), ncol=2)
N = all.sim.vals$n[seedVal]
## Now generate theta here
theta = MASS::mvrnorm(n = N, mu = muVals, Sigma = varCovMat)
reps1 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta)

## Grab the MIRT estimates here
model <- "
  F1 = 1-YYZ
  F2 = BBH-FNF
  COV = F1*F2
"
model <- gsub(x = model, pattern = "YYZ", replacement = all.sim.vals$nItem[seedVal])
model <- gsub(x = model, pattern = "BBH", replacement = all.sim.vals$nItem[seedVal]+1)
model <- gsub(x = model, pattern = "FNF", replacement = ncol(reps1$mplusMat))
item.type.rep <- c(rep("2PL", all.sim.vals$nItem[seedVal]), rep("graded", all.sim.vals$nItem[seedVal]))
sv1 <- mirt(data.frame(reps1$mplusMat), model = model, itemtype = item.type.rep)

loadDatFromMplus <- function(trueVals){
  ## Load Adon's custom scripts
  source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
  ## First prep the output data directory
  all.params <- cbind(1:ncol(trueVals$responses),trueVals$varCovMat[1,2], trueVals$a_z, trueVals$b_z,
                          trueVals$a, trueVals$b)
  ## Now do the true names
  prepVec <- c("item","rhoTrue","true_z_discrim","true_z_diff","true_grm_discrim")
  grmDiffVec <- (dim(all.params)[2] - length(prepVec))
  grmEstVec <- paste("true_grm_diff_", 1:grmDiffVec, sep='')
  prepVec2 <- c(prepVec, grmEstVec)
  colnames(all.params) <- c(prepVec2)
  all.params <- data.frame(all.params)
  all.params
  return(all.params)
}

vals1 <- loadDatFromMplus(reps1)

## Now go through each of these and add the estimated test reliability from alpha and omega
vals_all <- NULL
for(i in 1){
  ## First get the values
  vals_loop <- get(paste("vals", i, sep=''))
  #vals_loopb <- get(paste("vals", i, "b",sep=''))
  mod_loopM <- get(paste("sv",i,sep=''))
  ## Now merge these
  true.cols <- grep(pattern = "true", names(vals_loop), ignore.case = TRUE, value = TRUE)
  rep_loop <- get(paste("reps", i, sep=''))
  ## Now estimate alpha & omega
  test_dat <- as.data.frame(rep_loop$responses)
  #rel.alpha <- psych::alpha(test_dat)
  rel.alpha <- psych::alpha(polychoric(test_dat)$rho)
  rel.ome <- psych::omega(test_dat, poly=TRUE, nfactors = 3, plot = FALSE)
  rel.uni <- psych::unidim(test_dat, cor="poly")
  vals_loop$omega_h <- rel.ome$omega_h
  vals_loop$alpha <- rel.alpha$total$raw_alpha
  vals_loop$omega_t <- rel.ome$omega.tot
  vals_loop$omega_h <- rel.ome$omega_h
  vals_loop$G_six <- rel.ome$G6
  cor.mat <- psych::polychoric(rep_loop$responses)
  split.half.rel <- suppressWarnings(psych::guttman(r = cor.mat$rho))
  vals_loop$lambda1Rel = split.half.rel$lambda.1
  vals_loop$lambda2Rel = split.half.rel$lambda.2
  vals_loop$lambda3Rel = split.half.rel$lambda.3
  vals_loop$lambda4Rel = split.half.rel$lambda.4
  vals_loop$lambda5Rel = split.half.rel$lambda.5
  vals_loop$lambda6Rel = split.half.rel$lambda.6
  vals_loop$alpheFromOme <- rel.ome$alpha
  vals_loop$unidim <- rel.uni$uni["u"]
  
  ##  Now grab the true reliability values
  a = vals_loop$true_grm_discrim
  b = data.matrix(data.frame(vals_loop[,grep(pattern = "true_grm_diff", x = names(vals_loop))]))
  a_z = vals_loop$true_z_discrim
  b_z = vals_loop$true_z_diff
  vals_loop$trueHurdleRel <- hurdInfo(theta.grid = expand.grid(seq(-3, 3, .2), seq(-3, 3, .2)), a = a, b = b, a_z = a_z, b_z = b_z, muVals = muVals, rhoVal = rho)$out.rel

  ## Now organize the MIRT values here
  mirt.coef <- coef(mod_loopM, IRTpars=TRUE)
  ## First grab all of the a vals
  sus.vals <- grep(pattern = "Sus", x = names(mirt.coef))
  a_z <- unlist(lapply(mirt.coef[sus.vals], function(x) x[1]))
  b_z <- unlist(lapply(mirt.coef[sus.vals], function(x) x[3]))
  sev.vals <- grep(pattern = "Sev", x = names(mirt.coef))
  a <- unlist(lapply(mirt.coef[sev.vals], function(x) x[2]))
  b <- lapply(mirt.coef[sev.vals], function(x) x[-c(1:2)])
  ## First check the length of all b estimates
  ## if any are shorter -- add a large value to these b estimates
  b_check <- table(unlist(lapply(b, length)))
  if(length(b_check)>1){
    ## Find all difficulty values w/o the correct number of values
    length.vec <- unlist(lapply(b, length))
    correct.val <- max(length.vec)
    ## Grab all values w/ length != 4
    index <- which(length.vec != correct.val)
    for(k in index){
      ident.val <- paste("Sev_", k, sep='')
      length.to.add <- correct.val- length(b[[ident.val]])
      add.vals <- seq(4, 8, length.out = length.to.add)
      b[[ident.val]] <- c(b[[ident.val]], add.vals)
    }
  }
  b <- t(bind_rows(b))
  #b <- t(apply(b, 1, function(x) sort(x, decreasing = FALSE)))
  rhoEst <- unique(mirt.coef$GroupPars["par","COV_21"])
  vals_loopM <- bind_cols(a_z, b_z,a, b)
  colnames(vals_loopM)[1:3] <- c("est_z_discrim", "est_z_diff", "est_grm_discrim")
  colnames(vals_loopM)[-c(1:3)] <- paste("est_grm_diff_", 1:dim(vals_loopM[,-c(1:3)])[2], sep='')
  vals_loopM$estHurdleRel <- hurdInfo(theta.grid = expand.grid(seq(-3, 3, .2), seq(-3, 3, .2)), a = a, b = b, a_z = a_z, b_z = b_z, muVals = muVals, rhoVal = rhoEst)$out.rel
  vals_loopM$item <- 1:nrow(vals_loopM)
  ## Now add the mirt values
  vals_loop <- merge(vals_loop, vals_loopM, by=c("item"), suffixes = c("", "_MIRT"))
  
  ## Now do a basic grm model
  mod <- mirt::mirt(data.frame(rep_loop$responses), 1, itemtype = "graded")
  vals_loop$grmRel <-  1 / (1 + (1 / weighted.mean(testinfo(mod, Theta=seq(-5, 5, .1)), dnorm(seq(-5, 5, .1)))))
  ## Now do the same for only the >0 values
  iso.col <- grep(pattern = "Sev", x = colnames(rep_loop$mplusMat))
  mod <- mirt::mirt(data.frame(rep_loop$mplusMat[,iso.col]), 1, itemtype = "graded")
  vals_loop$grmRelExludeZero <-  1 / (1 + (1 / weighted.mean(testinfo(mod, Theta=seq(-5, 5, .1)), dnorm(seq(-5, 5, .1)))))
  vals_all <- dplyr::bind_rows(vals_all, vals_loop)
}

## Now save these
out.list <- list(allParams = vals_all)
saveRDS(out.list, file=out.file)

## Now clean up mplus dir
output.dir <- paste("./data/hurdleCollapse/seedVal_", seedVal, "_min", sep='')
system(paste("rm -rf", output.dir))
