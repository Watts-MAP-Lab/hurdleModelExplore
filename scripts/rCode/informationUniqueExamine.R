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

##### --declare-sim-params-------
## Sim params will need to be modified at a later time point
sim.param.nitems <- c(5,10)
sim.param.ncates <- c(7)
sim.param.discri <- c(1.5, 3)
sim.param.2pl.spread <- c(.5, 1)
sim.param.sample <- c(5000, 20000)
sim.param.faccor <- c(.1,.4,.8)
sim.param.difgrmF <- c(-2,-.5)
sim.param.difgrmC <- c(1,2)
sim.param.nreps <- c(3, 7)
sim.param.discri2 <- c(1.5, 3)
sim.iter <- 1:100
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.2pl.spread,sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.param.nreps, sim.param.discri2,sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCategory", "dicrim", "diffSpread", "n", "facCor", "difGrmF","difGrmC" ,"nCat","grmDiscrim","iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)
seedVal <- as.integer(commandArgs(TRUE))
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
## First create the data -- this will start with the maximum full dataset, 9 total respone categories, full range of difficulty parameters
## This will also showcase where I need to streamline code with custom functions
a = rep(all.sim.vals$dicrim[seedVal], all.sim.vals$nItem[seedVal])
b = genDiffGRM(num_items = all.sim.vals$nItem[seedVal], num_categories = all.sim.vals$nCat[seedVal], min = all.sim.vals$difGrmF[seedVal], max = all.sim.vals$difGrmC[seedVal], rnorm_var = .2)
a_z = rep(all.sim.vals$grmDiscrim[seedVal], all.sim.vals$nItem[seedVal])
## Need to generate 4 separate b_z levels
b_z1 = runif(all.sim.vals$nItem[seedVal], min = -.5, max=-.5+all.sim.vals$diffSpread[seedVal])
b_z2 = runif(all.sim.vals$nItem[seedVal], min = 0, max=0+all.sim.vals$diffSpread[seedVal])
b_z3 = runif(all.sim.vals$nItem[seedVal], min = 1, max=1+all.sim.vals$diffSpread[seedVal])
b_z4 = runif(all.sim.vals$nItem[seedVal], min = 1, max=1+all.sim.vals$diffSpread[seedVal])
b_z5 = runif(all.sim.vals$nItem[seedVal], min = 1.5, max=1.5+all.sim.vals$diffSpread[seedVal])
b_z6 = runif(all.sim.vals$nItem[seedVal], min = 2, max=2+all.sim.vals$diffSpread[seedVal])

muVals = c(0,0)
rho <- all.sim.vals$facCor[seedVal]
varCovMat = matrix(c(1,rho,rho,1), ncol=2)
N = all.sim.vals$n[seedVal]
## Now generate theta here
theta = MASS::mvrnorm(n = N, mu = muVals, Sigma = varCovMat)
reps1 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps2 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z2, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps3 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z3, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps4 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z4, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps5 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z5, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps6 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z6, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps1a = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps2a = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z2, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps3a = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z3, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps4a = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z4, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps5a = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z5, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps6a = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z6, muVals = muVals, varCovMat = varCovMat, theta = theta)

s1 = reps1$mplusMat
s2 = reps2$mplusMat
s3 = reps3$mplusMat
s4 = reps4$mplusMat
s5 = reps5$mplusMat
s6 = reps6$mplusMat

## 2. Create the mplus model formulas
## Create a function which will return the mplus model strings based on the column names
## of the input data matrix
retMplusHurdleInfo <- function(in.mat){
  ## First grab all of the presence variables
  binary.vars <- grep(pattern = "Sus_", x = colnames(in.mat), value = TRUE)
  sev.var <- grep(pattern = "Sev_", x = colnames(in.mat), value = TRUE)
  all.vars <- c(binary.vars, sev.var)
  p1 <- paste0("Presence BY ",paste0(binary.vars, collapse = " \n"))
  p1 <- gsub(pattern = "Sus_1 ", replacement = "Sus_1* ", x = p1)
  p2 <- paste0("Severity BY ",paste0(sev.var, collapse = " \n"))
  p2 <- gsub(pattern = "Sev_1 ", replacement = "Sev_1* ", x = p2)
  p3 <- paste0("Usevariables = ", paste0(all.vars, collapse = "\n"))
  p4 <- paste0("Categorical = ", paste0(all.vars, collapse = "\n"))
  outMplus <- MplusAutomation::mplusObject(TITLE = "testIRTree Model",
                                       MODEL = paste0(p1,"; \n",
                                                      p2, "; \n
                                         Presence@1; Severity@1; \n
                                         [Presence@0]; [Severity@0];"),
                                       rdata = in.mat,
                                       usevariables = all.vars,
                                       VARIABLE = paste0(p3, "; \n",p4,";"),
                                       ANALYSIS = "  Estimator = ML; \n Link = Logit; \n Integration = GAUSSHERMITE(15);")
  return(outMplus)
}
d1 = retMplusHurdleInfo(s1)
d2 = retMplusHurdleInfo(s2)
d3 = retMplusHurdleInfo(s3)
d4 = retMplusHurdleInfo(s4)
d5 = retMplusHurdleInfo(s5)
d6 = retMplusHurdleInfo(s6)

## 3. Estimate all mplus models
## Now create a function which will estimate these models
## Data org prep
modelPrep <- function(in.mat, seedVal, min=TRUE){
  ## First prep the output data directory
  randNum <- sample(1:100000000,size=1)
  output.dir <- paste("./data/hurdleCollapse/seedVal_", seedVal, sep='')
  if(min){
    output.dir <- paste(output.dir, "_min", sep='')
  }
  if(!dir.exists(output.dir)){
    dir.create(output.dir)
  }
  # Grab the max file val here
  max.val <- max(in.mat$rdata, na.rm=TRUE)
  output.file <- paste(output.dir, "/fileVal_", randNum,".inp", sep='')
  output.data <- paste(output.dir, "/fileVal_", randNum, ".txt", sep='') 
  ## Now prep the mplus values
  out <- mplusModeler(in.mat, dataout = output.data, modelout = output.file, run = 1)
  return(out)
}

## Now prep a list of all of these values for mclapply
strat.one3 <- list(d1, d2, d3, d4, d5, d6)
strat.one4 <- parallel::mclapply(strat.one3, FUN = function(x) modelPrep(x, seedVal = seedVal), mc.cores = length(strat.one3))
#strat.one4 <- lapply(strat.one3, FUN = function(x) modelPrep(x, seedVal = seedVal))

## 4. Organize all MPlus output
## Create a function which will grab and organize all mplus output
loadDatFromMplus <- function(seedVal, trueVals){
  ## Load Adon's custom scripts
  source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
  ## First prep the output data directory
  mod <- seedVal$results
  if(identical(mod$parameters, list())){
    all.params <- "modFail"
  }else{
   ## Now obtain estimated a, b, a_z, and b_z parameters
   a_z <- mod$parameters$unstandardized
   a_z <- a_z[grep(pattern = "PRESENCE.BY", x = a_z$paramHeader),]
   b_z <- mod$parameters$unstandardized
   b_z <- b_z[grep(pattern = "Thresholds", x = b_z$paramHeader),]
   b_z$itemVal <- strSplitMatrixReturn(b_z$param, "\\$")[,1]
   b_z <- b_z[which(b_z$itemVal %in% a_z$param),]
   b_z$estStand <- b_z$est / a_z$est
   a <- mod$parameters$unstandardized
   a <- a[grep(pattern = "SEVERITY.BY", x = a$paramHeader),]
   b <- mod$parameters$unstandardized
   b <- b[grep(pattern = "Thresholds", x = b$paramHeader),]
   b$itemVal <- strSplitMatrixReturn(b$param, "\\$")[,1]
   b$threshVal <- strSplitMatrixReturn(b$param, "\\$")[,2]
   b <- b[which(!b$itemVal %in% a_z$param),]
   ## Now translate the b into proper format
   b_org <- pivot_wider(b[,c("est", "itemVal", "threshVal")], id_cols = "itemVal", names_from = "threshVal", values_from = "est")[,-1]
   b_org <- data.frame(b_org)
   ## Now convert these difficulty metrics into the standard normal scale
   b_orgStand <- apply(b_org, 2, function(x) x / a$est)
   rho <- mod$parameters$unstandardized$est[which(mod$parameters$unstandardized$paramHeader=="SEVERITY.WITH")]
   ## Prep the output
   all.params <- cbind(1:length(a$est),rho,a_z$est, b_z$est, a$est, b_orgStand)
   all.paramsTrue <- cbind(trueVals$varCovMat[1,2], trueVals$a_z, trueVals$b_z,
                           trueVals$a, trueVals$b)
   prepVec <- c("item","rho","est_z_discrim","est_z_diff","est_grm_discrim")
   grmDiffVec <- (dim(all.params)[2] - length(prepVec))
   grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
   prepVec1 <- c(prepVec, grmEstVec)
   ## Now do the true names
   prepVec <- c("rhoTrue","true_z_discrim","true_z_diff","true_grm_discrim")
   grmDiffVec <- (dim(all.paramsTrue)[2] - length(prepVec))
   grmEstVec <- paste("true_grm_diff_", 1:grmDiffVec, sep='')
   prepVec2 <- c(prepVec, grmEstVec)
   all.params <- cbind(all.params, all.paramsTrue)
   colnames(all.params) <- c(prepVec1, prepVec2)
   all.params <- data.frame(all.params)
   all.params
  }
  return(all.params)
}

vals1 <- loadDatFromMplus(strat.one4[[1]], reps1)
vals2 <- loadDatFromMplus(strat.one4[[2]], reps2)
vals3 <- loadDatFromMplus(strat.one4[[3]], reps3)
vals4 <- loadDatFromMplus(strat.one4[[4]], reps4)
vals5 <- loadDatFromMplus(strat.one4[[5]], reps5)
vals6 <- loadDatFromMplus(strat.one4[[6]], reps6)

## Now go through each of these and add the estimated test reliability from alpha and omega
vals_all <- NULL
for(i in 1:6){
  ## First get the values
  vals_loop <- get(paste("vals", i, sep=''))
  rep_loop <- get(paste("reps", i, sep=''))
  ## Now estimate alpha & omega
  test_dat <- as.data.frame(rep_loop$responses)
  rel.alpha <- psych::alpha(test_dat)
  rel.ome <- psych::omega(test_dat, poly=TRUE, nfactors = 1)
  rel.ome2 <- psych::omega(test_dat, poly=TRUE, nfactors = 2)
  vals_loop$omega_h <- rel.ome2$omega_h
  vals_loop$alpha <- rel.alpha$total$raw_alpha
  vals_loop$omega_t <- rel.ome$omega.tot
  vals_loop$G_six <- rel.ome$G6
  vals_loop$alpheFromOme <- rel.ome$alpha
  ##  Now grab the true reliability values
  a = vals_loop$true_z_discrim
  b = vals_loop[,grep(pattern = "true_grm_diff", x = names(vals_loop))]
  a_z = vals_loop$true_z_discrim
  b_z = vals_loop$true_z_diff
  vals_loop$trueHurdleRel <- estHurdleRel(rep_loop,a = a, b = b,a_z = a_z, b_z = b_z)
  ## Do the same for estimated here -- need to protect against NA values in diff parameters
  a = vals_loop$est_z_discrim
  b = vals_loop[,grep(pattern = "est_grm_diff", x = names(vals_loop))]
  if(sum(is.na(b))>0){
    b[is.na(b)] <- 3
  }
  b = t(apply(b, 1, sort))
  a_z = vals_loop$est_z_discrim
  b_z = vals_loop$est_z_diff
  vals_loop$estHurdleRel <- estHurdleRel(rep_loop,a = a, b = b,a_z = a_z, b_z = b_z)
  vals_all <- dplyr::bind_rows(vals_all, vals_loop)
}

## Now do test re test?
all_trt <- list()
for(i in 1:6){
  rep_loop <- get(paste("reps", i, sep=''))
  rep_loopa <- get(paste("reps", i, "a",sep=''))
  ## Now estimate alpha & omega
  test_dat <- as.data.frame(rep_loop$responses)
  ## add ID and time
  test_dat$id <- 1:nrow(test_dat)
  test_dat$time <- 1
  test_data <- as.data.frame(rep_loopa$responses)
  ## add id & time
  test_data$id <- 1:nrow(test_dat)
  test_data$time <- 2
  for_trt <- bind_rows(test_dat, test_data)
  ## retest
  trtVal <- psych::testRetest(t1 = for_trt, keys = names(test_dat)[1:5],id="id", time="time",lmer=FALSE)
  all_trt[[i]] <- trtVal
}

## Now save these
out.list <- list(allParams = vals_all,
                 all_trt)
saveRDS(out.list, file=out.file)