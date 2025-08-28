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
library("mgcv")

##### --declare-sim-params-------
## Sim params will need to be modified at a later time point
source("./scripts/rCode/simParam.r")
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
add.val.2pl <- 1.5
add.val.grm <- 1.5
a = runif(n = all.sim.vals$nItem[seedVal], min = all.sim.vals$grmDiscrim[seedVal], all.sim.vals$grmDiscrim[seedVal] + add.val.grm)
b = genDiffGRM(num_items = all.sim.vals$nItem[seedVal], num_categories = all.sim.vals$nCat[seedVal], min = all.sim.vals$difGrmF[seedVal], max = all.sim.vals$difGrmF[seedVal]+2.5, rnorm_var = .3)
a_z = runif(n = all.sim.vals$nItem[seedVal], min = all.sim.vals$discrim2pl[seedVal], all.sim.vals$discrim2pl[seedVal] + add.val.2pl)
## Need to generate 4 separate b_z levels
b_z1 = runif(all.sim.vals$nItem[seedVal], min = all.sim.vals$dif2PL[seedVal], max=all.sim.vals$dif2PL[seedVal]+2)
muVals = c(0,0)
rho <- all.sim.vals$facCor[seedVal]
varCovMat = matrix(c(1,rho,rho,1), ncol=2)
N = all.sim.vals$n[seedVal]
## Now generate theta here
#theta = MASS::mvrnorm(n = N, mu = muVals, Sigma = varCovMat)
theta = MASS::mvrnorm(n = N, mu = muVals, Sigma = varCovMat)
#reps1 = system.time(simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta))
#reps1f = system.time(simulate_hurdle_responses_fast(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta, qpoints = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1))))
reps1 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta)
#reps1 = simulate_hurdle_responses_fast(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta, qpoints = expand.grid(seq(-6, 6, .2), seq(-6, 6, .2)))

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
  #source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
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
  # rel.alpha <- psych::alpha(polychoric(test_dat)$rho)
  # rel.ome <- psych::omega(test_dat, poly=TRUE, nfactors = 3, plot = FALSE)
  # rel.uni <- psych::unidim(test_dat, cor="poly")
  cor.mat <- psych::polychoric(rep_loop$responses)
  cor.mat1 <- cor.mat$rho
  rel.all <- psych::reliability(cor.mat$rho, n.obs = 15000, nfactors=3,plot=FALSE)
  vals_loop$omega_h <- rel.all$result.df[,"omega_h"]
  vals_loop$alpha <- rel.all$result.df[,"alpha"]
  vals_loop$omega_t <- rel.all$result.df[,"omega.tot"]
  vals_loop$Uni <- rel.all$result.df[,"Uni"]
  vals_loop$tau <- rel.all$result.df[,"tau"]
  vals_loop$cong <- rel.all$result.df[,"cong"]
  vals_loop$CFI <- rel.all$result.df[,"CFI"]
  vals_loop$ECV <- rel.all$result.df[,"ECV"]
  vals_loop$Beta <- rel.all$result.df[,"Beta"]
  vals_loop$EVR <- rel.all$result.df[,"EVR"]
  vals_loop$MAP <- rel.all$result.df[,"MAP"]
  vals_loop$skewVal <- as.numeric(psych::describe(rowSums(rep_loop$responses))["skew"])
  split.half.rel <- suppressWarnings(psych::guttman(r = cor.mat$rho))
  vals_loop$lambda1Rel = split.half.rel$lambda.1
  vals_loop$lambda2Rel = split.half.rel$lambda.2
  vals_loop$lambda3Rel = split.half.rel$lambda.3
  vals_loop$lambda4Rel = split.half.rel$lambda.4
  vals_loop$lambda5Rel = split.half.rel$lambda.5
  vals_loop$lambda6Rel = split.half.rel$lambda.6
  omega.sem.val <- omegaSem(m = cor.mat$rho, nfactors=1, plot = FALSE, n.obs = 15000)
  vals_loop$singleFactorOmegaT <- omega.sem.val$omegaSem$omega.tot

  ##  Now grab the true reliability values
  a = vals_loop$true_grm_discrim
  b = data.matrix(data.frame(vals_loop[,grep(pattern = "true_grm_diff", x = names(vals_loop))]))
  a_z = vals_loop$true_z_discrim
  b_z = vals_loop$true_z_diff
  vals_loop$trueHurdleRel <- hurdInfo(theta.grid = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)), a = a, b = b, a_z = a_z, b_z = b_z, muVals = muVals, rhoVal = rho)$out.rel
  ## NOw obtain the true reliability using the observed variances
  true.score.var <- var(reps1$theta$theta2)
  error.var <- var(reps1$theta$eapSev - reps1$theta$theta2)
  vals_loop$trueRel <- true.score.var / (true.score.var + error.var)
  ## Now try to grab the true rel from these data using the GAM method?
  tmp.dat <- reps1$theta
  rel.sev <- mgcv::gam(eapSev ~ s(theta1, k = 20) + s(theta2, k = 20), data = tmp.dat, REML=TRUE)
  rel.sus <- mgcv::gam(eapSus ~ s(theta1) + s(theta2), data = tmp.dat)
  prms.sus <- mgcv::gam(eapSus ~ s(theta1), data = tmp.dat)
  prms.sev <- mgcv::gam(eapSev ~ s(theta2), data = tmp.dat)
  vals_loop$trueRelSevGam <- summary(rel.sev)$r.sq
  vals_loop$trueRelSusGam <- summary(rel.sus)$r.sq
  vals_loop$truePrmseSevGam <- summary(prms.sev)$r.sq
  vals_loop$truePrmseSusGam <- summary(prms.sus)$r.sq
  
  ## Now compare hurd info obtained from MIRT
  ## Estimate test info using the expand grid with factor correlations
  ## Isolate those that sum to 90
  t2 <- mirt::testinfo(sv1, Theta = expand.grid(seq(-5, 5, .2), seq(-5, 5, .2)), degrees = c(0,90))
  vals_loop$mirtHurdleInfo0_90 <-  1 / (1 + (1 / weighted.mean(testinfo(sv1, Theta=expand.grid(seq(-5, 5, .1), seq(-5, 5, .1)), degrees = c(0,90)), dmnorm(mu = muVals, sigma = varCovMat, x = expand.grid(seq(-5, 5, .1), seq(-5, 5, .1))))))
  vals_loop$mirtHurdleInfo45_45 <-  1 / (1 + (1 / weighted.mean(testinfo(sv1, Theta=expand.grid(seq(-5, 5, .1), seq(-5, 5, .1)), degrees = c(45,45)), dmnorm(mu = muVals, sigma = varCovMat, x = expand.grid(seq(-5, 5, .1), seq(-5, 5, .1))))))
  vals_loop$mirtHurdleInfo90_0 <-  1 / (1 + (1 / weighted.mean(testinfo(sv1, Theta=expand.grid(seq(-5, 5, .1), seq(-5, 5, .1)), degrees = c(90,0)), dmnorm(mu = muVals, sigma = varCovMat, x = expand.grid(seq(-5, 5, .1), seq(-5, 5, .1))))))
  
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
  vals_loopM <- bind_cols(rhoEst,a_z, b_z,a, b)
  colnames(vals_loopM)[1:4] <- c("rhoEst","est_z_discrim", "est_z_diff", "est_grm_discrim")
  colnames(vals_loopM)[-c(1:4)] <- paste("est_grm_diff_", 1:dim(vals_loopM[,-c(1:4)])[2], sep='')
  vals_loopM$estHurdleRel <- hurdInfo(theta.grid = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)), a = a, b = b, a_z = a_z, b_z = b_z, muVals = muVals, rhoVal = rhoEst)$out.rel
  vals_loopM$item <- 1:nrow(vals_loopM)
  ## Now add the mirt values
  vals_loop <- merge(vals_loop, vals_loopM, by=c("item"), suffixes = c("", "_MIRT"))
  ## Now obtain the Reliability as assessed via GAM & simulation
  # simVals <- simulate_hurdle_responses(a, b, a_z, b_z, muVals = muVals, varCovMat = matrix(c(1,rhoEst,rhoEst,1), ncol=2), N = 100000)
  # 
  # ## Now perform the gam regressions
  # tmp.dat <- simVals$theta
  # rel.sev <- mgcv::gam(eapSev ~ s(theta1) + s(theta2), data = tmp.dat)
  # rel.sus <- mgcv::gam(eapSus ~ s(theta1) + s(theta2), data = tmp.dat)
  # prms.sus <- mgcv::gam(eapSus ~ s(theta1), data = tmp.dat)
  # prms.sev <- mgcv::gam(eapSev ~ s(theta2), data = tmp.dat)
  # vals_loop$estRelSevGam <- summary(rel.sev)$r.sq
  # vals_loop$estRelSusGam <- summary(rel.sus)$r.sq
  # vals_loop$estPrmseSevGam <- summary(prms.sev)$r.sq
  # vals_loop$estPrmseSusGam <- summary(prms.sus)$r.sq
  
  ## Now do a basic grm model
  mod <- mirt::mirt(data.frame(rep_loop$responses), 1, itemtype = "graded")
  vals_loop$grmRel <-  1 / (1 + (1 / weighted.mean(testinfo(mod, Theta=seq(-6, 6, .1)), dnorm(seq(-6, 6, .1)))))
  # mod <- mirt::mirt(data.frame(rep_loop$responses[-rm.index,]), 1, itemtype = "graded")
  # vals_loop$grmRel_NoZ <-  1 / (1 + (1 / weighted.mean(testinfo(mod, Theta=seq(-5, 5, .1)), dnorm(seq(-5, 5, .1)))))
  ## Now do the same for only the >0 values
  iso.col <- grep(pattern = "Sev", x = colnames(rep_loop$mplusMat))
  mod.rm <- mirt::mirt(data.frame(rep_loop$mplusMat[,iso.col]), 1, itemtype = "graded")
  vals_loop$grmRel_rmZeroOption <-  1 / (1 + (1 / weighted.mean(testinfo(mod.rm, Theta=seq(-6, 6, .1)), dnorm(seq(-6, 6, .1)))))
  vals_all <- dplyr::bind_rows(vals_all, vals_loop)
}

## Now prep the response patterns and theta values
#reps1$theta$responsePattern <- apply(reps1$responses, 1, function(x) paste(x, collapse = "-"))
reps1$theta$rowSums <- rowSums(reps1$responses)
reps1$theta$grmFacScore <- fscores(mod)
reps1$theta$grmFacScoreRM <- fscores(mod.rm)
reps1$theta$factorScores <- psych::fa.poly(reps1$responses, 1)$scores$scores[,1]
reps1$theta$estSus <- fscores(sv1)[,1]
reps1$theta$estSev <- fscores(sv1)[,2]
colnames(reps1$theta) <- c("trueSus", "trueSev", "eapTrueSusParam", "eapTrueSevParam","rowSum", "estGRM", "estGRMrm", "est1Fac", "estSus", "estSev")

## Now save these
out.list <- list(allParams = vals_all, thetaVals = reps1$theta)
saveRDS(out.list, file=out.file)

## Now clean up mplus dir
output.dir <- paste("./data/hurdleCollapse/seedVal_", seedVal, "_min", sep='')
system(paste("rm -rf", output.dir))
