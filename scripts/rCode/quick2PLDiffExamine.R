##### --load-library-------------
source("./scripts/rCode/hurdleFunctions.r")
suppressPackageStartupMessages(library("tidyverse"))
library("numDeriv")
suppressPackageStartupMessages(library(MplusAutomation))

##### --declare-sim-params-------
## Sim params will need to be modified at a later time point
sim.param.nitems <- c(5,8,11)
sim.param.ncates <- c(7)
sim.param.discri <- c(2)
sim.param.difficF <- c(-1,0)
sim.param.difficC <- c(1,2)
sim.param.sample <- c(10000)
sim.param.faccor <- c(.1,.4,.8)
sim.param.difgrmF <- c(-2,-.5)
sim.param.difgrmC <- c(1,2)
sim.iter <- 1:100
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.difficF, sim.param.difficC,sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCategory", "dicrim", "diff2PLF","diff2PLC", "n", "facCor", "difGrmF","difGrmC" ,"iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)
#seedVal <- as.integer(commandArgs(TRUE))
seedVal <- 1
set.seed(seedVal)


a = rep(all.sim.vals$dicrim[seedVal], all.sim.vals$nItem[seedVal])
b = genDiffGRM(num_items = all.sim.vals$nItem[seedVal], num_categories = 3, min = all.sim.vals$difGrmF[seedVal], max = all.sim.vals$difGrmC[seedVal], rnorm_var = .2)
a_z = rep(all.sim.vals$dicrim[seedVal], all.sim.vals$nItem[seedVal])
b_z = runif(all.sim.vals$nItem[seedVal], min = all.sim.vals$diff2PLF[seedVal], max=all.sim.vals$diff2PLC[seedVal])
muVals = c(0,0)
rho <- all.sim.vals$facCor[seedVal]
varCovMat = matrix(c(1,rho,rho,1), ncol=2)
N = all.sim.vals$n[seedVal]
reps = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z, muVals = muVals, varCovMat = varCovMat, N = N)


## I need to simulate data while modifying the difficulty parameters of the 2PL portion of the model
## so I need to make a group of b_z parameters
b_z1 = runif(all.sim.vals$nItem[seedVal], min = -2, max=0)
b_z2 = runif(all.sim.vals$nItem[seedVal], min = -1, max=1)
b_z3 = runif(all.sim.vals$nItem[seedVal], min = 0, max=1)
b_z4 = runif(all.sim.vals$nItem[seedVal], min = 1, max=2)
b_z5 = runif(all.sim.vals$nItem[seedVal], min = 2, max=3)
b_z6 = runif(all.sim.vals$nItem[seedVal], min = 3, max=4)


## NOw simulate a massive data set
theta = MASS::mvrnorm(n = 10000, Sigma = varCovMat, mu=c(0,0))
reps1 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps2 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z2, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps3 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z3, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps4 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z4, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps5 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z5, muVals = muVals, varCovMat = varCovMat, theta = theta)
reps6 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z6, muVals = muVals, varCovMat = varCovMat, theta = theta)


## Now run thse through the mplus process
createMplusHurdleMatrix <- function(in.dat){
  ## This function will be used to create the pseudo items
  ## for the hurdle model, specifically a susceptibility & severity
  ## scoring framework
  ## Now make the binary indicators first
  matrix.mplus1<- in.dat
  matrix.mplus1[which(matrix.mplus1[,1:ncol(in.dat)]>0)] <- 1
  ## Now do severity indicators here
  matrix.mplus2 <- in.dat
  matrix.mplus2[which(matrix.mplus2==0)] <- NA
  matrix.mplus <- cbind(matrix.mplus1, matrix.mplus2)
  ## Now make the column names
  col.names.sus <- paste("Sus_", 1:ncol(matrix.mplus1), sep='')
  col.names.sev <- paste("Sev_", 1:ncol(matrix.mplus2), sep='')
  matrix.mplus <- data.frame(matrix.mplus)
  colnames(matrix.mplus) <- c(col.names.sus, col.names.sev)
  return(matrix.mplus)
}
s1 = createMplusHurdleMatrix(reps1$responses)
s2 = createMplusHurdleMatrix(reps2$responses)
s3 = createMplusHurdleMatrix(reps3$responses)
s4 = createMplusHurdleMatrix(reps4$responses)
s5 = createMplusHurdleMatrix(reps5$responses)
s6 = createMplusHurdleMatrix(reps6$responses)


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


modelPrep <- function(in.mat, seedVal, min=TRUE){
  ## First prep the output data directory
  output.dir <- paste("./data/hurdleCollapse/seedVal_", seedVal, sep='')
  if(min){
    output.dir <- paste(output.dir, "_min", sep='')
  }
  if(!dir.exists(output.dir)){
    dir.create(output.dir)
  }
  # Grab the max file val here
  max.val <- max(in.mat$rdata, na.rm=TRUE)
  output.file <- paste(output.dir, "/fileVal_", max.val, ".inp", sep='')
  output.data <- paste(output.dir, "/fileVal_", max.val, ".txt", sep='') 
  ## Now prep the mplus values
  out <- mplusModeler(in.mat, dataout = output.data, modelout = output.file, run = 1)
  return(out)
}

z1 = modelPrep(d1, 99999999, min=FALSE)
z2 = modelPrep(d2, 99999998, min=FALSE)
z3 = modelPrep(d3, 99999997, min=FALSE)
z4 = modelPrep(d4, 99999996, min=FALSE)
z5 = modelPrep(d5, 99999995, min=FALSE)
z6 = modelPrep(d5, 99999994, min=FALSE)

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

j1 <- loadDatFromMplus(z1, reps1)
j2 <- loadDatFromMplus(z2, reps2)
j3 <- loadDatFromMplus(z3, reps3)
j4 <- loadDatFromMplus(z4, reps4)
j5 <- loadDatFromMplus(z5, reps5)
j6 <- loadDatFromMplus(z6, reps6)

## Now plot these results
hist(rowSums(reps1$responses))
hist(rowSums(reps2$responses))
hist(rowSums(reps3$responses))
hist(rowSums(reps4$responses))
hist(rowSums(reps5$responses))
hist(rowSums(reps6$responses))
plot(data.frame(diffID = 1:6, rhoEst = c(j1$rho[1],j2$rho[1],j3$rho[1],j4$rho[1],j5$rho[1],j6$rho[1])))