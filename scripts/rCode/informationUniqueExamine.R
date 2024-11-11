## This script will be used to run the winosrization hurdle model explore
## It is going to require simulation a zero-inflated test response pattern
## and changing response patterns by collapsing response categories into one another
## This will be performed using three separate strategies:
##  1. Winsorization: Taking the highest response category and collapsing it into the next highest
##  2. Rec Winsor: Collapsing the lowest resposne category into the 0 response
##  3. Total information preservation: Identify the difficulty parameter that by removing would preserev the total information; total information will be taken as the integral of the test informaiton function

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
seedVal <- as.integer(commandArgs(TRUE))
set.seed(seedVal)
## --run-test-sim ------------
## I will run through a single test sim here and see how I can best, and most efficiently run this
## I really need to focus on the memory leakage problem the R has, especially when running large loops
## First create the data -- this will start with the maximum full dataset, 9 total respone categories, full range of difficulty parameters
## This will also showcase where I need to streamline code with custom functions
a = rep(all.sim.vals$dicrim[seedVal], all.sim.vals$nItem[seedVal])
b = genDiffGRM(num_items = all.sim.vals$nItem[seedVal], num_categories = 9, min = all.sim.vals$difGrmF[seedVal], max = all.sim.vals$difGrmC[seedVal], rnorm_var = .2)
a_z = rep(all.sim.vals$dicrim[seedVal], all.sim.vals$nItem[seedVal])
b_z = runif(all.sim.vals$nItem[seedVal], min = all.sim.vals$diff2PLF[seedVal], max=all.sim.vals$diff2PLC[seedVal])
muVals = c(0,0)
rho <- all.sim.vals$facCor[seedVal]
varCovMat = matrix(c(1,rho,rho,1), ncol=2)
N = all.sim.vals$n[seedVal]
reps = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z, muVals = muVals, varCovMat = varCovMat, N = N)
## Now prep all of the winsor examines
retWinsorMat <- function(roof, mat){
  mat[mat>roof] <- roof
  return(mat)
}
retRevWinsorMat <- function(floor, mat){
  mat[mat<floor] <- 0
  ## Now change the total number of response categories
  non.zero.index <- which(mat>0)
  min.val <- min(mat[non.zero.index])
  diff.val <- min.val-1
  mat[non.zero.index] <- mat[non.zero.index]-diff.val 
  return(mat)
}

## This function will require some serious work to produce
## The problem is that item response endorsements aren't that easy to calculate within the IRTree framewrok?
## So while the end goal should be the sum of the derivatives for every response pattern
## This will take some work to actually produce
maximumInfoStrat <- function(numReps, mat, a, b, a_z, b_z){
  ## First get the prob of endorsements at theta levels
  ## I need to go ahead and make a list of functions 
  ## Which I can then take the derivative of to grab information at specific points
  ## SO I will loop through each of the response categories and create a function
  ## For probability of endorsment at all theta levels
  list.of.funcs <- list()
  num.categories <- dim(b)[2]+1
  num.items <- length(a)
  theta.vals1 <- seq(length)
  all.information <- matrix(ncol=2, )
  for(l in 1:nu.items){
    ## First create the probability of any endorsement
    fac.one <- prob.endorse.ret(a_z[l], b_z[l])
    ## Now go through the response categories and grab the function and estimate the derivative
    for(k in 1:num.categories){
      if(k==1){
        fac.two <- prob.endorse.grm.ret(a[k], b[k], min=TRUE)
      }else if(k>1 & k < num.categories){
        k1 <- k-1
        fac.two <- prob.endorse.grm.ret(a[k], b[k], b1 = b[k1])
      }else{
        fac.two <- prob.endorse.grm.ret(a[k],b[k],max=TRUE)
      }
      fac.thr <- function(x){fac.one(x) * fac.two(x)}
      ## Now grab the derivative of these values from -8:8
      theta.vals1 <- expand.grid(seq(-8,8,.2), 0)
      # grab the derivative
      
    }
  }
}

## Create a function which will return a function with the a * b values for endorsement
prob.endorse.ret <- function(a, b){
  out.func <- function(x) exp(a*(x-b))/(1+exp(a*(x-b)))
  return(out.func)
}
prob.endorse.grm.ret <- function(a,b, min=FALSE, max=FALSE, b1=NULL){
  if(min){
    out.func <- function(x){1-exp(a*(x-b))/(1+exp(a*(x-b)))}
  }
  if(max){
    out.func <- function(x){exp(a*(x - b))/(1+exp(a*(x - b)))}
  }
  else if (!min & !max){
    out.func <- function(x){exp(a*(x - b1))/(1+exp(a*(x - b1))) - exp(a*(x - b))/(1+exp(a*(x - b)))}
  }
  return(out.func)
}

strat.one <- lapply(3:7, function(x) retWinsorMat(x, reps$responses))
strat.two <- lapply(3:7, function(x) retRevWinsorMat(x, reps$responses))
## Now estimate the mplus models across all of these models
## This means I am going to have to create another function which will perform
## all of the prep for the MPlus models, this includes:
## 1. Create the data matrices
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
strat.one2 <- lapply(strat.one, createMplusHurdleMatrix)
strat.two2 <- lapply(strat.two, createMplusHurdleMatrix)

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
strat.one3 <- lapply(strat.one2, retMplusHurdleInfo)
strat.two3 <- lapply(strat.two2, retMplusHurdleInfo)

## 3. Estimate all mplus models
## Now create a function which will estimate these models
## Data org prep
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
strat.one4 <- parallel::mclapply(strat.one3, FUN = function(x) modelPrep(x, seedVal = seedVal), mc.cores = length(strat.one3))
#strat.one4 <- lapply(strat.one3, FUN = function(x) modelPrep(x, seedVal = seedVal))
strat.two4 <- parallel::mclapply(strat.two3, FUN = function(x) modelPrep(x, seedVal = seedVal, min=FALSE), mc.cores = length(strat.two3))
#strat.two4 <- lapply(strat.two3, FUN = function(x) modelPrep(x, seedVal = seedVal, min=FALSE))

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

vals1 <- lapply(strat.one4, function(x) loadDatFromMplus(x, reps))
vals2 <- lapply(strat.two4, function(x) loadDatFromMplus(x, reps))

## Now save these
out.list <- list(vals1, vals2)
out.file = paste("./data/hurdleCollapse/seedVal_", seedVal, ".RDS", sep='')
saveRDS(out.list, file=out.file)

