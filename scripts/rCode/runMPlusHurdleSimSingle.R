## First load packages
library(MplusAutomation)
library(foreach)
library(doParallel)
library(tidyverse)
source("./scripts/rCode/hurdleFunctions.r")
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
## Now run through a test case
## first sim data
## This script will read a single input which will be the row for the simulation to run
## First generate all of the simulation permutations
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


orig.dir <- getwd()

#setup parallel backend to use many processors
## Now split this across parallel loops
seedVal <- as.integer(commandArgs(TRUE))
out.directory <- paste(orig.dir, '/data/simdir_',seedVal ,sep='')
dir.create(out.directory)
setwd(out.directory)
theta1 <- seq(-8,8,0.2) # Quadrature points for theta1
theta2 <- seq(-8,8,0.2) # Quadrature points for theta1
theta <- expand.grid(theta1, theta2)
names(theta) <- c("theta1", "theta2")
qpoints <- theta
for(i in 1:nrow(all.sim.vals)){
  set.seed(seedVal)
  ## Check if files exist
  out.file <- paste('fileVal_', i,'.inp', sep='')
  out.file2 <- paste('fileVal_', i,'.out', sep='')
  out.file3 <- paste('fileVal_', i,'.RDS', sep='')
  dat.file <- paste('fileVal_', i,'.txt', sep='')
  if(file.exists(out.file3) & file.exists(dat.file)){
    print(paste(out.file3, " exists"))
  }else{
    ## Now generate the data
    row.1 <- c(1,all.sim.vals[i,6], all.sim.vals[i,7])
    row.2 <- c(all.sim.vals[i,6], 1, all.sim.vals[i,8])
    row.3 <- c(all.sim.vals[i,7], all.sim.vals[i,8], 1)
    varCovMat <- matrix(c(row.1, row.2, row.3), nrow = 3, byrow=T)
    set.seed(all.sim.vals[i,10])
    ## Now simulate the data
    out.data1 <- simulate_hurdle_responses(a = rep(all.sim.vals[i,3], all.sim.vals[i,1]),b_z = seq(all.sim.vals[i,4], 1, length.out=all.sim.vals[i,1]),a_z = rep(all.sim.vals[i,3], all.sim.vals[i,1]),
                                           b = genDiffGRM(num_items = all.sim.vals[i,1], num_categories = all.sim.vals[i,2], min = all.sim.vals[i,9], max = 1), muVals = c(0,0,0), varCovMat, all.sim.vals[i,5])
    for.mplus.dat <- out.data1$responses
    ## Now make the matrix for these data
    matrix.mplus1 <- matrix(NA, nrow = nrow(out.data1$responses), ncol=ncol(out.data1$responses))
    matrix.mplus2 <- matrix(NA, nrow = nrow(out.data1$responses), ncol=ncol(out.data1$responses))
    
    ## Now make the binary indicators first
    matrix.mplus1[,1:ncol(out.data1$responses)] <- out.data1$responses
    matrix.mplus1[which(matrix.mplus1[,1:ncol(out.data1$responses)]>0)] <- 1
    ## Now do severity indicators here
    matrix.mplus2 <- out.data1$responses
    matrix.mplus2[which(matrix.mplus2==0)] <- NA
    
    ## Now combine these
    matrix.mplus <- cbind(matrix.mplus1, matrix.mplus2)
    colnames(matrix.mplus) <- paste("X", 1:ncol(matrix.mplus), sep='')
    first.dim <- ncol(matrix.mplus1)
    second.dim1 <- first.dim+1
    second.dim2 <- first.dim*2
    binary.vars <- colnames(matrix.mplus)[1:first.dim]
    sev.var <- colnames(matrix.mplus)[second.dim1:second.dim2]
    all.vars <- colnames(matrix.mplus)
    p1 <- paste0("Presence BY ",paste0(binary.vars, collapse = " "))
    p1 <- gsub(pattern = "X1 ", replacement = "X1* ", x = p1)
    p2 <- paste0("Severity BY ",paste0(sev.var, collapse = " "))
    p2 <- gsub(pattern = sev.var[1], replacement = paste(sev.var[1], "*", sep=''), x = p2)
    p3 <- paste0("Usevariables = ", paste0(all.vars, collapse = " "))
    p4 <- paste0("Categorical = ", paste0(all.vars, collapse = " "))
    matrix.mplus <- data.frame(matrix.mplus)
    test <- MplusAutomation::mplusObject(TITLE = "testIRTree Model",
                                         MODEL = paste0(p1,"; \n",
                                                        p2, "; \n
                                         Presence@1; Severity@1; \n
                                         [Presence@0]; [Severity@0];"),
                                         rdata = matrix.mplus,
                                         usevariables = all.vars,
                                         VARIABLE = paste0(p3, "; \n",p4,";"),
                                         ANALYSIS = "  Estimator = ML; \n Link = Logit; \n Integration = GAUSSHERMITE(15);")
    mplusModeler(test, dataout=dat.file, modelout = out.file)
    runModels(target = out.file)
    mod <- readModels(target = out.file2)
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
    ## Attach these to data out
    out.data1$a <- cbind(out.data1$a, a$est)
    out.data1$b <- cbind(out.data1$b, b_orgStand)
    out.data1$a_z <- cbind(out.data1$a_z, a_z$est)
    out.data1$b_z <- cbind(out.data1$b_z, b_z$estStand)
    out.data1$estRho <- rho
    ## Now add the theta estimates
    prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
    itemtrace <- trace.line.pts(a$est, b_orgStand, a_z$est, b_z$estStand, theta)
    out.data1$theta$eap2PL_Hurdle <- NA
    out.data1$theta$se2PL_Hurdle <- NA
    out.data1$theta$eapGRM_Hurdle <- NA
    out.data1$theta$seGRM_Hurdle <- NA
    for(respondent in 1:nrow(out.data1$responses)) {
      pattern <- out.data1$responses[respondent,]
      lhood <- score(pattern, itemtrace, qpoints)
      
      eap2PL_Hurdle <- sum(lhood*prior*qpoints$theta1)/sum(lhood*prior)
      se2PL_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta1-eap2PL_Hurdle)^2)/sum(lhood*prior))
      out.data1$theta$eap2PL_Hurdle[respondent] <- eap2PL_Hurdle
      out.data1$theta$se2PL_Hurdle[respondent] <- se2PL_Hurdle
      
      eapGRM_Hurdle <- sum(lhood*prior*qpoints$theta2)/sum(lhood*prior)
      seGRM_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta2-eapGRM_Hurdle)^2)/sum(lhood*prior))
      out.data1$theta$eapGRM_Hurdle[respondent] <- eapGRM_Hurdle
      out.data1$theta$seGRM_Hurdle[respondent] <- seGRM_Hurdle
    }
    saveRDS(out.data1, file = out.file3)
    ## Now rm all named objects
    rm(row.1, row.2, row.3, varCovMat, out.data1, for.mplus.dat, matrix.mplus1, matrix.mplus2, matrix.mplus, first.dim, second.dim1, second.dim2, binary.vars, sev.var, all.vars, p1, p2, p3, p4, test, a_z, b_z, a, b, b_org, b_orgStand, out.data1, prior, itemtrace, pattern, lhood, eap2PL_Hurdle, se2PL_Hurdle, eapGRM_Hurdle, seGRM_Hurdle)
    gc()
  }
}
setwd(orig.dir)
print(paste("Done: ", seedVal))
