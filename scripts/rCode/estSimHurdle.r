## This script will be used to estimate the hurdle models on the ACCRE cluster
## This script will try to reduce the total amount of disk usage requried
## It will perfrom 4 discrete tasks:
## 1. Simulate the data required to esitmate the models
## 2. Esimtate the models by calling the julia script
## 3. Return true and esitmated parameters
## 4. Also return the esitmated theta parameters and their true counter parts

## Source and load libraries
library(tidyverse)
source("./scripts/rCode/hurdleFunctions.r")

## Get the command calls
i <- commandArgs(TRUE)

## This script will read a single input which will be the row for the simulation to run
## FIrst generate all of the simulation permutations
sim.param.nitems <- c(5,9)
sim.param.ncates <- c(3,5,7)
sim.param.discri <- c(3)
sim.param.diffic <- c(-1,1)
sim.param.sample <- c(10000)
sim.param.faccor <- c(.2, .4, .7)
sim.param.critco <- c(.2, .4, .7)
sim.param.difgrm <- c(-1, 0)
covarPat <- c("a","b","c")
sim.perms <- 1
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.diffic, sim.param.sample, sim.param.faccor, 
                            sim.param.critco, sim.param.difgrm, covarPat,sim.perms)

## Grab a random integer for the data files
random.int <- sample(1:2^30, 1)

## Now generate the data
row.1 <- c(1,all.sim.vals[i,6], 0)
row.2 <- c(all.sim.vals[i,6], 1, all.sim.vals[i,7])
row.3 <- c(0, all.sim.vals[i,7], 1)
varCovMat1 <- matrix(c(row.1, row.2, row.3), nrow = 3, byrow=T)
row.1 <- c(1,all.sim.vals[i,6], all.sim.vals[i,7])
row.2 <- c(all.sim.vals[i,6], 1, 0)
row.3 <- c(all.sim.vals[i,7], 0,1)
varCovMat2 <- matrix(c(row.1, row.2, row.3), nrow = 3, byrow=T)
row.1 <- c(1,all.sim.vals[i,6], all.sim.vals[i,7])
row.2 <- c(all.sim.vals[i,6], 1, all.sim.vals[i,7])
row.3 <- c(all.sim.vals[i,7], all.sim.vals[i,7], 1)
varCovMat3 <- matrix(c(row.1, row.2, row.3), nrow = 3, byrow=T)    
set.seed(all.sim.vals[i,10])
if(all.sim.vals[i,9]=="a"){varCovMat <- varCovMat1}
if(all.sim.vals[i,9]=="b"){varCovMat <- varCovMat2}
if(all.sim.vals[i,9]=="c"){varCovMat <- varCovMat3}

## Now simulate the data
out.data1 <- simulate_hurdle_responses(a = rep(all.sim.vals[i,3], all.sim.vals[i,1]),b_z = seq(all.sim.vals[i,4], 2, length.out=all.sim.vals[i,1]),a_z = rep(all.sim.vals[i,3], all.sim.vals[i,1]), 
    b = genDiffGRM(num_items = all.sim.vals[i,1], num_categories = all.sim.vals[i,2], min = all.sim.vals[i,8]), muVals = c(0,0,0), varCovMat, all.sim.vals[i,5])
#out.data2 <- simulate_hurdle_responses(a = rep(all.sim.vals[i,3], all.sim.vals[i,1]),b_z = seq(all.sim.vals[i,4], 2, length.out=all.sim.vals[i,1]),a_z = rep(all.sim.vals[i,3], all.sim.vals[i,1]), 
#    b = genDiffGRM(num_items = all.sim.vals[i,1], num_categories = all.sim.vals[i,2], min = all.sim.vals[i,8]), muVals = c(0,0,0), varCovMat2, all.sim.vals[i,5])
#out.data3 <- simulate_hurdle_responses(a = rep(all.sim.vals[i,3], all.sim.vals[i,1]),b_z = seq(all.sim.vals[i,4], 2, length.out=all.sim.vals[i,1]),a_z = rep(all.sim.vals[i,3], all.sim.vals[i,1]), 
#    b = genDiffGRM(num_items = all.sim.vals[i,1], num_categories = all.sim.vals[i,2], min = all.sim.vals[i,8]), muVals = c(0,0,0), varCovMat3, all.sim.vals[i,5])

## Now esimtate models here


tmp.file <- paste("./data/", random.int, "_tabs.csv", sep='')
write.csv(out.data1$tabs, tmp.file, quote=F, row.names=F)
out.directory <- paste('./data/simdir_',all.sim.vals[i,10] ,sep='')
out.file <- paste(out.directory,'/fileVal_', i,'.RData', sep='')

val1 <- system(paste("julia ./scripts/juliaCode/mHurdleFlex.jl",tmp.file), intern = TRUE)
#write.csv(out.data2$tabs, tmp.file, quote=F, row.names=F)
#val2 <- system(paste("julia ./scripts/juliaCode/mHurdleFlex.jl",tmp.file), intern = TRUE)
#write.csv(out.data3$tabs, tmp.file, quote=F, row.names=F)
#val3 <- system(paste("julia ./scripts/juliaCode/mHurdleFlex.jl",tmp.file), intern = TRUE)
#system(paste("rm ", tmp.file))
## Now gather parameter estimates
mod.est1 <- return_Mod_Params(val1, out.data1)
#mod.est2 <- return_Mod_Params(val2, out.data2)
#mod.est3 <- return_Mod_Params(val3, out.data3)

# Scale Scores
theta1 <- seq(-6,6,0.25) # Quadrature points for theta1
theta2 <- seq(-6,6,0.25) # Quadrature points for theta1
theta <- expand.grid(theta1, theta2)
names(theta) <- c("theta1", "theta2")
a_grm <- mod.est1$irtGRMParam$estDis
b_grm <- mod.est1$irtGRMParam[,grep(x = colnames(mod.est1$irtGRMParam), pattern="estDif")]
a_2PL <- mod.est1$irt2PLParam$estDis
b_2PL <- mod.est1$irt2PLParam$estDif

itemtrace2PL <- trace.line.pts.2PL(a_2PL, b_2PL, theta)
itemtraceGRM <- trace.line.pts.grm(a_grm, b_grm, theta)
itemtrace <- trace.line.pts(a_grm, b_grm, a_2PL, b_2PL, theta)

rho <- mod.est1$rho$estRho[1]

qpoints <- theta
prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
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

## Now attach estimated parameters
out.data1$estGRMVals <- mod.est1$irtGRMParam
out.data1$est2plVals <- mod.est1$irt2PLParam
out.data1$estRho <- mod.est1$rho

## Now save output
out.directory <- paste('./data/simdir_',all.sim.vals[i,10] ,sep='')
out.file <- paste(out.directory,'/fileVal_', i,'.RData', sep='')
dir.create(out.directory, showWarnings = FALSE)
saveRDS(out.data1, file=out.file)