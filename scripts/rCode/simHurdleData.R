## Clean slate
rm(list=ls())

## These functions will be used to create data for a hurdle model estimation

## First create a funciton which will simulate a 2PL model given:
## 1. Discrimination parameters --> vector of discrim parameters
## 2. Difficulty parameters --> vector of difficulty parameters
## 3. NUmber of participants --> integer
## 4. Theta patterns --> vector of \mu  & \sigma for normal dist

library(tidyverse)

sim2PLMod <- function(a, b, N, thetaP){
    theta <- rnorm(N, thetaP[1], thetaP[2])

    ## Now create output response data frame
    responses = matrix(NA, nrow = N, ncol = length(a))
    for(i in 1:N){
        for(j in 1:length(a)){
            p = 1 / (1 + exp(-a[j] * (theta[i] - b[j])))
            #p = 1-p
            responses[i,j] = rbinom(1, 1, p)
        }
    }
    ## Now write the output
    output.df = data.frame(responses, theta)
    ## NOw return this
    return(output.df) 
}

## Now create a function which will sample from a multivaraite 2PL model
## Similar to the 2 PL model, but now the discrimination parameters will be multivaraite
## As well as the patterns from the normal distritbuion which will be used to sample data from
sim2PLMultiMod <- function(muVals = c(0,0), varCovMat = diag(2), a = matrix(rnorm(20), ncol = 2), b = rnorm(10), N = 3000){
    ## First create the theta values
    theta <- MASS::mvrnorm(N, mu = muVals, Sigma = varCovMat)

    ## NOw create response pattern output matrix
    responses = matrix(NA, nrow = N, ncol = length(b))
    for(i in 1:N){
        for(j in 1:length(b)){
            p = sum(a[j,] * theta[i,]) - b[j]
            p = 1 / (1 + exp(-p))
            p= 1 - p
            responses[i,j] = rbinom(1, 1, p)
        }
    }
    ## Now write the output
    output.df = data.frame(responses, theta)
    ## NOw return this
    return(output.df)
}


## Quick helper function
genDiffGRM <- function(num_items= 20, num_categories=5, min=0, max=2){
    diffs <- t(apply(matrix(runif(num_items*(num_categories-1), min = min, max = max), num_items), 1, cumsum))
    diffs <- -(diffs - rowMeans(diffs))
    d <- diffs + rnorm(num_items)
    d <- t(apply(d, 1, sort))
    return(d)
}

## Create a function which will return values for a GRM model
# n: number of examinees
# m: number of items
# a: a vector of discrimination parameters (length m)
# b: a list of threshold parameters (each list element contains thresholds for that item)
# theta: latent trait levels for examinees (if NULL, they will be drawn from N(0, 1))
simulate_grm_data <- function(n, m, a, b, theta = NULL) {
  
  
  # Initialize response matrix
  responses <- matrix(NA, nrow = n, ncol = m)
  
  # Loop through items
  for (j in 1:m) {
    aj <- a[j]  # Discrimination parameter for item j
    bj <- b[j,]  # Threshold parameters for item j
    num_categories <- length(bj) + 1  # Number of response categories for item j
    
    # Loop through examinees
    for (i in 1:n) {
        ## Create response probs here
        resp_probs = itemtraceGRM(a=aj, b=bj, theta= )
      
    }
  }
  
  return(list(responses = responses, theta = theta))
}


## Create a function here which will return all of the probs from a GRM model
itemtraceGRM <- function(a, b, theta){
    num_categories = length(b) + 1
    out_prob <- rep(NA, num_categories)
    for(k in 1:num_categories){
        if(k==1){
            out_prob[k] <- 1 - exp(a*(theta - b[k])) / (1 + exp(a*(theta - b[k])))
        }
        if(k == num_categories){
            out_prob[k] <- exp(a*(theta - b[k-1]))/(1+exp(a*(theta - b[k-1])))
        }
        else if (k != 1 & k != num_categories){
            out_prob[k] <- (exp(a*(theta - b[k-1]))/(1+exp(a*(theta - b[k-1])))) - (exp(a*(theta - b[k]))/(1+exp(a*(theta - b[k]))))
        }
    }
    return(out_prob)
}

#tmp <-  simulate_grm_data(3000, 20, a = rnorm(10,mean = 2), b = genDiffGRM(num_items = 10, num_categories = 3), theta=rnorm(3000))
#vals <- rowSums(tmp$responses)
#plot(vals, tmp$theta)
#cor(vals, tmp$theta)

## Now I need to simulate the full hurdle model here
# a = 2PL discrim
# b = 2 PL diff
# a_z = GRM discrim vals (vector)
# b_z = GRM diff vals (matrix)
# mu = multivaraite means (vector)
# sigma = var cov mat (matrix)
# N = sample size
# k = categories
simulate_hurdle_responses <- function(a, b, a_z, b_z, muVals, varCovMat, N){
    ## First obtain the theta values
    theta <- MASS::mvrnorm(N, mu = muVals, Sigma = varCovMat)
    numItems = length(a)
    num_categories = dim(b)[2] + 2
    ## First go through an obtain all of the repsone probabilities for the 2PL protion of the model
    ## This is the susceptibility factor
    # Prep a matrix to store all of the prob of endorsments here
    responseProb2PL = matrix(NA, nrow = N, ncol = length(b_z))
    for(i in 1:N){
        for(j in 1:length(a)){
            #p = 1 / (1 + exp(-a[j] * (theta[i,1] - b[j])))
            p_z = 1 - exp(a_z[j] * (theta[i, 1] - b_z[j])) / (1 + exp(a_z[j] * (theta[i, 1] - b_z[j])))
            p_1 = 1 - p_z
            responseProb2PL[i,j] = p_1
        }
    }
    ## Now go through and obtain the response prob from the GRM portion
    initVector <- rep(NA, (N*numItems*(num_categories-1)))
    responseProbGRM = array(initVector, c(N, numItems, (num_categories-1)))
    ## Now go through and 4estimate response probs from the GRM
    for (j in 1:numItems) {
        aj <- a[j]  # Discrimination parameter for item j
        bj <- b[j,]  # Threshold parameters for item j
        # Loop through examinees
        for (i in 1:N) {
            category_probs <- itemtraceGRM(aj, bj, theta[i,2])  
            # Generate response based on category probabilities
            for(k in 1:num_categories-1){
                responseProbGRM[i, j,k] <- category_probs[k]
            }
        }
    }

    ## Now I need to turn this into response categories
    responses = matrix(NA, nrow = N, ncol = numItems)
    for(i in 1:N){
        for( j in 1:numItems){
            ## Create a vector of probabilities
            tmp_probs <- rep(NA, num_categories)
            for(k in 1:num_categories){
                if(k == 1){
                    tmp_probs[k] <- 1 - responseProb2PL[i,j]
                }else{
                    tmp_probs[k] <- ((responseProb2PL[i,j])) * responseProbGRM[i,j,k-1]
                }
            }
            responses[i,j] <- sample(1:num_categories, 1,prob = tmp_probs)-1
        }
    }
    ## Now caluclate the cross tabs
    responses2 <- data.frame(responses)
    crosstab <- responses2 %>%
        group_by(across(everything())) %>%
        summarise(Count = n(), .groups = "drop") %>%
        ungroup()
    ## Now return everything
    out.data = list(responses = responses, theta = theta, tabs = crosstab,
    a = a, b=b, a_z = a_z, b_z = b_z, muVals = muVals, varCovMat = varCovMat)
    return(out.data)
}

## Now take a make a function which will take the output of the julia function and return the parameters for the
## 2PL model, GRM model, the varCov mat, as well as theta estimates
return_Mod_Params <- function(juliaOutput, dataIn){
    ## Initialize all of the output
    ## First create the 2PL data frame with discrim and diff values
    nitems = dim(dataIn$responses)[2]
    irt2PLParam <- matrix(NA, nrow=nitems, ncol=2)
    ## Now do the GRM portion
    ncatgrm <- diff(range(dataIn$responses))
    ## First idenitfy the total number of difficult params needed
    irtGRMParam <- matrix(NA, nrow=dim(dataIn$responses)[2], ncol=ncatgrm)

    ## Clean up the Julia output here
    list_index <- length(juliaOutput)
    out_params <- juliaOutput[list_index]
    out_params <- strsplit(out_params, " ")[[1]]
    out_params <- gsub(pattern = "]", replacement="", x = out_params)
    out_params <- gsub(pattern = "\\[|\\]", replacement="", x = out_params)
    out_params <- gsub(pattern = ",", replacement="", x = out_params)
    out_params <- as.numeric(out_params)

    nitemsgrm = dim(dataIn)[2]
    paramGRM = nitems*ncatgrm
    param2PL = nitems*2
    nParmsPerItemGRM = ncatgrm
    nParmsPerItem2PL = 2
    nitemsgrm = nitems
    for (j in 1:nitems) {
        irtGRMParam[j,1] <- abs(out_params[(j-1)*nParmsPerItemGRM + 1])
        irt2PLParam[j,1] <- abs(out_params[nitemsgrm*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 1])
        irt2PLParam[j,2] <- out_params[nitemsgrm*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 2]
        for (k in 1:(ncatgrm-1)){
            irtGRMParam[j,k+1] <- out_params[(j-1)*nParmsPerItemGRM + 1 + k]
        }
    }
    ## Now sort the GRM params
    irtGRMParam[,2:dim(irtGRMParam)[2]] <- t(apply(irtGRMParam[,2:dim(irtGRMParam)[2]], 1, function(x) sort(x, decreasing = FALSE)))

    ## Now exponentiate the rho estimate
    rho_nexp = out_params[length(out_params)]
    rho = exp(rho_nexp) / (1 + exp(rho_nexp))
    rho_true = unique(dataIn$varCovMat[row(dataIn$varCovMat)!=col(dataIn$varCovMat)])

    ## Now create a data frame with the sim and the true params
    dataframeGRM = data.frame(trueDis = dataIn$a, estDis = irtGRMParam[,1], trueDif = dataIn$b, estDis = irtGRMParam[,-1],modelSource = "GRM")
    dataframe2PL = data.frame(trueDis = dataIn$a_z, estDis = irt2PLParam[,1], trueDif = dataIn$b_z, estDif = irt2PLParam[,2], modelSource = "2PL")
    dataframeRho = data.frame(trueRho = rho_true, estRho = rho)

    ## Now return all of these estimates
    out_vals = list(irt2PLParam = dataframe2PL, irtGRMParam = dataframeGRM, rho = dataframeRho)
    return(out_vals)
}

varCovMat <- matrix(c(1,.2,.2,1), nrow=2)
varCovMat <- matrix(c(c(1, .2, 0),c(.2, 1, .4), c(0, .4, 1)), nrow=3, byrow = TRUE)
n <- 5000
num.items <- 3
out.datA <- simulate_hurdle_responses(a = rep(3, num.items),b_z = seq(-1, 1, length.out=num.items),a_z = rep(3, num.items), b = genDiffGRM(num_items = num.items, num_categories = 3), muVals = c(0,0,0), varCovMat, 5000)
out.dat <- out.datA$responses
write.csv(out.dat, file="./data/testSimDat.csv", quote=F, row.names=F)
out.datT <- out.datA$tabs
write.csv(out.datT, file="./data/testSimDat2.csv", quote=F, row.names=F)
val <- system("julia ./scripts/juliaCode/mHurdleFlex.jl ./data/testSimDat.csv ./data/testSimDat2.csv", intern = TRUE)
checkVals <- return_Mod_Params(val, out.datA)

## Now create al of the simulated datasets here
sim.param.nitems <- c(5,9)
sim.param.ncates <- c(3,7)
sim.param.discri <- c(3)
sim.param.diffic <- c(-1,1)
sim.param.sample <- c(1000)
sim.param.faccor <- c(.2, .4, .7)
sim.param.critco <- c(.2, .4, .7)
sim.param.difgrm <- c(-1, 0)
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.diffic, sim.param.sample, sim.param.faccor, 
                            sim.param.critco, sim.param.difgrm)
sim.datasets <- 100
for(i in 1:nrow(all.sim.vals)){
    ## Frist create the varcovar matrix
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
    ## Now make the datasets here
    for(j in 1:sim.datasets){
        set.seed(j)
        cat("\r",paste("Outer loop: ", i, "Inner loop: ", j))
        out.directory <- paste('./data/simdir_', i,sep='')
        out.file <- paste(out.directory,'/fileVal_', j,'.RData', sep='')
        dir.create(out.directory, showWarnings = FALSE)
        ## Now simulate the data
        out.data1 <- simulate_hurdle_responses(a = rep(all.sim.vals[i,3], all.sim.vals[i,1]),b_z = seq(all.sim.vals[i,4], 2, length.out=all.sim.vals[i,1]),a_z = rep(all.sim.vals[i,3], all.sim.vals[i,1]), 
        b = genDiffGRM(num_items = all.sim.vals[i,1], num_categories = all.sim.vals[i,2], min = all.sim.vals[i,8]), muVals = c(0,0,0), varCovMat1, all.sim.vals[i,5])
        out.data2 <- simulate_hurdle_responses(a = rep(all.sim.vals[i,3], all.sim.vals[i,1]),b_z = seq(all.sim.vals[i,4], 2, length.out=all.sim.vals[i,1]),a_z = rep(all.sim.vals[i,3], all.sim.vals[i,1]), 
        b = genDiffGRM(num_items = all.sim.vals[i,1], num_categories = all.sim.vals[i,2], min = all.sim.vals[i,8]), muVals = c(0,0,0), varCovMat2, all.sim.vals[i,5])
        out.data3 <- simulate_hurdle_responses(a = rep(all.sim.vals[i,3], all.sim.vals[i,1]),b_z = seq(all.sim.vals[i,4], 2, length.out=all.sim.vals[i,1]),a_z = rep(all.sim.vals[i,3], all.sim.vals[i,1]), 
        b = genDiffGRM(num_items = all.sim.vals[i,1], num_categories = all.sim.vals[i,2], min = all.sim.vals[i,8]), muVals = c(0,0,0), varCovMat3, all.sim.vals[i,5])
        out.vals <- list(dat1 = out.data1, dat2 = out.data2, dat3 = out.data3)
        ## Now write the data
        saveRDS(out.vals, file=out.file)
    }
}
print("Done")

