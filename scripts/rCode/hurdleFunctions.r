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

genDiffGRM <- function(num_items=20, num_categories=5, min=0, max=2, rnorm_var = .5){
  ## Figure out how many difficulty params are needed
  num_dif = num_categories-1
  ## Now first create a vector of equidistant difficulty parameters for every item
  equi.vec <- seq(min, max, length.out = num_dif)
  ## Assign difficulty parameters here
  out.mat <- matrix(equi.vec, nrow=num_items, ncol=num_dif, byrow = TRUE)
  ## Add some error
  out.mat <- out.mat + rnorm(num_items * num_dif, sd = rnorm_var)
  out.mat <- t(apply(out.mat, 1, sort))
  return(out.mat)
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
    ## First go through an obtain all of the response probabilities for the 2PL portion of the model
    ## This is the susceptibility factor
    # Prep a matrix to store all of the prob of endorsements here
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
    ## Now prep the dataframe for MPlus hurdle extension
    ## Now make the matrix for these data
    matrix.mplus1 <- matrix(NA, nrow = nrow(responses), ncol=ncol(responses))
    matrix.mplus2 <- matrix(NA, nrow = nrow(responses), ncol=ncol(responses))
    
    ## Now make the binary indicators first
    matrix.mplus1<- responses
    matrix.mplus1[which(matrix.mplus1[,1:ncol(responses)]>0)] <- 1
    ## Now do severity indicators here
    matrix.mplus2 <- responses
    matrix.mplus2[which(matrix.mplus2==0)] <- NA
    matrix.mplus <- cbind(matrix.mplus1, matrix.mplus2)
    
    ## Now return everything
    out.data = list(responses = responses, theta = data.frame(theta), tabs = crosstab,
    a = a, b=b, a_z = a_z, b_z = b_z, muVals = muVals, varCovMat = varCovMat, mplusMat = matrix.mplus)
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
    ## First identify the total number of difficult params needed
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
    dataframeGRM = data.frame(trueDis = dataIn$a, estDis = irtGRMParam[,1], trueDif = dataIn$b, estDif = irtGRMParam[,-1],modelSource = "GRM")
    dataframe2PL = data.frame(trueDis = dataIn$a_z, estDis = irt2PLParam[,1], trueDif = dataIn$b_z, estDif = irt2PLParam[,2], modelSource = "2PL")
    dataframeRho = data.frame(trueRho = rho_true, estRho = rho)

    ## Now return all of these estimates
    out_vals = list(irt2PLParam = dataframe2PL, irtGRMParam = dataframeGRM, rho = dataframeRho)
    return(out_vals)
}

## Now do all factors used to estimate theta values here
trace.line.pts.2PL <- function(a_z, b_z, theta)	{
  nitems = length(a_z)
  itemtrace2PL <- list()
  probPos <- matrix(NA, nrow = nitems, ncol = length(theta[,1]))
  probNeg <- matrix(NA, nrow = nitems, ncol = length(theta[,1]))
  for (j in 1:nitems){
    probNeg[j,] <- 1 - exp(a_z[j]*(theta[,1] - b_z[j]))/(1+exp(a_z[j]*(theta[,1] - b_z[j])))
    probPos[j,] <- exp(a_z[j]*(theta[,1] - b_z[j]))/(1+exp(a_z[j]*(theta[,1] - b_z[j])))
  }
  itemtrace2PL[[1]] <- probNeg
  itemtrace2PL[[2]] <- probPos
  return(itemtrace2PL)
}

# GRM trace lines
trace.line.pts.grm <- function(a, b, theta)	{
  nitems <- length(a)
  n_categories <- dim(b)[2] + 1
  itemtraceGRM <- rep(list(matrix(NA, nrow=nitems, ncol=dim(theta)[1])), n_categories)
  for (j in 1:nitems){
    for (k in 1:n_categories){
      if (k == 1){
        itemtraceGRM[[1]][j,] <- 1 - exp(a[j]*(theta[,2] - b[j,k]))/(1+exp(a[j]*(theta[,2] - b[j,k])))
      }
      if (k > 1 & k < n_categories){
        itemtraceGRM[[k]][j,] <- exp(a[j]*(theta[,2] - b[j,k-1]))/(1+exp(a[j]*(theta[,2] - b[j,k-1]))) - exp(a[j]*(theta[,2] - b[j,k]))/(1+exp(a[j]*(theta[,2] - b[j,k])))
      }
      if (k == n_categories){
        itemtraceGRM[[k]][j,] <- exp(a[j]*(theta[,2] - b[j,k-1]))/(1+exp(a[j]*(theta[,2] - b[j,k-1])))
      }
    }
  }
  return(itemtraceGRM)
}

trace.line.pts <- function(a, b, a_z, b_z, theta){
  n_categories = dim(b)[2] + 2
  nitems = length(a)
  itemtrace = rep(list(matrix(NA, nrow=nitems, ncol=dim(theta)[1])), n_categories)
  for (j in 1:nitems){
    for(k in 0:n_categories-1){
        if(k == 0){
            itemtrace[[k+1]][j,] = trace.line.pts.2PL(a_z, b_z, theta)[[1]][j,]
        }
        else if(k > 0){
            itemtrace[[k+1]][j,] = trace.line.pts.2PL(a_z, b_z, theta)[[2]][j,] * trace.line.pts.grm(a, b, theta)[[k]][j,]
        }
    }
  }
  return(itemtrace)
}

score <- function(response_pattern, itemtrace, qPoints){
    lhood <- rep(1, dim(qPoints)[1])
    nitems <- dim(itemtrace[[1]])[1]
    for(item in 1:nitems){
        answerval = response_pattern[item]
        indexval = as.integer(answerval + 1)
        lhood <- lhood*itemtrace[[indexval]][item,]
    }
    return(lhood)
}