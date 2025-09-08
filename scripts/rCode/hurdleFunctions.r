## Add a function which will efficiently identify row matches
## This code was taken from: https://github.com/tagteam/prodlim/blob/master/R/row.match.R
## many thanks to the tagteam group
row.match<-  function(x, table, nomatch=NA){
    if (inherits(table,"matrix")) table <- as.data.frame(table)
    if (is.null(dim(x))) x <- as.data.frame(matrix(x,nrow=1))
    cx <- do.call("paste",c(x[,,drop=FALSE],sep="\r"))
    ct <- do.call("paste",c(table[,,drop=FALSE],sep="\r"))
    match(cx,ct,nomatch=nomatch)
  }


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
    ## Now return this
    return(output.df) 
}

## Quick helper function
genDiffGRM <- function(num_items= 20, num_categories=5, min=0, max=2){
    diffs <- t(apply(matrix(runif(num_items*(num_categories-1), min = min, max = max), num_items), 1, cumsum))
    diffs <- -(diffs - rowMeans(diffs))
    d <- diffs# + rnorm(num_items, sd = .2)
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
simulate_hurdle_responses <- function(a, b, a_z, b_z, muVals, varCovMat, N=NULL, theta=NULL){
  ## First obtain the theta values
  if(!is.null(N) & is.null(theta)){  
    theta <- MASS::mvrnorm(N, mu = muVals, Sigma = varCovMat)
  }
  if(!is.null(theta) & is.null(N)){
    N = dim(theta)[1]
  }
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

  ## Now calculate the cross tabs
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
  ## Now make the column names
  col.names.sus <- paste("Sus_", 1:ncol(matrix.mplus1), sep='')
  col.names.sev <- paste("Sev_", 1:ncol(matrix.mplus2), sep='')
  matrix.mplus <- data.frame(matrix.mplus)
  colnames(matrix.mplus) <- c(col.names.sus, col.names.sev)
  ## Now add the true theta estimate using the real hurdle parameters
  ## Create any global vars needed
  qpoints <- expand.grid(seq(-6, 6, .2), seq(-6, 6, .2))
  prior <- mvtnorm::dmvnorm(qpoints, muVals, varCovMat)
  itemtrace = trace.line.pts(a, b, a_z, b_z, qpoints)
  theta = data.frame(theta1 = theta[,1],
                     theta2 = theta[,2])
  theta$eapSus <- NA
  theta$eapSev <- NA
  theta$mapSus <- NA
  theta$mapSev <- NA
  ## Now loop through every response option as specified in the crosstab
  ## Make a progress bar
  pb <- txtProgressBar(min = 0, max = nrow(crosstab), style = 3) # style = 3 for a full bar
  for(i in 1:nrow(crosstab)){
    pattern <- crosstab[i,1:ncol(crosstab)-1] ## pattern to score
    lhood <- score(pattern, itemtrace, qpoints) ## lik of response 
    denom <- lhood*prior
    sum.denom <- sum(denom)
    eap2PL_Hurdle <- sum(denom*qpoints$Var1)/sum.denom ## take sum over all values for 2PL model
    eapGRM_Hurdle <- sum(denom*qpoints$Var2)/sum.denom ## same for GRM model
    ## Now do the MAP estimate as well
    map2PL_Hurdle <- qpoints[which.max(denom),1] ## take the median value here
    mapGRM_Hurdle <- qpoints[which.max(denom),2] ## same for GRM
    ## Now put these estimates into the theta data frame
    ## First identify which rows in responses2 have these same response patterns
    #index <- apply(responses2, 1, function(a) apply(pattern, 1, function(b) all(a==b))) ## for some reason this line of code is very slow
    index <- which(row.match(responses2, pattern, nomatch = 0)==1)
    ## Now assign values
    theta$eapSus[index] <- eap2PL_Hurdle
    theta$eapSev[index] <- eapGRM_Hurdle
    theta$mapSus[index] <- map2PL_Hurdle
    theta$mapSev[index] <- mapGRM_Hurdle
    setTxtProgressBar(pb, i)
  }
  ## Now return everything
  out.data = list(responses = responses, theta = theta, tabs = crosstab,
  a = a, b=b, a_z = a_z, b_z = b_z, muVals = muVals, varCovMat = varCovMat, mplusMat = matrix.mplus)
  return(out.data)
}

simulate_hurdle_responses_fast <- function(
    a, b, a_z, b_z, muVals, varCovMat,
    N = NULL, theta = NULL,
    # scoring grid cache (optional)
    qpoints = NULL, prior = NULL, itemtrace = NULL,
    # control
    parallel = FALSE, progress = FALSE, grid_seq = seq(-6, 6, by = 0.2)
){
  # Dependencies
  requireNamespace("MASS")
  requireNamespace("mvtnorm")
  requireNamespace("data.table")
  requireNamespace("matrixStats")
  if (parallel) requireNamespace("future.apply")
  
  # 1) Latent traits
  if (!is.null(N) && is.null(theta)) {
    theta <- MASS::mvrnorm(N, mu = muVals, Sigma = varCovMat)
  }
  if (!is.null(theta) && is.null(N)) {
    N <- nrow(theta)
  }
  if (is.null(N) || is.null(theta)) stop("Provide either N or theta.")
  
  theta <- as.matrix(theta)
  if (ncol(theta) != 2) stop("theta must be N x 2 (susceptibility, severity).")
  
  # 2) Basic dims
  J <- length(a)
  if (length(a_z) != J || length(b_z) != J) stop("a_z and b_z must be length J.")
  if (!is.matrix(b) || nrow(b) != J) stop("b must be a J x m matrix.")
  m <- ncol(b)               # thresholds per item
  num_categories <- m + 2    # 0 for non-endorse + (m+1) positive categories
  
  # 3) Vectorized susceptibility (2PL) probabilities: P(endorse)
  #    plogis supports vectorized location/scale; build N x J matrices
  T1 <- matrix(theta[, 1], nrow = N, ncol = J)
  Bz <- matrix(rep(b_z, each = N), nrow = N, ncol = J)
  Az <- matrix(rep(a_z, each = N), nrow = N, ncol = J)
  p_endorse <- plogis(q = T1, location = Bz, scale = Az)  # N x J
  
  # 4) Draw responses per item without building a big 3D array
  responses <- matrix(NA_integer_, nrow = N, ncol = J)
  
  # Handy local function for fast categorical sampling from row-wise probabilities
  # probs: N x K, rows sum to 1
  .sample_rowwise_cat <- function(probs) {
    # cumulative probs
    cp <- matrixStats::rowCumsums(probs)
    # one uniform per row
    u <- stats::runif(nrow(cp))
    # first category where cp >= u
    idx <- max.col(cp >= u, ties.method = "first")
    # idx in 1..K; we return 0..K-1
    idx - 1L
  }
  
  # 5) GRM per item (vectorized inside item), combine with hurdle, sample
  T2 <- theta[, 2]
  for (j in seq_len(J)) {
    aj <- a[j]
    bj <- b[j, ]  # length m
    # cumulative probabilities for graded response model
    # C_k = P(Y >= k | theta2) = plogis(theta2, location = bj[k], scale = aj)
    C <- plogis(outer(T2, bj, FUN = function(x, y) x), location = matrix(bj, nrow = N, ncol = m, byrow = TRUE), scale = aj)
    # category probs (positive categories 1..m+1)
    # cat1 = 1 - C1; cat_k = C_{k-1} - C_k (2..m); cat_{m+1} = C_m
    if (m == 1L) {
      P_pos <- cbind(1 - C[, 1], C[, 1, drop = TRUE])
    } else {
      mid <- C[, 1:(m - 1), drop = FALSE] - C[, 2:m, drop = FALSE]
      P_pos <- cbind(1 - C[, 1, drop = TRUE], mid, C[, m, drop = TRUE])
    }
    # hurdle combine with susceptibility
    # full categories: 0 (non-endorse) + 1..(m+1) positive
    probs <- cbind(1 - p_endorse[, j], p_endorse[, j] * P_pos)  # N x (m+2)
    # numeric safety: small negative due to floating errors
    probs[probs < 0] <- 0
    rs <- rowSums(probs)
    # normalize just in case of rounding; mostly redundant but safe
    probs <- probs / rs
    # sample responses for item j
    responses[, j] <- .sample_rowwise_cat(probs)
  }
  
  # 6) Cross-tab of patterns and mapping using data.table
  DT <- data.table::as.data.table(responses)
  # Count unique patterns
  tabs <- DT[, .N, by = names(DT)]
  data.table::setnames(tabs, "N", "Count")
  # Attach pattern id to tabs
  tabs[, pattern_id := .I]
  # Map each row to its pattern_id via keyed join
  data.table::setkeyv(tabs, names(DT))
  DT[, pattern_id := tabs[DT, on = names(DT), x.pattern_id]]
  pattern_id <- DT[["pattern_id"]]
  
  # 7) Mplus matrices (vectorized)
  matrix.mplus1 <- (responses > 0L) + 0L
  matrix.mplus2 <- responses
  matrix.mplus2[matrix.mplus2 == 0L] <- NA_integer_
  matrix.mplus <- data.frame(
    matrix.mplus1,
    matrix.mplus2,
    check.names = FALSE
  )
  colnames(matrix.mplus) <- c(
    paste0("Sus_", seq_len(J)),
    paste0("Sev_", seq_len(J))
  )
  
  # 8) Scoring setup (EAP/MAP via grid), cached if provided
  if (is.null(qpoints)) {
    qpoints <- expand.grid(Var1 = grid_seq, Var2 = grid_seq)
  } else {
    # ensure names
    if (!all(c("Var1", "Var2") %in% names(qpoints))) {
      colnames(qpoints) <- c("Var1", "Var2")
    }
  }
  if (is.null(prior)) {
    prior <- mvtnorm::dmvnorm(as.matrix(qpoints), mean = muVals, sigma = varCovMat)
  }
  if (is.null(itemtrace)) {
    # user-supplied function must exist
    if (!exists("trace.line.pts", mode = "function")) {
      stop("trace.line.pts function not found. Provide itemtrace or define the function.")
    }
    itemtrace <- trace.line.pts(a, b, a_z, b_z, qpoints)
  }
  
  # Prepare result holders for theta + EAP/MAP
  theta_df <- data.frame(
    theta1 = theta[, 1],
    theta2 = theta[, 2]
  )
  theta_df$eapSus <- NA_real_
  theta_df$eapSev <- NA_real_
  theta_df$mapSus <- NA_real_
  theta_df$mapSev <- NA_real_
  
  # 9) Build representative patterns and index sets once
  item_cols <- names(DT)[seq_len(J)]
  # Split row indices by pattern_id
  idx_by_pat <- split(seq_len(N), pattern_id)
  # Representative pattern per id as integer vector
  reps <- lapply(idx_by_pat, function(idx) as.integer(responses[idx[1], ]))
  
  # 10) Score patterns: sequential or parallel
  score_one <- function(p) {
    # p is integer vector length J
    # Expected: score(pattern, itemtrace, qpoints) returns likelihood over grid (vector length nrow(qpoints))
    lhood <- score(p, itemtrace, qpoints)
    post <- lhood * prior
    Z <- sum(post)

    # EAPs
    eap1 <- sum(post * qpoints$Var1) / Z
    eap2 <- sum(post * qpoints$Var2) / Z
    # MAP (grid argmax)
    k <- which.max(post)
    map1 <- qpoints$Var1[k]
    map2 <- qpoints$Var2[k]
    c(eap1 = eap1, eap2 = eap2, map1 = map1, map2 = map2)
  }
  
  if (parallel) {
    # parallel over patterns
    # user must set a plan() beforehand; otherwise future.apply uses sequential
    stats_list <- future.apply::future_lapply(reps, score_one, future.seed = TRUE)
  } else {
    if (progress) {
      pb <- utils::txtProgressBar(min = 0, max = length(reps), style = 3)
    }
    stats_list <- vector("list", length(reps))
    for (i in seq_along(reps)) {
      stats_list[[i]] <- score_one(reps[[i]])
      if (progress) utils::setTxtProgressBar(pb, i)
    }
    if (progress) close(pb)
  }
  
  # 11) Assign scores back to all rows per pattern
  for (i in seq_along(idx_by_pat)) {
    idx <- idx_by_pat[[i]]
    s <- stats_list[[i]]
    theta_df$eapSus[idx] <- s[["eap1"]]
    theta_df$eapSev[idx] <- s[["eap2"]]
    theta_df$mapSus[idx] <- s[["map1"]]
    theta_df$mapSev[idx] <- s[["map2"]]
  }
  
  # 12) Return
  list(
    responses = responses,
    theta = theta_df,
    tabs = as.data.frame(tabs),
    a = a, b = b, a_z = a_z, b_z = b_z,
    muVals = muVals, varCovMat = varCovMat,
    mplusMat = matrix.mplus,
    qpoints = qpoints, prior = prior, itemtrace = itemtrace
  )
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

# Compute 2PL (two-parameter logistic) item trace lines for a dichotomous indicator
# given a latent trait dimension (theta[,1]).
#
# Inputs:
# - a_z: vector of discriminations (slopes) for the 2PL indicator, length = nitems
# - b_z: vector of difficulties (thresholds) for the 2PL indicator, length = nitems
# - theta: matrix of person parameters; uses column 1 (theta[,1]) for this 2PL
#
# Output:
# - A list of length 2:
#   [[1]]: matrix of P(Z = 0 | theta) for each item (rows) across persons (cols)
#   [[2]]: matrix of P(Z = 1 | theta) for each item (rows) across persons (cols)
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

# Compute GRM (graded response model) category trace lines given a latent trait
# dimension (theta[,2]).
#
# In the GRM, the category probabilities are derived from cumulative (graded)
# logistic functions. For item j with m categories, there are m-1 thresholds b[j,k].
#
# Inputs:
# - a: vector of item discriminations (slopes), length = nitems
# - b: matrix of ordered thresholds, dimension = nitems x (m-1)
#      thresholds increase with k for proper GRM behavior
# - theta: matrix of person parameters; uses column 2 (theta[,2]) for this GRM
#
# Output:
# - A list of length m (number of categories).
#   Each element is an nitems x nperson matrix of P(Y = k | theta) for category k.
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

# Combine 2PL (binary) and GRM (polytomous) trace lines into a single set of
# category trace lines.
#
# Interpretation:
# - Assume a two-part model with a binary indicator Z ~ 2PL (using theta[,1])
#   and, conditional on Z=1, a polytomous response Y ~ GRM (using theta[,2]).
# - The resulting categories are:
#   k = 0: the "structural zero" or non-endorsement category tied to Z=0
#   k = 1..m: the m GRM categories (where m = number of GRM categories),
#             each weighted by P(Z=1) from the 2PL.
#
# Inputs:
# - a: vector (length nitems) of GRM discriminations
# - b: matrix (nitems x (m-1)) of GRM thresholds; m = dim(b)[2] + 1 categories
# - a_z: vector (length nitems) of 2PL discriminations
# - b_z: vector (length nitems) of 2PL difficulties
# - theta: matrix of person parameters; uses:
#          theta[,1] for the 2PL, theta[,2] for the GRM
#
# Output:
# - A list of length (m + 1) = dim(b)[2] + 2.
#   Each element is an nitems x nperson matrix.
#   Indexing:
#     [[1]]         = P(k=0 | theta) = P(Z=0 | theta) from 2PL
#     [[k+1]] (k>0) = P(k | theta)   = P(Z=1 | theta) * P_GRM(Y=k | theta)
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

dmnorm <- function(mu, sigma, x){
  k <- ncol(sigma)
  x <- t(x)
  dmn <- exp((-1/2)*diag(t(x-mu)%*%solve(sigma)%*%(x- 
                                                     mu)))/sqrt(((2*pi)^k)*det(sigma))  
  dmn
}

trace.line.pts.grm.expr <- function(a, b, var = "t") {
  nitems <- length(a)
  K <- ncol(b) + 1
  t_sym <- as.name(var)
  
  # Each category k contains an nitems x 1 matrix of expressions (one per item)
  expr_list  <- vector("list", K)
  dexpr_list <- vector("list", K)
  
  for (k in 1:K) {
    mat_expr  <- matrix(vector("list", nitems), nrow = nitems, ncol = 1)
    mat_dexpr <- matrix(vector("list", nitems), nrow = nitems, ncol = 1)
    
    for (j in 1:nitems) {
      if (k == 1) {
        # P(Y=1) = 1 - s(a*(t - b1))
        expr <- substitute(
          1 - exp(aj*(tt - bj1)) / (1 + exp(aj*(tt - bj1))),
          list(aj = a[j], bj1 = b[j, 1], tt = t_sym)
        )
      } else if (k < K) {
        # P(Y=k) = s(a*(t - b_{k-1})) - s(a*(t - b_k))
        expr <- substitute(
          exp(aj*(tt - bjm1))/(1 + exp(aj*(tt - bjm1))) -
            exp(aj*(tt - bj ))/(1 + exp(aj*(tt - bj ))),
          list(aj = a[j], bjm1 = b[j, k - 1], bj = b[j, k], tt = t_sym)
        )
      } else {
        # P(Y=K) = s(a*(t - b_{K-1}))
        expr <- substitute(
          exp(aj*(tt - bjm1)) / (1 + exp(aj*(tt - bjm1))),
          list(aj = a[j], bjm1 = b[j, K - 1], tt = t_sym)
        )
      }
      
      mat_expr[j, 1]  <- list(expr)
      # Derivative wrt t (symbol var)
      mat_dexpr[j, 1] <- list(D(expr, var))
    }
    
    expr_list[[k]]  <- mat_expr
    dexpr_list[[k]] <- mat_dexpr
  }
  
  list(expr = expr_list, dexpr = dexpr_list, var = var)
}

# Helper: evaluate all expressions at numeric t values (vector) and return an array
# Result dims: (nitems) x (length(t)) x (K categories)
eval.trace.expr <- function(expr_obj, tvals) {
  expr_list <- expr_obj$expr
  var <- expr_obj$var
  K <- length(expr_list)
  nitems <- nrow(expr_list[[1]])
  
  out <- array(NA_real_, dim = c(nitems, length(tvals), K))
  for (k in 1:K) {
    for (j in 1:nitems) {
      f <- expr_list[[k]][j, 1][[1]]
      out[j, , k] <- vapply(tvals, function(tv) eval(f, list2env(setNames(list(tv), var))), numeric(1))
    }
  }
  out
}

# Helper: evaluate all derivative expressions at numeric t values (vector) and return an array
# Result dims: (nitems) x (length(t)) x (K categories)
eval.dtrace.expr <- function(expr_obj, tvals) {
  dexpr_list <- expr_obj$dexpr
  var <- expr_obj$var
  K <- length(dexpr_list)
  nitems <- nrow(dexpr_list[[1]])
  
  out <- array(NA_real_, dim = c(nitems, length(tvals), K))
  for (k in 1:K) {
    for (j in 1:nitems) {
      f <- dexpr_list[[k]][j, 1][[1]]
      out[j, , k] <- vapply(tvals, function(tv) eval(f, list2env(setNames(list(tv), var))), numeric(1))
    }
  }
  out
}

# Compute test and item information for the original hurdle GRM approach:
# - Information for the GRM part is computed on the GRM probabilities directly,
#   then scaled by P(Z=1) from the 2PL (hurdle) component.
#
# Model structure:
# - Z ~ 2PL on theta[,1]; if Z=1, Y ~ GRM on theta[,2].
#
# Steps:
# 1) Compute P(Z=1) from 2PL and P_GRM(Y=k) from GRM.
# 2) Compute numerical derivatives d/dtheta of GRM category probabilities via
# 3) Compute category information: (dP_k/dtheta)^2 / P_k (for each GRM category k).
# 4) Sum category information across GRM categories per item/person.
# 5) Multiply item-level GRM information by P(Z=1) to get hurdle information.
# 6) Sum across items to get test information at each theta point.
# 7) Compute reliability using a bivariate normal weighting.
#
hurdInfo_PLbyGRM <- function(theta.grid = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)),
                       a, b, a_z, b_z,
                       muVals = c(0,0), rhoVal = .2,
                       eps = 1e-5) {
  
  # 1) Probabilities
  # P(Z=1) from 2PL (uses theta[,1]); matrix: nitems x nperson
  theta2pl_probs_one <- trace.line.pts.2PL(a_z, b_z, theta.grid)[[2]]
  
  # GRM category probabilities list (length m); each element: nitems x nperson
  thetagrm_probs_two <- trace.line.pts.grm(a, b, theta.grid)
  
  # 2) Numerical derivative setup (perturb theta in both dims by +/- eps)
  grm_minus <- theta.grid - eps
  grm_plus  <- theta.grid + eps
  
  # GRM probabilities at perturbed grids (used for central differences)
  thetagrm_probs_twoMin <- trace.line.pts.grm(a, b, grm_minus)
  thetagrm_probs_twoMax <- trace.line.pts.grm(a, b, grm_plus)
  
  # 3) Central difference derivatives per GRM category
  new.deriv <- thetagrm_probs_two
  for (i in 1:length(out.info)) {
    # dP_k/dtheta ≈ [P_k(theta + eps) - P_k(theta - eps)] / (2*eps)
    new.deriv[[i]] <- (thetagrm_probs_twoMax[[i]] - thetagrm_probs_twoMin[[i]]) / (2 * eps)
  }
  
  # 4) Square derivatives for information calculation: (dP_k/dtheta)^2
  new.derivsq <- new.deriv
  for (i in 1:length(new.derivsq)) {
    new.derivsq[[i]] <- new.derivsq[[i]]^2
  }
  
  # 5) Per-category information: I_k(theta) = (dP_k/dtheta)^2 / P_k(theta)
  out.info <- thetagrm_probs_two
  for (i in 1:length(out.info)) {
    out.info[[i]] <- new.derivsq[[i]] / thetagrm_probs_two[[i]]
  }
  
  # 6) Sum across GRM categories (per item-person matrix)
  out.info.add <- out.info[[1]]
  for (i in 2:length(out.info)) {
    out.info.add <- out.info.add + out.info[[i]]
  }
  
  # 7) Apply hurdle weighting: multiply by P(Z=1) item-wise
  # Result remains nitems x nperson
  out.info.add.mult <- out.info.add
  for (i in 1:nrow(out.info.add)) {
    out.info.add.mult[i, ] <- out.info.add[i, ] * theta2pl_probs_one[i, ]
  }
  
  # 8) Aggregate to test information per theta (sum over items)
  out.info.sum <- colSums(out.info.add.mult)
  
  # 9) Reliability via weighted mean of inverse information under BVN prior
  out.info.sum.inv <- out.info.sum^-1
  varCovMat <- matrix(c(1, rhoVal, rhoVal, 1), ncol = 2)
  weights <- mvtnorm::dmvnorm(theta.grid, mean = muVals, sigma = varCovMat)
  out.rel <- 1 / (1 + weighted.mean(out.info.sum.inv, weights))
  
  # Package outputs
  out.list <- list(
    out.rel   = out.rel,        # scalar reliability
    test.info = out.info.sum,   # vector of test information at each theta point
    item.info = out.info        # list of per-category item information matrices
  )
  return(out.list)
}

## Now try the hurdle function as just the mean of the 2PL info and the GRM info
## I think this will makes things much easier to use?
estHurdleRel <- function(simVals, a, b, a_z, b_z, thetaVals = expand.grid(seq(-7, 7, .1),seq(-7, 7, .1)), muVals = c(0,0),rhoVal = .2 ){
  ## Estimate pre mirt model
  est.data <- simVals$mplusMat[,grep(pattern = "Sev", x = names(simVals$mplusMat))]
  est.data2 <- simVals$mplusMat[,grep(pattern = "Sus", x = names(simVals$mplusMat))]
  if(length(unique(unlist(lapply(apply(est.data, 2,table), length))))>1){
    ## Insert some artificial response into the data
    ## First identify which column has the issue
    unique.resp.vals <- lapply(apply(est.data, 2,table), length)
    ## Now idenitfy columns
    inject.vals <- which.min(unique.resp.vals)
    n.cats <- dim(b)[2]+1
    for(i in inject.vals){
      est.data[,i] <- sample(1:n.cats, size = nrow(est.data), replace=TRUE)
    }
  }
  sv <- suppressWarnings(mirt(data=est.data, 1,itemtype= 'graded', pars = 'values'))
  ## NOw organize true and est vals
  slopeInt <- matrix(0, length(a), dim(b)[2] + 1)
  ## Now make a dataframe of all of the discrim and diff values
  input.vals <- data.frame(cbind(a, b))
  for(i in 1:length(a)){
    input.vector <- unlist(as.vector(input.vals[i,]))
    slopeInt[i, ] <- traditional2mirt(input.vector, cls='graded', ncat=dim(b)[2]+1)
  }
  ## Now replace sv with these mirt values
  sv$value[sv$name == 'a1'] <- slopeInt[,1]
  ## Now do this for the rest of the difficulty values
  d.values <- grep(pattern="^d", x = names(table(sv$name)), value = TRUE)
  index.val <- 2
  for(i in d.values){
    sv$value[sv$name == i] <- slopeInt[,index.val]
    index.val <- index.val + 1
  }
  sv$est <- FALSE
  mod <- suppressWarnings(mirt(est.data, 1, pars=sv))
  ## Now obtain all of the information values
  vals.2pl <- trace.line.pts.2PL(a_z = a_z, b_z = b_z, theta = thetaVals)
  item.info.2pl <- vals.2pl[[1]] * vals.2pl[[2]] * (a_z^2)
  test.info.2pl <- apply(item.info.2pl, 2, sum)
  ## NOw make sure this is equivalent to the MIRT item info
  tmp <- mirt(est.data2, 1, pars='values', itemtype = "2PL")
  slopeInt <- matrix(0, length(a), 2)
  ## Now make a dataframe of all of the discrim and diff values
  input.vals <- data.frame(cbind(a_z, b_z))
  for(i in 1:length(a)){
    input.vector <- unlist(as.vector(input.vals[i,]))
    slopeInt[i, ] <- traditional2mirt(c(input.vector,0, 1), cls='2PL')[1:2]
  }
  tmp$value[tmp$name == 'a1'] <- slopeInt[,1]
  tmp$value[tmp$name == 'd'] <- slopeInt[,2]
  tmp$est <- FALSE
  mod2 <- suppressWarnings(mirt(est.data2, 1, pars=tmp))
  vals.2pl <- testinfo(mod2, Theta = thetaVals[,1],individual=TRUE)
  ## Now obtain the grm information here
  vals.grm <- testinfo(mod, Theta = thetaVals[,2],individual=TRUE)
  ## Now add these and take the mean
  test.info.combined <- apply(vals.2pl + vals.grm, 1, sum)
  test.info.combined <- test.info.combined / 2
  ## Now multiply by prob of theta value
  varCovMat = matrix(c(1,rhoVal,rhoVal,1), ncol=2)
  #test.info.weight <- test.info * dmnorm(muVals, varCovMat, thetaVals)
  ## Now estimate the rel with these info functions
  out.rel <- 1 / (1+(1/weighted.mean(test.info.combined, dmnorm(muVals, varCovMat, thetaVals))))
  return(out.rel)
}


estHurdleRelGen <- function(a, b, a_z, b_z, thetaVals = expand.grid(seq(-7, 7, .1),seq(-7, 7, .1)), muVals = c(0,0),rhoVal = .2 ){
  ## Estimate pre mirt model
  est.data <- simVals$mplusMat[,grep(pattern = "Sev", x = names(simVals$mplusMat))]
  est.data2 <- simVals$mplusMat[,grep(pattern = "Sus", x = names(simVals$mplusMat))]
  if(length(unique(unlist(lapply(apply(est.data, 2,table), length))))>1){
    ## Insert some artificial respones into the data
    ## FIrst idenitfy which column has the issue
    unique.resp.vals <- lapply(apply(est.data, 2,table), length)
    ## Now idenitfy columns
    inject.vals <- which.min(unique.resp.vals)
    n.cats <- dim(b)[2]+1
    for(i in inject.vals){
      est.data[,i] <- sample(1:n.cats, size = nrow(est.data), replace=TRUE)
    }
  }
  sv <- suppressWarnings(mirt(data=est.data, 1,itemtype= 'graded', pars = 'values'))
  ## NOw organize true and est vals
  slopeInt <- matrix(0, length(a), dim(b)[2] + 1)
  ## Now make a dataframe of all of the discrim and diff values
  input.vals <- data.frame(cbind(a, b))
  for(i in 1:length(a)){
    input.vector <- unlist(as.vector(input.vals[i,]))
    slopeInt[i, ] <- traditional2mirt(input.vector, cls='graded', ncat=dim(b)[2]+1)
  }
  ## Now replace sv with these mirt values
  sv$value[sv$name == 'a1'] <- slopeInt[,1]
  ## Now do this for the rest of the difficulty values
  d.values <- grep(pattern="^d", x = names(table(sv$name)), value = TRUE)
  index.val <- 2
  for(i in d.values){
    sv$value[sv$name == i] <- slopeInt[,index.val]
    index.val <- index.val + 1
  }
  sv$est <- FALSE
  mod <- suppressWarnings(mirt(est.data, 1, pars=sv))
  ## Now obtain all of the information values
  vals.2pl <- trace.line.pts.2PL(a_z = a_z, b_z = b_z, theta = thetaVals)
  item.info.2pl <- vals.2pl[[1]] * vals.2pl[[2]] * (a_z^2)
  test.info.2pl <- apply(item.info.2pl, 2, sum)
  ## NOw make sure this is equivalent to the MIRT item info
  tmp <- mirt(est.data2, 1, pars='values', itemtype = "2PL")
  slopeInt <- matrix(0, length(a), 2)
  ## Now make a dataframe of all of the discrim and diff values
  input.vals <- data.frame(cbind(a_z, b_z))
  for(i in 1:length(a)){
    input.vector <- unlist(as.vector(input.vals[i,]))
    slopeInt[i, ] <- traditional2mirt(c(input.vector,0, 1), cls='2PL')[1:2]
  }
  tmp$value[tmp$name == 'a1'] <- slopeInt[,1]
  tmp$value[tmp$name == 'd'] <- slopeInt[,2]
  tmp$est <- FALSE
  mod2 <- suppressWarnings(mirt(est.data2, 1, pars=tmp))
  vals.2pl <- testinfo(mod2, Theta = thetaVals[,1],individual=TRUE)
  ## Now obtain the grm information here
  vals.grm <- testinfo(mod, Theta = thetaVals[,2],individual=TRUE)
  ## Now add these and take the mean
  test.info.combined <- apply(vals.2pl + vals.grm, 1, sum)
  test.info.combined <- test.info.combined / 2
  ## Now multiply by prob of theta value
  varCovMat = matrix(c(1,rhoVal,rhoVal,1), ncol=2)
  #test.info.weight <- test.info * dmnorm(muVals, varCovMat, thetaVals)
  ## Now estimate the rel with these info functions
  out.rel <- 1 / (1 + (1 / weighted.mean(test.info.combined, dmnorm(muVals, varCovMat, thetaVals))))
  return(out.rel)
}


# Compute test and item information for a hurdle GRM model over a 2D theta grid,
# and summarize reliability by weighting information with a bivariate normal prior.
#
# Model structure:
# - Hurdle mechanism: Z ~ 2PL on theta[,1]; if Z>0, response Y ~ GRM on theta[,2].
# - Category probabilities for the hurdle model are:
#   k=0: P(Z=0)
#   k=1..m: P(Z=1) * P_GRM(Y=k)
#
# Information:
# - Uses the GRM item information approximation:
#   I(theta) = sum_k ( (dP_k/dtheta)^2 / P_k ), for k over observable categories.
# - Derivatives are obtained via a symmetric finite difference on theta.
#
# Inputs:
# - theta.grid: nperson x 2 matrix of ability points to evaluate (default is a grid)
# - a, b: GRM parameters (a: vector of slopes; b: matrix of thresholds)
# - a_z, b_z: 2PL parameters for the hurdle component
# - muVals: mean vector for weighting distribution (bivariate normal)
# - rhoVal: correlation for weighting distribution (variance fixed at 1 per dim)
# - h: step size for numerical derivative (not used explicitly
#
# Outputs (list):
# - out.rel: scalar reliability estimate = 1 / (1 + weighted.mean(1 / I_test(theta)))
# - test.info: vector (length = nrow(theta.grid)) of test information at each grid point
# - item.info: list of length m (GRM categories) with per-category item information matrices

hurdInfo <- function(theta.grid = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)),
                     a, b, a_z, b_z,
                     muVals = c(0, 0),
                     rhoVal = .2,
                     eps = 1e-5,
                     irtMod = NULL) {
  
  # 1) Get component probabilities
  # P(Z=1) from 2PL (uses theta[,1]); take the "positive" category [[2]]
  theta2pl_probs_one <- trace.line.pts.2PL(a_z, b_z, theta.grid)[[2]]
  
  # GRM category probabilities P_GRM(Y=k) (uses theta[,2]); list of length m
  thetagrm_probs_two <- trace.line.pts.grm(a, b, theta.grid)
  
  # 2) Numerical derivatives via symmetric finite differences for the hurdle probabilities
  # Create perturbed theta grids (small step applied to both dimensions)
  grm_minus <- theta.grid
  grm_minus[,2] <- grm_minus[,2] - eps
  grm_plus  <- theta.grid
  grm_plus[,2] <- theta.grid[,2] + eps
  
  # Hurdle model category probabilities at theta - eps, theta + eps, and theta (base)
  # Each is a list of length (m + 1): [[1]] = k=0 (Z=0), [[2..m+1]] = k=1..m
  thetagrm_probs_twoMin  <- trace.line.pts(a, b, a_z, b_z, grm_minus)
  thetagrm_probs_twoMax  <- trace.line.pts(a, b, a_z, b_z, grm_plus)
  thetagrm_probs_base    <- trace.line.pts(a, b, a_z, b_z, theta.grid)
  
  # Initialize derivative list to match GRM categories (length m)
  # We'll compute d/dtheta for the observable hurdle categories k = 1..m
  new.deriv <- thetagrm_probs_two
  
  # Map GRM categories (index 1..m) to hurdle list indices (2..m+1)
  index <- 1
  for (i in 2:length(thetagrm_probs_twoMax)) {
    # Central difference derivative for category k = i-1:
    # dP_k/dtheta ≈ [P_k(theta + eps) - P_k(theta - eps)] / (2*eps)
    new.deriv[[index]] <- (thetagrm_probs_twoMax[[i]] - thetagrm_probs_twoMin[[i]]) / (2 * eps)
    index <- index + 1
  }
  
  # Square the derivatives for information calculation
  new.derivsq <- new.deriv
  for (i in 1:length(new.derivsq)) {
    new.derivsq[[i]] <- new.derivsq[[i]]^2
  }
  
  # 3) Compute per-category item information:
  # I_k(theta) = (dP_k/dtheta)^2 / P_k(theta), for k = 1..m (exclude k=0)
  out.info <- thetagrm_probs_two
  for (i in 1:length(out.info)) {
    plus.index <- i + 1  # shift to match hurdle categories (skip k=0 at [[1]])
    out.info[[i]] <- new.derivsq[[i]] / thetagrm_probs_base[[plus.index]]
  }
  
  # 4) Sum information across categories to get item-level test information at each theta
  out.info.add <- out.info[[1]]
  for (i in 2:length(out.info)) {
    out.info.add <- out.info.add + out.info[[i]]
  }
  
  # 5) Aggregate to test information per grid point (sum over items = row-wise sum, then colSums)
  # Here, out.info.add is an nitems x nperson matrix; summing over items yields a length nperson vector
  out.info.sum <- colSums(out.info.add)
  
  # 6) Convert information to error variance and compute reliability
  out.info.sum.inv <- out.info.sum^-1
  
  # Bivariate normal weighting over theta.grid with mean muVals and correlation rhoVal
  varCovMat <- matrix(c(1, rhoVal, rhoVal, 1), ncol = 2)
  
  # Reliability approximation:
  # rel = 1 / (1 + E_theta[1 / I_test(theta)]), expectation taken under N(mu, Sigma)
  weights <- mvtnorm::dmvnorm(theta.grid, muVals, varCovMat)
  out.rel <- 1 / (1 + weighted.mean(out.info.sum.inv, weights))
  
  # Package outputs
  out.list <- list(
    out.rel   = out.rel,        # scalar reliability
    test.info = out.info.sum,   # vector of test information at each theta point
    item.info = out.info        # list of per-category item information matrices
  )
  return(out.list)
}
