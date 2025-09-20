# Generate GRM (graded response model) difficulty/threshold parameters for items.
#
# For a GRM with m categories, each item has m-1 ordered thresholds.
# This function:
# 1) Creates a base set of m-1 equidistant thresholds between [min, max]
# 2) Replicates them across items
# 3) Adds random noise
# 4) Sorts thresholds within each item to ensure increasing order
#
# Inputs:
# - num_items: number of items to generate
# - num_categories: number of response categories (m); thresholds per item = m - 1
# - min, max: range over which the base thresholds are spaced
# - rnorm_var: standard deviation of Gaussian noise added to thresholds
#
# Output:
# - A matrix of size num_items x (num_categories - 1)
#   Each row contains the ordered thresholds for one item.
genDiffGRM <- function(num_items = 20, num_categories = 5, min = 0, max = 2, rnorm_var = 0.5) {
  # Number of thresholds per item (m - 1)
  num_dif <- num_categories - 1
  # Base equidistant thresholds between min and max (length = m - 1)
  equi.vec <- seq(min, max, length.out = num_dif)
  # Replicate the base thresholds for all items (each row = one item)
  out.mat <- matrix(equi.vec, nrow = num_items, ncol = num_dif, byrow = TRUE)
  # Add Gaussian noise to each threshold to introduce item variability
  out.mat <- out.mat + rnorm(num_items * num_dif, sd = rnorm_var)
  # Ensure thresholds are ordered within each item (ascending), as required by GRM
  out.mat <- t(apply(out.mat, 1, sort))
  # Return item-by-threshold matrix
  return(out.mat)
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

# trace.line.pts.2PL <- function(a_z, b_z, theta) {
#   D <- 1.702
#   nitems = length(a_z)
#   itemtrace2PL <- list()
#   probPos <- matrix(NA, nrow = nitems, ncol = length(theta[,1]))
#   probNeg <- matrix(NA, nrow = nitems, ncol = length(theta[,1]))
#   for (j in 1:nitems){
#     probNeg[j,] <- 1 - pnorm((a_z[j]/D) * (theta[,1] - b_z[j]))
#     probPos[j,] <- pnorm((a_z[j]/D) * (theta[,1] - b_z[j]))
#   }
#   itemtrace2PL[[1]] <- probNeg
#   itemtrace2PL[[2]] <- probPos
#   return(itemtrace2PL)
# }

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

# trace.line.pts.grm <- function(a, b, theta) {
#   D <- 1.702
#   nitems <- length(a)
#   n_categories <- dim(b)[2] + 1
#   itemtraceGRM <- rep(list(matrix(NA, nrow=nitems, ncol=dim(theta)[1])), n_categories)
#   for (j in 1:nitems){
#     for (k in 1:n_categories){
#       if (k == 1){
#         itemtraceGRM[[1]][j,] <- 1 - pnorm((a[j]/D) * (theta[,2] - b[j,k]))
#       }
#       if (k > 1 & k < n_categories){
#         itemtraceGRM[[k]][j,] <- pnorm((a[j]/D) * (theta[,2] - b[j,k-1])) -
#           pnorm((a[j]/D) * (theta[,2] - b[j,k]))
#       }
#       if (k == n_categories){
#         itemtraceGRM[[k]][j,] <- pnorm((a[j]/D) * (theta[,2] - b[j,k-1]))
#       }
#     }
#   }
#   return(itemtraceGRM)
# }

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
  for (i in 1:length(new.deriv)) {
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
  out.info.add.multf <- out.info.add
  
  for (i in 1:nrow(out.info.add)) {
    out.info.add.mult[i, ] <- out.info.add[i, ] * theta2pl_probs_one[i, ]
    out.info.add.multf[i, ] <- theta2pl_probs_one[i, ]
  }
  
  # 8) Aggregate to test information per theta (sum over items)
  out.info.sum <- colSums(out.info.add.mult)
  out.info.sumf <- colSums(out.info.add.multf)

  # 9) Reliability via weighted mean of inverse information under BVN prior
  out.info.sum.inv <- out.info.sum^-1
  varCovMat <- matrix(c(1, rhoVal, rhoVal, 1), ncol = 2)
  weights <- mvtnorm::dmvnorm(theta.grid, mean = muVals, sigma = varCovMat)
  out.rel <- 1 / (1 + weighted.mean(out.info.sum.inv, weights, na.rm=TRUE))
  
  # Package outputs
  out.list <- list(
    out.rel   = out.rel,        # scalar reliability
    test.info = out.info.sum,   # vector of test information at each theta point
    item.info = out.info        # list of per-category item information matrices
  )
  return(out.list)
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

# Compute category response probabilities for a GRM (Graded Response Model)
# given a discrimination parameter (a), a vector of thresholds (b), and a latent trait (theta).
#
# Inputs:
#   a     : numeric scalar, item discrimination parameter
#   b     : numeric vector of length (K-1), ordered category thresholds
#   theta : numeric scalar, latent trait value
#
# Details:
#   - For a GRM with K categories, there are K-1 thresholds in b.
#   - The model uses cumulative category functions:
#       P(Y >= k | theta) = logistic( a * (theta - b_{k-1}) ) for k = 2..K
#     with the convention P(Y >= 1 | theta) = 1.
#   - Category probabilities are computed by differences of cumulative probabilities:
#       P(Y = 1) = 1 - P(Y >= 2)
#       P(Y = k) = P(Y >= k) - P(Y >= k+1) for k = 2..K-1
#       P(Y = K) = P(Y >= K)
#   - This function returns a length-K vector of category probabilities that should sum to 1.
#   - The code uses explicit exp()/(1+exp()) rather than plogis() for the logistic,
#     which is equivalent but can be less numerically stable for large |a*(theta - b)|.
#
# Output:
#   A numeric vector of length K = length(b) + 1 with probabilities for categories 1..K.
itemtraceGRM <- function(a, b, theta){
  # Number of observed categories (K): thresholds + 1
  num_categories = length(b) + 1
  
  # Initialize output vector
  out_prob <- rep(NA, num_categories)
  
  # Loop over categories k = 1..K
  for(k in 1:num_categories){
    if(k==1){
      # First category probability:
      # P(Y=1) = 1 - P(Y >= 2) = 1 - logistic( a*(theta - b1) )
      out_prob[k] <- 1 - exp(a*(theta - b[k])) / (1 + exp(a*(theta - b[k])))
    }
    if(k == num_categories){
      # Last category probability:
      # P(Y=K) = P(Y >= K) = logistic( a*(theta - b_{K-1}) )
      out_prob[k] <- exp(a*(theta - b[k-1]))/(1+exp(a*(theta - b[k-1])))
    }
    else if (k != 1 & k != num_categories){
      # Middle categories: difference of cumulative logits
      # P(Y=k) = P(Y >= k) - P(Y >= k+1)
      #        = logistic( a*(theta - b_{k-1}) ) - logistic( a*(theta - b_k) )
      out_prob[k] <- (exp(a*(theta - b[k-1]))/(1+exp(a*(theta - b[k-1])))) - (exp(a*(theta - b[k]))/(1+exp(a*(theta - b[k]))))
    }
  }
  
  # Returns length-K probability vector summing approximately to 1
  return(out_prob)
}

# Compute the likelihood of a response pattern over a grid of quadrature points
# using precomputed item trace probabilities for a hurdle/IRT model.
#
# Inputs:
#   response_pattern : numeric vector of length nitems, containing observed category
#                      responses for each item. Assumes categories are coded 0..(K-1).
#   itemtrace        : list of length K (number of categories), where each element is an
#                      nitems x nq matrix of probabilities:
#                        itemtrace[[k+1]][i, q] = P(Y_i = k | theta_q)
#                      for item i and quadrature point q. Note the +1 shift to map
#                      category k (0-based) to list index k+1 (1-based in R).
#   qPoints          : data.frame or matrix of size nq x ndim (not used directly here
#                      except for its number of rows to size the likelihood vector).
#
# Details:
#   - The function assumes local independence given theta: the likelihood at each
#     quadrature point is the product over items of the item-category probabilities.
#   - For each item, it finds the observed category, maps it to the appropriate
#     slice of the itemtrace list, and multiplies the corresponding probability
#     vector across quadrature points into the running likelihood.
#
# Output:
#   lhood : numeric vector of length nq, where each entry is the likelihood of the
#           response pattern at the corresponding quadrature point.
score <- function(response_pattern, itemtrace, qPoints){
  # Initialize likelihood across all quadrature points to 1 (neutral element for multiplication)
  lhood <- rep(1, dim(qPoints)[1])
  
  # Number of items is the number of rows in any of the itemtrace matrices
  nitems <- dim(itemtrace[[1]])[1]
  
  # Loop over items and accumulate the product of probabilities
  for(item in 1:nitems){
    # Observed category response for this item (assumed 0..K-1)
    answerval = response_pattern[item]
    
    # Convert 0-based category to 1-based list index
    indexval = as.integer(answerval + 1)
    
    # Multiply the likelihood by P(Y_item = observed_category | theta_q) for all q
    lhood <- lhood*itemtrace[[indexval]][item,]
  }
  
  # Return the likelihood vector over quadrature points
  return(lhood)
}

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


simulate_hurdle_responses <- function(a, b, a_z, b_z, muVals, varCovMat, N=NULL, theta=NULL){
  require(dplyr)
  require(MASS)
  require(mvtnorm)
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


# Set seed
set.seed(16)
## Now generate some known values
## Fixed a_z across all sim
n.item <- 15
a_z <- rep(3, n.item) ## Generate n.item 2PL discrimination parameters, all identical 
a_z <- a_z
b_z <- runif(n.item, min=-2, max = 2) ## Generate n.item 2PL discrimination parameters randomly sampled between -2 & 0 
a <- rep(3, n.item) ## Generate n.item GRM discrimination parameters, all identical 
a <- a
b <- genDiffGRM(num_items = n.item, num_categories = 3, min = -2, max = 2, rnorm_var = 0.7) ## This function call generates k-1 GRM difficulty locaitons, where k == the number of GRM response categories 
## Fixed factor cor
rho <- .4
## Latent means == 0
muVals <- c(0,0)
## Now obtain the reliability across these
# hurdInfo_derivHurdleGRM = hurdInfo(a = a, b = b, a_z = a_z, b_z = b_z, rhoVal = rho, muVals = muVals) ## out.rel == 0.8596858
hurdInfo_grmBy2PL = hurdInfo_PLbyGRM(a = a, b = b, a_z = a_z, b_z = b_z, rhoVal = rho, muVals = muVals) ## out.rel == 0.7610786

## Now estimate the rel across various ranges of 2PL difficulty
# out.plot <- data.frame(twoPL = seq(-9, 9, 1), out.rel = NA)
# for(i in 1:nrow(out.plot)){
#   b_z <- runif(n.item, min=out.plot$twoPL[i], max = out.plot$twoPL[i]+2) ## Generate n.item 2PL discrimination parameters randomly sampled between -2 & 0 
#   out.plot$out.rel[i] = hurdInfo_PLbyGRM(theta.grid = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)), a = a, b = b, a_z = a_z, b_z = b_z, rhoVal = rho, muVals = muVals)$out.rel
# }
# plot(out.plot)


## Now examine the variance within a single theta pattern
## theta will be the mean for both Sus & Sev
theta <- matrix(0, nrow=100000, ncol = 2)
## Sim responses
n.item = 7 ## number of items
a = rep(2, n.item) ## GRM discrim
b = genDiffGRM(num_items = n.item, num_categories = 4, min = -3, max = 1, rnorm_var = .3) ## GRM difficulty
a_z = rep(2, n.item) ## 2PL discrim
b_z = runif(n.item, min = 1, max = 3) ## 2PL difficulty
varCovMat = matrix(c(1,rho,rho,1), ncol=2) ## true covar matrix for theta
## This next line geenrates the responses using the fixed theta values
out_reps = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z, muVals = muVals, varCovMat = varCovMat, theta = theta)
var(out_reps$theta) ## Here is the variance of the theta, EAP & MAP theta estimates
## Now obtain the test info for these items
theta.grid = expand.grid(seq(-6, 6, .2), seq(-6, 6, .2)) 
out.info <- hurdInfo_PLbyGRM(a = a, b = b, a_z = a_z, b_z = b_z, rhoVal = rho, muVals = muVals, theta.grid = theta.grid)
out.info$test.info[which(theta.grid[,1]==0 & theta.grid[,2]==0)]
## Now estimate the test info derived S.E.E.
1 / sqrt(out.info$test.info[which(theta.grid[,1]==0 & theta.grid[,2]==0)])
