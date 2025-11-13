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

## Quick helper function
genDiffGRM <- function(num_items= 20, num_categories=5, min=0, max=2){
    diffs <- t(apply(matrix(runif(num_items*(num_categories-1), min = min, max = max), num_items), 1, cumsum))
    diffs <- -(diffs - rowMeans(diffs))
    d <- diffs# + rnorm(num_items, sd = .2)
    d <- t(apply(d, 1, function(x) sort(x, decreasing = TRUE)))
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
 out.mat <- t(apply(out.mat, 1, function(x) sort(x, decreasing = TRUE)))
 return(out.mat)
}

genDiffGRM <- function(num_items = 20, num_categories = 5, min = 0, max = 2, min_gap = 0.15) {
  stopifnot(num_categories >= 2)
  num_dif <- num_categories - 1
  
  # Generate as in your original function
  diffs <- t(apply(matrix(runif(num_items * num_dif, min = min, max = max), num_items), 1, cumsum))
  diffs <- -(diffs - rowMeans(diffs))
  
  # Enforce minimum spacing between adjacent thresholds
  out <- matrix(NA_real_, nrow = num_items, ncol = num_dif)
  for (i in seq_len(num_items)) {
    inc <- sort(diffs[i, ], decreasing = FALSE)
    gaps <- diff(inc)
    cur_min_gap <- min(gaps)
    
    if (is.finite(cur_min_gap) && cur_min_gap < min_gap) {
      # Scale around the mean so that all gaps are increased proportionally
      m <- mean(inc)
      scale_factor <- min_gap / cur_min_gap
      inc <- (inc - m) * scale_factor + m
    }
    
    out[i, ] <- sort(inc, decreasing = TRUE)
  }
  
  out
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


## Create a function which will simulate HURDLE responses
## this function is using the MIRT package and will require estimates to be
## in the slope intercept form that MIRT uses
sim_mirt_hurdle <- function(a, d, a_z, d_z, rho, THETA = NULL, 
                            N = 15000, 
                            theta_grid = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)),
                            estModel = FALSE) {
  
  require(mirt)
  require(dplyr)
  
  # Generate bivariate normal theta values if not provided
  muVals <- rep(0, 2)  # Mean of 0 for both dimensions
  varCovMat <- matrix(c(1, rho, rho, 1), ncol = 2)  # Correlation matrix
  
  # Set up item parameter matrices for mirt
  # Total items = 2 * number of hurdle items (2PL + corresponding graded response)
  n_item <- length(a_z) * 2
  n_diff <- ncol(d) + 1  # Number of difficulty parameters per item
  
  # Initialize discrimination parameter matrix (items x factors)
  a_mat <- matrix(0, nrow = n_item, ncol = 2)
  n_2pl <- n_item / 2  # Number of 2PL (hurdle) items
  grm_vec <- seq(n_2pl + 1, n_item, 1)  # Indices for graded response items
  
  # Assign discrimination parameters
  a_mat[1:n_2pl, 1] <- a_z        # 2PL items load on Factor 1 (susceptibility)
  a_mat[grm_vec, 2] <- a          # Graded items load on Factor 2 (severity)
  
  # Initialize intercept parameter matrix
  d_mat <- matrix(NA, nrow = n_item, ncol = n_diff)
  d_mat[1:n_2pl, 1] <- d_z        # 2PL difficulty parameters
  d_mat[grm_vec, 2:ncol(d_mat)] <- d  # Graded response difficulty parameters
  
  # Define item types for mirt
  item_type <- c(rep("2PL", n_2pl), rep("graded", n_2pl))
  
  # Simulate item responses
  sim_vals <- simdata(a = a_mat, d = d_mat, #sigma = varCovMat, 
                      Theta = THETA, itemtype = item_type)

  theta_grid_mat <- as.matrix(theta_grid)
  
  # Apply hurdle mechanism: set graded responses to NA when 2PL item = 0
  # This creates the hurdle structure where you must "pass" the 2PL to answer graded items
  pl2_dat <- sim_vals[, 1:n_2pl]           # 2PL responses (hurdle component)
  grm_dat <- sim_vals[, grm_vec]           # Graded responses (severity component)
  grm_dat[pl2_dat == 0] <- NA              # Missing if didn't endorse hurdle item
  sim_vals_w_na <- bind_cols(pl2_dat, grm_dat)
  ## NOw make the observed test score responses
  obs_dat = grm_dat
  obs_dat[is.na(obs_dat)] <- 0
  obs_dat = obs_dat + pl2_dat
  
  ## Ensure all response options have at least 1 endorsement
  if(length(unique(apply(sim_vals_w_na[,grm_vec], 2, function(x) length(table(x)))) )>1){
    ## Flag error status
    stop("Difficulty intercepts provided null endorsment pattern")
  }
  
  # Create mirt model specification string
  # F1 = susceptibility factor, F2 = severity factor, COV = their correlation
  model <- "
  F1 = 1-YYZ
  F2 = BBH-FNF
  COV = F1*F2
  "
  # Replace placeholders with actual item numbers
  model <- gsub(x = model, pattern = "YYZ", replacement = n_2pl)
  model <- gsub(x = model, pattern = "BBH", replacement = (n_2pl + 1))
  model <- gsub(x = model, pattern = "FNF", replacement = ncol(sim_vals))
  
  # Create "true" model with known parameters (for comparison/validation)
  pars_df <- mirt(sim_vals_w_na, model = model, itemtype = item_type, pars = 'values')
  
  # Fix all parameters to their true values
  # 2PL discrimination parameters (Factor 1)
  pars_df$value[pars_df$name == "a1"][1:n_2pl] <- a_z
  pars_df$est[pars_df$name == "a1"][1:n_2pl] <- FALSE
  
  # Graded response discrimination parameters (Factor 2)  
  pars_df$value[pars_df$name == "a2"][grm_vec] <- a
  pars_df$est[pars_df$name == "a2"][grm_vec] <- FALSE
  
  # 2PL difficulty parameters
  pars_df$value[pars_df$name == "d"] <- d_z
  pars_df$est[pars_df$name == "d"] <- FALSE
  
  # Graded response difficulty parameters
  ## First make sure every response has the proper number of rows
  grm_df_vals = paste("d", 1:ncol(d), sep='')
  count_vals = table(pars_df$name)
  for (i in 1:ncol(d)) {
    d_val <- paste("d", i, sep = '')
    pars_df$value[pars_df$name == d_val] <- d[, i]
    pars_df$est[pars_df$name == d_val] <- FALSE
  }
  ## Check for equality
  if(length(unique(count_vals[grm_df_vals]))>1){
    ## First isolate the rows that are stable
    stable_pre = pars_df[which(pars_df$class=="dich"),]
    stable_post = pars_df[which(pars_df$class=="GroupPars"),]
    ## Now grab the rows to be modified
    unstable_mod = pars_df[which(pars_df$class=="graded"),]
    ## Find the best candidate template
    full_candidate = c("a1", "a2", grm_df_vals)
    ## Find the lowest value with the full length parameter set
    full_candidate = names(which(table(unstable_mod$item) == length(full_candidate))[1])
    ## Now grab the prototype item set
    prototype = unstable_mod[which(unstable_mod$item==full_candidate),]
    ## First find the proper number of rows
    prop_nrow = nrow(prototype)
    ## Repeat the prototype the total number of GRM items
    rep_count = length(unique(unstable_mod$item))
    # Repeat the data frame n times
    prototype <- prototype[rep(seq_len(prop_nrow), times = rep_count), ]
    ## Now fix the item names
    prototype$item = rep(unique(names(table(unstable_mod$item))), each = prop_nrow)
    ## Now fix the slope & intercept values
    grm_df_vals = paste("d", 1:ncol(d), sep='')
    for (i in 1:ncol(d)) {
      d_val <- paste("d", i, sep = '')
      prototype$value[prototype$name == d_val] <- d[, i]
      prototype$est[prototype$name == d_val] <- FALSE
    }
    ## Now do slopes
    # Graded response discrimination parameters (Factor 2)  
    prototype$value[prototype$name == "a2"] <- a
    ## Now recombine everything and fix parnums
    pars_df <- bind_rows(stable_pre, prototype, stable_post)
    pars_df$parnum = 1:nrow(pars_df)
  }
  
  # Factor correlation
  pars_df$value[pars_df$name == "COV_21"] <- rho
  pars_df$est[pars_df$name == "COV_21"] <- FALSE
  
  # Fit model with fixed parameters
  sv2 <- mirt(sim_vals_w_na, model = model, itemtype = item_type, pars = pars_df)
  
  # Get factor scores from true parameter model
  fscores_truE <- mirt::fscores(sv2, method = "MAP")
  
  # Repeat trace line and information calculations for true model
  prob_2pl <- probtrace(sv2, Theta = theta_grid_mat)
  
  # Extract 2PL endorsement probabilities
  prob_2pl_cor <- grep(x = colnames(prob_2pl), pattern = ".P.1")
  prob_2pl_cor <- prob_2pl_cor[prob_2pl_cor <= n_item]
  prob_2pl_cor <- prob_2pl[, prob_2pl_cor]
  
  # Item information for graded items (Factor 2)
  item_info <- testinfo(sv2, Theta = theta_grid_mat, degrees = c(90, 0), 
                        individual = TRUE)[, grm_vec]
  test_info_2pl <- testinfo(sv2, Theta = theta_grid_mat, degrees = c(0, 90), 
                            individual = FALSE)
  # Hurdle-adjusted item information
  if (dim(theta_grid_mat)[1]==1){
    item_info<-t(item_info)
    prob_2pl_cor<-t(prob_2pl_cor)
  }
  item_info_hurd <- item_info*prob_2pl_cor
  
  test_info_hurd_tru <- rowSums(item_info_hurd)  # Total test information
  # Only true model results
  theta_vals <- data.frame(THETA, fscores_truE)
  colnames(theta_vals) <- c("trueSus", "trueSev", "truP_susMAP", "truP_sevMAP")
  
  test_info_hurd_tru_inv <- test_info_hurd_tru^-1
  
  out_vals <- list(
    theta_vals = theta_vals, 
    tru_mod = sv2,
    tru_hurd_test_info = test_info_hurd_tru,
    tru_hurd_item_info = item_info_hurd,#_tru,
    prob_endorse_tru = prob_2pl_cor,#_tru,
    test_info_2pl = test_info_2pl,
    grm_info = item_info,
    responses = sim_vals, 
    responses_na = sim_vals_w_na,
    obs_responses = obs_dat
    #      out_rel = out.rel
  )
  #  }
  
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

# Return test information matrix I(theta) at specified theta(s)
# Uses caches if provided to avoid repeated calls
ret_test_info_mat <- function(sv2 = NULL,
                              theta_grid = expand.grid(seq(-6, 6, 0.5), seq(-6, 6, 0.5)),
                              caches = NULL) {
  theta_grid_mat <- as.matrix(theta_grid)
  all_coef <- coef(sv2, simplify = TRUE)$items
  n_item <- nrow(all_coef)
  n_2pl <- n_item / 2
  grm_vec <- seq(n_2pl + 1L, n_item)
  
  if (!is.null(caches) && !is.null(caches$theta_grid) &&
      nrow(caches$theta_grid) == nrow(theta_grid_mat) &&
      all(caches$theta_grid == theta_grid_mat)) {
    prob_2pl <- caches$probtrace
    item_info <- caches$item_info_sev[, grm_vec, drop = FALSE]
    test_info_2pl <- as.numeric(caches$test_info_sus)
    # subset endorsement P(1) columns
    idx_2pl_p1 <- grep(x = colnames(prob_2pl), pattern = ".P.1")
    idx_2pl_p1 <- idx_2pl_p1[idx_2pl_p1 <= n_item]
    prob_2pl_cor <- prob_2pl[, idx_2pl_p1, drop = FALSE]
  } else {
    prob_2pl <- probtrace(sv2, Theta = theta_grid_mat)
    idx_2pl_p1 <- grep(x = colnames(prob_2pl), pattern = ".P.1")
    idx_2pl_p1 <- idx_2pl_p1[idx_2pl_p1 <= n_item]
    prob_2pl_cor <- prob_2pl[, idx_2pl_p1, drop = FALSE]
    
    item_info <- testinfo(sv2, Theta = theta_grid_mat,
                          degrees = c(90, 0), individual = TRUE)[, grm_vec, drop = FALSE]
    test_info_2pl <- as.numeric(testinfo(sv2, Theta = theta_grid_mat,
                                         degrees = c(0, 90), individual = FALSE))
  }
  
  if (nrow(theta_grid_mat) == 1L) {
    item_info <- t(item_info)
    prob_2pl_cor <- t(prob_2pl_cor)
    item_info_hurd <- item_info * prob_2pl_cor
    test_info_sev <- sum(item_info_hurd)
  } else {
    item_info_hurd <- item_info * prob_2pl_cor
    test_info_sev <- rowSums(item_info_hurd)
  }
  
  mat <- cbind(theta_grid_mat, I11 = test_info_2pl, I22 = test_info_sev)
  nm <- apply(mat[, 1:2, drop = FALSE], 1, function(x) paste(x[1], x[2], sep = "_"))
  lst <- lapply(seq_len(nrow(mat)), function(i) diag(c(mat[i, "I11"], mat[i, "I22"])))
  names(lst) <- nm
  list(I_list = lst, I_df = as.data.frame(mat))
}

# A(theta) = I(theta) + Sigma_prior^{-1}
A_mat <- function(mirt_model, Sigma_prior = diag(2), theta) {
  I_out <- ret_test_info_mat(sv2 = mirt_model, theta_grid = matrix(theta, nrow = 1))
  I_mat <- I_out$I_list[[1]]  # 2x2 diagonal
  Q <- solve(Sigma_prior)
  A <- I_mat + Q
  A
}

# Bias function (vector): solve((I + Sigma_prior %*% I(theta))^{-1}) %*% theta
bias_mat <- function(Sigma_prior = diag(2), mirt_model, theta) {
  I_out <- ret_test_info_mat(sv2 = mirt_model, theta_grid = matrix(theta, nrow = 1))
  I_mat <- I_out$I_list[[1]]  # 2x2 diagonal
  solve(diag(2) + Sigma_prior %*% I_mat) %*% theta
}

non_lin_reg = function(mirt_model, Sigma_prior, theta) {
  A_mat = A_mat(mirt_model, Sigma_prior, theta)
  I_out <- ret_test_info_mat(sv2 = mirt_model, theta_grid = matrix(theta, nrow = 1))
  I_mat <- I_out$I_list[[1]]  # 2x2 diagonal
  A_mat_inv = solve(A_mat) ## inv A mat
  A_mat_inv %*% I_mat %*% theta
}

# Gauss–Hermite expectation under Z ~ N(0, Sigma_prior), 2D
gh_expectation_2d <- function(f, Sigma_prior = diag(2), order = 15) {
  gh <- fastGHQuad::gaussHermiteData(order)
  x <- gh$x
  w <- gh$w
  grid <- expand.grid(x1 = x, x2 = x)
  w2 <- as.numeric(outer(w, w))
  
  L <- chol(Sigma_prior)
  Xmat <- as.matrix(grid)             # (order^2 x 2)
  Zmat <- sqrt(2) * (Xmat %*% t(L))   # (order^2 x 2)
  
  vals <- apply(Zmat, 1L, f)
  sum(vals * w2) / pi                 # d=2; divide by pi
}

# rectangle integration
E_g_integral <- function(sv2,
                         Sigma_prior = diag(2),
                         bounds = cbind(c(-5, -5), c(5, 5))) {
  f <- function(z) {
    # bias_mat returns a 2x1 matrix; extract the second component
    g2 <- as.numeric(bias_mat(Sigma_prior = Sigma_prior, mirt_model = sv2, theta = z)[2, 1])
    dens <- mvtnorm::dmvnorm(z, mean = c(0,0), sigma = Sigma_prior)
    g2 * dens
  }
  hcubature(f, lowerLimit = bounds[, 1], upperLimit = bounds[, 2])$integral
}

# rectangle integration
E_g2_integral <- function(sv2,
                          Sigma_prior = diag(2),
                          bounds = cbind(c(-5, -5), c(5, 5))) {
  mean_val <- E_g_integral(sv2, Sigma_prior, bounds= bounds)
  f <- function(z, mean_val) {
    theta <- c(z[1], z[2])
    g <- bias_mat(Sigma_prior = Sigma_prior, mirt_model = sv2, theta = theta)[2,]
    dens <- mvtnorm::dmvnorm(z, mean = c(0,0), sigma = Sigma_prior)
    (g - mean_val)^2 * dens
  }
  hcubature(f, mean_val = mean_val,lowerLimit = bounds[, 1], upperLimit = bounds[, 2])$integral
}

# --- I(theta) using your information function ---
get_I_theta <- function(sv2, theta) {
  # theta is numeric length-2
  I_out <- ret_test_info_mat(sv2 = sv2, theta_grid = matrix(theta, nrow = 1))
  I_out$I_list[[1]]  # 2x2 diag: diag(I11, I22)
}

# --- A(theta), g(theta), and V_e(theta) per Hao ---
A_of_theta <- function(I_theta, Sigma_prior) {
  I_theta + solve(Sigma_prior)
}

g_of_theta <- function(sv2, Sigma_prior, theta) {
  I_theta <- get_I_theta(sv2, theta)
  A_theta <- A_of_theta(I_theta, Sigma_prior)
  # g(theta) = A^{-1} I theta
  solve(A_theta, I_theta %*% matrix(theta, ncol = 1))
}

Ve_of_theta <- function(sv2, Sigma_prior, theta) {
  I_theta <- get_I_theta(sv2, theta)
  A_theta <- A_of_theta(I_theta, Sigma_prior)
  Ainv <- solve(A_theta)
  # V_e(theta) = A^{-1} I A^{-1}
  Ainv %*% I_theta %*% Ainv
}

# --- Reliability via GH using info and GH function ---
severity_reliability_GH <- function(sv2,
                                    Sigma_prior = diag(2),  # scoring prior inside A(theta)
                                    order = 15) {
  # g_sev(theta): severity component of g(theta)
  g_sev_fun <- function(theta_vec) {
    as.numeric(g_of_theta(sv2, Sigma_prior, theta_vec)[2, 1])
  }
  # v_sev(theta): conditional error variance for severity (2,2 element of V_e)
  v_sev_fun <- function(theta_vec) {
    as.numeric(Ve_of_theta(sv2, Sigma_prior, theta_vec)[2, 2])
  }
  
  Eg    <- gh_expectation_2d(g_sev_fun, Sigma_prior, order)                    # E[g_sev]
  Vtrue <- gh_expectation_2d(function(z) (g_sev_fun(z) - Eg)^2, Sigma_prior, order)  # Var(g_sev)
  Eerr  <- gh_expectation_2d(v_sev_fun, Sigma_prior, order)                    # E[V_e(θ)[2,2]]
  
  rho <- Vtrue / (Vtrue + Eerr)
  list(rho = rho, V_true = Vtrue, E_err = Eerr, Eg = Eg)
}

# Build mirt model string once
make_model <- function(n_2pl, n_item_total) {
  sprintf("F1 = 1-%d\nF2 = %d-%d\nCOV = F1*F2", n_2pl, n_2pl + 1L, n_item_total)
}

# Convert traditional parameters to mirt (vectorized wrapper)
to_mirt <- function(a_grm, b_grm, a_2pl, b_2pl, n_cat) {
  pl2_conv <- apply(cbind(a_2pl, b_2pl, 0, 0), 1,
                    function(x) traditional2mirt(x, cls = "2PL", ncat = 2))
  grm_conv <- apply(cbind(a_grm, b_grm), 1,
                    function(x) traditional2mirt(x, cls = "graded", ncat = n_cat))
  a_z_m <- pl2_conv["a1", ]
  d_z_m <- pl2_conv["d", ]
  a_m   <- grm_conv["a1", ]
  d_m   <- t(apply(grm_conv[-1, , drop = FALSE], 2, function(x) sort(x, decreasing = TRUE)))
  list(a_m = a_m, d_m = d_m, a_z_m = a_z_m, d_z_m = d_z_m)
}

# Unidimensional GRM reliability (better empirically): 1 − mean(SE^2) / var(scores)
# Uses EAP or MAP SE if available; falls back to testinfo-based approximation if SE not available.
grm_reliability_emp <- function(mod, method = "EAP") {
  fs <- try(fscores(mod, method = method, full.scores = TRUE), silent = TRUE)
  se_attr <- try(attr(fs, "SE"), silent = TRUE)
  if (!inherits(fs, "try-error") && !inherits(se_attr, "try-error") && !is.null(se_attr)) {
    scores <- as.numeric(fs)
    ses    <- as.numeric(se_attr)
    return(1 - mean(ses^2, na.rm = TRUE) / var(scores, na.rm = TRUE))
  } else {
    # Fallback: integrate 1 / I(theta) (MLE-style approx). Not ideal for MAP.
    Theta <- seq(-6, 6, by = 0.1)
    info  <- as.numeric(testinfo(mod, Theta = Theta))
    return(1 / (1 + weighted.mean(1 / info, dnorm(Theta))))
  }
}

## Function to estimate reliability from non linear regression
severity_reliability <- function(sv2, Sigma_prior = diag(2), bounds = cbind(c(-5, -5), c(5, 5))) {
  Eg  <- E_g_integral(sv2, Sigma_prior, bounds)
  Eg2 <- E_g2_integral(sv2, Sigma_prior, bounds)
  Rel = 1 - Eg2
  list(Eg = Eg, Eg2= Eg2, Rel = Rel)
}
