## Rewrite all of this wround MIRT
library(mirt)
## Create a function which will return MIRT models, simulated responses & hurdle information functions
## when given all hurdle parameters
sim_mirt_hurdle <- function(a, b, a_z, b_z, rho, THETA = NULL, 
                            N = 15000, theta_grid = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)),
                            estModel = FALSE) {
  
  require(mirt)
  require(dplyr)
  
  # Generate bivariate normal theta values if not provided
  # THETA represents latent traits: column 1 = susceptibility, column 2 = severity
  if (!is.null(N) & is.null(THETA)) {
    muVals <- rep(0, 2)  # Mean of 0 for both dimensions
    varCovMat <- matrix(c(1, rho, rho, 1), ncol = 2)  # Correlation matrix
    THETA <- MASS::mvrnorm(N, mu = muVals, Sigma = varCovMat)
  }
  
  # Set up item parameter matrices for mirt
  # Total items = 2 * number of hurdle items (2PL + corresponding graded response)
  n_item <- length(a_z) * 2
  n_diff <- ncol(b) + 1  # Number of difficulty parameters per item
  
  # Initialize discrimination parameter matrix (items x factors)
  a_mat <- matrix(0, nrow = n_item, ncol = 2)
  n_2pl <- n_item / 2  # Number of 2PL (hurdle) items
  grm_vec <- seq(n_2pl + 1, n_item, 1)  # Indices for graded response items
  
  # Assign discrimination parameters
  a_mat[1:n_2pl, 1] <- a_z        # 2PL items load on Factor 1 (susceptibility)
  a_mat[grm_vec, 2] <- a          # Graded items load on Factor 2 (severity)
  
  # Initialize difficulty parameter matrix
  b_mat <- matrix(NA, nrow = n_item, ncol = n_diff)
  b_mat[1:n_2pl, 1] <- b_z        # 2PL difficulty parameters
  b_mat[grm_vec, 2:ncol(b_mat)] <- b  # Graded response difficulty parameters
  
  # Define item types for mirt
  item_type <- c(rep("2PL", n_2pl), rep("graded", n_2pl))
  
  # Simulate item responses
  sim_vals <- simdata(a = a_mat, d = b_mat, sigma = varCovMat, 
                      Theta = THETA, itemtype = item_type)
  theta_grid_mat <- as.matrix(theta_grid)
  
  # Apply hurdle mechanism: set graded responses to NA when 2PL item = 0
  # This creates the hurdle structure where you must "pass" the 2PL to answer graded items
  pl2_dat <- sim_vals[, 1:n_2pl]           # 2PL responses (hurdle component)
  grm_dat <- sim_vals[, grm_vec]           # Graded responses (severity component)
  grm_dat[pl2_dat == 0] <- NA              # Missing if didn't endorse hurdle item
  sim_vals_w_na <- bind_cols(pl2_dat, grm_dat)
  
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
  
  # Estimate mirt model if requested
  if (estModel) {
    sv1 <- mirt(sim_vals_w_na, model = model, itemtype = item_type)
    fscores_estE <- mirt::fscores(sv1, method = "EAP")  # Expected A Posteriori scores

    # Calculate trace lines and information for estimated model
    prob_2pl <- probtrace(sv1, Theta = theta_grid_mat)
    
    # Extract probability of endorsement for 2PL items (P.1 = prob of response = 1)
    prob_2pl_cor <- grep(x = colnames(prob_2pl), pattern = ".P.1")
    prob_2pl_cor <- prob_2pl_cor[which(prob_2pl_cor <= n_item)]
    prob_2pl_cor <- prob_2pl[, prob_2pl_cor]
    
    # Get item information for graded items with respect to Factor 2 (severity)
    # degrees = c(90,0) means 90° for F1, 0° for F2 (i.e., information wrt F2 only)
    item_info <- testinfo(sv1, Theta = theta_grid_mat, degrees = c(90, 0), 
                          individual = TRUE)[, grm_vec]
    
    # Adjust item information by probability of endorsing hurdle item
    # This gives "effective" information accounting for hurdle mechanism
    item_info_hurd <- item_info
    for (i in 1:ncol(item_info_hurd)) {
      item_info_hurd[, i] <- item_info[, i] * prob_2pl_cor[, i]
    }
    
    # Store estimated model results
    prob_2pl_est <- prob_2pl
    item_info_est <- item_info
    item_info_hurd_est <- item_info_hurd
  }
  
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
  pars_df$value[pars_df$name == "d"] <- b_z
  pars_df$est[pars_df$name == "d"] <- FALSE
  
  # Graded response difficulty parameters
  for (i in 1:ncol(b)) {
    d_val <- paste("d", i, sep = '')
    pars_df$value[pars_df$name == d_val] <- b[, i]
    pars_df$est[pars_df$name == d_val] <- FALSE
  }
  
  # Factor correlation
  pars_df$value[pars_df$name == "COV_21"] <- rho
  pars_df$est[pars_df$name == "COV_21"] <- FALSE
  
  # Fit model with fixed parameters
  sv2 <- mirt(sim_vals_w_na, model = model, itemtype = item_type, pars = pars_df)
  
  # Get factor scores from true parameter model
  fscores_truE <- mirt::fscores(sv2, method = "EAP")

  # Repeat trace line and information calculations for true model
  prob_2pl <- probtrace(sv2, Theta = theta_grid_mat)
  
  # Extract 2PL endorsement probabilities
  prob_2pl_cor <- grep(x = colnames(prob_2pl), pattern = ".P.1")
  prob_2pl_cor <- prob_2pl_cor[which(prob_2pl_cor <= n_item)]
  prob_2pl_cor <- prob_2pl[, prob_2pl_cor]
  
  # Item information for graded items (Factor 2)
  item_info <- testinfo(sv2, Theta = theta_grid_mat, degrees = c(90, 0), 
                        individual = TRUE)[, grm_vec]
  
  # Hurdle-adjusted item information
  item_info_hurd <- item_info
  for (i in 1:ncol(item_info_hurd)) {
    item_info_hurd[, i] <- item_info[, i] * prob_2pl_cor[, i]
  }
  
  # Store true model results
  prob_2pl_tru <- prob_2pl
  item_info_tru <- item_info
  item_info_hurd_tru <- item_info_hurd
  test_info_hurd_tru <- rowSums(item_info_hurd)  # Total test information
  
  # Organize output based on whether model was estimated
  if (estModel) {
    theta_vals <- data.frame(THETA, fscores_estE, fscores_truE)
    colnames(theta_vals) <- c("trueSus", "trueSev", "estP_susEAP", "estP_sevEAP", 
                              "truP_susEAP", "truP_sevEAP")
    
    out_vals <- list(
      theta_vals = theta_vals, 
      est_mod = sv1, 
      tru_mod = sv2,
      est_hurd_trace = item_info_hurd_est, 
      tru_hurd_trace = item_info_hurd_tru,
      prob_endorse_est = prob_2pl_est, 
      prob_endorse_tru = prob_2pl_tru,
      grm_est_trace = item_info_est, 
      grm_tru_trace = item_info_tru,
      responses = sim_vals, 
      responses_na = sim_vals_w_na
    )
  } else {
    # Only true model results
    theta_vals <- data.frame(THETA, fscores_truE)
    colnames(theta_vals) <- c("trueSus", "trueSev", "truP_susEAP", "truP_sevEAP")
    
    # Calculate marginal reliability using test information
    muVals <- rep(0, 2)  # Should be defined earlier for consistency
    varCovMat <- matrix(c(1, rho, rho, 1), ncol = 2)  # Should reuse earlier definition
    weights <- mvtnorm::dmvnorm(theta_grid, muVals, varCovMat)
    test_info_hurd_tru_inv <- test_info_hurd_tru^-1
    out.rel <- 1 / (1 + weighted.mean(test_info_hurd_tru_inv, weights))
    
    out_vals <- list(
      theta_vals = theta_vals, 
      tru_mod = sv2,
      tru_hurd_test_info = test_info_hurd_tru,
      tru_hurd_trace = item_info_hurd_tru,
      prob_endorse_tru = prob_2pl_tru,
      grm_tru_trace = item_info_tru,
      responses = sim_vals, 
      responses_na = sim_vals_w_na, 
      out_rel = out.rel
    )
  }
  
  return(out_vals)
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

## Now run this function
## First generate all of the inputs
n_items = 14
a = runif(n_items, min = .8, max = 2)
b = genDiffGRM(num_items = n_items, num_categories = 3, min = -2, max = 1)
a_z = runif(n_items, min = .8, max = 2)
b_z = runif(n_items, min = 0, max = 2)
N = 15000
muVals = rep(0, 2)
rho = .4
varCovMat = matrix(c(1,rho,rho,1), ncol=2)
THETA = MASS::mvrnorm(N, mu = muVals, Sigma = varCovMat)
## Convert from trad to mirt here
pl2_conv <- apply(cbind(a_z, b_z, 0,0), 1, function(x)traditional2mirt(x, cls="2PL", ncat=2) )
grm_conv <- apply(cbind(a, b), 1, function(x)traditional2mirt(x, cls="graded", ncat=3) )
## Now reassign
a_z <- pl2_conv["a1",]
b_z <- pl2_conv["d",]
a <- grm_conv["a1",]
b <- cbind(grm_conv[3,], grm_conv[2,])
sim_vals = sim_mirt_hurdle(a = a, b = b, a_z = a_z, b_z = b_z, rho = .4, THETA=THETA)
var(sim_vals$theta_vals)
true.score.var <- var(sim_vals$theta_vals$trueSev)
error.var <- var(sim_vals$theta_vals$truP_sevEAP - sim_vals$theta_vals$trueSev)
trueRel <- true.score.var / (true.score.var + error.var)
sim_vals$out_rel

## Now examine estimated test rel versus the 2PL difficulty range
out.plot <- data.frame(twoPL = seq(-4, 1, 1), out.rel = NA)
for(i in 1:nrow(out.plot)){
  a_z <- rep(3, n_items)
  b_z <- runif(n_items, min=out.plot$twoPL[i], max = out.plot$twoPL[i]+2) ## Generate n.item 2PL discrimination parameters randomly sampled between -2 & 0
  pl2_conv <- apply(cbind(a_z, b_z, 0,0), 1, function(x)traditional2mirt(x, cls="2PL", ncat=2) )
  a_z <- pl2_conv["a1",]
  b_z <- pl2_conv["d",]
  out.plot$out.rel[i] = sim_mirt_hurdle(a = a, b = b, a_z = a_z, b_z = b_z, rho = rho)$out_rel
}
plot(out.plot)

## Now generate these responses for mean theta values
sus.vals <- c(-1, 0, 1)
sev.vals <- c(-1, 0, 1)
all.iter <- expand.grid(sus.vals, sev.vals)
all.iter <- data.frame(all.iter)
all.iter$simVar <- NA
all.iter$testInfoVar <- NA
b_z <- runif(n = n_items, min= 0, max = 2)
a_z <- rep(3, n_items)
b_z <- runif(n_items, min=-1, max = 1) ## Generate n.item 2PL discrimination parameters randomly sampled between -2 & 0
pl2_conv <- apply(cbind(a_z, b_z, 0,0), 1, function(x)traditional2mirt(x, cls="2PL", ncat=2) )
a_z <- pl2_conv["a1",]
b_z <- pl2_conv["d",]
theta_grid = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1))
for(i in 1:nrow(all.iter)){
  THETA = matrix(c(all.iter[i,1],all.iter[i,2]), nrow=300000, ncol = 2, byrow = TRUE)
  ## Now generate the response values
  sim_vals = sim_mirt_hurdle(a = a, b = b, a_z = a_z, b_z = b_z, rho = .4, THETA=THETA, theta_grid = theta_grid)
  ident_var = which(theta_grid[,1]==all.iter[i,1] & theta_grid[,2]==all.iter[i,2])
  ## Examine simulated variance
  all.iter$simVar[i] <- var(sim_vals$theta_vals)[4,4]
  ## Now obtain variance from test information function
  all.iter$testInfoVar[i] <- sim_vals$tru_hurd_test_info[ident_var]^-1
}
## Now plot this
library(ggplot2)
for.plot <- reshape2::melt(all.iter, id.vars=c("Var1", "Var2"))
ggplot(for.plot, aes(x=variable, y = value)) +
  geom_bar(stat="identity", position="dodge") +
  facet_grid(Var1 ~ Var2) +
  theme_bw() +
  xlab("Var Source") +
  ylab("Variance")
