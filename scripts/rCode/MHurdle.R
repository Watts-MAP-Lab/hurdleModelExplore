library(mvtnorm)
library(plyr)

data_out <- read.csv("./data/ScaleWithOutcomesOSF.csv")
names(data_out)

response <- data_out[ , c(2:15)] # Select 14 item responses

tabs_all <- apply(response, 2, table)
tabs_all <- apply(tabs_all, 2, prop.table)
tabs_all

patterns <- ddply(response, ~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14, summarize, n=length(x1)) # Get frequency of each response pattern
patterns <- ddply(response, ~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14, summarize, n=length(x1)) # Get frequency of each response pattern

tmp <- as.formula(paste("~", paste(colnames(response), collapse = "+")))
patterns <- ddply(response, tmp, summarize, n=length(x1)) # Get frequency of each response pattern


patterns <- as.data.frame(patterns)
patterns <- patterns[order(patterns$x1, patterns$x2, patterns$x3, patterns$x4, patterns$x5, patterns$x6, patterns$x7,
                           patterns$x8, patterns$x9, patterns$x10, patterns$x11, patterns$x12, patterns$x13, patterns$x14),]
colnames(patterns) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "r")
x1 <- patterns$x1
x2 <- patterns$x2
x3 <- patterns$x3
x4 <- patterns$x4
x5 <- patterns$x5
x6 <- patterns$x6
x7 <- patterns$x7
x8 <- patterns$x8
x9 <- patterns$x9
x10 <- patterns$x10
x11 <- patterns$x11
x12 <- patterns$x12
x13 <- patterns$x13
x14 <- patterns$x14

x <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14)
r <- patterns$r

nitems <- ncol(x)
n <- sum(r)
data <- data.frame(I(x), r) # Dataframe with all item responses and frequencies

theta1 <- seq(-6,6,0.25) # Quadrature points for theta1
theta2 <- seq(-6,6,0.25) # Quadrature points for theta1
theta <- expand.grid(theta1, theta2)
names(theta) <- c("theta1", "theta2")

###################################
#### Functions for trace lines ####
###################################

nitemsgrm <- nitems
ncatgrm <- max(x)
ncat2PL <- 2
mat <- matrix(0, nrow=nitems, ncol=length(theta$theta1))
itemtraceGRM <- lapply(seq_len(ncatgrm), function(X) mat)
itemtrace2PL <- lapply(seq_len(ncat2PL), function(X) mat)
itemtrace <- lapply(seq_len(ncatgrm+1), function(X) mat)

# 2PL trace lines

trace.line.pts.2PL <- function(a_z, b_z, theta)	{
  for (j in 1:nitems){
    for (k in 0:3){
      if (k == 0){
        itemtrace2PL[[1]][j,] <- 1 - exp(a_z[j,]*(theta[,1] - b_z[j,1]))/(1+exp(a_z[j,]*(theta[,1] - b_z[j,1])))
      }
      else{
        itemtrace2PL[[2]][j,] <- exp(a_z[j,]*(theta[,1] - b_z[j,1]))/(1+exp(a_z[j,]*(theta[,1] - b_z[j,1])))
      }
    }
  }
  return(itemtrace2PL)
}

# GRM trace lines
trace.line.pts.grm <- function(a, b, theta)	{
  for (j in 1:nitems){
    for (k in 1:3){
      if (k == 1){
        itemtraceGRM[[1]][j,] <- 1 - exp(a[j,]*(theta[,2] - b[j,k]))/(1+exp(a[j,]*(theta[,2] - b[j,k])))
      }
      if (k == 2){
        itemtraceGRM[[2]][j,] <- exp(a[j,]*(theta[,2] - b[j,k-1]))/(1+exp(a[j,]*(theta[,2] - b[j,k-1]))) - exp(a[j,]*(theta[,2] - b[j,k]))/(1+exp(a[j,]*(theta[,2] - b[j,k])))
      }
      if (k == 3){
        itemtraceGRM[[3]][j,] <- exp(a[j,]*(theta[,2] - b[j,k-1]))/(1+exp(a[j,]*(theta[,2] - b[j,k-1])))
      }
    }
  }
  return(itemtraceGRM)
}

trace.line.pts <- function(a, b, a_z, b_z, theta){
  for (j in 1:nitems){
    for (k in 0:3){
      if (k == 0){
        itemtrace[[k+1]][j,] <- (1-trace.line.pts.2PL(a_z, b_z, theta)[[2]][j,])
      }
      if (k == 1){
        itemtrace[[k+1]][j,] <- (trace.line.pts.2PL(a_z, b_z, theta)[[2]][j,])*trace.line.pts.grm(a, b, theta)[[1]][j,]
      }
      if (k == 2){
        itemtrace[[k+1]][j,] <- (trace.line.pts.2PL(a_z, b_z, theta)[[2]][j,])*trace.line.pts.grm(a, b, theta)[[2]][j,]
      }
      if (k == 3){
        itemtrace[[k+1]][j,] <- (trace.line.pts.2PL(a_z, b_z, theta)[[2]][j,])*trace.line.pts.grm(a, b, theta)[[3]][j,]
      }
    }
  }
  return(itemtrace)
}

#############################
#### Likelihood function ####
#############################

ll.grm.ip <- function(p,testdata,theta) {
  nParmsPerItemGRM <- ncatgrm
  nParmsPerItem2PL <- ncat2PL
  a <- matrix(c(rep(-1, nitems)), nitems, 1)
  b <- matrix(c(rep(-1, nitems*(ncatgrm-1))), nitems, ncatgrm-1)
  a_z <- matrix(c(rep(-1, nitems)), nitems, 1)
  b_z <- matrix(c(rep(-1, nitems)), nitems, 1)
  rho_exp <- 0.5
  for (j in 1:nitemsgrm) {
    a[j,] <- p[(j-1)*nParmsPerItemGRM + 1]
    a_z[j,] <- p[nitems*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 1]
    b_z[j,] <- p[nitems*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 2]
    for (k in 1:(ncatgrm-1)){
      b[j,k] <- p[(j-1)*nParmsPerItemGRM + 1 + k]
    }
  }
  rho_exp <- p[nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL + 1]
  rho <- exp(rho_exp)/(1+exp(rho_exp))
  itemtrace <- trace.line.pts(a, b, a_z, b_z, theta)
  expected <- rep(0,length(testdata$r))
  for (i in 1:length(testdata$r)) {
    posterior <- dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
    for (item in 1:ncol(testdata$x)) {
      x <- I(testdata$x[i,item])
      posterior <- posterior*itemtrace[[x+1]][item,] 
    }
    expected[[i]] <- sum(posterior) 
  }
  l <- (-1)*(sum(testdata$r*log(expected)))
  print(l)
  return(l)
}

nParmsPerItemGRM <- ncatgrm
nParmsPerItem2PL <- ncat2PL

a <- matrix(1, nitems, 1)
a_z <- matrix(1, nitems, 1)
b <- matrix(c(1,2), nitems, ncatgrm-1, byrow=T)
b_z <- matrix(-1, nitems, 1, byrow=T)
rho_exp <- 0.4
p <- rep(0, nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL + 1)

for (j in 1:nitems) {
  p[(j-1)*nParmsPerItemGRM + 1] <- a[j,]
  p[nitemsgrm*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 1] <- a_z[j,]
  p[nitemsgrm*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 2] <- b_z[j,]
  for (k in 1:(ncatgrm-1)){
    p[(j-1)*nParmsPerItemGRM + 1 + k] <- b[j,k]
  }
}
p[nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL + 1] <- rho_exp

system.time(ll.grm.ip(p, data, theta))

# Starting values
system.time(MHurdle <- optim(par = p, fn = ll.grm.ip, testdata = data, hessian = FALSE,
                 theta = theta, control = list(trace = 5, maxit = 350,000, abstol = 0.001), method = "L-BFGS-B"))

# Scale Scores

Dep <- MHurdle$par
Dep <- matrix(Dep, nrow = length(p), ncol = 1)

a_grm <- Dep[seq(1, nitems*nParmsPerItemGRM, nParmsPerItemGRM)] # Select a parameters from GRM
a_grm <- as.data.frame(a_grm, nrow = nitems)

b_grm <- matrix(0, nitems, ncatgrm - 1) # Select b parameters from GRM
for (i in 1:nitems){
  b_grm[i,] <- Dep[c(ncatgrm*i - 1, ncatgrm*i),]
}
b_grm <- as.data.frame(b_grm, nrow = nitems, ncol = ncatgrm - 1)

a_2PL <- Dep[seq(nitems*nParmsPerItemGRM + 1, nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL, 2)] # Select a parameters from 2PL
a_2PL <- as.data.frame(a_2PL, nrow = nitems)

b_2PL <- Dep[seq(nitems*nParmsPerItemGRM + 2, nitems*nParmsPerItemGRM + nitems*nParmsPerItem2PL, 2)] # Select b parameters from 2PL
b_2PL <- as.data.frame(b_2PL, nrow = nitems)

rho <- exp(Dep[c(length(p)),])/(1+exp(Dep[c(length(p)),]))

itemtrace2PL <- trace.line.pts.2PL(a_2PL, b_2PL, theta)
itemtraceGRM <- trace.line.pts.grm(a_grm, b_grm, theta)
itemtrace <- trace.line.pts(a_grm, b_grm, a_2PL, b_2PL, theta)

qpoints <- theta
prior <- dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))

score <- function(pattern) {
  lhood <- rep(1,length(qpoints$theta1))
  for (item in 1:nitems){
    if(pattern[item]==0) lhood <- lhood*itemtrace[[1]][item,]
    if(pattern[item]==1) lhood <- lhood*itemtrace[[2]][item,]
    if(pattern[item]==2) lhood <- lhood*itemtrace[[3]][item,]
    if(pattern[item]==3) lhood <- lhood*itemtrace[[4]][item,]
  }
  lhood
}

for(respondent in 1:nrow(data_out)) {
  pattern <- response[respondent,grep("x[[:digit:]]", names(response))]
  lhood <- score(pattern)
  
  eap2PL_Hurdle <- sum(lhood*prior*qpoints$theta1)/sum(lhood*prior)
  se2PL_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta1-eap2PL_Hurdle)^2)/sum(lhood*prior))
  data_out$eap2PL_Hurdle[respondent] <- eap2PL_Hurdle
  data_out$se2PL_Hurdle[respondent] <- se2PL_Hurdle
  
  eapGRM_Hurdle <- sum(lhood*prior*qpoints$theta2)/sum(lhood*prior)
  seGRM_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta2-eapGRM_Hurdle)^2)/sum(lhood*prior))
  data_out$eapGRM_Hurdle[respondent] <- eapGRM_Hurdle
  data_out$seGRM_Hurdle[respondent] <- seGRM_Hurdle
}

# GLM for sample outcome variable

mod <- glm(data_out$D_MDE30 ~ data_out$eap2PL_Hurdle + data_out$eapGRM_Hurdle, family = "binomial")
summary(mod)
exp(coef(mod))
