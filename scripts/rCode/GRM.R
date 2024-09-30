library(mvtnorm)
library(plyr)

data_out <- read.csv("./data/ScaleWithOutcomesOSF.csv")
names(data_out)

response <- data_out[ , c(2:15)] # Select 14 item responses

tabs_all <- apply(response, 2, table)
tabs_all <- apply(tabs_all, 2, prop.table)
tabs_all

patterns <- ddply(response, ~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14, summarize, n=length(x1)) # Get frequency of each response pattern
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

theta <- seq(-6,6,0.25) # Quadrature points for theta
theta

###################################
#### Functions for trace lines ####
###################################

nitemsgrm <- nitems
ncatgrm <- max(x) + 1
mat <- matrix(0, nrow=nitems, ncol=length(theta))
itemtrace <- lapply(seq_len(ncatgrm+1), function(X) mat)


# GRM trace lines

trace.line.pts.grm <- function(a, b, theta)	{
  for (j in 1:nitems){
    for (k in 0:(ncatgrm-1)){
      if (k == 0){
        itemtrace[[k+1]][j,] <- 1 - exp(a[j,]*(theta - b[j,k+1]))/(1+exp(a[j,]*(theta - b[j,k+1])))
      }
      if (k == 1){
        itemtrace[[k+1]][j,] <- exp(a[j,]*(theta - b[j,k]))/(1+exp(a[j,]*(theta - b[j,k]))) - exp(a[j,]*(theta - b[j,k+1]))/(1+exp(a[j,]*(theta - b[j,k+1])))
      }
      if (k == 2){
        itemtrace[[k+1]][j,] <- exp(a[j,]*(theta - b[j,k]))/(1+exp(a[j,]*(theta - b[j,k]))) - exp(a[j,]*(theta - b[j,k+1]))/(1+exp(a[j,]*(theta - b[j,k+1])))
      }
      if (k == 3){
        itemtrace[[k+1]][j,] <- exp(a[j,]*(theta - b[j,k]))/(1+exp(a[j,]*(theta - b[j,k])))
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
  a <- matrix(c(rep(-1, nitems)), nitems, 1)
  b <- matrix(c(rep(-1, nitems*(ncatgrm-1))), nitems, ncatgrm-1)
  for (j in 1:nitems) {
    a[j,] <- p[(j-1)*nParmsPerItemGRM + 1]
    for (k in 1:(ncatgrm-1)){
      b[j,k] <- p[(j-1)*nParmsPerItemGRM + 1 + k]
    }
  }
  itemtrace <- trace.line.pts.grm(a, b, theta)
  expected <- rep(0,length(testdata$r))
  for (i in 1:length(testdata$r)) {
    posterior <- dnorm(theta, 0, 1)
    for (item in 1:ncol(testdata$x)) {
      x <- I(testdata$x[i,item])   
      posterior <- posterior*itemtrace[[x+1]][item,] 
    }
    expected[[i]] <- sum(posterior) 
  }
  l <- (-1)*(sum(testdata$r*log(expected)))
}

# Starting values

nParmsPerItem <- ncatgrm

a <- matrix(1, nitems, 1)
b <- matrix(c(1,2,3), nitems, ncatgrm-1, byrow=T)
p <- rep(0, nitems*nParmsPerItem)

for (j in 1:nitemsgrm) {
  p[(j-1)*nParmsPerItem + 1] <- a[j,]
  for (k in 1:(ncatgrm-1)){
    p[(j-1)*nParmsPerItem + 1 + k] <- b[j,k]
  }
}

GRM <- optim(par = p, fn = ll.grm.ip, testdata = data, hessian = TRUE,
                 theta = theta, control = list(trace = TRUE, maxit = 5, abstol = 0.001), method = "BFGS")

# Scale Scores

Dep <- GRM$par
Dep <- matrix(Dep, nrow = length(p), ncol = 1)

a_grm <- Dep[seq(1, nitems*nParmsPerItemGRM, nParmsPerItemGRM)] # Select a parameters from GRM
a_grm <- as.data.frame(a_grm, nrow = nitems)

b_grm <- matrix(0, nitems, ncatgrm - 1) # Select b parameters from GRM
for (i in 1:nitems){
  b_grm[i,] <- Dep[c(ncatgrm*i - 2, ncatgrm*i - 1, ncatgrm*i),]
}
b_grm <- as.data.frame(b_grm, nrow = nitems, ncol = ncatgrm - 1)

itemtrace <- trace.line.pts.grm(a_grm, b_grm, theta)

qpoints <- theta
qpoints <- data.frame(theta)
names(qpoints) <- "theta"
prior <- dnorm(theta, 0, 1)

score <- function(pattern) {
  lhood <- rep(1,length(qpoints$theta))
  for (item in 1:nitems){
    if(pattern[item]==0) lhood <- lhood*itemtrace[[1]][item,]
    if(pattern[item]==1) lhood <- lhood*itemtrace[[2]][item,]
    if(pattern[item]==2) lhood <- lhood*itemtrace[[3]][item,]
    if(pattern[item]==3) lhood <- lhood*itemtrace[[4]][item,]
  }
  lhood
}

for(respondent in 1:nrow(data_out)) {
  pattern <- response[respondent, grep("x[[:digit:]]", names(response))]
  lhood <- score(pattern)
  eapGRM <- sum(lhood*prior*qpoints$theta)/sum(lhood*prior)
  seGRM <- sqrt(sum(lhood*prior*(qpoints$theta-eapGRM)^2)/sum(lhood*prior))
  data_out$eapGRM[respondent] <- eapGRM
  data_out$seGRM[respondent] <- seGRM
}

# GLM for sample outcome variable

mod <- glm(data_out$D_MDE30 ~ data_out$eapGRM, family = "binomial")
summary(mod)
exp(coef(mod))