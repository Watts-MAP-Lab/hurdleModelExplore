## This file will be used to try to obtain test retest using baysesian methods
##### --load-library-------------
source("./scripts/rCode/hurdleFunctions.r")
suppressPackageStartupMessages(library("tidyverse"))
library("numDeriv")
suppressPackageStartupMessages(library(MplusAutomation))
library("mirt")
library("psych")
library("mgcv")
library("rstan")

##### --declare-sim-params-------
## Sim params will need to be modified at a later time point
source("./scripts/rCode/simParam.r")
seedVal <- as.integer(commandArgs(TRUE))
seedVal <- 129
set.seed(seedVal)

## --run-test-sim ------------
## I will run through a single test sim here and see how I can best, and most efficiently run this
## I really need to focus on the memory leakage problem the R has, especially when running large loops
## First create the data -- this will start with the maximum full dataset, 9 total response categories, full range of difficulty parameters
## This will also showcase where I need to streamline code with custom functions
## Run a spread check for discrim values
add.val.2pl <- 1.5
add.val.grm <- 1.5
a = runif(n = all.sim.vals$nItem[seedVal], min = all.sim.vals$grmDiscrim[seedVal], all.sim.vals$grmDiscrim[seedVal] + add.val.grm)
b = genDiffGRM(num_items = all.sim.vals$nItem[seedVal], num_categories = all.sim.vals$nCat[seedVal], min = all.sim.vals$difGrmF[seedVal], max = all.sim.vals$difGrmF[seedVal]+2.5, rnorm_var = .3)
a_z = runif(n = all.sim.vals$nItem[seedVal], min = all.sim.vals$discrim2pl[seedVal], all.sim.vals$discrim2pl[seedVal] + add.val.2pl)
## Need to generate 4 separate b_z levels
b_z1 = runif(all.sim.vals$nItem[seedVal], min = all.sim.vals$dif2PL[seedVal], max=all.sim.vals$dif2PL[seedVal]+2)
muVals = c(0,0)
rho <- all.sim.vals$facCor[seedVal]
varCovMat = matrix(c(1,rho,rho,1), ncol=2)
N = all.sim.vals$n[seedVal]
## Now generate theta here
#theta = MASS::mvrnorm(n = N, mu = muVals, Sigma = varCovMat)
theta = MASS::mvrnorm(n = N, mu = muVals, Sigma = varCovMat)
#reps1 = system.time(simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta))
#reps1f = system.time(simulate_hurdle_responses_fast(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta, qpoints = expand.grid(seq(-6, 6, .1), seq(-6, 6, .1))))
reps1 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta)
for.stan.2pl <- reps1$mplusMat[, 1:7]
## Now estimate the stan model
stan_data = list(
   n = nrow(for.stan.2pl),
   p = 7,
   Y = for.stan.2pl
 )

#stan_est = stan(file = "./scripts/stanCode/stanHurdle2PL.stan", data = stan_data, chains = 2, cores = 2, warmup = 1000, iter = 2000)
for.stan.grm <- reps1$mplusMat[, 8:14]
## Now estimate the stan model
# int<lower=2, upper=8> K; //number of categories
# int <lower=0> n_student;
# int <lower=1> n_item;
# int<lower=1,upper=K> Y[n_student,n_item];
stan_data = list(
  K = 3,
  n_student = nrow(for.stan.grm),
  n_item = 7,
  Y = for.stan.grm
)

#stan_grm = stan(file = "./scripts/stanCode/stanHurdleGRM.stan", data = stan_data, chains = 2, cores = 2, warmup = 30, iter = 60)


## NOw organize the data fro stan
prep_grm <- reps1$mplusMat[, 8:14]
prep_grm$id = 1:nrow(prep_grm)
## Now melt
prep_grm <- reshape2::melt(prep_grm, id.var="id")
## rm na values
prep_grm <- prep_grm[complete.cases(prep_grm),]
prep_grm$variable <- as.numeric(factor(prep_grm$variable))
prep_grm$id <- as.numeric(factor(prep_grm$id))

## Now prep the data
# int<lower=0> N; // number of subjects
# int<lower=0> p; // number of items -- equivalent for 2PL and GRM portions
# int y[n,p]; // response vector -- binary
# // Now do the GRM data down here -- responses should be vectorized to accomodate missing data
# int k; // number of grm categories excluding 0
# int<lower=1> NG; // number of grm responses
# int<lower=1,upper=NH> jj[NG];  // respondent for observation n
# int<lower=1,upper=NH> kk[NG];  // question for observation n
# int<lower=0,upper=k> y[NG];   // correctness for observation n
stan_data = list(
  N = nrow(for.stan.2pl),
  p = 7,
  Y = for.stan.2pl,
  k = 3,
  NG = nrow(prep_grm),
  jj = prep_grm$id,
  kk = prep_grm$variable,
  y = prep_grm$value
)
stan_hurd = stan(file = "./scripts/stanCode/stanHurdleHURDLE.stan", data = stan_data, chains = 2, cores = 2, warmup = 1000, iter = 2000)

## Compare the estimated theta values
sum.vals <- summary(stan_hurd)
iso.theta = sum.vals$summary[1:15000,"mean"]
cor(reps1$theta$theta1, iso.theta)
iso.theta2 = sum.vals$summary[15001:30000,"mean"]
cor(reps1$theta$theta2, iso.theta2)
cor(reps1$theta$eapSev, iso.theta2)
plot(reps1$theta$theta2, iso.theta2)
plot(reps1$theta$eapSev, iso.theta2)

iso.2pl.discrim <- sum.vals$summary[grep(pattern = "alpha_pl", x = rownames(sum.vals$summary))[1:7],"mean"]
plot(iso.2pl.discrim, reps1$a_z)
cor(iso.2pl.discrim, reps1$a_z)
plot(iso.2pl.discrim-reps1$a_z)
iso.grm.discrim <- sum.vals$summary[grep(pattern = "alpha_grm", x = rownames(sum.vals$summary))[1:7],"mean"]
plot(iso.grm.discrim, reps1$a)
cor(iso.grm.discrim, reps1$a)
plot(iso.grm.discrim-reps1$a)
iso.grm.diff <- sum.vals$summary[grep(pattern = "beta_pl", x = rownames(sum.vals$summary))[1:7],"mean"]
plot(iso.grm.diff, reps1$b_z)
cor(iso.grm.diff, reps1$b_z)
plot(iso.grm.diff-reps1$b_z)

## Now fit a GAM comparing the standard deviation of the theta 2 estimates versus the real theta value
iso.theta2 = sum.vals$summary[15001:30000,"mean"]
iso.theta2.see = sum.vals$summary[15001:30000,"sd"]
weighted.mean(iso.theta2.see, w = dmnorm(mu = c(0,0), sigma = reps1$varCovMat, x = cbind(iso.theta, iso.theta2)))


model = gam(iso.theta2.see ~ s(iso.theta2, k = 3) + s(iso.theta))

## Now try to plot this?
# Create sequences for x1 and x2
x1_seq <- seq(-3, 3, .1)
x2_seq <- seq(-3, 3, .1)

# Create a data frame for predictions
pred_grid <- expand.grid(iso.theta2 = x1_seq, iso.theta = x2_seq)
pred_grid$predicted_y <- predict(model, newdata = pred_grid)
z_matrix <- matrix(pred_grid$predicted_y, nrow = length(x1_seq), ncol = length(x2_seq))
library(plotly)
plot_ly(x = x1_seq, y = x2_seq, z = z_matrix, type = "surface") %>%
  layout(title = "GAM Predicted Surface",
         scene = list(xaxis = list(title = "x1"),
                      yaxis = list(title = "x2"),
                      zaxis = list(title = "Predicted Y")))

visreg(model, "iso.theta2", gg=TRUE) + xlab("Estimated Severity") + ylab("Standard Deviation") + theme_bw()
true.score.var <- var(reps1$theta$theta2)
error.var <- var(reps1$theta$eapSev - reps1$theta$theta2)
trueRel <- true.score.var / (true.score.var + error.var)
