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
set.seed(16)
n_items = 7
a = runif(n_items, min = .8, max = 2)
b = genDiffGRM(num_items = n_items, num_categories = 3, min = -2, max = 1)
a_z = runif(n_items, min = .8, max = 2)
b_z = runif(n_items, min = 0, max = 2)
N = 30000
muVals = rep(0, 2)
rho = .4
varCovMat = matrix(c(1,rho,rho,1), ncol=2)
THETA = MASS::mvrnorm(N, mu = muVals, Sigma = varCovMat)

## NOw construct a THETA with fixed values
# THETA <-expand.grid(seq(-4, 4, .25), seq(-4, 4, .25))
# THETA <- as.matrix(THETA)
# THETA <- do.call("rbind", replicate(5, THETA, simplify = FALSE))

## Convert from trad to mirt here
pl2_conv <- apply(cbind(a_z, b_z, 0,0), 1, function(x)traditional2mirt(x, cls="2PL", ncat=2) )
grm_conv <- apply(cbind(a, b), 1, function(x)traditional2mirt(x, cls="graded", ncat=3) )
## Now reassign
a_z <- pl2_conv["a1",]
b_z <- pl2_conv["d",]
a <- grm_conv["a1",]
b <- cbind(grm_conv[3,], grm_conv[2,])
reps1 = sim_mirt_hurdle(a = a, b = b, a_z = a_z, b_z = b_z, rho = .4, THETA=THETA)

## Now create the compositie matrix
composite.irt <- reps1$responses[,1:(ncol(reps1$responses)/2)]
composite.irt2 <- reps1$responses_na[,(ncol(reps1$responses)/2 + 1) : ncol(reps1$responses)] 
#composite.irt2 <- composite.irt2 + 1
composite.irt2[is.na(composite.irt2)] <- 0
composite.irt <- composite.irt + composite.irt2
cross_tabs <- composite.irt2 %>%
  group_by(across(everything())) %>%
  summarise(Count = n(), .groups = "drop") %>%
  ungroup()

## Use this one for the subject specific theta estimates
# stan_data = list(
#   N = nrow(composite.irt2),
#   p = ncol(composite.irt2),
#   y = composite.irt2,
#   k = length(unique(table(composite.irt))),
#   Zero = c(0,0)
# )
# stan_hurd2 = stan(file = "./scripts/stanCode/stanHurdleHURDLE_noRecode.stan", data = stan_data, chains = 1, cores = 1, warmup = 1000, iter = 2000)
# saveRDS(stan_hurd2, file="./data/tmpBayesHurd.RDS")
stan_hurd2 <- readRDS("./data/tmpBayesHurd.RDS")

## Use this one for response pattern specific estimates
## There is a modest performance improvement with the response pattern, but the speed improves
## quickly as the sample size increases
# stan_data = list(
#   M = nrow(cross_tabs),
#   p = ncol(composite.irt2),
#   y_pat = cross_tabs[,1:ncol(composite.irt2)],
#   n_pat = cross_tabs$Count,
#   k = length(unique(table(composite.irt))),
#   Zero = c(0,0)
# )
# stan_hurd2 = stan(file = "./scripts/stanCode/stanHurdleHURDLE_noRecodeRepPat.stan", data = stan_data, chains = 1, cores = 1, warmup = 1000, iter = 2000)

## Compare the estimated theta values
N = nrow(THETA)
sum.vals <- summary(stan_hurd2)
bayes_theta = matrix(sum.vals$summary[1:(N*2),"mean"], ncol = 2, byrow = TRUE)
cor(reps1$theta$trueSus, bayes_theta[,1])
cor(reps1$theta$trueSev, bayes_theta[,2])
cor(reps1$theta$truP_susEAP, bayes_theta[,2])
cor(reps1$theta$truP_sevEAP, bayes_theta[,1])
plot(reps1$theta$trueSus, x=bayes_theta[,2])
plot(reps1$theta$trueSev, x=bayes_theta[,1])
plot(reps1$theta$truP_sevEAP, bayes_theta[,1])
plot(reps1$theta$truP_susEAP, bayes_theta[,2])
colnames(bayes_theta) = c("bayes_Sus", "bayes_Sev")
reps1$theta_vals <- bind_cols(reps1$theta_vals, bayes_theta)
cor(reps1$theta)

iso.2pl.discrim <- sum.vals$summary[grep(pattern = "alpha_pl", x = rownames(sum.vals$summary))[1:ncol(composite.irt)],"mean"]
plot(exp(iso.2pl.discrim), a_z)
cor(iso.2pl.discrim, a_z)
plot(iso.2pl.discrim-a_z)
plot(exp(iso.2pl.discrim)-a_z)
iso.grm.discrim <- sum.vals$summary[grep(pattern = "log_alpha_grm", x = rownames(sum.vals$summary))[1:n_items],"mean"]
plot(iso.grm.discrim, log(a))
cor(iso.grm.discrim, log(a))
plot(iso.grm.discrim-log(a))
plot(exp(iso.grm.discrim)-a)
iso.grm.diff <- sum.vals$summary[grep(pattern = "beta_pl", x = rownames(sum.vals$summary))[1:n_items],"mean"]
plot(iso.grm.diff, b_z)
cor(iso.grm.diff, b_z)
plot(iso.grm.diff-reps1$b_z)

## Now fit a GAM comparing the standard deviation of the theta 2 estimates versus the real theta value
iso.theta2 = bayes_theta
iso.theta2.see = matrix(sum.vals$summary[1:(N*2),"sd"], ncol = 2, byrow = TRUE)
weighted.mean(iso.theta2.see[,2], w = dmnorm(mu = c(0,0), sigma =var(bayes_theta), x = iso.theta2))
## combine and run the GAM models
gam.data = data.frame(bayes_theta, iso.theta2.see)
true.score.var <- var(reps1$theta$trueSev)
error.var <- var(reps1$theta$truP_sevEAP - reps1$theta$trueSev)
trueRel <- true.score.var / (true.score.var + error.var)
1 / (1 + weighted.mean(iso.theta2.see[,2], w = dmnorm(mu = c(0,0), sigma =var(bayes_theta), x = bayes_theta)))

## Plot these 3-d values
library(plotly)
data = data.frame(x = gam.data[,1], y = gam.data[,2], se = gam.data[,4])
data = data.frame(x = reps1$theta_vals$trueSus, y = reps1$theta_vals$trueSev, var = gam.data[,4])
fig1 <- plot_ly(data, x = ~x, y = ~y, z = ~var) %>% add_markers()%>% 
  layout(xaxis = list(title = "Sus Var"),
         yaxis = list(title = "Sev Var"),
         zaxis = list(title = "Theta Variance"))
# data = data.frame(x = gam.data$bayes_Sus, y = gam.data$bayes_Sev, se = gam.data[,4])
# fig2 <- plot_ly(data, x = ~x, y = ~y, z = ~se) %>% add_markers()

## Now recreate this with the estimated test information function
# data = data.frame(x = gam.data[,1], y = gam.data[,2], se = gam.data[,4])
# data = data.frame(x = reps1$theta_vals$trueSus, y = reps1$theta_vals$trueSev, se = reps1$tru_hurd_test_info)
data = data.frame(expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)))
data$var = reps1$tru_hurd_test_info^-1
## NOw reduce data to values wihtin 
colnames(data) <- c("x", "y", "var")
data <- data[which(data$x>-4 & data$x<4 & data$y>-4 & data$y<4),]
fig <- plot_ly(data, x = ~x, y = ~y, z = ~var) %>% add_markers() %>% 
  layout(xaxis = list(title = "Sus Var"),
         yaxis = list(title = "Sev Var"),
         zaxis = list(title = "Theta Variance"))


## Now put these next to each other
subplot(fig1, fig)%>% 
  layout(xaxis = list(title = "Sus Var"),
         yaxis = list(title = "Sev Var"),
         zaxis = list(title = "Theta Variance"),
         xaxis2 = list(title = "Sus Var"),
         yaxis2 = list(title = "Sev Var"),
         zaxis2 = list(title = "Theta Variance"))


data = data.frame(x = gam.data[,1], y = gam.data[,2], var = gam.data[,4])
fig1 <- plot_ly(data, x = ~x, y = ~y, z = ~var) %>% add_markers()%>% 
  layout(xaxis = list(title = "Sus Var"),
         yaxis = list(title = "Sev Var"),
         zaxis = list(title = "Theta Variance"))

## Now merge the theoretical var versus the observed variance values
dataTarg = data.frame(x = gam.data[,1], y = gam.data[,2], var = gam.data[,4])
## Now I need to round these theta estimates so they can match w/ the theoretical values
# expand.grid(seq(-6, 6, .1), seq(-6, 6, .1))
## I need to round to the closest .1 value
dataTarg[,1] <- round(dataTarg[,1], 1)
dataTarg[,2] <- round(dataTarg[,2], 1)
dataSource = data.frame(expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)))
dataSource$var = reps1$tru_hurd_test_info^-1
colnames(dataSource) <- c("x", "y", "var")
## Now merge
data.merge <- merge(dataTarg, dataSource, by=c("x", "y"), suffixes = c("_Bayes", "_Theoretical"))
## Now plot these
ggplot(data.merge, aes(x=var_Bayes, y=var_Theoretical)) +
  geom_point() +
  theme_bw()
  
data = data.frame(expand.grid(seq(-6, 6, .1), seq(-6, 6, .1)))
data$var = reps1$tru_hurd_test_info^-1
## NOw reduce data to values wihtin 
colnames(data) <- c("x", "y", "var")
data <- data[which(data$x>-1.1 & data$x<2.5 & data$y>-.6 & data$y<2.6),]
fig <- plot_ly(data, x = ~x, y = ~y, z = ~var) %>% add_markers() %>% 
  layout(xaxis = list(title = "Sus Var"),
         yaxis = list(title = "Sev Var"),
         zaxis = list(title = "Theta Variance"))
