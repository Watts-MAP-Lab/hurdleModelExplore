## Load library
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
library(ggplot2)
library(visreg)

# ## This will be used to read all hurdle collapse values
# in.dat1 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[1-9].RDS", full.names = TRUE)
# in.dat2 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[0-9][0-9].RDS", full.names = TRUE)
# in.dat3 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[1][0-4][0-4].RDS", full.names = TRUE)
# 
# ## Combine these
# in.dat <- c(in.dat1, in.dat2, in.dat3)

in.dat <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = ".RDS", full.names = TRUE)
#   
# ## Collapse all of the estimated models
# all.collapse <- NULL
# # all.expand <- NULL
# for(i in in.dat){
#   tmp <- readRDS(i)
#   ## Identify the seed val
#   seedVal1 <- basename(i)
#   seedVal <- strsplit(strsplit(basename(i), "_")[[1]][2], ".RDS")[[1]][1]
#   ## Now print this
#   #print(c(seedVal1, seedVal))
#   ## Grab the coefficeints
#   tmp.dat <- tmp$allParams
#   ## Now identify the lower 2PL diff value
#   tmp.dat$plFloor <- rep(1:4, each = dim(tmp.dat)[1] / 4)
#   tmp.dat$seedVal <- seedVal
#   all.collapse <- dplyr::bind_rows(all.collapse, tmp.dat)
# }
# saveRDS(all.collapse, file = "./data/allCollapse.RDS")
# saveRDS(all.expand, file = "./data/allExpand.RDS")
all.collapse <- readRDS("./data/allCollapse.RDS")
#all.expand <- readRDS("./data/allExpand.RDS")

## Grab the rho error
all.collapse$error <- all.collapse$rhoTrue - all.collapse$rho
## Now run an anova
## Attach real vals
sim.param.nitems <- c(8,16)
sim.param.ncates <- c(3,5,7)
sim.param.discri <- c(1.2,2.4)
sim.param.2pl.spread <- c(2)
sim.param.sample <- c(10000)
sim.param.faccor <- c(.3,.8)
sim.param.difgrmF <- c(-2,-.5)
sim.param.difgrmC <- c(2)
sim.param.discri2 <- c(1.2, 2.4)
sim.iter <- 1:50
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.2pl.spread,sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.param.discri2,sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCat", "discrim", "diffSpread", "n", "facCor", "difGrmF","difGrmC","grmDiscrim","iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)
all.dat.collapse <- merge(all.collapse, all.sim.vals, by.x = "seedVal", by.y="seed")

## Turn these into factors
for(i in c("nItem", "nCat", "discrim", "diffSpread", "facCor","n", "rhoTrue","difGrmF","difGrmC" ,"iter", "plFloor", "grmDiscrim")){
  all.dat.collapse[,i] <- factor(all.dat.collapse[,i])
}
all.dat.collapse$plFloor <- factor(all.dat.collapse$plFloor, levels=c(1,2,3,4))
all.dat.collapse$sampSize <- all.dat.collapse$n

## Now i need to isolate unique permutations
all.dat.collapse$uniqueIdent <- paste(all.dat.collapse$seedVal, all.dat.collapse$iter,all.dat.collapse$plFloor)
all.dat.collapse <- all.dat.collapse[!duplicated(all.dat.collapse$uniqueIdent),]

## Now run anova
# mod.one <- lm(error ~ (nItem + nCat + diffSpread + facCor + n + difGrmF + difGrmC + plFloor)^4, data = all.dat.collapse)
# mod.one.res <- car::Anova(mod.one)
# effectsize::eta_squared(mod.one.res)
#alias(mod.one)
#car::Anova(mod.one)

## Now do the same for expand
## Now do this for the reliability values
# mod.one <- lm(alpha ~ (nItem + nCat + diffSpread + facCor + n + difGrmF + difGrmC + plFloor)^4, data = all.dat.collapse)
# mod.one.res <- car::Anova(mod.one)
# effectsize::eta_squared(mod.one.res)
## Now make all error values here
# all.dat.collapse$alpha <- all.dat.collapse$alpheFromOme
# all.dat.collapse$alpha <- all.dat.collapse$alpha - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$alpha <- all.dat.collapse$alpheFromOme - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$G_six <- all.dat.collapse$G_six - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$omega_t <- all.dat.collapse$omega_t - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$omega_h <- all.dat.collapse$omega_h - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$estHurdleRel <- all.dat.collapse$estHurdleRel - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$estHurdleRel2 <- all.dat.collapse$estHurdleRel2 - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$estHurdleRelB <- all.dat.collapse$estHurdleRelB - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$estHurdleRel2B <- all.dat.collapse$estHurdleRel2B - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$unidim <- all.dat.collapse$unidim - all.dat.collapse$trueHurdleRel2
# all.dat.collapse$grmRel <- all.dat.collapse$grmRel - all.dat.collapse$trueHurdleRel2
## Now do the same with the root mean sqaured error^2)
all.dat.collapse$alpha <-            sqrt((all.dat.collapse$alpheFromOme)^2)
all.dat.collapse$alpha <-            sqrt((all.dat.collapse$alpha - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$alpha <-            sqrt((all.dat.collapse$alpheFromOme - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$G_six <-            sqrt((all.dat.collapse$G_six - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$omega_t <-          sqrt((all.dat.collapse$omega_t - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$omega_h <-          sqrt((all.dat.collapse$omega_h - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$estHurdleRel <-     sqrt((all.dat.collapse$estHurdleRel - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$estHurdleRel2 <-    sqrt((all.dat.collapse$estHurdleRel2 - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$estHurdleRelB <-    sqrt((all.dat.collapse$estHurdleRelB - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$estHurdleRel2B <-   sqrt((all.dat.collapse$estHurdleRel2B - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$unidim <-           sqrt((all.dat.collapse$unidim - all.dat.collapse$trueHurdleRel)^2)
all.dat.collapse$grmRel <-           sqrt((all.dat.collapse$grmRel - all.dat.collapse$trueHurdleRel)^2)




all.me.vals <- c("nItem","facCor","difGrmF","nCat", "discrim","grmDiscrim")
d1 <- summarySE(data = all.dat.collapse, measurevar = "alpha", groupvars = c("plFloor"))
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "alpha", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1 <- ggplot(d1, aes(x=facLevel, y=alpha)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = alpha - se, ymax = alpha + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")

## Now do this for the omega values too
d1 <- summarySE(data = all.dat.collapse, measurevar = "omega_h", groupvars = c("plFloor"))
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "omega_h", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1a <- ggplot(d1, aes(x=facLevel, y=omega_h)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = omega_h - se, ymax = omega_h + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free") 

d1 <- summarySE(data = all.dat.collapse, measurevar = "omega_t", groupvars = c("plFloor"))
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "omega_t", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1b <- ggplot(d1, aes(x=facLevel, y=omega_t)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = omega_t - se, ymax = omega_t + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free") 
  
d1 <- summarySE(data = all.dat.collapse, measurevar = "G_six", groupvars = c("plFloor"))
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "G_six", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1c <- ggplot(d1, aes(x=facLevel, y=G_six)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = G_six - se, ymax = G_six + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")
## Now do our hurdle info examine
d1 <- summarySE(data = all.dat.collapse, measurevar = "grmRel", groupvars = c("plFloor"), na.rm = TRUE)
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "grmRel", groupvars = c(i), na.rm=TRUE)
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1d <- ggplot(d1, aes(x=facLevel, y=grmRel)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = grmRel - se, ymax = grmRel + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")
## Now do our hurdle info examine
d1 <- summarySE(data = all.dat.collapse, measurevar = "estHurdleRel2B", groupvars = c("plFloor"))
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "estHurdleRel2B", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1e <- ggplot(d1, aes(x=facLevel, y=estHurdleRel2B)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = estHurdleRel2B - se, ymax = estHurdleRel2B + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")

## ANd here is the true info from a hurdle model
## Now do our hurdle info examine
d1 <- summarySE(data = all.dat.collapse, measurevar = "trueHurdleRel", groupvars = c("plFloor"))
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "trueHurdleRel", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1G <- ggplot(d1, aes(x=facLevel, y=trueHurdleRel)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = trueHurdleRel - se, ymax = trueHurdleRel + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")

## Now grab the other two forms 
d1 <- summarySE(data = all.dat.collapse, measurevar = "estHurdleRel", groupvars = c("plFloor"))
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "estHurdleRel", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1f <- ggplot(d1, aes(x=facLevel, y=estHurdleRel)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = estHurdleRel - se, ymax = estHurdleRel + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")

## Now grab the other two forms 
d1 <- summarySE(data = all.dat.collapse, measurevar = "unidim", groupvars = c("plFloor"))
d1$colVal <- "plFloor"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "unidim", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1g <- ggplot(d1, aes(x=facLevel, y=unidim)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = unidim - se, ymax = unidim + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")


all.one <- ggpubr::ggarrange(p1, p1c, p1b, p1a, p1d,p1e,p1f,p1g,ncol=4, nrow = 2, labels = c("A", "B", "C", "D", "E","F"))

ggsave(filename = "./mainEffectExplore2.png", plot = all.one, width = 16, height = 8, dpi = 300, units = "in")
ggsave(filename = "~/Downloads/mainEffectExplore2.png", plot = all.one, width = 16, height = 8, dpi = 300, units = "in")

ggpubr::ggarrange(all.one, p1G)

## Now make a histogram of the amount of zero inflation across various 2pl difficulty patterns
# source("./scripts/rCode/hurdleFunctions.r")
# suppressPackageStartupMessages(library("tidyverse"))
# library("numDeriv")
# library("mirt")
# seedVal <- 1
# a = rep(all.sim.vals$discrim[seedVal], all.sim.vals$nItem[seedVal])
# b = genDiffGRM(num_items = all.sim.vals$nItem[seedVal], num_categories = all.sim.vals$nCat[seedVal], min = all.sim.vals$difGrmF[seedVal], max = all.sim.vals$difGrmC[seedVal], rnorm_var = .3)
# a_z = rep(all.sim.vals$grmDiscrim[seedVal], all.sim.vals$nItem[seedVal])
# ## Need to generate 4 separate b_z levels
# b_z1 = runif(all.sim.vals$nItem[seedVal], min = -2, max=-2+all.sim.vals$diffSpread[seedVal])
# b_z2 = runif(all.sim.vals$nItem[seedVal], min = -1, max=-1+all.sim.vals$diffSpread[seedVal])
# b_z3 = runif(all.sim.vals$nItem[seedVal], min = 0, max=0+all.sim.vals$diffSpread[seedVal])
# b_z4 = runif(all.sim.vals$nItem[seedVal], min = 1, max=1+all.sim.vals$diffSpread[seedVal])
# muVals = c(0,0)
# rho <- all.sim.vals$facCor[seedVal]
# varCovMat = matrix(c(1,rho,rho,1), ncol=2)
# N = all.sim.vals$n[seedVal]
# ## Now generate theta here
# theta = MASS::mvrnorm(n = N, mu = muVals, Sigma = varCovMat)
# library(dplyr)
# reps1 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z1, muVals = muVals, varCovMat = varCovMat, theta = theta)
# reps2 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z2, muVals = muVals, varCovMat = varCovMat, theta = theta)
# reps3 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z3, muVals = muVals, varCovMat = varCovMat, theta = theta)
# reps4 = simulate_hurdle_responses(a = a, b = b, a_z = a_z, b_z = b_z4, muVals = muVals, varCovMat = varCovMat, theta = theta)
# 

## Now look into 2pl diff error
all.collapse$diffError2pl <- all.collapse$est_z_diff - all.collapse$true_z_diff
all.collapse$diffError2plB <- all.collapse$est_z_diff_Bayes - all.collapse$true_z_diff

## Now plot this error by 0-inflation level
d1 <- summarySE(data = all.collapse, measurevar = "diffError2pl", groupvars = c("plFloor"))
p1g <- ggplot(d1, aes(x=plFloor, y=diffError2pl)) +
  #geom_bar(position="dodge", stat="identity") +
  geom_violin(data = all.collapse, mapping = aes(x=plFloor, y=diffError2pl, group = cut_width(plFloor, 1)), scale = "width") +
  theme_bw()



# ## Now look into the two way interaction between grmDiscrim & discrim acros all of these
# tmp <- cbind(rowSums(reps1$responses),rowSums(reps2$responses),rowSums(reps3$responses),rowSums(reps4$responses))
# for.plot <- reshape2::melt(tmp)
# ggplot(for.plot, aes(x=value)) +
#   geom_histogram(bins = 20) +
#   facet_wrap(Var2~.)
plot1 <- summarySE(data = all.dat.collapse, measurevar = "alpha", groupvars = c("discrim", "grmDiscrim"))
p1 <- ggplot(plot1, aes(x=discrim, y=alpha, group = grmDiscrim, fill=grmDiscrim)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = alpha - se, ymax = alpha + se),position="dodge") +
  coord_cartesian(ylim=c(.65, .9))
plot2 <- summarySE(data = all.dat.collapse, measurevar = "omega_h", groupvars = c("discrim", "grmDiscrim"))
p2 <- ggplot(plot2, aes(x=discrim, y=omega_h, group = grmDiscrim, fill=grmDiscrim)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = omega_h - se, ymax = omega_h + se),position="dodge") +
  coord_cartesian(ylim=c(.8, .95))
plot3 <- summarySE(data = all.dat.collapse, measurevar = "omega_t", groupvars = c("discrim", "grmDiscrim"))
p3 <- ggplot(plot3, aes(x=discrim, y=omega_t, group = grmDiscrim, fill=grmDiscrim)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = omega_t - se, ymax = omega_t + se),position="dodge") +
  coord_cartesian(ylim=c(.8, .95))
plot4 <- summarySE(data = all.dat.collapse, measurevar = "G_six", groupvars = c("discrim", "grmDiscrim"))
p4 <- ggplot(plot4, aes(x=discrim, y=G_six, group = grmDiscrim, fill=grmDiscrim)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = G_six - se, ymax = G_six + se),position="dodge") +
  coord_cartesian(ylim=c(.8, .95))
ggpubr::ggarrange(p1, p2, p3, p4, ncol=2, nrow = 2, labels = c("A", "B", "C", "D"))

## Now model each of these outcomes
mod.1 <- lm(alpha ~ (nItem + diffSpread + facCor + difGrmF + difGrmC + difGrmF + nCat + sampSize + discrim + grmDiscrim)^3, data = all.dat.collapse)
mod.2 <- lm(omega_t ~ (nItem + diffSpread + facCor + difGrmF + difGrmC + difGrmF + nCat + sampSize + discrim + grmDiscrim)^3, data = all.dat.collapse)
mod.3 <- lm(omega_h ~ (nItem + diffSpread + facCor + difGrmF + difGrmC + difGrmF + nCat + sampSize + discrim + grmDiscrim)^3, data = all.dat.collapse)
mod.4 <- lm(G_six ~ (nItem + diffSpread + facCor + difGrmF + difGrmC + difGrmF + nCat + sampSize + discrim + grmDiscrim)^3, data = all.dat.collapse)

## Now run through the anovas here
anova.one <- car::Anova(mod.1)
anova.two <- car::Anova(mod.2)
anova.thr <- car::Anova(mod.3)
anova.fou <- car::Anova(mod.4)

## Now effect sizes
eff.one <- effectsize::cohens_f(anova.one)
eff.two <- effectsize::cohens_f(anova.two)
eff.thr <- effectsize::cohens_f(anova.thr)
eff.fou <- effectsize::cohens_f(anova.fou)

## Order
