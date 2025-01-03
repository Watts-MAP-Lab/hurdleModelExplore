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
all.collapse <- NULL
# all.expand <- NULL
for(i in in.dat){
  tmp <- readRDS(i)
  ## Idenitfy the seed val
  seedVal1 <- basename(i)
  seedVal <- strsplit(strsplit(basename(i), "_")[[1]][2], ".RDS")[[1]][1]
  ## Now print this
  #print(c(seedVal1, seedVal))
  ## Grab the coefficeints
  tmp.dat <- tmp$allParams
  ## Now idenitfy the lower 2PL diff value
  tmp.dat$plFloor <- rep(1:6, each = dim(tmp.dat)[1] / 6)
  tmp.dat$seedVal <- seedVal
  all.collapse <- dplyr::bind_rows(all.collapse, tmp.dat)
}
# saveRDS(all.collapse, file = "./data/allCollapse.RDS")
# saveRDS(all.expand, file = "./data/allExpand.RDS")
#all.collapse <- readRDS("./data/allCollapse.RDS")
#all.expand <- readRDS("./data/allExpand.RDS")

## Grab the rho error
all.collapse$error <- all.collapse$rhoTrue - all.collapse$rho

## Now run an anova
## Attach real vals
sim.param.nitems <- c(5,10)
sim.param.ncates <- c(7)
sim.param.discri <- c(1.5, 3)
sim.param.2pl.spread <- c(.5, 1)
sim.param.sample <- c(5000, 20000)
sim.param.faccor <- c(.1,.4,.8)
sim.param.difgrmF <- c(-2,-.5)
sim.param.difgrmC <- c(1,2)
sim.param.nreps <- c(3, 7)
sim.param.discri2 <- c(1.5, 3)
sim.iter <- 1:100
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.2pl.spread,sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.param.nreps, sim.param.discri2,sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCategory", "discrim", "diffSpread", "n", "facCor", "difGrmF","difGrmC" ,"nCat","grmDiscrim","iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)
all.dat.collapse <- merge(all.collapse, all.sim.vals, by.x = "seedVal", by.y="seed")
## Turn these into factors
for(i in c("nItem", "nCat", "discrim", "diffSpread", "facCor","n", "rhoTrue","difGrmF","difGrmC" ,"iter", "plFloor", "grmDiscrim")){
  all.dat.collapse[,i] <- factor(all.dat.collapse[,i])
}
all.dat.collapse$plFloor <- factor(all.dat.collapse$plFloor, levels=c(1,2,3,4,5,6))
all.dat.collapse$sampSize <- all.dat.collapse$n

## Now i need to isolate unique permutations
all.dat.collapse$uniqueIdent <- paste(all.dat.collapse$seedVal, all.dat.collapse$item, all.dat.collapse$nCol)
all.dat.collapse <- all.dat.collapse[!duplicated(all.dat.collapse),]

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

all.me.vals <- c("nItem","diffSpread","facCor","difGrmF","difGrmC","nCat", "n", "discrim","grmDiscrim")
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
  facet_wrap(colVal ~ ., scales = "free") +
  coord_cartesian(ylim=c(.55, .87))

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
  facet_wrap(colVal ~ ., scales = "free") +
  coord_cartesian(ylim=c(.8, .95))

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
  facet_wrap(colVal ~ ., scales = "free") +
  coord_cartesian(ylim=c(.8, .95))
  
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
  facet_wrap(colVal ~ ., scales = "free") +
  coord_cartesian(ylim=c(.75, .95))

ggpubr::ggarrange(p1, p1c, p1b, p1a, ncol=2, nrow = 2, labels = c("A", "B", "C", "D"))

## Now look into the two way interaction between grmDiscrim & discrim acros all of these
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
