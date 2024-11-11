## Load library
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
library(ggplot2)
library(visreg)
## This will be used to read all hurdle collapse values
in.dat1 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[1-9].RDS", full.names = TRUE)
in.dat2 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[0-9][0-9].RDS", full.names = TRUE)
in.dat3 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[1][0-4][0-4].RDS", full.names = TRUE)

## Combine these
in.dat <- c(in.dat1, in.dat2, in.dat3)
in.dat <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = ".RDS", full.names = TRUE)
  
## Collapse all of the estimated models
all.collapse <- NULL
all.expand <- NULL
for(i in in.dat){
  tmp <- readRDS(i)
  ## Idenitfy the seed val
  seedVal1 <- basename(i)
  seedVal <- strsplit(strsplit(basename(i), "_")[[1]][2], ".RDS")[[1]][1]
  ## Now print this
  #print(c(seedVal1, seedVal))
  ## Remove all indices where modFail
  tmp.1 <- tmp[[1]]
  tmp.1 <- lapply(tmp.1, function(x) cbind(x, ncol(x)))
  tmp.1 <- dplyr::bind_rows(tmp.1)
  tmp.1$seedVal <- seedVal
  all.collapse <- dplyr::bind_rows(all.collapse, tmp.1)
  tmp.2 <- tmp[[2]]
  tmp.2 <- lapply(tmp.2, function(x) cbind(x, ncol(x)))
  tmp.2 <- dplyr::bind_rows(tmp.2)
  tmp.2$seedVal <- seedVal
  all.expand <- dplyr::bind_rows(all.expand, tmp.2)
  
}
## Grab the rho error
all.collapse$error <- all.collapse$rhoTrue - all.collapse$rho
plot(all.collapse$error, all.collapse$`ncol(x)`)
cor(all.collapse$error, all.collapse$`ncol(x)`)
all.expand$error <- all.expand$rhoTrue - all.expand$rho
plot(all.expand$error, all.expand$`ncol(x)`)
cor(all.expand$error, all.expand$`ncol(x)`)

## Now run an anova
## Attach real vals
sim.param.nitems <- c(5,8,11)
sim.param.ncates <- c(7)
sim.param.discri <- c(2)
sim.param.difficF <- c(-1,0)
sim.param.difficC <- c(1,2)
sim.param.sample <- c(10000)
sim.param.faccor <- c(.1,.4,.8)
sim.param.difgrmF <- c(-2,-.5)
sim.param.difgrmC <- c(1,2)
sim.iter <- 1:100
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.difficF, sim.param.difficC,sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCategory", "discrim", "diff2PLF","diff2PLC", "n", "facCor", "difGrmF","difGrmC" ,"iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)
all.dat.collapse <- merge(all.collapse, all.sim.vals, by.x = "seedVal", by.y="seed")
## Turn these into factors
for(i in c("nItem", "nCategory", "discrim", "diff2PLF","diff2PLC", "n", "difGrmF","difGrmC" ,"iter", "ncol(x)", "rhoTrue")){
  all.dat.collapse[,i] <- factor(all.dat.collapse[,i])
}
all.dat.collapse$nCol <- all.dat.collapse$`ncol(x)`
## Now i need to isolate unique permutations
all.dat.collapse$uniqueIdent <- paste(all.dat.collapse$seedVal, all.dat.collapse$item, all.dat.collapse$nCol)
all.dat.collapse <- all.dat.collapse[!duplicated(all.dat.collapse),]

## Now run anova
mod.one <- lm(error ~ (nCol + nItem + diff2PLF + diff2PLC + facCor + difGrmF + difGrmC)^3, data = all.dat.collapse)
mod.one.res <- car::Anova(mod.one)
effectsize::eta_squared(mod.one.res)
#alias(mod.one)
#car::Anova(mod.one)
library(visreg)
## Now ggplot these me
all.me.vals <- c("nCol", "nItem", "diff2PLF","diff2PLC", "facCor", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.collapse, measurevar = "rho", groupvars = c("nCol", "rhoTrue"))
d1$colVal <- "nCol"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals[-1]){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "rho", groupvars = c(i, "rhoTrue"))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1 <- ggplot(d1, aes(x=facLevel, y=rho, group=rhoTrue, fill=rhoTrue)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")
## Now do the same for expand
all.dat.expand <- merge(all.expand, all.sim.vals, by.x = "seedVal", by.y="seed")
for(i in c("nItem", "nCategory", "discrim", "diff2PLF","diff2PLC", "n", "facCor", "difGrmF","difGrmC" ,"iter", "ncol(x)", "rhoTrue")){
  all.dat.expand[,i] <- factor(all.dat.expand[,i])
}
all.dat.expand$nCol <- all.dat.expand$`ncol(x)`
all.dat.expand$nCol <- factor(all.dat.expand$nCol)
all.dat.expand$uniqueIdent <- paste(all.dat.expand$seedVal, all.dat.expand$item, all.dat.expand$nCol)
all.dat.expand <- all.dat.expand[!duplicated(all.dat.expand),]
mod.two <- lm(error ~ (nCol + nItem + diff2PLF + diff2PLC + difGrmF + difGrmC + rhoTrue)^4, data = all.dat.expand)
aov.res <- car::Anova(mod.two)
ef.sizes <- effectsize::cohens_f(aov.res)
effectsize::eta_squared(aov.res)
# visreg(mod.two, "nCol")
# visreg(mod.two, "nCol", "rhoTrue")
#visreg(mod.one, "nCol")

all.me.vals <- c("nCol", "nItem", "diff2PLF","diff2PLC", "facCor", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = "nCol")
d1$colVal <- "nCol"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals[-1]){
  d2 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = i)
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
## Plot this
p2 <- ggplot(d1, aes(x=facLevel, y=rho)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se)) +
  facet_wrap(colVal ~ ., scales = "free")

## Now do all of these with the true Rho
all.me.vals <- c("nCol", "nItem", "diff2PLF","diff2PLC", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c("nCol", "rhoTrue"))
d1$colVal <- "nCol"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals[-1]){
  d2 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c(i, "rhoTrue"))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p2 <- ggplot(d1, aes(x=facLevel, y=rho, group=rhoTrue, fill=rhoTrue)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")

## Now do these all with ncol
all.me.vals <- c("nItem", "diff2PLF","diff2PLC", "facCor", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c("rhoTrue", "nCol"))
d1$colVal <- "nCol"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals[-1]){
  d2 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c(i, "nCol"))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
ggplot(d1, aes(x=facLevel, y=rho, group=nCol, fill=nCol)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")

## Now do all two-way interactions with the nCOl varaible
all.me.vals <- c("nItem", "diff2PLF","diff2PLC", "facCor", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c("rhoTrue", "nCol"))
d1$colVal <- "rhoTrue"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals[-1]){
  d2 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c(i, "nCol"))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
ggplot(d1, aes(x=facLevel, y=rho, group=nCol, fill=nCol)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")
