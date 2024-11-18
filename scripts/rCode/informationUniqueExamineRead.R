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
# in.dat <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = ".RDS", full.names = TRUE)
#   
# ## Collapse all of the estimated models
# all.collapse <- NULL
# all.expand <- NULL
# for(i in in.dat){
#   tmp <- readRDS(i)
#   ## Idenitfy the seed val
#   seedVal1 <- basename(i)
#   seedVal <- strsplit(strsplit(basename(i), "_")[[1]][2], ".RDS")[[1]][1]
#   ## Now print this
#   #print(c(seedVal1, seedVal))
#   ## Remove all indices where modFail
#   tmp.1 <- tmp[[1]]
#   tmp.1 <- lapply(tmp.1, function(x) cbind(x, ncol(x)))
#   tmp.1 <- dplyr::bind_rows(tmp.1)
#   tmp.1$seedVal <- seedVal
#   all.collapse <- dplyr::bind_rows(all.collapse, tmp.1)
#   tmp.2 <- tmp[[2]]
#   tmp.2 <- lapply(tmp.2, function(x) cbind(x, ncol(x)))
#   tmp.2 <- dplyr::bind_rows(tmp.2)
#   tmp.2$seedVal <- seedVal
#   all.expand <- dplyr::bind_rows(all.expand, tmp.2)
# }
# saveRDS(all.collapse, file = "./data/allCollapse.RDS")
# saveRDS(all.expand, file = "./data/allExpand.RDS")
all.collapse <- readRDS("./data/allCollapse.RDS")
all.expand <- readRDS("./data/allExpand.RDS")

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
all.collapse$nResp <- all.collapse$`ncol(x)` - 16
all.expand$nResp <- all.expand$`ncol(x)` - 16
all.dat.collapse <- merge(all.collapse, all.sim.vals, by.x = "seedVal", by.y="seed")
## Turn these into factors
for(i in c("nItem", "nCategory", "discrim", "diff2PLF","diff2PLC", "n", "difGrmF","difGrmC" ,"iter", "nResp", "rhoTrue")){
  all.dat.collapse[,i] <- factor(all.dat.collapse[,i])
}
## Now i need to isolate unique permutations
all.dat.collapse$uniqueIdent <- paste(all.dat.collapse$seedVal, all.dat.collapse$item, all.dat.collapse$nCol)
all.dat.collapse <- all.dat.collapse[!duplicated(all.dat.collapse),]

## Now run anova
mod.one <- lm(error ~ (nResp + nItem + diff2PLF + diff2PLC + facCor + difGrmF + difGrmC)^4, data = all.dat.collapse)
mod.one.res <- car::Anova(mod.one)
effectsize::eta_squared(mod.one.res)
#alias(mod.one)
#car::Anova(mod.one)
library(visreg)
## Now ggplot these me
all.me.vals <- c("nResp", "nItem", "diff2PLF","diff2PLC", "difGrmF","difGrmC", "rhoTrue")
d1 <- summarySE(data = all.dat.collapse, measurevar = "rho", groupvars = c("nResp"))
d1$colVal <- "nCol"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals[-1]){
  d2 <- summarySE(data = all.dat.collapse, measurevar = "rho", groupvars = c(i))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p1 <- ggplot(d1, aes(x=facLevel, y=rho)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")
## Now do the two-way interaction with the true rho values
all.me.vals <- c("nItem", "diff2PLF","diff2PLC", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.collapse, measurevar = "rho", groupvars = c("nResp", "rhoTrue"))
d1$colVal <- "nResp"
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
for(i in c("nItem", "nCategory", "nResp","discrim", "diff2PLF","diff2PLC", "n", "difGrmF","difGrmC" ,"iter", "ncol(x)", "rhoTrue")){
  all.dat.expand[,i] <- factor(all.dat.expand[,i])
}
all.dat.expand$uniqueIdent <- paste(all.dat.expand$seedVal, all.dat.expand$item, all.dat.expand$nCol)
all.dat.expand <- all.dat.expand[!duplicated(all.dat.expand),]
mod.two <- lm(error ~ (nResp + nItem + diff2PLF + diff2PLC + difGrmF + difGrmC + rhoTrue)^4, data = all.dat.expand)
aov.res <- car::Anova(mod.two)
ef.sizes <- effectsize::cohens_f(aov.res)
effectsize::eta_squared(aov.res)
# visreg(mod.two, "nCol")
# visreg(mod.two, "nCol", "rhoTrue")
#visreg(mod.one, "nCol")
all.me.vals <- c("nItem", "diff2PLF","diff2PLC", "rhoTrue", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = "nResp")
d1$colVal <- "nResp"
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
all.me.vals <- c("nResp", "nItem", "diff2PLF","diff2PLC", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c("nResp", "rhoTrue"))
d1$colVal <- "nResp"
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

## Now do the same for the nResp
all.me.vals <- c("nItem", "diff2PLF","diff2PLC", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c("rhoTrue", "nResp"))
d1$colVal <- "rhoTrue"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c(i, "nResp"))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p2 <- ggplot(d1, aes(x=facLevel, y=rho, group=nResp, fill=nResp)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")

## Now all threeway interaction between rhoTrue & nResp
all.me.vals <- c("diff2PLF","diff2PLC", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c("nItem", "rhoTrue", "nResp"))
d1$colVal <- "nItem"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c(i, "rhoTrue","nResp"))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p3 <- ggplot(d1, aes(x=facLevel, y=rho, group=nResp, fill=nResp)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se),position="dodge") +
  facet_grid(rhoTrue ~ colVal, scales = "free")

d1 <- summarySE(data = all.dat.expand, measurevar = "rho", groupvars = c("difGrmF", "difGrmC", "rhoTrue"))
p3 <- ggplot(d1, aes(x=difGrmF, y=rho, group=rhoTrue, fill=rhoTrue)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = rho - se, ymax = rho + se),position="dodge") +
  facet_grid(. ~ difGrmC, scales = "free")

## Now do 2PL difficulty estimation error?
## First create the error terms
all.collapse$diff2PLError <- all.collapse$true_z_diff - all.collapse$est_z_diff
all.expand$diff2PLError <- all.expand$est_z_diff - all.expand$true_z_diff
## Now model these?
all.dat.collapse <- merge(all.collapse, all.sim.vals, by.x = "seedVal", by.y="seed")
## Turn these into factors
for(i in c("nItem", "nCategory", "discrim", "diff2PLF","diff2PLC", "n", "difGrmF","difGrmC" ,"iter", "nResp", "rhoTrue")){
  all.dat.collapse[,i] <- factor(all.dat.collapse[,i])
}
mod.thr <- lm(diff2PLError ~ (nResp + nItem + diff2PLF + diff2PLC + facCor + difGrmF + difGrmC)^4, data = all.dat.collapse)
mod.res <- car::Anova(mod.thr)
effectsize::eta_squared(mod.res)
## Now do the expand
all.dat.expand <- merge(all.expand, all.sim.vals, by.x = "seedVal", by.y="seed")
for(i in c("nItem", "nCategory", "nResp","discrim", "diff2PLF","diff2PLC", "n", "difGrmF","difGrmC" ,"iter", "ncol(x)", "rhoTrue")){
  all.dat.expand[,i] <- factor(all.dat.expand[,i])
}
mod.fou <- lm(diff2PLError ~ (nResp + nItem + diff2PLF + diff2PLC + facCor + difGrmF + difGrmC)^4, data = all.dat.expand)
mod.res <- car::Anova(mod.fou)
effectsize::eta_squared(mod.res)[order(effectsize::eta_squared(mod.res)[,2], decreasing = TRUE),]

all.me.vals <- c("nItem", "diff2PLF","diff2PLC", "rhoTrue", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "diff2PLError", groupvars = "nResp")
d1$colVal <- "nResp"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.expand, measurevar = "diff2PLError", groupvars = i)
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p2 <- ggplot(d1, aes(x=facLevel, y=diff2PLError)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = diff2PLError - se, ymax = diff2PLError + se)) +
  facet_wrap(colVal ~ ., scales = "free") +
  ylab("Error (Est - True)")

## Now do the same for the nResp
all.me.vals <- c("nItem", "diff2PLF","diff2PLC", "difGrmF","difGrmC")
d1 <- summarySE(data = all.dat.expand, measurevar = "diff2PLError", groupvars = c("rhoTrue", "nResp"))
d1$colVal <- "rhoTrue"
colnames(d1)[1] <- "facLevel"
for(i in all.me.vals){
  d2 <- summarySE(data = all.dat.expand, measurevar = "diff2PLError", groupvars = c(i, "nResp"))
  d2$colVal <- i
  colnames(d2)[1] <- "facLevel"
  d1 <- dplyr::bind_rows(d1, d2)
}
p2 <- ggplot(d1, aes(x=facLevel, y=diff2PLError, group=nResp, fill=nResp)) +
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = diff2PLError - se, ymax = diff2PLError + se),position="dodge") +
  facet_wrap(colVal ~ ., scales = "free")
