## Load library
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
library(ggplot2)
library(visreg)
# library("knitr")
# library("kableExtra")

## CP data over
system("rsync -hvrPt --ignore-existing rosena@login.accre.vu:/home/rosena/hurdleModelExplore/data/hurdleCollapse/*RDS ./data/hurdleCollapse/")

# ## This will be used to read all hurdle collapse values
# in.dat1 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[1-9].RDS", full.names = TRUE)
# in.dat2 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[0-9][0-9].RDS", full.names = TRUE)
# in.dat3 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[1][0-4][0-4].RDS", full.names = TRUE)
#
# ## Combine these
# in.dat <- c(in.dat1, in.dat2, in.dat3)

in.dat <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = ".RDS", full.names = TRUE)
#
# Collapse all of the estimated models
# all.collapse <- NULL
# # # all.expand <- NULL
# for(i in in.dat){
#  tmp <- readRDS(i)
#  ## Identify the seed val
#  seedVal1 <- basename(i)
#  seedVal <- strsplit(strsplit(basename(i), "_")[[1]][2], ".RDS")[[1]][1]
#  ## Now print this
#  #print(c(seedVal1, seedVal))
#  ## Grab the coefficients
#  tmp.dat <- tmp$allParams
#  ## Now identify the lower 2PL diff value
#  #tmp.dat$plFloor <- rep(1:4, each = dim(tmp.dat)[1] / 4)
#  tmp.dat$seedVal <- seedVal
#  all.collapse <- dplyr::bind_rows(all.collapse, tmp.dat)
# }
# saveRDS(all.collapse, file = "./data/allCollapse.RDS")
#saveRDS(all.expand, file = "./data/allExpand.RDS")
all.collapse <- readRDS("./data/allCollapse.RDS")
#all.expand <- readRDS("./data/allExpand.RDS")

## Grab the rho error
#all.collapse$error <- all.collapse$rhoTrue - all.collapse$rho
## Now run an anova
## Attach real vals
source("./scripts/rCode/simParam.r")
all.dat.collapse <- merge(all.collapse, all.sim.vals, by.x = "seedVal", by.y="seed")

## Turn these into factors
for(i in c("nItem", "nCat", "grmDiscrim", "diffSpread", "facCor","n", "rhoTrue","difGrmF","difGrmC" ,"iter", "discrim2pl", "dif2PL")){
  all.dat.collapse[,i] <- factor(all.dat.collapse[,i])
}
all.dat.collapse$sampSize <- all.dat.collapse$n

# Now deal with multiple values within each dataset
all.dat.collapse <- all.dat.collapse[!duplicated(all.dat.collapse$seedVal),]
## Now remove middle 2PL difficulty value

#all.dat.collapse <- all.dat.collapse[-which(all.dat.collapse$dif2PL=="-1"),]

## NOw rm impossible values
#all.dat.collapse$alpha[which(all.dat.collapse$alpha<0)] <- NA
#all.dat.collapse$estHurdleRel_MIRT[which(all.dat.collapse$estHurdleRel_MIRT<.1)] <- NA

mean(sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo0_90)^2))
mean(sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo90_0)^2))
mean(sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo45_45)^2))
mean(sqrt((all.dat.collapse$trueRel-all.dat.collapse$omega_t)^2))
mean(sqrt((all.dat.collapse$trueRel-all.dat.collapse$alpha)^2))
mean(sqrt((all.dat.collapse$trueRel-all.dat.collapse$grmRel)^2))
mean(sqrt((all.dat.collapse$trueRel-all.dat.collapse$estHurdleRel)^2))

median(sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo0_90)^2))
median(sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo90_0)^2))
median(sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo45_45)^2))
median(sqrt((all.dat.collapse$trueRel-all.dat.collapse$grmRel)^2))
median(sqrt((all.dat.collapse$trueRel-all.dat.collapse$estHurdleRel)^2))

all.dat.collapse$rseHurd090 <- sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo0_90)^2)
all.dat.collapse$rseHurd900 <- sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo90_0)^2)
all.dat.collapse$rseHurd4545 <- sqrt((all.dat.collapse$trueRel-all.dat.collapse$mirtHurdleInfo45_45)^2)
all.dat.collapse$rseGRMRel <- sqrt((all.dat.collapse$trueRel-all.dat.collapse$grmRel)^2)
all.dat.collapse$rseHurdRel <- sqrt((all.dat.collapse$trueRel-all.dat.collapse$estHurdleRel)^2)
all.dat.collapse$rseOmega <- sqrt((all.dat.collapse$trueRel-all.dat.collapse$omega_t)^2)
all.dat.collapse$rseAlpha <- sqrt((all.dat.collapse$trueRel-all.dat.collapse$alpha)^2)

## Now prep the anova models for all interactions up to four ways
# first declare all outcomes
all.out <- c("alpha", "omega_t", "omega_h", "lambda1Rel", "lambda2Rel","lambda3Rel","alpheFromOme",
             "lambda4Rel", "lambda5Rel", "lambda6Rel", "estHurdleRel", "grmRelExludeZero", "grmRel")
all.out <- c("omega_t", "omega_h", "lambda6Rel", "estHurdleRel", "grmRelExludeZero", "grmRel","alpha")
all.out <- c("singleFactorOmegaT_NoZ", "singleFactorOmegaT","estHurdleRel", "grmRelExludeZero", "grmRel","alpha", "lambda6Rel")
all.out <- c("singleFactorOmegaT", "omega_t", "grmRel", "estHurdleRel","grmRel_rmZeroOption", "alpha") 
all.out <- c("singleFactorOmegaT", "grmRel", "estHurdleRel", "alpha", "trueRel", "rseGRMRel", "rseHurdRel", "rseOmega", "rseAlpha") 
pred.vals <- c("nItem","nCat","facCor","difGrmF","dif2PL", "grmDiscrim", "discrim2pl")
#pred.vals <- c("dif2PL")

## NOw prep all of the models with these data
all.anova.mods <- list()
all.ef.vals <- list()
all.ef.vals2 <- list()
for(i in 1:length(all.out)){
  ## first declare the model
  model.tmp <- as.formula(paste(all.out[i], "~(", paste(pred.vals, collapse = "+"), ")^4"))
  tmp.anova <- lm(model.tmp, data = all.dat.collapse)
  all.anova.mods[[i]] <- tmp.anova
  ## Now estimate the f values
  f.vals <- car::Anova(tmp.anova)
  eff.size.vals <- effectsize::omega_squared(f.vals, partial = TRUE)
  eff.size.vals2 <- effectsize::cohens_f(tmp.anova, partial = TRUE)
  all.ef.vals[[i]] <- eff.size.vals
  all.ef.vals2[[i]] <- eff.size.vals2
}

## Now do one with the skew values
## Now combine these values
all.ef.vals.iso <- lapply(all.ef.vals, function(x) x$Omega2)
#all.ef.vals.iso <- lapply(all.ef.vals, function(x) x$Eta2)

## combine these
all.ef.vals.iso <- dplyr::bind_cols(all.ef.vals[[1]]$Parameter,all.ef.vals.iso)
colnames(all.ef.vals.iso) <- c("Parameter", all.out)

## Now prep the main ef table
main.ef.table <- all.ef.vals.iso[1:7,]
## Print the main effect table here
# library(dplyr)
# main.ef.table.round <- main.ef.table %>% dplyr::mutate(across(where(is.numeric), round, 2))
# main.ef.table.round <- main.ef.table.round[,c("Parameter", "alpha", "omega_t", "G_six","grmRel", "grmRelExludeZero", "estHurdleRel")]
# table_out <- kable(main.ef.table.round, "latex")
# cat(table_out, file="./reports/mainEFTable.tex")
# 
# ## Now do the two-way interaction vals here
# two.way.ef.table <- all.ef.vals.iso[8:53,]
# two.way.ef.table.3 <- two.way.ef.table[which(apply(two.way.ef.table[2:6], 1, function(x) sum(x > .01))>1),]
# two.way.ef.table.3 <- two.way.ef.table.3 %>% dplyr::mutate(across(where(is.numeric), round, 2))
# two.way.ef.table.3 <- two.way.ef.table.3[,c("Parameter", "alpha", "omega_t", "grmRel", "grmRelExludeZero", "estHurdleRel")]
# table_out <- kable(two.way.ef.table.3, "latex")
# cat(table_out, file="./reports/interactionEFTable.tex")


## Alright, I am making an executive decision here -- I am only interested in the following reliability metrics:
## 1. Omega values Total
## 2. Cronbach's alpha -- not from the omega function
## 3. I will do the average split half reliability -- not sure the full reasoning for this one
##    1. I am going to go with G6 --> because: considers the amount of variance in each item that can be accounted for the linear regression of all of the other items (the squared multiple correlation or smc), or more precisely, the variance of the errors, 
## 4. The three IRT derived metrics include:
##    1. estHurdleRel_MIRT --> The focus of this manuscript
##    2. grmRelExcludeZero --> GRM removing all 0 responses
##    3. grmRel --> GRM with all possible responses

## Get our fully crossed sample sizes here?
ful.cross <- summarySE(data = all.dat.collapse, measurevar = c("alpha"), groupvars = c("nItem", "facCor","difGrmF", "nCat", "discrim2pl", "grmDiscrim", "dif2PL"), na.rm=TRUE)
ful.cross$N

## Get mean alpha values across 2PL diff condition
summarySE(data = all.dat.collapse, measurevar = c("alpha"), groupvars = c("dif2PL"), na.rm=TRUE)

## Explore main effect of 2PL here for pi grm & pi h
summarySE(data = all.dat.collapse, measurevar = c("estHurdleRel"), groupvars = c("dif2PL"), na.rm=TRUE)
summarySE(data = all.dat.collapse, measurevar = c("grmRel"), groupvars = c("dif2PL"), na.rm=TRUE)
summarySE(data = all.dat.collapse, measurevar = c("grmRel_rmZeroOption"), groupvars = c("dif2PL"), na.rm=TRUE)
summarySE(data = all.dat.collapse, measurevar = c("alpha"), groupvars = c("dif2PL"), na.rm=TRUE)
summarySE(data = all.dat.collapse, measurevar = c("singleFactorOmegaT"), groupvars = c("dif2PL"), na.rm=TRUE)
summarySE(data = all.dat.collapse, measurevar = c("grmRel_rmZeroOption"), groupvars = c("facCor"), na.rm=TRUE)

## Now compare all permutations of rmse estimates
ful.cross1 <- summarySE(data = all.dat.collapse, measurevar = c("rseGRMRel"), groupvars = c("nItem", "facCor","difGrmF", "nCat", "discrim2pl", "grmDiscrim", "dif2PL"), na.rm=TRUE)
ful.cross2 <- summarySE(data = all.dat.collapse, measurevar = c("rseHurdRel"), groupvars = c("nItem", "facCor","difGrmF", "nCat", "discrim2pl", "grmDiscrim", "dif2PL"), na.rm=TRUE)


## Prep the data for figure 1 -- which will be examining the main effect of plFLoor across
## all reliability metrics
iso.vars <- c(c("nItem", "nCat", "discrim2pl", "diffSpread", "n", "facCor", "difGrmF","difGrmC","grmDiscrim","dif2PL","iter"), all.out, "rhoTrue", "seedVal", "sampSize")
## Now revalue levels of the iso vars here
all.dat.collapse$nItem <-  plyr::revalue(all.dat.collapse[,"nItem"], replace = c("5" = "5", "10" = "10"))
all.dat.collapse$nCat <-   plyr::revalue(all.dat.collapse[,"nCat"], replace = c("3" = "3", "6" = "6"))
all.dat.collapse$discrim2pl <-  plyr::revalue(all.dat.collapse[,"discrim2pl"], replace = c("0.8" = "[0.8:2.0]", "1.8" = "[1.8:3.0]"))
all.dat.collapse$facCor <-  plyr::revalue(all.dat.collapse[,"facCor"], replace = c("0.2" = "0.2", "0.8" = "0.8"))
all.dat.collapse$difGrmF <-  plyr::revalue(all.dat.collapse[,"difGrmF"], replace = c("-2" = "[-3:-1]", "0.5" = "[-0.5:2.5]"))
all.dat.collapse$grmDiscrim<-  plyr::revalue(all.dat.collapse[,"grmDiscrim"], replace = c("1.8" = "[0.8:2.0]", "4.5" = "[1.8:3.0]"))
all.dat.collapse$dif2PL <-  plyr::revalue(all.dat.collapse[,"dif2PL"], replace = c("-2" = "[-2:0]", "1" = "[1:3]"))

## Now change the name of the outcome variables
iso.dat <- reshape2::melt(all.dat.collapse[,c(iso.vars)], id.vars=c("seedVal","rhoTrue","nItem","nCat","facCor","difGrmF","dif2PL", "grmDiscrim","n","discrim2pl","iter","difGrmC", "sampSize","diffSpread"))
## Now change the name of the outcome variables
iso.dat$outcome <- "Raw"
iso.dat$outcome[which(iso.dat$variable %in% c("rseGRMRel", "rseHurdRel", "rseOmega", "rseAlpha"))] <- "RMSE"
iso.dat$variable <- plyr::revalue(iso.dat$variable, c("estHurdleRel" = "Hurdle", "grmRel" = "GRM", "singleFactorOmegaT" = "Omega", "alpha" = "Alpha",  "trueRel" = "True",
                                                      "rseGRMRel" = "GRM", "rseHurdRel" = "Hurdle", "rseOmega" = "Omega", "rseAlpha" = "Alpha"))
iso.dat$variable <- factor(iso.dat$variable, levels = c("Alpha", "Omega", "GRM", "Hurdle", "True"))

## Plot the main effect of plFloor
## First rm impossible values
d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("dif2PL", "variable", "outcome"), na.rm=TRUE)
#d1 <- d1[which(d1$variable %in% c("alpha", "omega_h", "omega_t", "lambda5Rel", "estHurdleRel_MIRT", "grmRel", "grmRelExludeZero")),]
#d1 <- d1[which(d1$variable %in% all.out),]
## Now prep the effect size values her 
# ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[5,-1])), ef.val = t(all.ef.vals.iso[5,-1]))
# ef.vals$roundVal = round(ef.vals$ef.val, 2)
# ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
#dd1 <- d1[which(d1$variable %in% c("alphaFromOme", "alpha","lambda4Rel","omega_h", "omega_t", "lambda5Rel", "estHurdleRel_MIRT", "grmRel", "grmRelExludeZero")),]
## Now plot this
plFloorME1 <- ggplot(d1[which(d1$outcome=="Raw"),], aes(x=dif2PL, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0.5, .9)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  ylab("Estimated Reliability") +
  xlab("2PL Difficulty")
plFloorME2 <- ggplot(d1[which(d1$outcome=="RMSE"),], aes(x=dif2PL, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, .3)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  ylab("RMSE") +
  xlab("2PL Difficulty")
plFloorME <- ggpubr::ggarrange(plFloorME1, plFloorME2)

## Now make one of these for each main effect value, but only the RMSE values
## This is factor correlation
iso.dat <- reshape2::melt(all.dat.collapse[,c(iso.vars)], id.vars=c("seedVal","rhoTrue","nItem","nCat","facCor","difGrmF","dif2PL", "grmDiscrim","n","discrim2pl","iter","difGrmC", "sampSize","diffSpread"))
## Now change the name of the outcome variables
iso.dat$outcome <- "Raw"
iso.dat$outcome[which(iso.dat$variable %in% c("rseGRMRel", "rseHurdRel", "rseOmega", "rseAlpha"))] <- "RMSE"
iso.dat$variable <- plyr::revalue(iso.dat$variable, c("estHurdleRel" = "Hurdle", "grmRel" = "GRM", "singleFactorOmegaT" = "Omega", "alpha" = "Alpha",  "trueRel" = "True",
                                                      "rseGRMRel" = "GRM", "rseHurdRel" = "Hurdle", "rseOmega" = "Omega", "rseAlpha" = "Alpha"))
iso.dat$variable <- factor(iso.dat$variable, levels = c("Alpha", "Omega", "GRM", "Hurdle", "True"))
iso.dat <- iso.dat[-which(iso.dat$outcome=="Raw"),]
all.out <- c("Alpha", "Omega", "GRM", "Hurdle")
d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("facCor", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[3,-1])), ef.val = t(all.ef.vals.iso[3,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
facCorME <- ggplot(d1, aes(x=facCor, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("RMSE") +
  xlab("Factor Correlation")# +
  #geom_label(data = ef.vals, aes(label = labelT), x = 1, y=.97)
## This is number of items
d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nItem", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[1,-1])), ef.val = t(all.ef.vals.iso[1,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
nItemsME <- ggplot(d1, aes(x=nItem, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, .2)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("") +
  xlab("Number of Items")# +
 # geom_label(data = ef.vals, aes(label = labelT), x = 1, y=.97)

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nCat", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[2,-1])), ef.val = t(all.ef.vals.iso[2,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
nCatME <- ggplot(d1, aes(x=nCat, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, .2)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("") +
  xlab("Number of Response Options")# +
  #geom_label(data = ef.vals, aes(label = labelT), x = 1, y=.97)

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("difGrmF", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[4,-1])), ef.val = t(all.ef.vals.iso[4,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
difGRMME <- ggplot(d1, aes(x=difGrmF, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 0.25)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("Reliability") +
  xlab("GRM Difficulty") #+
 # geom_label(data = ef.vals, aes(label = labelT), x = 1, y=.97)

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("discrim2pl", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[7,-1])), ef.val = t(all.ef.vals.iso[7,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
discrim2PLMME <- ggplot(d1, aes(x=discrim2pl, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 0.25)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("") +
  xlab("2PL Discrim") #+
  #geom_label(data = ef.vals, aes(label = labelT), x = 1, y=.97)

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("grmDiscrim", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[6,-1])), ef.val = t(all.ef.vals.iso[6,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
discrimGRMMME <- ggplot(d1, aes(x=grmDiscrim, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 0.25)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("") +
  xlab("GRM Discrim") #+
  #geom_label(data = ef.vals, aes(label = labelT), x = 1, y=.97)

## Now plot all of these secondary me explores
plFloorME
ggpubr::ggarrange(facCorME, nItemsME, nCatME, difGRMME, discrim2PLMME, discrimGRMMME)
plot.all.me <- ggpubr::ggarrange(facCorME, nItemsME, nCatME, difGRMME, discrim2PLMME, discrimGRMMME, labels = "AUTO")

## Now print each of these
ggsave(filename = "./reports/figure4_mainEffect2PLDiff.png", plot = plFloorME, dpi = 300, width = 12, height = 8, bg = "white")
ggsave(filename = "./reports/figure5_allOtherME.png", plot = plot.all.me, dpi = 300, width = 24, height = 14, bg = "white")
ggsave(filename = "./texOut/figures/figure4_mainEffect2PLDiff.png", plot = plFloorME, dpi = 300, width = 12, height = 8, bg = "white")
ggsave(filename = "./texOut/figures/figure5_allOtherME.png", plot = plot.all.me, dpi = 300, width = 24, height = 14, bg = "white")

## Write a csv with all of the effect sizes
write.csv(x = all.ef.vals.iso, file = "~/Desktop/efVals.csv", quote=F, row.names = FALSE)

## Now do the largest two-way interaction?
d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nItem", "grmDiscrim","discrim2pl","variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
## Now plot this
plFloorIE <- ggplot(d1, aes(x=discrim2pl, y=value, group=nItem, fill=nItem, color=nItem)) +
  #geom_bar(position="dodge", stat="identity") +
  geom_line(size = 2, position=position_dodge(width=1)) +
  theme_minimal() +
  facet_grid(grmDiscrim ~ variable) +
  coord_cartesian(ylim = c(0, .5)) +
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),position=position_dodge(width =1), size = 2) +
  ylab("Reliability") +
  scale_fill_grey() +
  scale_color_grey() +
  xlab("2PL Difficulty") +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  labs(fill='Number of Items', color="Number of Items") 

## Now do an ANCOVA with SKEW predicting the reliability metrics
mod <- lm(alpha ~ skewVal * dif2PL, data = all.dat.collapse)

full.plot <- ggpubr::ggarrange(plFloorME, ggpubr::ggarrange(facCorME, nItemsME, nCatME, difGRMME, discrim2PLMME, discrimGRMMME), ncol = 2, labels = "AUTO")
full.plot

# ## Now do the two way interaction for the discrim params
d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nItem", "dif2PL","variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
## Now plot this
plFloorIE1 <- ggplot(d1, aes(x=dif2PL, y=value, group=nItem, fill=nItem, color=nItem)) +
  #geom_bar(position="dodge", stat="identity") +
  geom_line(size = 2, position=position_dodge(width=1)) +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, .4)) +
  geom_errorbar(aes(ymin = ifelse(value - se < 0, 0, value - se), ymax = value + se),position=position_dodge(width =1), size = 2) +
  ylab("Reliability") +
  scale_fill_grey() +
  scale_color_grey() +
  xlab("2PL Difficulty") +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  labs(fill='Number of Items', color="Number of Items") 

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nItem", "grmDiscrim","variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
## Now plot this
plFloorIE1 <- ggplot(d1, aes(x=grmDiscrim, y=value, group=nItem, fill=nItem, color=nItem)) +
  #geom_bar(position="dodge", stat="identity") +
  geom_line(size = 2, position=position_dodge(width=1)) +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, .4)) +
  geom_errorbar(aes(ymin = ifelse(value - se < 0, 0, value - se), ymax = value + se),position=position_dodge(width =1), size = 2) +
  ylab("RMSE") +
  scale_fill_grey() +
  scale_color_grey() +
  xlab("GRM Discrimination") +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  labs(fill='Number of Items', color="Number of Items") 

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nItem", "discrim2pl","variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
## Now plot this
plFloorIE2 <- ggplot(d1, aes(x=discrim2pl, y=value, group=nItem, fill=nItem, color=nItem)) +
  #geom_bar(position="dodge", stat="identity") +
  geom_line(size = 2, position=position_dodge(width=1)) +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, .4)) +
  geom_errorbar(aes(ymin = ifelse(value - sd < 0, 0, value - sd), ymax = value + sd),position=position_dodge(width =1), size = 2) +
  ylab("") +
  scale_fill_grey() +
  scale_color_grey() +
  xlab("2PL Discrimination") +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  labs(fill='Number of Items', color="Number of Items") 

plFloorIE <- ggpubr::ggarrange(plFloorIE1, plFloorIE2, legend = "bottom", common.legend = TRUE, labels = "AUTO")

ggsave(filename = "./reports/figure6_twoWay1.png", plot = plFloorIE, dpi = 300, width = 12, height = 8)
ggsave(filename = "./texOut/figures/figure6_twoWay1.png", plot = plFloorIE, dpi = 300, width = 12, height = 8)

## Now do the two-way interaction between the 2PL & GRM discrimination
d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("discrim2pl", "grmDiscrim","variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
pldiscrimIE2 <- ggplot(d1, aes(x=discrim2pl, y=value, group=grmDiscrim, fill=grmDiscrim, color=grmDiscrim)) +
  #geom_bar(position="dodge", stat="identity") +
  geom_line(size = 2, position=position_dodge(width=1)) +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, .4)) +
  geom_errorbar(aes(ymin = ifelse(value - sd < 0, 0, value - sd), ymax = value + sd),position=position_dodge(width =1), size = 2) +
  ylab("Reliability") +
  scale_fill_grey() +
  scale_color_grey() +
  xlab("2PL Discrimination") +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  labs(fill='GRM Discrimination', color="GRM Discrimination") 

## Now do this three-way interaction


## Now plot the main effects for the raw units here
iso.dat <- reshape2::melt(all.dat.collapse[,c(iso.vars)], id.vars=c("seedVal","rhoTrue","nItem","nCat","facCor","difGrmF","dif2PL", "grmDiscrim","n","discrim2pl","iter","difGrmC", "sampSize","diffSpread"))
## Now change the name of the outcome variables
iso.dat$outcome <- "Raw"
iso.dat$outcome[which(iso.dat$variable %in% c("rseGRMRel", "rseHurdRel", "rseOmega", "rseAlpha"))] <- "RMSE"
iso.dat$variable <- plyr::revalue(iso.dat$variable, c("estHurdleRel" = "Hurdle", "grmRel" = "GRM", "singleFactorOmegaT" = "Omega", "alpha" = "Alpha",  "trueRel" = "True",
                                                      "rseGRMRel" = "GRM", "rseHurdRel" = "Hurdle", "rseOmega" = "Omega", "rseAlpha" = "Alpha"))
iso.dat$variable <- factor(iso.dat$variable, levels = c("Alpha", "Omega", "GRM", "Hurdle", "True"))
iso.dat <- iso.dat[which(iso.dat$outcome=="Raw"),]
all.out <- c("Alpha", "Omega", "GRM", "Hurdle", "True")
d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("facCor", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[3,-1])), ef.val = t(all.ef.vals.iso[3,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
facCorME <- ggplot(d1, aes(x=facCor, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("Reliability") +
  xlab("Factor Correlation")# +

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nItem", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
nItemME <- ggplot(d1, aes(x=nItem, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("Reliability") +
  xlab("n Items")# +

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nCat", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
nCatME <- ggplot(d1, aes(x=nCat, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("") +
  xlab("Number of Response Options")# +

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("difGrmF", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[4,-1])), ef.val = t(all.ef.vals.iso[4,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
difGRMME <- ggplot(d1, aes(x=difGrmF, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("Reliability") +
  xlab("GRM Difficulty") #+
# geom_label(data = ef.vals, aes(label = labelT), x = 1, y=.97)

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("discrim2pl", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[7,-1])), ef.val = t(all.ef.vals.iso[7,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
discrim2PLMME <- ggplot(d1, aes(x=discrim2pl, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("") +
  xlab("2PL Discrim") #+

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("grmDiscrim", "variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
ef.vals <- data.frame(variable = rownames(t(all.ef.vals.iso[6,-1])), ef.val = t(all.ef.vals.iso[6,-1]))
ef.vals$roundVal = round(ef.vals$ef.val, 2)
ef.vals$labelT = paste("es = ", ef.vals$roundVal, sep='')
discrimGRMMME <- ggplot(d1, aes(x=grmDiscrim, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),position="dodge", width = 0, size = 3) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), axis.text.y = element_text(size = 11)) +
  ylab("") +
  xlab("GRM Discrim") #+
all.raw.me <- ggpubr::ggarrange(facCorME, nItemME, nCatME, difGRMME, discrim2PLMME, discrimGRMMME, labels = "AUTO")
ggsave(filename = "./reports/figure5_allOtherME_Supplement.png", plot = all.raw.me, dpi = 300, width = 24, height = 14, bg = "white")
ggsave(filename = "./texOut/figures/figure5_allOtherME_Supplement.png", plot = all.raw.me, dpi = 300, width = 24, height = 14, bg = "white")

## Now do the two two-way interactions down here
d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nItem", "grmDiscrim","variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
## Now plot this
plFloorIE1 <- ggplot(d1, aes(x=grmDiscrim, y=value, group=nItem, fill=nItem, color=nItem)) +
  #geom_bar(position="dodge", stat="identity") +
  geom_line(size = 2, position=position_dodge(width=1)) +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(.4, 1)) +
  geom_errorbar(aes(ymin = ifelse(value - sd < 0, 0, value - sd), ymax = value + sd),position=position_dodge(width =1), size = 2) +
  ylab("Reliability") +
  scale_fill_grey() +
  scale_color_grey() +
  xlab("GRM Discrimination") +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  labs(fill='Number of Items', color="Number of Items") 

d1 <- summarySE(iso.dat, measurevar = "value", groupvars = c("nItem", "discrim2pl","variable"), na.rm=TRUE)
d1 <- d1[which(d1$variable %in% all.out),]
## Now plot this
plFloorIE2 <- ggplot(d1, aes(x=discrim2pl, y=value, group=nItem, fill=nItem, color=nItem)) +
  #geom_bar(position="dodge", stat="identity") +
  geom_line(size = 2, position=position_dodge(width=1)) +
  theme_minimal() +
  facet_grid(. ~ variable) +
  coord_cartesian(ylim = c(.4, 1)) +
  geom_errorbar(aes(ymin = ifelse(value - sd < 0, 0, value - sd), ymax = value + sd),position=position_dodge(width =1), size = 2) +
  ylab("") +
  scale_fill_grey() +
  scale_color_grey() +
  xlab("2PL Discrimination") +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black")) +
  labs(fill='Number of Items', color="Number of Items") 

plFloorIE <- ggpubr::ggarrange(plFloorIE1, plFloorIE2, legend = "bottom", common.legend = TRUE, labels = "AUTO")

ggsave(filename = "./reports/figure6_Int_Supplement.png", plot = plFloorIE, dpi = 300, width = 12, height = 8, bg = "white")
ggsave(filename = "./texOut/figures/figure6_Int_Supplement.png", plot = plFloorIE, dpi = 300, width = 12, height = 8, bg = "white")

## Now plot the correlation of reliability metrics for my information here
cor.map <- all.dat.collapse[,iso.vars[-c(1,2,23:34)]]
cor.map <- cor.map[complete.cases(cor.map),]
library(GGally)
#ggpairs(all.dat.collapse[,c("alpha","omega_h", "omega_t", "lambda6Rel", "estHurdleRel_MIRT", "grmRel", "grmRelExludeZero")])
#out <- ggpairs(all.dat.collapse[,c("alpha","omega_h", "omega_t", "lambda5Rel", "estHurdleRel", "grmRel", "grmRelExludeZero", "dif2PL")], ggplot2::aes(colour = dif2PL))
#ggsave(filename = "./reports/reliabilityCorMat.png", plot=out, width=20, height = 16, units="in", dpi=300)
#ggsave(filename = "~/Downloads/reliabilityCorMat.png", plot=out, width=20, height = 16, units="in", dpi=300)
out <- ggpairs(all.dat.collapse[,c("alpha", "estHurdleRel", "omega_t","singleFactorOmegaT","grmRel" ,"skewVal","dif2PL")], ggplot2::aes(colour = dif2PL))

## Now run some quick interaction plots
library(visreg)
all.visreg <- list()
all.ancova.mods <- list()
for(i in 1:length(all.out)){
  ## first declare the model
  #model.tmp <- as.formula(paste(all.out[i], "~(", paste(pred.vals, collapse = "+"), ")^4"))
  model.tmp <- as.formula(paste(all.out[i], "~(skewVal * dif2PL)"))
  tmp.anova <- lm(model.tmp, data = all.dat.collapse)
  all.ancova.mods[[i]] <- tmp.anova
  ## Now estimate the f values
  p1 <- visreg(tmp.anova, "skewVal", "dif2PL", overlay=TRUE, gg=TRUE, alpha = 1) +coord_cartesian(xlim=c(-1, 6), ylim=c(0.6, 1)) +
    theme_bw() +
    scale_fill_grey() +
    scale_color_grey()
  all.visreg[[i]] <- p1
}
# out.plot <- ggpubr::ggarrange(all.visreg[[1]], all.visreg[[2]], all.visreg[[3]],
#                   all.visreg[[4]], all.visreg[[5]], all.visreg[[6]],legend = "bottom")
# ggsave(filename = "./reports/figureANCOR_twoWay1.png", plot = out.plot, dpi = 300, width = 9, height = 8)

## Now try to do this all in a single figure
all.out <- c("rseHurdRel", "rseGRMRel", "rseOmega", "rseAlpha") 
dat.tmp <- all.dat.collapse[,c("skewVal", "dif2PL", "iter", all.out)]
dat.tmp <- reshape2::melt(dat.tmp, id.vars=c("skewVal", "dif2PL", "iter"))
dat.tmp$variable <- plyr::revalue(dat.tmp$variable, c("rseHurdRel" = "Hurdle", "rseGRMRel" = "GRM", "rseOmega" = "Omega", "rseAlpha" = "Alpha"))
dat.tmp$variable <- factor(dat.tmp$variable, levels = c("Alpha", "Omega", "GRM", "Hurdle"))

#all.out <- c("Alpha", "Omega", "GRM", "Hurdle")

tmp2 <- ggplot(data = dat.tmp, aes(x=skewVal, y=value, group = dif2PL, fill = dif2PL, color = dif2PL)) +
  geom_smooth(method="lm", se = FALSE, size=2) +
  #geom_point() +
  theme_minimal() +
  scale_fill_grey() +
  scale_color_grey() + 
  facet_grid(. ~ variable) +
  xlab("Skew") +
  ylab("RMSE") +
  coord_cartesian(xlim=c(-1, 6), ylim = c(0, .4)) +
  guides(fill=guide_legend(title="2PL Difficulty"), color=guide_legend(title="2PL Difficulty")) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), legend.key.width = unit(1, "in"))

all.out <- c("estHurdleRel", "grmRel", "singleFactorOmegaT", "alpha", "trueRel") 
dat.tmp <- all.dat.collapse[,c("skewVal", "dif2PL", "iter", all.out)]
dat.tmp <- reshape2::melt(dat.tmp, id.vars=c("skewVal", "dif2PL", "iter"))
dat.tmp$variable <- plyr::revalue(dat.tmp$variable, c("estHurdleRel" = "Hurdle", "grmRel" = "GRM", "singleFactorOmegaT" = "Omega", "alpha" = "Alpha", "trueRel" = "True"))
dat.tmp$variable <- factor(dat.tmp$variable, levels = c("Alpha", "Omega", "GRM", "Hurdle", "True"))

tmp1 <- ggplot(data = dat.tmp, aes(x=skewVal, y=value, group = dif2PL, fill = dif2PL, color = dif2PL)) +
  geom_smooth(method="lm", se = FALSE, size=2) +
  #geom_point() +
  theme_minimal() +
  scale_fill_grey() +
  scale_color_grey() + 
  facet_grid(. ~ variable) +
  xlab("Skew") +
  ylab("Reliability") +
  coord_cartesian(xlim=c(-1, 6), ylim = c(0.5,1)) +
  guides(fill=guide_legend(title="2PL Difficulty"), color=guide_legend(title="2PL Difficulty")) +
  theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"), legend.key.width = unit(1, "in"))
tmp <- ggpubr::ggarrange(tmp1, tmp2, labels = "AUTO", common.legend = TRUE, legend = "bottom")

ggsave(filename = "./reports/figureANCOR_twoWay1.png", plot=tmp, width = 20, height = 8, dpi= 400, units="in", bg="white")
ggsave(filename = "./texOut/figures/figureANCOR_twoWay1.png", plot=tmp, width = 20, height = 8, dpi= 300, units="in", bg="white")

#out.plot <- ggpubr::ggarrange(p1, p2, p3, p4, p5, p6,legend = "bottom", nrow = 2, ncol=3)
#ggsave(filename = "./reports/figureANCOR_twoWay1.png", plot=out.plot, width = 6, height = 5, dpi=300, units="in")

## Print the main effect table here
library(dplyr)
library(kableExtra)
all.out <- c("Alpha", "Omega", "GRM", "Hurdle")
main.ef.table.round <- main.ef.table %>% dplyr::mutate(across(where(is.numeric), round, 2))
main.ef.table.round <- main.ef.table.round[,c("Parameter","rseAlpha", "rseOmega","rseGRMRel", "rseHurdRel")]
table_out <- kable(main.ef.table.round, "latex")
cat(table_out, file="./reports/mainEFTable.tex")

## Now do the two-way interaction vals here
two.way.ef.table <- all.ef.vals.iso[8:98,]
two.way.ef.table.3 <- two.way.ef.table[which(apply(two.way.ef.table[2:5], 1, function(x) sum(x > .01))>1),]
two.way.ef.table.3 <- two.way.ef.table.3 %>% dplyr::mutate(across(where(is.numeric), round, 2))
two.way.ef.table.3 <- two.way.ef.table.3[,c("Parameter","rseAlpha", "rseOmega","rseGRMRel", "rseHurdRel")]
table_out <- kable(two.way.ef.table.3, "latex")
cat(table_out, file="./reports/interactionEFTable.tex")
