## This script will be used to read the AUD 87 data
in.dat <- read.csv("./data/aud87/AUD-87_combined.csv")

## Load library(s)
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
source("./scripts/rCode/hurdleFunctions.r")
library(mirt)
library(tidyverse)
library(MplusAutomation)

## Organize the data for julia job submission
## Make sure to rm WD1 & WD2 as these are only found in half of the data
names(in.dat)
question.names <- names(in.dat)[3:89]
all.cats <- unique(strSplitMatrixReturn(question.names, "_")[,1])
out.dat <- list()
all.rel.dat <- NULL
iter.index <- 1
## Remove questions 1 & 2 form the WD subscales
for(i in all.cats){
  ## Isolate questions
  iso.dat <- grep(pattern = i, x = names(in.dat)[1:89])
  if(i == "WD"){
    iso.dat <- iso.dat[-c(1:2)]
  }
  iso.dat <- in.dat[,iso.dat]
  ## rm the na rows
  iso.dat <- iso.dat[complete.cases(iso.dat),]
  ## Now write these to file
  crosstab <- iso.dat %>%
    group_by(across(everything())) %>%
    summarise(Count = n(), .groups = "drop") %>%
    ungroup()
  out.file <- paste("./data/aud87/", i, "_audTab.csv", sep='')
  write.csv(crosstab, out.file, quote=F, row.names=F)
  out.dat[[iter.index]] <- iso.dat
  iter.index <- iter.index + 1
  ## Now prep the dat5a for mplus
  mat.one <- iso.dat
  mat.one[mat.one>1] <- 1
  mat.two <- iso.dat
  mat.two[mat.two==0] <- NA
  #mat.two <- mat.two-l
  ## Now prep the mplus code
  all.dat <- cbind(mat.one, mat.two)
  colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
  ## Now prep the mplus file
  first.dim <- ncol(mat.one)
  second.dim1 <- first.dim+1
  second.dim2 <- first.dim*2
  binary.vars <- colnames(all.dat)[1:first.dim]
  sev.var <- colnames(all.dat)[second.dim1:second.dim2]
  all.vars <- colnames(all.dat)
  model <- "
  F1 = 1-YYZ
  F2 = BBH-FNF
  COV = F1*F2
  "
  model <- gsub(x = model, pattern = "YYZ", replacement = first.dim)
  model <- gsub(x = model, pattern = "BBH", replacement = second.dim1)
  model <- gsub(x = model, pattern = "FNF", replacement = second.dim2)
  item.type.rep <- c(rep("2PL", ncol(mat.one)), rep("graded", ncol(mat.one)))
  sv1 <- mirt(data.frame(all.dat), model = model, itemtype = item.type.rep)
  mirt.coef <- coef(sv1, IRTpars=TRUE)
  ## First grab all of the a vals
  sus.vals <- 1:first.dim
  a_z <- unlist(lapply(mirt.coef[sus.vals], function(x) x[1]))
  b_z <- unlist(lapply(mirt.coef[sus.vals], function(x) x[3]))
  sev.vals <- grep(pattern = "Sev", x = names(mirt.coef))
  sev.vals <- second.dim1:second.dim2
  a <- unlist(lapply(mirt.coef[sev.vals], function(x) x[2]))
  b <- lapply(mirt.coef[sev.vals], function(x) x[-c(1:2)])
  b <- t(dplyr::bind_rows(b))
  #b <- t(apply(b, 1, function(x) sort(x, decreasing = FALSE)))
  rhoEst <- unique(mirt.coef$GroupPars["par","COV_21"])
  ## Now prep these all for output
  all.params <- cbind(i, 1:length(a),rhoEst,a_z, b_z, a, b)
  prepVec <- c("subscale","item","rho","est_z_discrim","est_z_diff","est_grm_discrim")
  grmDiffVec <- (dim(all.params)[2] - length(prepVec))
  grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
  prepVec <- c(prepVec, grmEstVec)
  all.params <- data.frame(all.params)
  colnames(all.params) <- prepVec
  ## Grab the reliability metrics here
  ## Now estimate the reliability of these data
  hurdle.rel <- hurdInfo(a = a, b = b, a_z = a_z, b_z = b_z, muVals = c(0,0), rhoVal = rhoEst, theta.grid = expand.grid(seq(-5, 5, .1), seq(-5, 5, .1)))
  ## Now do the GRM rel here
  mod <- mirt::mirt(data.frame(iso.dat), 1, itemtype = "graded")
  grmRel <-  1 / (1 + (1 / weighted.mean(testinfo(mod, Theta=seq(-5, 5, .1)), dnorm(seq(-5, 5, .1)))))
  ## Now do all of the split half assessments
  cor.mat <- psych::polychoric(iso.dat)
  split.half.rel <- psych::guttman(r = cor.mat$rho)
  omega.rel <- psych::omega(m = iso.dat, poly = TRUE, plot = FALSE, fm = "pc")
  cor.mat1 <- cor.mat$rho
  rel.all <- psych::reliability(cor.mat$rho, n.obs = nrow(all.dat), nfactors=3,plot=FALSE)
  omega_h <- rel.all$result.df[,"omega_h"]
  alpha <- rel.all$result.df[,"alpha"]
  omega_t <- rel.all$result.df[,"omega.tot"]
  Uni <- rel.all$result.df[,"Uni"]
  tau <- rel.all$result.df[,"tau"]
  cong <- rel.all$result.df[,"cong"]
  CFI <- rel.all$result.df[,"CFI"]
  ECV <- rel.all$result.df[,"ECV"]
  Beta <- rel.all$result.df[,"Beta"]
  EVR <- rel.all$result.df[,"EVR"]
  MAP <- rel.all$result.df[,"MAP"]
  skewVal <- as.numeric(psych::describe(rowSums(iso.dat))["skew"])
  kurtVal <- as.numeric(psych::describe(rowSums(iso.dat))["kurtosis"])
  
  split.half.rel <- suppressWarnings(psych::guttman(r = cor.mat$rho))
  lambda1Rel = split.half.rel$lambda.1
  lambda2Rel = split.half.rel$lambda.2
  lambda3Rel = split.half.rel$lambda.3
  lambda4Rel = split.half.rel$lambda.4
  lambda5Rel = split.half.rel$lambda.5
  lambda6Rel = split.half.rel$lambda.6
  omega.sem.val <- psych::omegaSem(m = cor.mat$rho, nfactors=1, plot = FALSE, n.obs = 15000)
  singleFactorOmegaT <- omega.sem.val$omegaSem$omega.tot
  ## Now do the same for only the >0 values
  iso.col <- grep(pattern = "Sev", x = colnames(rep_loop$mplusMat))
  mod <- mirt::mirt(data.frame(all.dat[,sev.vals]), 1, itemtype = "graded")
  grmRel_rmZeroOption <-  1 / (1 + (1 / weighted.mean(testinfo(mod, Theta=seq(-5, 5, .1)), dnorm(seq(-5, 5, .1)))))
  
  out.rel.df <- data.frame(hurdleRel = hurdle.rel$out.rel, grmRel = grmRel,omegaTRel = omega.rel$omega.tot, 
                           alpha = alpha, grmRel_rmZeroOption = grmRel_rmZeroOption, singleFactorOmegaT = singleFactorOmegaT,
                           scaleSkew = skewVal, scaleKurt = kurtVal)
  
  out.rel.df$scale <- i
  all.rel.dat <- bind_rows(all.rel.dat, out.rel.df)
}

## Now plot these
## Now grab the score values for all of these so I can make a histogram for each of these
m2 <- lapply(out.dat, function(x) table(rowSums(x)))
names(m2) <- names(all.cats)

## Now plot these
plot.dat1 <- reshape2::melt(all.rel.dat, id.vars=c("scale"))
#plot.dat1 <- plot.dat1[which(plot.dat1$variable %in% c("hurdleRel", "grmRel", "lambda3Rel", "omegaTRel")),]
plot.dat1 <- plot.dat1[-which(plot.dat1$variable %in% c("scaleSkew", "omegaTRel", "grmRel_rmZeroOption", "scaleKurt")),]
plot.dat1$variable <- plyr::revalue(plot.dat1$variable, c("alpha" = "Alpha", "hurdleRel" = "Hurdle", "grmRel" = "GRM", "singleFactorOmegaT" = "Omega"))
plot.dat1$scale <- factor(plot.dat1$scale, levels = unique(plot.dat1$scale[order(plot.dat1$scale)]))
plot.dat1$variable <- factor(plot.dat1$variable, levels =c("Alpha", "Omega", "GRM", "Hurdle"))

p1 <- ggplot(plot.dat1, aes(x=variable, y = value)) +
  geom_bar(stat="identity") +
  facet_wrap(scale ~ .) +
  #theme(axis.text.x = element_text(angle=35, vjust=.25)) +
  coord_cartesian(ylim=c(.7, 1)) +
  ylab("Rel") +
  xlab("Metric") +
  theme_bw()# +
  #theme(axis.text.x = element_text(angle=90, vjust=.9))
plot.dat2 <- reshape2::melt(m2, id.vars=c("scale"))
for(i in 1:length(all.cats)){
  rep.val <- all.cats[i]
  plot.dat2[which(plot.dat2[,3]==i),3] <- rep.val
}
plot.dat2$scale <- factor(plot.dat2$L1, levels = unique(plot.dat1$scale[order(plot.dat1$scale)]))
p2 <- ggplot(plot.dat2, aes(x=Var1, y = value)) +
  geom_bar(stat="identity") +
  facet_wrap(L1 ~ .) +
  theme_bw() + 
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks = 0:10) +
  xlab("Sum score") +
  ylab("n")
library(ggpubr)
out.plot <- ggarrange(p2, p1, labels = "AUTO")
ggsave(filename = "./reports/figure7_NESARCData.png", plot = out.plot, dpi=300, width = 13, height=7, units="in")

## Now make a loop to pair all of these plots together within domain
fig.list <- list()
fig.iter <- 1
for(i in sort(unique(plot.dat1$scale))){
  ## Isolate all data across pairs
  vals.rel <- plot.dat1[which(plot.dat1$scale == i),]
  vals.ss <- plot.dat2[which(plot.dat2$L1==i),]
  ## Now make a plot of these two together
  p1 <- ggplot(vals.rel, aes(x=variable, y = value)) +
    geom_bar(stat="identity") +
    #theme(axis.text.x = element_text(angle=35, vjust=.5)) +
    coord_cartesian(ylim=c(.65, 1)) +
    scale_y_continuous(expand = c(0,0)) +
    ylab("") +
    xlab("") +
    theme_minimal() +
    theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"),panel.grid.major.x = element_blank())
  p2 <- ggplot(vals.ss, aes(x=Var1, y = value)) +
    geom_bar(stat="identity") +
    theme_minimal() + 
    coord_cartesian(xlim=c(0,10), ylim=c(0,1250)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks = 0:10) +
    xlab("") +
    ylab("") +
    theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"),panel.grid.major.x = element_blank())
  out.p <- ggarrange(p2, p1, ncol = 2)
  out.p <- annotate_figure(out.p, fig.lab = i, fig.lab.pos = "top.left", fig.lab.face = "bold")
  fig.list[[fig.iter]] <- out.p
  fig.iter <- fig.iter + 1
}

out.plot<- ggarrange(fig.list[[1]], fig.list[[2]], fig.list[[3]], fig.list[[4]], fig.list[[5]],
                     fig.list[[6]], fig.list[[7]], fig.list[[8]], fig.list[[9]], fig.list[[10]],
                     fig.list[[11]],ncol = 3, nrow=4)
ggsave(filename = "./reports/figure7_NESARCData.png", plot = out.plot, dpi=300, width = 24, height=9, units="in")




## Now write this as a latex object
library(dplyr)
library(kableExtra)
out.table <- data.frame(Subdomain=all.cats, 
                        totalQuestion = unlist(lapply(out.dat, function(x) dim(x)[2])),
                        meanVals = unlist(lapply(out.dat, function(x) psych::describe(rowSums(x))$mean)),
                        sdVals = unlist(lapply(out.dat, function(x) psych::describe(rowSums(x))$sd)),
                        skewVals = unlist(lapply(out.dat, function(x) psych::describe(rowSums(x))$skew)))
out.table$skewVals <- round(out.table$skewVals, 2)
out.table$meanVals <- round(out.table$meanVals, 2)
out.table$sdVals <- round(out.table$sdVals, 2)

table_out <- kable(out.table, "latex")
cat(table_out, file="./reports/audOut.tex")

out.table <- data.frame(Subdomain=all.cats, 
                        totalQuestion = unlist(lapply(out.dat, function(x) dim(x)[2])),
                        meanVals = unlist(lapply(out.dat, function(x) psych::describe(rowSums(x))$mean)),
                        sdVals = unlist(lapply(out.dat, function(x) psych::describe(rowSums(x))$sd)),
                        skewVals = unlist(lapply(out.dat, function(x) psych::describe(rowSums(x))$skew)))
descrip.table <- data.frame(Subdomain = c("LL","CD","PP","CR","RI","GU","HZ","TL","SI","WD","TS"),
                            Description = c("Larger/longer amount than intended",
                                            "Cut down activities to drink",
                                            "Physical/psychological harm",
                                            "Craving",
                                            "Role interference",
                                            "Gave up activities to drink",
                                            "Hazardous use",
                                            "Tolerance",
                                            "Social/interpersonal harm",
                                            "Withdrawal",
                                            "Time Spent"))
out.table <- merge(descrip.table, out.table)
out.table2 <- merge(out.table, all.rel.dat, by.x=c("Subdomain", "skewVals"), by.y=c("scale", "scaleSkew"))
## Now order the outcome
out.table2 <- out.table2[,c("Subdomain", "Description", "totalQuestion", "meanVals", "sdVals", "skewVals","hurdleRel", "grmRel", "singleFactorOmegaT", "alpha")]
out.table2$skewVals <- round(out.table2$skewVals, 2)
out.table2$meanVals <- round(out.table2$meanVals, 2)
out.table2$sdVals <- round(out.table2$sdVals, 2)
out.table2$hurdleRel <- round(out.table2$hurdleRel, 2)
out.table2$grmRel <- round(out.table2$grmRel, 2)
out.table2$singleFactorOmegaT <- round(out.table2$singleFactorOmegaT, 2)
out.table2$alpha <- round(out.table2$alpha, 2)

table_out <- kable(out.table2, "latex")
cat(table_out, file="./reports/audOut.tex")


## Now prep the table with the Hurdle versus GRM model estimates for the RI subscale
iso.dat <- grep(pattern = "RI", x = names(in.dat)[1:89])
iso.dat <- in.dat[,iso.dat]
## rm the na rows
iso.dat <- iso.dat[complete.cases(iso.dat),]
## Now write these to file
crosstab <- iso.dat %>%
  group_by(across(everything())) %>%
  summarise(Count = n(), .groups = "drop") %>%
  ungroup()
## Now prep the dat5a for mplus
mat.one <- iso.dat
mat.one[mat.one>1] <- 1
mat.two <- iso.dat
mat.two[mat.two==0] <- NA
#mat.two <- mat.two-l
## Now prep the mplus code
all.dat <- cbind(mat.one, mat.two)
colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
## Now prep the mplus file
first.dim <- ncol(mat.one)
second.dim1 <- first.dim+1
second.dim2 <- first.dim*2
binary.vars <- colnames(all.dat)[1:first.dim]
sev.var <- colnames(all.dat)[second.dim1:second.dim2]
all.vars <- colnames(all.dat)
model <- "
  F1 = 1-YYZ
  F2 = BBH-FNF
  COV = F1*F2
  "
model <- gsub(x = model, pattern = "YYZ", replacement = first.dim)
model <- gsub(x = model, pattern = "BBH", replacement = second.dim1)
model <- gsub(x = model, pattern = "FNF", replacement = second.dim2)
item.type.rep <- c(rep("2PL", ncol(mat.one)), rep("graded", ncol(mat.one)))
sv1 <- mirt(data.frame(all.dat), model = model, itemtype = item.type.rep)
mirt.coef <- coef(sv1, IRTpars=TRUE)
## First grab all of the a vals
sus.vals <- 1:first.dim
a_z <- unlist(lapply(mirt.coef[sus.vals], function(x) x[1]))
b_z <- unlist(lapply(mirt.coef[sus.vals], function(x) x[3]))
sev.vals <- grep(pattern = "Sev", x = names(mirt.coef))
sev.vals <- second.dim1:second.dim2
a <- unlist(lapply(mirt.coef[sev.vals], function(x) x[2]))
b <- lapply(mirt.coef[sev.vals], function(x) x[-c(1:2)])
b <- t(dplyr::bind_rows(b))
#b <- t(apply(b, 1, function(x) sort(x, decreasing = FALSE)))
rhoEst <- unique(mirt.coef$GroupPars["par","COV_21"])
## Now prep these all for output
all.params <- cbind("RI", 1:length(a),rhoEst,a_z, b_z, a, b)
prepVec <- c("subscale","item","rho","est_z_discrim","est_z_diff","est_grm_discrim")
grmDiffVec <- (dim(all.params)[2] - length(prepVec))
grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
prepVec <- c(prepVec, grmEstVec)
all.params <- data.frame(all.params)
colnames(all.params) <- prepVec
mod <- mirt::mirt(data.frame(iso.dat), 1, itemtype = "graded")
uni_grm <- coef(mod, IRTpars=TRUE)
uni_coef <- rbind(uni_grm$RI_1, uni_grm$RI_2,uni_grm$RI_3,uni_grm$RI_4, uni_grm$RI_5, uni_grm$RI_6)
all.params <- cbind(all.params, uni_coef)


## Now write this data to latex
out.table2 <- all.params[,-1]
rownames(out.table2) <- NULL
out.table2$rho <- as.numeric(out.table2$rho)
out.table2$est_z_discrim <- as.numeric(out.table2$est_z_discrim)
out.table2$est_z_diff <- as.numeric(out.table2$est_z_diff)
out.table2$est_grm_discrim <- as.numeric(out.table2$est_grm_discrim)
out.table2$est_grm_diff_1 <- as.numeric(out.table2$est_grm_diff_1)
out.table2$est_grm_diff_2 <- as.numeric(out.table2$est_grm_diff_2)
out.table2$est_grm_diff_3 <- as.numeric(out.table2$est_grm_diff_3)

out.table2 <- out.table2 %>%
  mutate(across(where(is.numeric), round, digits = 2))
table_out <- kable(out.table2, "latex")
cat(table_out, file="./reports/suppTableRICoef.tex")
