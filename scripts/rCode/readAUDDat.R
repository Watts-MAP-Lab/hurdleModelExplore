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
question.names <- names(in.dat)[2:88]
all.cats <- unique(strSplitMatrixReturn(question.names, "_")[,1])

for(i in all.cats){
  ## Isolate questions
  iso.dat <- grep(pattern = i, x = names(in.dat)[1:88])
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
}

i <- all.cats[1]
iso.dat <- grep(pattern = i, x = names(in.dat)[1:88])
iso.dat <- in.dat[,iso.dat]
## rm the na rows
iso.dat <- iso.dat[complete.cases(iso.dat),]
## Now make the IRTree data needed
## Make two matricies same dim as the iso.dat
mat.one <- iso.dat
mat.two <- iso.dat
mat.one[mat.one>1] <- 1
mat.two[mat.two==0] <- NA
## Now bind
all.dat <- cbind(mat.one, mat.two)
colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
## Now prep the mplus file
first.dim <- ncol(mat.one)
second.dim1 <- first.dim+1
second.dim2 <- first.dim*2
binary.vars <- colnames(all.dat)[1:first.dim]
sev.var <- colnames(all.dat)[second.dim1:second.dim2]
all.vars <- colnames(all.dat)
p1 <- paste0("Presence BY ",paste0(binary.vars, collapse = " "))
p1 <- gsub(pattern = "X1", replacement = "X1*", x = p1)
p2 <- paste0("Severity BY ",paste0(sev.var, collapse = " "))
p2 <- gsub(pattern = sev.var[1], replacement = paste(sev.var[1], "*", sep=''), x = p2)
p3 <- paste0("Usevariables = ", paste0(all.vars, collapse = " "))
p4 <- paste0("Categorical = ", paste0(all.vars, collapse = " "))
all.dat <- data.frame(all.dat)
test <- MplusAutomation::mplusObject(TITLE = "testIRTree Model",
                                     MODEL = paste0(p1,"; \n",
                                                    p2, "; \n
                                           Presence@1; Severity@1; \n
                                           [Presence@0]; [Severity@0];"),
                                     rdata = all.dat,
                                     usevariables = all.vars,
                                     VARIABLE = paste0(p3, "; \n",p4,";"),
                                     ANALYSIS = "  Estimator = ML; \n Link = Logit; \n Integration = GAUSSHERMITE(15);")
mplusModeler(test, dataout="test.txt", modelout = "test.inp")
#runModels()

## Now collapse the upper most response category by 1
mat.one <- iso.dat
mat.two <- iso.dat
mat.one[mat.one>1] <- 1
mat.two[mat.two==0] <- NA
mat.two[mat.two==5] <- 4
## Now bind
all.dat <- cbind(mat.one, mat.two)
colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
## Now prep the mplus file
first.dim <- ncol(mat.one)
second.dim1 <- first.dim+1
second.dim2 <- first.dim*2
binary.vars <- colnames(all.dat)[1:first.dim]
sev.var <- colnames(all.dat)[second.dim1:second.dim2]
all.vars <- colnames(all.dat)
p1 <- paste0("Presence BY ",paste0(binary.vars, collapse = " "))
p1 <- gsub(pattern = "X1", replacement = "X1*", x = p1)
p2 <- paste0("Severity BY ",paste0(sev.var, collapse = " "))
p2 <- gsub(pattern = sev.var[1], replacement = paste(sev.var[1], "*", sep=''), x = p2)
p3 <- paste0("Usevariables = ", paste0(all.vars, collapse = " "))
p4 <- paste0("Categorical = ", paste0(all.vars, collapse = " "))
all.dat <- data.frame(all.dat)
test <- MplusAutomation::mplusObject(TITLE = "testIRTree Model",
                                     MODEL = paste0(p1,"; \n",
                                                    p2, "; \n
                                           Presence@1; Severity@1; \n
                                           [Presence@0]; [Severity@0];"),
                                     rdata = all.dat,
                                     usevariables = all.vars,
                                     VARIABLE = paste0(p3, "; \n",p4,";"),
                                     ANALYSIS = "  Estimator = ML; \n Link = Logit; \n Integration = GAUSSHERMITE(15);")
mplusModeler(test, dataout="test1.txt", modelout = "test1.inp")
#runModels()

## Now by two
mat.one <- iso.dat
mat.two <- iso.dat
mat.one[mat.one>1] <- 1
mat.two[mat.two==0] <- NA
mat.two[mat.two==5] <- 4
mat.two[mat.two==4] <- 3

## Now bind
all.dat <- cbind(mat.one, mat.two)
colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
## Now prep the mplus file
first.dim <- ncol(mat.one)
second.dim1 <- first.dim+1
second.dim2 <- first.dim*2
binary.vars <- colnames(all.dat)[1:first.dim]
sev.var <- colnames(all.dat)[second.dim1:second.dim2]
all.vars <- colnames(all.dat)
p1 <- paste0("Presence BY ",paste0(binary.vars, collapse = " "))
p1 <- gsub(pattern = "X1", replacement = "X1*", x = p1)
p2 <- paste0("Severity BY ",paste0(sev.var, collapse = " "))
p2 <- gsub(pattern = sev.var[1], replacement = paste(sev.var[1], "*", sep=''), x = p2)
p3 <- paste0("Usevariables = ", paste0(all.vars, collapse = " "))
p4 <- paste0("Categorical = ", paste0(all.vars, collapse = " "))
all.dat <- data.frame(all.dat)
test <- MplusAutomation::mplusObject(TITLE = "testIRTree Model",
                                     MODEL = paste0(p1,"; \n",
                                                    p2, "; \n
                                           Presence@1; Severity@1; \n
                                           [Presence@0]; [Severity@0];"),
                                     rdata = all.dat,
                                     usevariables = all.vars,
                                     VARIABLE = paste0(p3, "; \n",p4,";"),
                                     ANALYSIS = "  Estimator = ML; \n Link = Logit; \n Integration = GAUSSHERMITE(15);")
mplusModeler(test, dataout="test2.txt", modelout = "test2.inp")
#runModels()

## Now by three
mat.one <- iso.dat
mat.two <- iso.dat
mat.one[mat.one>1] <- 1
mat.two[mat.two==0] <- NA
mat.two[mat.two==5] <- 4
mat.two[mat.two==4] <- 3
mat.two[mat.two==3] <- 2

## Now bind
all.dat <- cbind(mat.one, mat.two)
colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
## Now prep the mplus file
first.dim <- ncol(mat.one)
second.dim1 <- first.dim+1
second.dim2 <- first.dim*2
binary.vars <- colnames(all.dat)[1:first.dim]
sev.var <- colnames(all.dat)[second.dim1:second.dim2]
all.vars <- colnames(all.dat)
p1 <- paste0("Presence BY ",paste0(binary.vars, collapse = " "))
p1 <- gsub(pattern = "X1", replacement = "X1*", x = p1)
p2 <- paste0("Severity BY ",paste0(sev.var, collapse = " "))
p2 <- gsub(pattern = sev.var[1], replacement = paste(sev.var[1], "*", sep=''), x = p2)
p3 <- paste0("Usevariables = ", paste0(all.vars, collapse = " "))
p4 <- paste0("Categorical = ", paste0(all.vars, collapse = " "))
all.dat <- data.frame(all.dat)
test <- MplusAutomation::mplusObject(TITLE = "testIRTree Model",
                                     MODEL = paste0(p1,"; \n",
                                                    p2, "; \n
                                           Presence@1; Severity@1; \n
                                           [Presence@0]; [Severity@0];"),
                                     rdata = all.dat,
                                     usevariables = all.vars,
                                     VARIABLE = paste0(p3, "; \n",p4,";"),
                                     ANALYSIS = "  Estimator = ML; \n Link = Logit; \n Integration = GAUSSHERMITE(15);")
mplusModeler(test, dataout="test3.txt", modelout = "test3.inp")
runModels()
mods <- readModels()
## Now try a julia version of this -- the negative covarainces are somewhat upsetting me?
## Perhaps because the 0 inflation is not that bad?
mod <- readModels(target = "test.out")
## Now obtain estimated a, b, a_z, and b_z parameters
a_z <- mod$parameters$unstandardized
a_z <- a_z[grep(pattern = "PRESENCE.BY", x = a_z$paramHeader),]
b_z <- mod$parameters$unstandardized
b_z <- b_z[grep(pattern = "Thresholds", x = b_z$paramHeader),]
b_z$itemVal <- strSplitMatrixReturn(b_z$param, "\\$")[,1]
b_z <- b_z[which(b_z$itemVal %in% a_z$param),]
b_z$estStand <- b_z$est / a_z$est
a <- mod$parameters$unstandardized
a <- a[grep(pattern = "SEVERITY.BY", x = a$paramHeader),]
b <- mod$parameters$unstandardized
b <- b[grep(pattern = "Thresholds", x = b$paramHeader),]
b$itemVal <- strSplitMatrixReturn(b$param, "\\$")[,1]
b$threshVal <- strSplitMatrixReturn(b$param, "\\$")[,2]
b <- b[which(!b$itemVal %in% a_z$param),]
## Now translate the b into proper format
b_org <- pivot_wider(b[,c("est", "itemVal", "threshVal")], id_cols = "itemVal", names_from = "threshVal", values_from = "est")[,-1]
b_org <- data.frame(b_org)
## Now convert these difficulty metrics into the standard normal scale
b_orgStand <- apply(b_org, 2, function(x) x / a$est)
rho <- mod$parameters$unstandardized$est[which(mod$parameters$unstandardized$paramHeader=="SEVERITY.WITH")]
## Attach these to data out
theta1 <- seq(-8,8,0.2) # Quadrature points for theta1
theta2 <- seq(-8,8,0.2) # Quadrature points for theta1
theta <- expand.grid(theta1, theta2)
names(theta) <- c("theta1", "theta2")
prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
itemtrace <- trace.line.pts(a$est, b_orgStand, a_z$est, b_z$estStand, theta)
iso.dat$eap2PL_Hurdle <- NA
iso.dat$se2PL_Hurdle <- NA
iso.dat$eapGRM_Hurdle <- NA
iso.dat$seGRM_Hurdle <- NA
for(respondent in 1:nrow(iso.dat)) {
  pattern <- iso.dat[respondent,1:6]
  qpoints <- theta
  lhood <- score(pattern, itemtrace, qpoints)
  
  eap2PL_Hurdle <- sum(lhood*prior*qpoints$theta1)/sum(lhood*prior)
  se2PL_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta1-eap2PL_Hurdle)^2)/sum(lhood*prior))
  iso.dat$eap2PL_Hurdle[respondent] <- eap2PL_Hurdle
  iso.dat$se2PL_Hurdle[respondent] <- se2PL_Hurdle
  
  eapGRM_Hurdle <- sum(lhood*prior*qpoints$theta2)/sum(lhood*prior)
  seGRM_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta2-eapGRM_Hurdle)^2)/sum(lhood*prior))
  iso.dat$eapGRM_Hurdle[respondent] <- eapGRM_Hurdle
  iso.dat$seGRM_Hurdle[respondent] <- seGRM_Hurdle
}
