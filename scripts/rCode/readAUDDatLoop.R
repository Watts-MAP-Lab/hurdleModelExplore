## This script will be used to read the AUD 87 data
in.dat <- read.csv("~/Documents/hurdleModelExplore/data/aud87/AUD-87_combined.csv")

## Load library(s)
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
source("~/Documents/hurdleModelExplore/scripts/rCode/hurdleFunctions.r")
library(mirt)
library(tidyverse)
library(MplusAutomation)
names(in.dat)
question.names <- names(in.dat)[3:89]
all.cats <- unique(strSplitMatrixReturn(question.names, "_")[,1])
i <- all.cats[1]
iso.dat <- grep(pattern = i, x = names(in.dat)[1:89])
iso.dat <- in.dat[,iso.dat]
## rm the na rows
iso.dat <- iso.dat[complete.cases(iso.dat),]
## Now make the IRTree data needed
## Make two matricies same dim as the iso.dat
mat.one <- iso.dat
mat.two <- iso.dat
mat.one[mat.one>1] <- 1
mat.two[mat.two==0] <- NA
setwd("./data/aud87/audMPlus/")
for(i in all.cats){
  iso.dat <- grep(pattern = i, x = names(in.dat)[1:89])
  iso.dat <- in.dat[,iso.dat]
  iso.dat <- iso.dat[complete.cases(iso.dat),]
  for(l in 2:6){
    mat.one <- iso.dat
    mat.one[mat.one>1] <- 1
    mat.two <- iso.dat
    mat.two[mat.two==0] <- NA
    mat.two[mat.two>=l] <- l
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
    p1 <- paste0("Presence BY ",paste0(binary.vars, collapse = " \n"))
    p1 <- gsub(pattern = "X1 ", replacement = "X1* ", x = p1)
    p2 <- paste0("Severity BY ",paste0(sev.var, collapse = "\n"))
    p2 <- gsub(pattern = sev.var[1], replacement = paste(sev.var[1], "*", sep=''), x = p2)
    p3 <- paste0("Usevariables = ", paste0(all.vars, collapse = "\n"))
    p4 <- paste0("Categorical = ", paste0(all.vars, collapse = "\n"))
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
    out.file.name <- paste(i,"_", l, ".inp", sep='')
    out.data.name <- paste(i,"_", l, ".txt", sep='')
    mplusModeler(test, dataout=out.data.name, modelout = out.file.name)
    runModels(target = out.file.name)
  }
}

## Now loop through these models and prep the data
tmp.dat <- list()
tmp.dat.scores <- list()
theta1 <- seq(-8,8,0.2) # Quadrature points for theta1
theta2 <- seq(-8,8,0.2) # Quadrature points for theta1
theta <- expand.grid(theta1, theta2)
names(theta) <- c("theta1", "theta2")
for(index in c(1:11)){
  i <- all.cats[index]
  all.dat <- NULL
  scores.dat <- list()
  for(l in 2:6){
    iso.dat <- grep(pattern = i, x = names(in.dat)[1:89])
    iso.dat <- in.dat[,iso.dat]
    iso.dat <- iso.dat[complete.cases(iso.dat),]
    mat.one <- iso.dat
    mat.one[mat.one>1] <- 1
    mat.two <- iso.dat
    mat.two[mat.two==0] <- NA
    mat.two[mat.two>=l] <- l
    #mat.two <- mat.two-l
    ## Now prep the mplus code
    all.dat <- cbind(mat.one, mat.two)
    colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
    in.file.name <- paste(i,"_", l, ".out", sep='')
    mod <- readModels(target=in.file.name)
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
    ## Now combine these all and attach them to output
    ## Fix column names here 
    all.params <- cbind(i, l, 1:length(a$est),rho,a_z$est, b_z$est, a$est, b_orgStand)
    prepVec <- c("subscale","capval","item","rho","est_z_discrim","est_z_diff","est_grm_discrim")
    grmDiffVec <- (dim(all.params)[2] - length(prepVec))
    grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
    prepVec <- c(prepVec, grmEstVec)
    all.params <- data.frame(all.params)
    colnames(all.params) <- prepVec
    if(l == 2){
      targ <- all.params
    }else{
      targ <- bind_rows(targ, all.params)
    }
    ## Now do the EAP factor score estimation here
    iso.dat$eap2PL_Hurdle <- NA
    iso.dat$se2PL_Hurdle <- NA
    iso.dat$eapGRM_Hurdle <- NA
    iso.dat$seGRM_Hurdle <- NA
    prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
    itemtrace <- trace.line.pts(a$est, b_orgStand, a_z$est, b_z$estStand, theta)
    for(respondent in 1:nrow(iso.dat)) {
      pattern <- iso.dat[respondent,1:dim(mat.one)[2]]
      pattern[pattern>l] <- l
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
    scores.dat[[l]] <- iso.dat
  }
  tmp.dat[[i]]$parameters <- targ
  tmp.dat[[i]]$scores <- iso.dat
  tmp.dat.scores[[i]] <- scores.dat
}
tmp.dat1 <- tmp.dat

## Now do the reverse windsor
for(i in all.cats){
  iso.dat <- grep(pattern = i, x = names(in.dat)[1:89])
  iso.dat <- in.dat[,iso.dat]
  iso.dat <- iso.dat[complete.cases(iso.dat),]
  for(l in 2:3){
    mat.one <- iso.dat
    mat.one[mat.one>1] <- 1
    mat.two <- iso.dat
    mat.two[mat.two==0] <- NA
    mat.two[mat.two<=l] <- 1
    ## Now reduce all l so the minimum is 1
    l.min <- l-1
    for(f in 1:ncol(mat.two)){
      index <- which(mat.two[,f] > l)
      mat.two[index,f] <- mat.two[index,f]-l.min
    }
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
    p1 <- paste0("Presence BY ",paste0(binary.vars, collapse = " \n"))
    p1 <- gsub(pattern = "X1 ", replacement = "X1* ", x = p1)
    p2 <- paste0("Severity BY ",paste0(sev.var, collapse = "\n"))
    p2 <- gsub(pattern = sev.var[1], replacement = paste(sev.var[1], "*", sep=''), x = p2)
    p3 <- paste0("Usevariables = ", paste0(all.vars, collapse = "\n"))
    p4 <- paste0("Categorical = ", paste0(all.vars, collapse = "\n"))
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
    out.file.name <- paste(i,"_", l, "_Min.inp", sep='')
    out.data.name <- paste(i,"_", l, "_Min.txt", sep='')
    mplusModeler(test, dataout=out.data.name, modelout = out.file.name)
    runModels(target = out.file.name)
  }
}
tmp.dat <- list()
tmp.dat.scores <- list()
theta1 <- seq(-8,8,0.2) # Quadrature points for theta1
theta2 <- seq(-8,8,0.2) # Quadrature points for theta1
theta <- expand.grid(theta1, theta2)
names(theta) <- c("theta1", "theta2")
for(index in c(1:11)){
  i <- all.cats[index]
  all.dat <- NULL
  scores.dat <- list()
  for(l in 2:3){
    iso.dat <- grep(pattern = i, x = names(in.dat)[1:89])
    iso.dat <- in.dat[,iso.dat]
    iso.dat <- iso.dat[complete.cases(iso.dat),]
    mat.two <- iso.dat
    mat.one <- iso.dat
    mat.one[mat.one>1] <- 1
    mat.two <- iso.dat
    mat.two[mat.two==0] <- NA
    mat.two[mat.two<l & mat.two>0] <- 1
    ## Now reduce all l so the minimum is 1
    mat.two <- mat.two - diff(c(1, min(mat.two, na.rm = TRUE)))
    #mat.two <- mat.two-l
    ## Now prep the mplus code
    all.dat <- cbind(mat.one, mat.two)
    colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
    in.file.name <- paste(i,"_", l, "_Min.out", sep='')
    mod <- readModels(target=in.file.name)
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
    ## Now combine these all and attach them to output
    ## Fix column names here 
    all.params <- cbind(i, l, 1:length(a$est),rho,a_z$est, b_z$est, a$est, b_orgStand)
    prepVec <- c("subscale","capval","item","rho","est_z_discrim","est_z_diff","est_grm_discrim")
    grmDiffVec <- (dim(all.params)[2] - length(prepVec))
    grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
    prepVec <- c(prepVec, grmEstVec)
    all.params <- data.frame(all.params)
    colnames(all.params) <- prepVec
    if(l == 2){
      targ <- all.params
    }else{
      targ <- bind_rows(targ, all.params)
    }
    ## Now do the EAP factor score estimation here
    iso.dat$eap2PL_Hurdle <- NA
    iso.dat$se2PL_Hurdle <- NA
    iso.dat$eapGRM_Hurdle <- NA
    iso.dat$seGRM_Hurdle <- NA
    prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
    itemtrace <- trace.line.pts(a$est, b_orgStand, a_z$est, b_z$estStand, theta)
  }
  tmp.dat[[i]]$parameters <- targ
  tmp.dat[[i]]$scores <- iso.dat
  tmp.dat.scores[[i]] <- scores.dat
}
tmp.dat2 <- tmp.dat

## Examine differences in parameter estimate values across method 1 & 2 here
iso.dat <- tmp.dat1[[1]]$parameters
iso.dat$method <- "depression"
iso.dat2 <- tmp.dat2[[1]]$parameters
iso.dat2$method <- "expansion"
## Combine these
iso.datA <- bind_rows(iso.dat, iso.dat2)
tmp <- merge(iso.dat, iso.dat2, by=c("subscale", "capval", "item"), suffixes = c("_depression", "_expansion"))
plot(tmp$est_z_diff_depression, tmp$est_z_diff_expansion)

## Now run some pure 1 dimensional GRM models
all.output <- NULL
for(index in c(1:11)){
  i <- all.cats[index]
  iso.dat <- grep(pattern = i, x = names(in.dat)[1:89])
  iso.dat <- in.dat[,iso.dat]
  iso.dat <- iso.dat[complete.cases(iso.dat),]
  mod.grm <- mirt::mirt(iso.dat, model=1, itemtype = "graded")
  ## Now get the AIC & BIC from the hurdle models for these data
  in.file.name <- paste(i,"_6.out", sep='')
  mod.hur <- readModels(target=in.file.name)
  ## Now prep this output
  output.summarie <- data.frame("SCALE" = i, "GRM_BIC" = mod.grm@Fit$BIC, "HUR_BIC" = mod.hur$summaries$BIC,
                                "GRM_AIC" = mod.grm@Fit$AIC, "HUR_AIC" = mod.hur$summaries$AIC)
  all.output <- bind_rows(all.output, output.summarie)
}
## Now plot these
all.output <- reshape2::melt(all.output, id.vars="SCALE")
all.output$TYPE <- strSplitMatrixReturn(charactersToSplit = all.output$variable, "_")[,2]
all.output$MOD <- strSplitMatrixReturn(charactersToSplit = all.output$variable, "_")[,1]
tmp <- ggplot(all.output, aes(x=MOD, y=value)) +
  geom_bar(stat="identity") +
  facet_grid(TYPE ~ SCALE, scales="free") +
  ylab("[A-B]IC") +
  coord_cartesian(ylim=c(5000, 30000)) +
  theme(axis.text.x = element_text(angle=35))
  
ggsave(filename = "~/Documents/hurdleModelExplore/data/GRMversusHurdleInformation.png",plot=tmp)
