## tabula rasa
rm(list=ls())

## First go ahead and load library(s) I might need
library(tidyverse)
library(doParallel)
library(mirt)
library(foreach)
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
source("./scripts/rCode/hurdleFunctions.r")


## Make a function to return the hurdle estimates
return_Mod_Params <- function(juliaOutput, dataIn){
    ## Initialize all of the output
    ## First create the 2PL data frame with discrim and diff values
    nitems = dim(dataIn)[2]
    irt2PLParam <- matrix(NA, nrow=nitems, ncol=2)
    ## Now do the GRM portion
    ncatgrm <- diff(range(dataIn))
    ## First idenitfy the total number of difficult params needed
    irtGRMParam <- matrix(NA, nrow=dim(dataIn)[2], ncol=ncatgrm)

    ## Clean up the Julia output here
    list_index <- length(juliaOutput)
    out_params <- juliaOutput[list_index]
    out_params <- strsplit(out_params, " ")[[1]]
    out_params <- gsub(pattern = "]", replacement="", x = out_params)
    out_params <- gsub(pattern = "\\[|\\]", replacement="", x = out_params)
    out_params <- gsub(pattern = ",", replacement="", x = out_params)
    out_params <- as.numeric(out_params)

    nitemsgrm = dim(dataIn)[2]
    paramGRM = nitems*ncatgrm
    param2PL = nitems*2
    nParmsPerItemGRM = ncatgrm
    nParmsPerItem2PL = 2
    nitemsgrm = nitems
    for (j in 1:nitems) {
        irtGRMParam[j,1] <- out_params[(j-1)*nParmsPerItemGRM + 1]
        irt2PLParam[j,1] <- out_params[nitemsgrm*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 1]
        irt2PLParam[j,2] <- out_params[nitemsgrm*nParmsPerItemGRM + (j-1)*nParmsPerItem2PL + 2]
        for (k in 1:(ncatgrm-1)){
            irtGRMParam[j,k+1] <- out_params[(j-1)*nParmsPerItemGRM + 1 + k]
        }
    }

    ## Now exponentiate the rho estimate
    rho_nexp = out_params[length(out_params)]
    rho = exp(rho_nexp) / (1 + exp(rho_nexp))


    ## Now reaturn these
    out_params = list(grmParams = irtGRMParam, irt2PLParams = irt2PLParam, rho = rho)
    return(out_params)
}

## This script will be used to examine the CBCL data
in_dat <- read.csv("./data/ABCD_CBCL_BL.csv")
in_dict <- readxl::read_xlsx("./data/CBCLitems.xlsx")

## AB == Aggressive Behavior
## AD == Anxious/Depressed
## AP == Attention Problems
## RB == Rule-Breaking Behavior
## SC == Somatic Complaints
## SP == Social Problems 
## TP == Thought Problems 
## WD == Withdrawn/Depressed
## Other == Other 


## Identify all possible scale values
unique_scales <- unique(in_dict$subscale)
all_item_position <- lapply(unique_scales, function(x) which(in_dict$subscale == x))
## Now make a list of lists with all of the subscales and which items map onto each
all_data <- list()
for(i in 1:length(all_item_position)){
  ## Grab the column values in the data
  index <- all_item_position[[i]]
  index <- index + 2
  ## Now grab the resp data
  rep_vals <- in_dat[,index]
  ## Drop any columns with less than 10 response endorsments
  rm.index <- which(apply(rep_vals, 2, function(x) length(which(x==2)))<15)
  if(length(rm.index > 0)){
    rep_vals <- rep_vals[,-c(rm.index)]
  }
  ## Now cross tab these
  crosstab <- rep_vals %>%
    group_by(across(everything())) %>%
    summarise(Count = n(), .groups = "drop") %>%
    ungroup()
  out_data <- list(rep_vals = rep_vals, crosstab = crosstab)
  all_data[[i]] <- out_data
}

## Now run one of these through the hurdle model in julia
cl <- makeCluster(9)
registerDoParallel(cl)
all.mods <- foreach(i = 1:length(all_item_position), .packages = c("MplusAutomation")) %dopar%{
  in_resp <- all_data[[i]]$rep_vals
  in_resp <- in_resp[complete.cases(in_resp),]
  in_tabs <- all_data[[i]]$crosstab
  in_tabs <- in_tabs[complete.cases(in_tabs),]
  ## Now prep the dat5a for mplus
  mat.one <- in_resp
  mat.one[mat.one>1] <- 1
  mat.two <- in_resp
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
  out.file.name <- paste("./data/cbclModels/",i, ".inp", sep='')
  out.data.name <- paste("./data/cbclModels/",i, ".txt", sep='')
  mplusModeler(test, dataout=out.data.name, modelout = out.file.name, run = TRUE)
  in.file.name <- paste("./data/cbclModels/",i, ".out", sep='')
  mod <- readModels(target = in.file.name)
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
  b_org <- tidyr::pivot_wider(b[,c("est", "itemVal", "threshVal")], id_cols = "itemVal", names_from = "threshVal", values_from = "est")[,-1]
  b_org <- data.frame(b_org)
  ## Now convert these difficulty metrics into the standard normal scale
  b_orgStand <- apply(b_org, 2, function(x) x / a$est)
  rho <- mod$parameters$unstandardized$est[which(mod$parameters$unstandardized$paramHeader=="SEVERITY.WITH")]
  ## Now prep these all for output
  all.params <- cbind(i, 1:length(a$est),rho,a_z$est, b_z$est, a$est, b_orgStand)
  prepVec <- c("subscale","item","rho","est_z_discrim","est_z_diff","est_grm_discrim")
  grmDiffVec <- (dim(all.params)[2] - length(prepVec))
  grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
  prepVec <- c(prepVec, grmEstVec)
  all.params <- data.frame(all.params)
  colnames(all.params) <- prepVec
  all.params
  ## Now do the factor scores here
  ## Now do the EAP factor score estimation here
  in_resp$eap2PL_Hurdle <- NA
  in_resp$se2PL_Hurdle <- NA
  in_resp$eapGRM_Hurdle <- NA
  in_resp$seGRM_Hurdle <- NA
  theta1 <- seq(-8,8,0.2) # Quadrature points for theta1
  theta2 <- seq(-8,8,0.2) # Quadrature points for theta1
  theta <- expand.grid(theta1, theta2)
  names(theta) <- c("theta1", "theta2")
  prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
  source("./scripts/rCode/hurdleFunctions.r")
  itemtrace <- trace.line.pts(a$est, b_orgStand, a_z$est, b_z$estStand, theta)
  for(respondent in 1:nrow(in_resp)) {
    pattern <- in_resp[respondent,1:dim(mat.one)[2]]
    qpoints <- theta
    lhood <- score(pattern, itemtrace, qpoints)
    
    eap2PL_Hurdle <- sum(lhood*prior*qpoints$theta1)/sum(lhood*prior)
    se2PL_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta1-eap2PL_Hurdle)^2)/sum(lhood*prior))
    in_resp$eap2PL_Hurdle[respondent] <- eap2PL_Hurdle
    in_resp$se2PL_Hurdle[respondent] <- se2PL_Hurdle
    
    eapGRM_Hurdle <- sum(lhood*prior*qpoints$theta2)/sum(lhood*prior)
    seGRM_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta2-eapGRM_Hurdle)^2)/sum(lhood*prior))
    in_resp$eapGRM_Hurdle[respondent] <- eapGRM_Hurdle
    in_resp$seGRM_Hurdle[respondent] <- seGRM_Hurdle
  }
  out.list <- list(params = all.params, factor_scores = in_resp)
  out.list
}
stopCluster(cl)
names(all.mods) <- unique_scales

## Now do the same using MIRT
cl <- makeCluster(9)
registerDoParallel(cl)
all.modsMIRT <- foreach(i = 1:length(all_item_position), .packages = c("mirt")) %dopar%{
  in_resp <- all_data[[i]]$rep_vals
  in_resp <- in_resp[complete.cases(in_resp),]
  in_tabs <- all_data[[i]]$crosstab
  in_tabs <- in_tabs[complete.cases(in_tabs),]
  ## Now prep the dat5a for mplus
  mat.one <- in_resp
  mat.one[mat.one>1] <- 1
  mat.two <- in_resp
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
  mod <- mirt::mirt(data.frame(in_resp), 1, itemtype = "graded")
  grmRel <-  1 / (1 + (1 / weighted.mean(testinfo(mod, Theta=seq(-5, 5, .1)), dnorm(seq(-5, 5, .1)))))
  ## Now do all of the split half assessments
  cor.mat <- psych::polychoric(in_resp)
  split.half.rel <- psych::guttman(r = cor.mat$rho)
  ## Make sure we have all of the other reliability estimates here too
  omega.rel <- psych::omega(m = in_resp, poly = TRUE, plot = FALSE, fm = "pc")
  rel.all <- psych::reliability(cor.mat$rho, n.obs = nrow(all.dat), nfactors=3,plot=FALSE)
  skewVal <- as.numeric(psych::describe(rowSums(in_resp))["skew"])
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
  #iso.col <- grep(pattern = "Sev", x = colnames(rep_loop$mplusMat))
  mod <- mirt::mirt(data.frame(all.dat[,sev.vals]), 1, itemtype = "graded")
  grmRel_rmZeroOption <-  1 / (1 + (1 / weighted.mean(testinfo(mod, Theta=seq(-5, 5, .1)), dnorm(seq(-5, 5, .1)))))
  
  out.rel.df <- data.frame(hurdleRel = hurdle.rel$out.rel, grmRel = grmRel,omegaTRel = omega.rel$omega.tot, 
                           alpha = rel.all$result.df[,"alpha"], grmRel_rmZeroOption = grmRel_rmZeroOption, singleFactorOmegaT = singleFactorOmegaT,
                           scaleSkew= skewVal)
  ## Now do the factor scores here
  ## Now do the EAP factor score estimation here
  in_resp$eap2PL_Hurdle <- NA
  in_resp$se2PL_Hurdle <- NA
  in_resp$eapGRM_Hurdle <- NA
  in_resp$seGRM_Hurdle <- NA
  theta1 <- seq(-6,6,0.2) # Quadrature points for theta1
  theta2 <- seq(-6,6,0.2) # Quadrature points for theta1
  theta <- expand.grid(theta1, theta2)
  names(theta) <- c("theta1", "theta2")
  prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rhoEst, rhoEst, 1), 2, 2, byrow = T))
  source("./scripts/rCode/hurdleFunctions.r")
  itemtrace <- trace.line.pts(a, b, a_z, b_z, theta)
  for(respondent in 1:nrow(in_resp)) {
    pattern <- in_resp[respondent,1:dim(mat.one)[2]]
    qpoints <- theta
    lhood <- score(pattern, itemtrace, qpoints)
    
    eap2PL_Hurdle <- sum(lhood*prior*qpoints$theta1)/sum(lhood*prior)
    se2PL_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta1-eap2PL_Hurdle)^2)/sum(lhood*prior))
    in_resp$eap2PL_Hurdle[respondent] <- eap2PL_Hurdle
    in_resp$se2PL_Hurdle[respondent] <- se2PL_Hurdle
    
    eapGRM_Hurdle <- sum(lhood*prior*qpoints$theta2)/sum(lhood*prior)
    seGRM_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta2-eapGRM_Hurdle)^2)/sum(lhood*prior))
    in_resp$eapGRM_Hurdle[respondent] <- eapGRM_Hurdle
    in_resp$seGRM_Hurdle[respondent] <- seGRM_Hurdle
  }
  out.list <- list(params = all.params, factor_scores = in_resp, rel_metrics = out.rel.df)
  out.list
}
stopCluster(cl)
names(all.modsMIRT) <- unique_scales

## Now merge mirt and mplus
m1 <- bind_rows(lapply(all.mods, function(x) x$params))
m2 <- bind_rows(lapply(all.modsMIRT, function(x) x$params))

merged.dat1 <- merge(m1, m2, by=c("subscale", "item"), suffixes = c("_mplus", "_mirt"))

## Now do the factor scores?
# m1 <- bind_rows(lapply(all.mods, function(x) x$factor_scores))
# m2 <- bind_rows(lapply(all.modsMIRT, function(x) x$factor_scores))
# 
# merged.dat1 <- merge(m1, m2, by=c("subscale", "item"), suffixes = c("_mplus", "_mirt"))


## Write these data for Ashley
saveRDS(object = all.mods, file = "./data/forAshleyCBCL.RDS")
saveRDS(object = all.modsMIRT, file = "./data/forAshleyCBCL.RDS")
#all.modsMIRT <- readRDS("./data/forAshleyCBCL.RDS")

## Examine the reliability estimates
m1 <- bind_rows(lapply(all.modsMIRT, function(x) x$rel_metrics))
m1$Scale <- names(all.modsMIRT)

## Now grab the score values for all of these so I can make a histogram for each of these
m2 <- lapply(all_data, function(x) table(rowSums(x$rep_vals)))
names(m2) <- names(all.modsMIRT)

## Isolate metrics of interest

## Now plot these
plot.dat1 <- reshape2::melt(m1, id.vars=c("Scale"))
#plot.dat1 <- plot.dat1[which(plot.dat1$variable %in% c("hurdleRel", "grmRel", "lambda3Rel", "omegaTRel")),]
plot.dat1 <- plot.dat1[-which(plot.dat1$variable %in% c("scaleSkew", "omegaTRel", "grmRel_rmZeroOption")),]
plot.dat1 <- plot.dat1[-which(plot.dat1$Scale %in% c("OTHER")),]

plot.dat1$variable <- plyr::revalue(plot.dat1$variable, c("alpha" = "Alpha", "singleFactorOmegaT" = "Omega", "hurdleRel" = "Hurdle", "grmRel" = "GRM", "singleFactorOmegaT" = "Omega"))

plot.dat1$Scale <- factor(plot.dat1$Scale, levels = unique(plot.dat1$Scale[order(plot.dat1$Scale)]))
plot.dat1$variable <- factor(plot.dat1$variable, levels =c("Alpha", "Omega", "GRM", "Hurdle"))


p1 <- ggplot(plot.dat1, aes(x=variable, y = value)) +
  geom_bar(stat="identity") +
  facet_wrap(Scale ~ .) +
  #theme(axis.text.x = element_text(angle=35, vjust=.5)) +
  coord_cartesian(ylim=c(.5, 1)) +
  ylab("Rel") +
  xlab("Metric") +
  theme_bw()
  #theme(axis.text.x = element_text(angle=35, vjust=.5))
plot.dat2 <- reshape2::melt(m2, id.vars=c("Scale"))
plot.dat2$L1 <- factor(plot.dat2$L1, levels = unique(plot.dat2$L1[order(plot.dat2$L1)]))
plot.dat2 <- plot.dat2[-which(plot.dat2$L1 %in% c("OTHER")),]
p2 <- ggplot(plot.dat2, aes(x=Var1, y = value)) +
  geom_bar(stat="identity") +
  facet_wrap(L1 ~ .) +
  theme_bw() + 
  coord_cartesian(xlim=c(0,10)) +
  scale_x_continuous(breaks = 0:10) +
  xlab("Sum score") +
  ylab("n")
library(ggpubr)
ggarrange(p2, p1)

## Now save this plot
out.plot <- ggarrange(p2, p1, labels = "AUTO")
ggsave(filename = "./reports/figure8_CBCLData.png", plot = out.plot, dpi=300, width = 13, height=7, units="in")

## Now make a loop to pair all of these plots together within domain
fig.list <- list()
fig.iter <- 1
for(i in sort(unique(plot.dat1$Scale))){
  ## Isolate all data across pairs
  vals.rel <- plot.dat1[which(plot.dat1$Scale == i),]
  vals.ss <- plot.dat2[which(plot.dat2$L1==i),]
  ## Now make a plot of these two together
  p1 <- ggplot(vals.rel, aes(x=variable, y = value)) +
    geom_bar(stat="identity") +
    #theme(axis.text.x = element_text(angle=35, vjust=.5)) +
    coord_cartesian(ylim=c(.6, 1)) +
    ylab("") +
    xlab("") +
    theme_minimal() +
    theme(text = element_text(color="black", size=16), axis.text = element_text(color="black"),panel.grid.major.x = element_blank())
  p2 <- ggplot(vals.ss, aes(x=Var1, y = value)) +
    geom_bar(stat="identity") +
    theme_minimal() + 
    coord_cartesian(xlim=c(0,10), ylim=c(0,6500)) +
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
          fig.list[[6]], fig.list[[7]], fig.list[[8]], ncol = 3, nrow=3)
ggsave(filename = "./reports/figure8_CBCLData.png", plot = out.plot, dpi=300, width = 20, height=9, units="in")
ggsave(filename = "./texOut/figures/figure8_CBCLData.png", plot = out.plot, dpi=300, width = 20, height=9, units="in")

## Now make a latex table for the descriptive stats
## Now write this as a latex object
library(dplyr)
library(kableExtra)
out.table <- data.frame(Subdomain=unique_scales, 
                        totalQuestion = unlist(lapply(all_data, function(x) dim(x$rep_vals)[2])),
                        meanVals = unlist(lapply(all_data, function(x) psych::describe(rowSums(x$rep_vals))$mean)),
                        sdVals = unlist(lapply(all_data, function(x) psych::describe(rowSums(x$rep_vals))$sd)),
                        skewVals = unlist(lapply(all_data, function(x) psych::describe(rowSums(x$rep_vals))$skew)))
out.table$skewVals <- round(out.table$skewVals, 2)
out.table$meanVals <- round(out.table$meanVals, 2)
out.table$sdVals <- round(out.table$sdVals, 2)

table_out <- kable(out.table, "latex")
cat(table_out, file="./reports/cbclOut.tex")

######## Now do updated table here
out.table <- data.frame(Subdomain=unique_scales, 
                        totalQuestion = unlist(lapply(all_data, function(x) dim(x$rep_vals)[2])),
                        meanVals = unlist(lapply(all_data, function(x) psych::describe(rowSums(x$rep_vals))$mean)),
                        sdVals = unlist(lapply(all_data, function(x) psych::describe(rowSums(x$rep_vals))$sd)),
                        skewVals = unlist(lapply(all_data, function(x) psych::describe(rowSums(x$rep_vals))$skew)))
descrip.table <- data.frame(Subdomain = c("AB","AD","AP","RB","SC","SP","TP","WD"),
                            Description = c("Aggressive Behavior",
                                            "Anxious/Depressed",
                                            "Attention Problems",
                                            "Rule-Breaking Behavior",
                                            "Somatic Complaints",
                                            "Social Problems",
                                            "Thought Problems",
                                            "Withdrawn/Depressed"))
out.table <- merge(descrip.table, out.table)
out.table2 <- merge(out.table, m1, by.x=c("Subdomain", "skewVals"), by.y=c("Scale", "scaleSkew"))
## Now order the outcome
out.table2 <- out.table2[,c("Subdomain", "Description", "totalQuestion", "meanVals", "sdVals", "skewVals","hurdleRel", "grmRel", "singleFactorOmegaT", "alpha")]
out.table2$skewVals <- round(out.table2$skewVals, 2)
out.table2$meanVals <- round(out.table2$meanVals, 2)
out.table2$sdVals <- round(out.table2$sdVals, 2)
out.table2$hurdleRel <- round(out.table2$hurdleRel, 2)
out.table2$grmRel <- round(out.table2$grmRel, 2)
out.table2$singleFactorOmegaT <- round(out.table2$singleFactorOmegaT, 2)
out.table2$singleFactorOmegaT <- round(out.table2$alpha, 2)
table_out <- kable(out.table2, "latex")
cat(table_out, file="./reports/cbclOut.tex")

## Now make a table of the parameter estimates for the most skewed dataset
## I am only going to do this for the category of interest: the WD scale
out.table3 <- all.modsMIRT$WD$params[,-c(1)]
colnames(out.table3) <- c("Item", "Factor Correlation", "2PL Discrimination", "2PL Difficulty", "GRM Discrimination","GRM DIfficulty")
rownames(out.table3) <- NULL
table_out <- kable(out.table3, "latex")
cat(table_out, file="./reports/cbclOutWDParam.tex")


## Quick plot of factor core
to.plot <- reshape2::melt(unlist(lapply(all.mods, function(x) unique(x$params$rho))))
to.plot$subscale <- rownames(to.plot)
p1 <- ggplot(to.plot, aes(x=subscale, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() + coord_cartesian(ylim=c(.45, .95))

to.plot <- reshape2::melt(unlist(lapply(all.modsMIRT, function(x) unique(x$params$rho))))
to.plot$subscale <- rownames(to.plot)
p2 <- ggplot(to.plot, aes(x=subscale, y=abs(value))) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() + coord_cartesian(ylim=c(.45, .95))


to.plot <- reshape2::melt(unlist(lapply(all.mods, function(x) unique(x$params$est_z_diff))))
to.plot$subscale <- rownames(to.plot)
to.plot$subscale <- gsub('[[:digit:]]+', '', to.plot$subscale)
ggplot(to.plot, aes(x=subscale, y=value)) +
  geom_violin() +
  theme_bw() + coord_cartesian(ylim=c(0,9))
to.plot <- reshape2::melt(unlist(lapply(all.modsMIRT, function(x) unique(x$params$est_z_diff))))
to.plot$subscale <- rownames(to.plot)
to.plot$subscale <- gsub('[[:digit:]]+', '', to.plot$subscale)
ggplot(to.plot, aes(x=subscale, y=value)) +
  geom_violin() +
  theme_bw() + coord_cartesian(ylim=c(0,9))


## Now look into the total E and total I variables
## This will take a lot of coding, most likely performed outside of R
## But I will prep the data here
## First grab the complete index for all I & E variables
all.e <- c("AB", "RB")
all_ext <- which(in_dict$subscale %in% all.e)
index <- all_ext
index <- index + 2
## Now grab the resp data
rep_vals <- in_dat[,index]
## Drop any columns with less than 10 response endorsements
rm.index <- which(apply(rep_vals, 2, function(x) length(which(x==2)))<15)
if(length(rm.index > 0)){
  rep_vals <- rep_vals[,-c(rm.index)]
}
rep_vals <- rep_vals[complete.cases(in_dat),]
mat.one <- rep_vals
mat.one[mat.one>1] <- 1
mat.two <- rep_vals
mat.two[mat.two==0] <- NA
#mat.two <- mat.two-l
## Now prep the mplus code
all.dat <- cbind(mat.one, mat.two)
colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
## Now print these data
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
out.file.name <- paste("./data/cbclModels/allExt.inp", sep='')
out.data.name <- paste("./data/cbclModels/allExt.txt", sep='')
mplusModeler(test, dataout=out.data.name, modelout = out.file.name, run = TRUE)
in.file.name <- paste("./data/cbclModels/allExt.out", sep='')
mod <- readModels(target = in.file.name)
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
b_org <- tidyr::pivot_wider(b[,c("est", "itemVal", "threshVal")], id_cols = "itemVal", names_from = "threshVal", values_from = "est")[,-1]
b_org <- data.frame(b_org)
## Now convert these difficulty metrics into the standard normal scale
b_orgStand <- apply(b_org, 2, function(x) x / a$est)
rho <- mod$parameters$unstandardized$est[which(mod$parameters$unstandardized$paramHeader=="SEVERITY.WITH")]
## Now prep these all for output
all.params <- cbind(i, 1:length(a$est),rho,a_z$est, b_z$est, a$est, b_orgStand)
prepVec <- c("subscale","item","rho","est_z_discrim","est_z_diff","est_grm_discrim")
grmDiffVec <- (dim(all.params)[2] - length(prepVec))
grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
prepVec <- c(prepVec, grmEstVec)
all.params <- data.frame(all.params)
colnames(all.params) <- prepVec
all.params
## Now do the factor scores here
## Now do the EAP factor score estimation here
all.dat$eap2PL_Hurdle <- NA
all.dat$se2PL_Hurdle <- NA
all.dat$eapGRM_Hurdle <- NA
all.dat$seGRM_Hurdle <- NA
theta1 <- seq(-8,8,0.2) # Quadrature points for theta1
theta2 <- seq(-8,8,0.2) # Quadrature points for theta1
theta <- expand.grid(theta1, theta2)
names(theta) <- c("theta1", "theta2")
prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
source("./scripts/rCode/hurdleFunctions.r")
itemtrace <- trace.line.pts(a$est, b_orgStand, a_z$est, b_z$estStand, theta)
for(respondent in 1:nrow(all.dat)) {
  pattern <- all.dat[respondent,1:dim(mat.one)[2]]
  qpoints <- theta
  lhood <- score(pattern, itemtrace, qpoints)
  
  eap2PL_Hurdle <- sum(lhood*prior*qpoints$theta1)/sum(lhood*prior)
  se2PL_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta1-eap2PL_Hurdle)^2)/sum(lhood*prior))
  all.dat$eap2PL_Hurdle[respondent] <- eap2PL_Hurdle
  all.dat$se2PL_Hurdle[respondent] <- se2PL_Hurdle
  
  eapGRM_Hurdle <- sum(lhood*prior*qpoints$theta2)/sum(lhood*prior)
  seGRM_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta2-eapGRM_Hurdle)^2)/sum(lhood*prior))
  all.dat$eapGRM_Hurdle[respondent] <- eapGRM_Hurdle
  all.dat$seGRM_Hurdle[respondent] <- seGRM_Hurdle
}
out.listE <- list(params = all.params, factor_scores = all.dat)

## Now do all Int
all.i <- c("AD", "WD", "SC")
all_int <- which(in_dict$subscale %in% all.i)
index <- all_int
index <- index + 2
## Now grab the resp data
rep_vals <- in_dat[,index]
## Drop any columns with less than 10 response endorsements
rm.index <- which(apply(rep_vals, 2, function(x) length(which(x==2)))<15)
if(length(rm.index > 0)){
  rep_vals <- rep_vals[,-c(rm.index)]
}
rep_vals <- rep_vals[complete.cases(in_dat),]
mat.one <- rep_vals
mat.one[mat.one>1] <- 1
mat.two <- rep_vals
mat.two[mat.two==0] <- NA
#mat.two <- mat.two-l
## Now prep the mplus code
all.dat <- cbind(mat.one, mat.two)
colnames(all.dat) <- paste("X", 1:ncol(all.dat), sep='')
## Now print these data
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
out.file.name <- paste("./data/cbclModels/allInt.inp", sep='')
out.data.name <- paste("./data/cbclModels/allInt.txt", sep='')
mplusModeler(test, dataout=out.data.name, modelout = out.file.name, run = TRUE)
in.file.name <- paste("./data/cbclModels/allInt.out", sep='')
mod <- readModels(target = in.file.name)
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
b_org <- tidyr::pivot_wider(b[,c("est", "itemVal", "threshVal")], id_cols = "itemVal", names_from = "threshVal", values_from = "est")[,-1]
b_org <- data.frame(b_org)
## Now convert these difficulty metrics into the standard normal scale
b_orgStand <- apply(b_org, 2, function(x) x / a$est)
rho <- mod$parameters$unstandardized$est[which(mod$parameters$unstandardized$paramHeader=="SEVERITY.WITH")]
## Now prep these all for output
all.params <- cbind(i, 1:length(a$est),rho,a_z$est, b_z$est, a$est, b_orgStand)
prepVec <- c("subscale","item","rho","est_z_discrim","est_z_diff","est_grm_discrim")
grmDiffVec <- (dim(all.params)[2] - length(prepVec))
grmEstVec <- paste("est_grm_diff_", 1:grmDiffVec, sep='')
prepVec <- c(prepVec, grmEstVec)
all.params <- data.frame(all.params)
colnames(all.params) <- prepVec
all.params
## Now do the factor scores here
## Now do the EAP factor score estimation here
all.dat$eap2PL_Hurdle <- NA
all.dat$se2PL_Hurdle <- NA
all.dat$eapGRM_Hurdle <- NA
all.dat$seGRM_Hurdle <- NA
theta1 <- seq(-8,8,0.2) # Quadrature points for theta1
theta2 <- seq(-8,8,0.2) # Quadrature points for theta1
theta <- expand.grid(theta1, theta2)
names(theta) <- c("theta1", "theta2")
prior <- mvtnorm::dmvnorm(theta, c(0, 0), matrix(c(1, rho, rho, 1), 2, 2, byrow = T))
source("./scripts/rCode/hurdleFunctions.r")
itemtrace <- trace.line.pts(a$est, b_orgStand, a_z$est, b_z$estStand, theta)
for(respondent in 1:nrow(all.dat)) {
  pattern <- all.dat[respondent,1:dim(mat.one)[2]]
  qpoints <- theta
  lhood <- score(pattern, itemtrace, qpoints)
  
  eap2PL_Hurdle <- sum(lhood*prior*qpoints$theta1)/sum(lhood*prior)
  se2PL_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta1-eap2PL_Hurdle)^2)/sum(lhood*prior))
  all.dat$eap2PL_Hurdle[respondent] <- eap2PL_Hurdle
  all.dat$se2PL_Hurdle[respondent] <- se2PL_Hurdle
  
  eapGRM_Hurdle <- sum(lhood*prior*qpoints$theta2)/sum(lhood*prior)
  seGRM_Hurdle <- sqrt(sum(lhood*prior*(qpoints$theta2-eapGRM_Hurdle)^2)/sum(lhood*prior))
  all.dat$eapGRM_Hurdle[respondent] <- eapGRM_Hurdle
  all.dat$seGRM_Hurdle[respondent] <- seGRM_Hurdle
}
out.listI <- list(params = all.params, factor_scores = all.dat)

## Now look at the pairs across E & I hurdle factors
all.vals <- bind_cols(out.listE$factor_scores, out.listI$factor_scores)
plot.vals <- all.vals[,c(53,55,119,121)]
colnames(plot.vals) <- c("pl2Ext", "grmExt", "pl2Int", "grmInt")
GGally::ggpairs(plot.vals)
## Look into E ~ I variables
lm(grmExt ~ grmInt, data = plot.vals)

## Now run a change point really quickly
model <- list(
  grmExt ~ 1 + grmInt,
  ~ grmInt
)
library(mcp)
mod <- mcp(model, data = plot.vals, cores = 3)
model <- list(
  pl2Ext ~ 1 + pl2Int,
  ~ pl2Int
)
mod2 <- mcp(model, data = plot.vals, cores = 3)
