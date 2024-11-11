## tabula rasa
rm(list=ls())

## First go ahead and load library(s) I might need
library(tidyverse)
library(doParallel)
library(mirt)
library(foreach)
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")

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

## Write these data for Ashley
saveRDS(object = all.mods, file = "./data/forAshleyCBCL.RDS")

## Quick plot of factor core
to.plot <- reshape2::melt(unlist(lapply(all.mods, function(x) unique(x$params$rho))))
to.plot$subscale <- rownames(to.plot)
ggplot(to.plot, aes(x=subscale, y=value)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() + coord_cartesian(ylim=c(.45, .95))
