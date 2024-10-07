## tabula rasa
rm(list=ls())

## First go ahead and load library(s) I might need
library(tidyverse)
library(doParallel)
library(mirt)
library(foreach)

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
  ## Now cross tab these
  crosstab <- rep_vals %>%
    group_by(across(everything())) %>%
    summarise(Count = n(), .groups = "drop") %>%
    ungroup()
  out_data <- list(rep_vals = rep_vals, crosstab = crosstab)
  all_data[[i]] <- out_data
}


## Now run one of these through the hurdle model in julia
cl <- makeCluster(3)
registerDoParallel(cl)
all.mods <- foreach(i = 1:length(all_item_position)) %dopar%{
  in_resp <- all_data[[i]]$rep_vals
  in_resp <- in_resp[complete.cases(in_resp),]
  in_tabs <- all_data[[i]]$crosstab
  in_tabs <- in_tabs[complete.cases(in_tabs),]
  ## Now write these
  out_fn1 <- paste("./data/CBCL_scale_",i,"_Resp.csv", sep='')
  out_fn2 <- paste("./data/CBCL_scale_",i,"_Tabs.csv", sep='')
  write.csv(x = in_resp, file = out_fn1, quote=F, row.names=F)
  write.csv(x = in_tabs, file = out_fn2, quote=F, row.names=F)
  ## Now run model in julia
  sys.command <- paste("julia ./scripts/juliaCode/mHurdleFlex.jl ", out_fn1, out_fn2, sep=' ')
  val <- system(sys.command, intern = TRUE)
  out_fn3 <- paste("./data/CBCL_scale_",i,"_Params.RData", sep='')
  saveRDS(val, file = out_fn3)
  val
}

## Now read all of the data from ACCRE here
cbcl_1_res = system("cat ./data/accreResults/r_output_66294209_1.txt", intern=TRUE)
cbcl_1_dat = read.csv("./data/CBCL_scale_1_Resp.csv")
cbcl_1_params = return_Mod_Params(cbcl_1_res, cbcl_1_dat)
#cbcl_2_res = system("cat ./data/accreResults/r_output_662", intern=TRUE)
#cbcl_3_res = system("cat ./data/accreResults/r_output_66294209_1.txt", intern=TRUE)
cbcl_4_res = system("cat ./data/accreResults/r_output_66275684_4.txt", intern=TRUE)
cbcl_4_dat = read.csv("./data/CBCL_scale_4_Resp.csv")
cbcl_4_params = return_Mod_Params(cbcl_4_res, cbcl_4_dat)
#cbcl_5_res = system("cat ./data/accreResults/r_output_66294209_1.txt", intern=TRUE)
#cbcl_6_res = system("cat ./data/accreResults/r_output_66294209_1.txt", intern=TRUE)
cbcl_7_res = system("cat ./data/accreResults/r_output_66294215_7.txt", intern=TRUE)
cbcl_7_dat = read.csv("./data/CBCL_scale_7_Resp.csv")
cbcl_7_params = return_Mod_Params(cbcl_7_res, cbcl_7_dat)
#cbcl_8_res = system("cat ./data/accreResults/r_output_66294209_1.txt", intern=TRUE)
#cbcl_9_res = system("cat ./data/accreResults/r_output_66294209_1.txt", intern=TRUE)

