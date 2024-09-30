## tabula rasa
rm(list=ls())

## First go ahead and load library(s) I might need
library(tidyverse)
library(doParallel)
library(mirt)
library(foreach)

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
cl <- makeCluster(4)
registerDoParallel(cl)
all.mods <- foreach(i = 1:length(all_item_position), .packages = c("mirt", "catIrt", "dplyr")) %dopar%{
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
  out_fn3 <- paste("./data/CBCL_scale_",i,"_Params.csv", sep='')
  saveRDS(val, file = out_fn3)
  val
}