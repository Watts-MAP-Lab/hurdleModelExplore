## Load library
source("~/GitHub/adroseHelperScripts/R/afgrHelpFunc.R")
library(ggplot2)
library(visreg)
# library("knitr")
# library("kableExtra")
library(psych)
library(purrr)
library(dplyr)


## CP data over
#system("rsync -hvrPt --ignore-existing rosena@login.accre.vu:/home/rosena/hurdleModelExplore/data/hurdleCollapse/*RDS ./data/hurdleCollapse/")
# scp arosen@sooner.oscer.ou.edu:/scratch/arosen/hurdleModelExplore/data/hurdleCollapse/* ./
# rsync -hvrPt --ignore-existing arosen@sooner.oscer.ou.edu:/scratch/arosen/hurdleModelExplore/data/hurdleCollapse/* ./data/hurdleCollapse/
# ## This will be used to read all hurdle collapse values
# in.dat1 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[1-9].RDS", full.names = TRUE)
# in.dat2 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[0-9][0-9].RDS", full.names = TRUE)
# in.dat3 <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = "*_[1][0-4][0-4].RDS", full.names = TRUE)
#
# ## Combine these
# in.dat <- c(in.dat1, in.dat2, in.dat3)
in.dat <- list.files(path="~/Documents/hurdleModelExplore/data/hurdleCollapse", pattern = ".RDS", full.names = TRUE)
all.collapse <- map_dfr(in.dat, function(file) {
  tmp <- readRDS(file)
  seedVal <- strsplit(strsplit(basename(file), "_")[[1]][2], ".RDS")[[1]][1]

  ## Estimate skew value here
  all_dat = tmp$est_mod@Data$data
  ## isolate 2pl and grm data
  n_items = ncol(all_dat)
  vec_2pl = seq(1, n_items/2)
  vec_grm = seq(max(vec_2pl)+1, n_items)
  grm_dat = all_dat[,vec_grm]
  pl2_dat = all_dat[,vec_2pl]
  obs_dat <- grm_dat
  obs_dat[is.na(obs_dat)] <- 0
  obs_dat <- obs_dat + pl2_dat
  skew_val_tmp = describe(rowSums(obs_dat))[,"skew"]
  tmp$allParams %>%
    mutate(seedVal = seedVal,skewVal = skew_val_tmp)
})
saveRDS(all.collapse, file = "./data/allCollapse.RDS")
#in.dat2 <- list.files(path = "~/Documents/extraDrive/Documents_Extra/hurdleModelExplore/data/hurdleCollapse/", pattern = ".RDS", full.names = TRUE)
# all.collapse <- map_dfr(in.dat2, function(file) {
#   tmp <- readRDS(file)
#   seedVal <- strsplit(strsplit(basename(file), "_")[[1]][2], ".RDS")[[1]][1]
#   
#   ## Estimate skew value here
#   all_dat = tmp$est_mod@Data$data
#   ## isolate 2pl and grm data
#   n_items = ncol(all_dat)
#   vec_2pl = seq(1, n_items/2)
#   vec_grm = seq(max(vec_2pl)+1, n_items)
#   grm_dat = all_dat[,vec_grm]
#   pl2_dat = all_dat[,vec_2pl]
#   obs_dat <- grm_dat
#   obs_dat[is.na(obs_dat)] <- 0
#   obs_dat <- obs_dat + pl2_dat
#   skew_val_tmp = describe(rowSums(obs_dat))[,"skew"]
#   tmp$allParams %>%
#     mutate(seedVal = seedVal,skewVal = skew_val_tmp)
# })
# saveRDS(all.collapse, file = "./data/allCollapse2.RDS")
#saveRDS(all.expand, file = "./data/allExpand.RDS")
all.collapse <- readRDS("./data/allCollapse.RDS")
all.collapse2 <- readRDS("./data/allCollapse2.RDS")

#all.collapseOld <- readRDS("~/Documents/extraDrive/Documents_Extra/hurdleModelExplore/data/")
#all.expand <- readRDS("./data/allExpand.RDS")

## Grab the rho error
#all.collapse$error <- all.collapse$rhoTrue - all.collapse$rho
## Now run an anova
## Attach real vals
source("./scripts/rCode/simParam.r")
all.dat.collapse <- merge(all.collapse, all.sim.vals, by.x = "seedVal", by.y="seed")
all.dat.collapse2 <- merge(all.collapse2, all.sim.vals, by.x = "seedVal", by.y="seed")
all.dat.collapse_new = all.dat.collapse
all.dat.collapse_old = all.dat.collapse2
all.dat.dat = merge(all.dat.collapse, all.dat.collapse2, 
                    by = names(all.dat.collapse)[c(1,2,32:43)],
                    suffixes = c("_new", "_old"))
#all.dat.collapse = all.dat.collapse_old

## Turn these into factors
for(i in c("nItem", "nCat", "grmDiscrim", "facCor","difGrmF","difGrmC" ,"iter", "discrim2pl", "dif2PL")){
  all.dat.collapse[,i] <- factor(all.dat.collapse[,i])
}
all.dat.collapse$sampSize <- all.dat.collapse$n

# Now deal with multiple values within each dataset
all.dat.collapse <- all.dat.collapse[!duplicated(all.dat.collapse$seedVal),]
## Now remove middle 2PL difficulty value

#all.dat.collapse <- all.dat.collapse[-which(all.dat.collapse$dif2PL=="-1"),]

## Now rm impossible values
#all.dat.collapse$alpha[which(all.dat.collapse$alpha<0)] <- NA
#all.dat.collapse$estHurdleRel_MIRT[which(all.dat.collapse$estHurdleRel_MIRT<.1)] <- NA
mean(sqrt((all.dat.collapse$trueGamRelSev-all.dat.collapse$omega.tot)^2))
mean(sqrt((all.dat.collapse$trueGamRelSev-all.dat.collapse$alpha)^2))
mean(sqrt((all.dat.collapse$trueGamRelSev-all.dat.collapse$grmRel)^2))
mean(sqrt((all.dat.collapse$trueGamRelSev-all.dat.collapse$estInfoRel)^2), na.rm=TRUE)
#mean(sqrt((all.dat.collapse$trueGamRelSev-all.dat.collapse$trueInfoRel)^2), na.rm=FALSE)

##### --DEFINE-VARIABLES-FOR-ANALYSIS------ 

# Define outcome variables (following your naming convention)
all.out <- c("alpha", "omega.tot", "grmRel", "estInfoRel", "trueInfoRel")

# Define predictor variables  
pred.vals <- c("nItem", "nCat", "facCor", "difGrmF", "dif2PL", "grmDiscrim", "discrim2pl")

# Calculate RMSE values
all.dat.collapse$rseAlpha <- sqrt((all.dat.collapse$alpha - all.dat.collapse$trueGamRelSev)^2)
all.dat.collapse$rseOmega <- sqrt((all.dat.collapse$omega.tot - all.dat.collapse$trueGamRelSev)^2)  
all.dat.collapse$rseGRMRel <- sqrt((all.dat.collapse$grmRel - all.dat.collapse$trueGamRelSev)^2)
all.dat.collapse$rseHurdRel <- sqrt((all.dat.collapse$estInfoRel - all.dat.collapse$trueGamRelSev)^2)
all.dat.collapse$rseHurdRelT <- sqrt((all.dat.collapse$trueInfoRel - all.dat.collapse$trueGamRelSev)^2)

## Check counts of each sim perm
count.check = summarySE(data = all.dat.collapse, measurevar = "rseAlpha",na.rm = TRUE,
                        groupvars = c("nItem", "nCat", "facCor", "difGrmF", "dif2PL", "grmDiscrim", "discrim2pl"))

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

## Now I need to try to prep the line graph Hao has described:
#More importantly, the presentation of the simulation result still revolves around AN(C)OVA. 
#I think the primary purpose of this paper is to propose a suitable formula of reliability for the huddle model 
#and show it can estimate the actual reliability properly, so a graph similar to Figure 5 should come first 
#and be the focus. Note the typo of y-label of the second row. 
#Instead of showing the main effects (each panel comparing two levels of a single factor), 
#the Figure should have 4 panels displaying the 4 crossed conditions between 2 factors, 
#and within each of the 4 panels is a line graph with different lines showing the four different measures' 
#RMSEs and x axis being the 32 crossed conditions between the remaining 5 factors (other than the 2 defining the 4 panels). 
#This way, all 128 simulation conditions are displayed, and the reader can see that the new measure has the smallest RMSE in 
#estimating the true reliability. Which two factors should define the four panels and how the combinations of the remaining 
#5 factors be arranged will be determined in a way to make the graph clearer. Maybe some try-and-errors are needed. 
#If a factor has very small effects (main and all relevant interactions), then it does not make much difference and its two levels can be averaged, 
#reducing the number of points on the graph. It is hard to tell from Figure 5 as it contains only the main effects. 
#There should also be a numerical tally showing how many of the 128 conditions see hurdle reliability with the smallest RMSE among the four. 
#Let me know if you need help here.
# First, create the long format data with RMSE values
library(tidyr)
library(dplyr)
# Your existing code up to creating long_rmse...
long_rmse <- all.dat.collapse %>%
  select(all_of(pred.vals), rseAlpha, rseOmega, rseGRMRel, rseHurdRel, rseHurdRelT, seedVal) %>%
  pivot_longer(cols = c(rseAlpha, rseOmega, rseGRMRel, rseHurdRel, rseHurdRelT),
               names_to = "method", 
               values_to = "rmse") %>%
  mutate(method = factor(method, 
                         levels = c("rseAlpha", "rseOmega", "rseGRMRel", "rseHurdRel", "rseHurdRelT"),
                         labels = c("Alpha", "Omega", "GRM Rel", "Hurdle Rel (Est)", "Hurdle Rel (True)")))

# Now this should work:
# Define factors
panel_factors <- c("discrim2pl", "dif2PL") 
x_axis_factors <- c("nCat", "facCor", "nItem", "difGrmF", "grmDiscrim")

# Your existing summarySE and plot_df creation...
plot_df <- summarySE(
  data = long_rmse,
  measurevar = "rmse",
  groupvars = c(panel_factors, x_axis_factors, "method")
)

plot_df$panel_label <- paste("2PL Dis:", plot_df$discrim2pl, "| 2PL Diff:", plot_df$dif2PL)

create_x_axis_with_order <- function(df, factor_order) {
  do.call(paste, c(df[factor_order], sep = " | "))
}

plot_df$x_axis_test = create_x_axis_with_order(plot_df, x_axis_factors)

# Filter plot_df for the methods we want
filtered_df <- plot_df[plot_df$method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel (Est)"), ]

# For each combination of panel_label and x_axis_test, find the method with lowest RMSE
# Modified to handle ties by giving preference to Hurdle Rel (Est)
win_results <- by(filtered_df, 
                  list(filtered_df$panel_label, filtered_df$x_axis_test), 
                  function(x) {
                    if(nrow(x) > 0) {
                      min_rmse <- min(x$rmse)
                      # Find all methods that have the minimum RMSE
                      min_indices <- which(x$rmse == min_rmse)
                      
                      # Check if "Hurdle Rel (Est)" is among the tied methods
                      hurdle_in_tie <- any(x$method[min_indices] == "Hurdle Rel (Est)")
                      
                      if(hurdle_in_tie) {
                        # If Hurdle Rel (Est) is tied for minimum, choose it
                        winning_method <- "Hurdle Rel (Est)"
                      } else {
                        # Otherwise, use the first minimum (original behavior)
                        winning_method <- as.character(x$method[min_indices[1]])
                      }
                      
                      data.frame(
                        panel_label = x$panel_label[1],
                        x_axis_test = x$x_axis_test[1], 
                        winning_method = winning_method
                      )
                    }
                  })

# Convert the by() result to a data frame
win_df <- do.call(rbind, win_results)
rownames(win_df) <- NULL

# Count wins per method within each panel
win_counts <- aggregate(x = rep(1, nrow(win_df)), 
                        by = list(panel_label = win_df$panel_label, 
                                  winning_method = win_df$winning_method), 
                        FUN = sum)
names(win_counts)[3] <- "wins"

# Create complete combinations of all panels and all methods
all_methods <- c("Alpha", "Omega", "GRM Rel", "Hurdle Rel (Est)")
all_panels <- unique(filtered_df$panel_label)
complete_combinations <- expand.grid(panel_label = all_panels, 
                                     winning_method = all_methods, 
                                     stringsAsFactors = FALSE)

# Merge with actual win counts, filling in 0s for missing combinations
win_counts_complete <- merge(complete_combinations, win_counts, 
                             by = c("panel_label", "winning_method"), 
                             all.x = TRUE)
win_counts_complete$wins[is.na(win_counts_complete$wins)] <- 0

# Create separate data frames for each method
alpha_data <- win_counts_complete[win_counts_complete$winning_method == "Alpha", ]
alpha_data$label <- paste0("α: ", alpha_data$wins)

omega_data <- win_counts_complete[win_counts_complete$winning_method == "Omega", ]
omega_data$label <- paste0("ω: ", omega_data$wins)

grm_data <- win_counts_complete[win_counts_complete$winning_method == "GRM Rel", ]
grm_data$label <- paste0("π: ", grm_data$wins)

hurdle_data <- win_counts_complete[win_counts_complete$winning_method == "Hurdle Rel (Est)", ]
hurdle_data$label <- paste0("πₕ: ", hurdle_data$wins)

# Plot with LARGER text sizes throughout
p1 <- plot_df %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel (Est)")) %>% 
  ggplot(aes(x = x_axis_test, y = rmse, group = method, color = method, linetype = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  # Alpha text (top) - INCREASED SIZE
  geom_text(data = alpha_data, 
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 1.1, 
            inherit.aes = FALSE, 
            size = 5, color = "#000000", family = "mono") +
  # Omega text (second) - INCREASED SIZE
  geom_text(data = omega_data, 
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 2.2, 
            inherit.aes = FALSE, 
            size = 5, color = "#404040", family = "mono") +
  # GRM text (third) - INCREASED SIZE
  geom_text(data = grm_data, 
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 3.3, 
            inherit.aes = FALSE, 
            size = 5, color = "#808080", family = "mono") +
  # Hurdle text (fourth) - INCREASED SIZE
  geom_text(data = hurdle_data, 
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.05, vjust = 4.4, 
            inherit.aes = FALSE, 
            size = 5, color = "#C0C0C0", family = "mono", face="bold") +
  facet_wrap(~ panel_label, ncol = 2, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          
                                "Omega" = "#404040",          
                                "GRM Rel" = "#808080",        
                                "Hurdle Rel (Est)" = "#C0C0C0")) +  
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel (Est)" = "solid"))+ 
  labs(title = NULL,
       x = NULL, 
       y = "RMSE",
       color = "Measure",
       linetype = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10, color = "black"),  # INCREASED
        axis.text.y = element_text(color = "black", size = 16),  # INCREASED
        axis.title = element_text(color = "black", size = 18),   # INCREASED
        plot.title = element_text(color = "black", size = 20),   # INCREASED
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),   # INCREASED
        legend.position = "bottom")

# Combine all annotation data into one label per panel
combined_annotation <- alpha_data %>%
  select(panel_label, alpha_label = label) %>%
  left_join(omega_data  %>% select(panel_label, omega_label  = label), by = "panel_label") %>%
  left_join(grm_data    %>% select(panel_label, grm_label    = label), by = "panel_label") %>%
  left_join(hurdle_data %>% select(panel_label, hurdle_label = label), by = "panel_label") %>%
  mutate(
    # Pad left column labels to same width so right column lines up
    alpha_label = formatC(alpha_label, width = max(nchar(alpha_label)), flag = "-"),
    omega_label = formatC(omega_label, width = max(nchar(omega_label)), flag = "-"),
    line1 = paste0(alpha_label, "  ", grm_label),
    line2 = paste0(omega_label, "  ", hurdle_label),
    label = paste(line1, line2, sep = "\n")
  )
# Left column: Alpha + Omega
left_annotation <- alpha_data %>%
  select(panel_label, alpha_label = label) %>%
  left_join(omega_data %>% select(panel_label, omega_label = label), by = "panel_label") %>%
  mutate(label = paste(alpha_label, omega_label, sep = "\n"))

# Right column: GRM + Hurdle
right_annotation <- grm_data %>%
  select(panel_label, grm_label = label) %>%
  left_join(hurdle_data %>% select(panel_label, hurdle_label = label), by = "panel_label") %>%
  mutate(label = paste(grm_label, hurdle_label, sep = "\n"))
plot_df$dis2PL = "[0.5:2.0]"
plot_df$dis2PL[which(plot_df$discrim2pl==2.3)] = "[2.0:3.5]"
plot_df$diff2PL = "[-2:0]"
plot_df$diff2PL[which(plot_df$dif2PL==0)] = "[0:2]"

p1 <- plot_df %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel (Est)")) %>% 
  ggplot(aes(x = x_axis_test, y = rmse, group = method, color = method, linetype = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  # geom_text(data = combined_annotation,
  #           aes(x = -Inf, y = Inf, label = label),
  #           hjust = -0.05, vjust = 1.1,
  #           inherit.aes = FALSE,
  #           size = 5, color = "#000000", family = "mono",
  #           lineheight = 1.6) +
  facet_grid(diff2PL~ dis2PL, scales = "free") + 
  scale_color_manual(values = c("Alpha"            = "#000000",
                                "Omega"            = "#404040",
                                "GRM Rel"          = "#808080",
                                "Hurdle Rel (Est)" = "#C0C0C0")) +
  scale_linetype_manual(values = c("Alpha"            = "solid",
                                   "Omega"            = "dashed",
                                   "GRM Rel"          = "dotted",
                                   "Hurdle Rel (Est)" = "solid")) +
  labs(title    = NULL,
       x        = NULL,
       y        = "RMSE",
       color    = "Measure",
       linetype = "Measure") +
  theme_minimal() +
  theme(
    axis.text.x   = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, color = "black"),
    axis.text.y   = element_text(color = "black", size = 16),
    axis.title    = element_text(color = "black", size = 18),
    plot.title    = element_text(color = "black", size = 20),
    legend.text   = element_text(color = "black", size = 18),
    legend.title  = element_text(color = "black", size = 20),
    strip.text    = element_text(color = "black", size = 16),
    legend.position = "bottom"
  )

print(p1)

ggsave(filename = "./reports/figure_mainEffect.png", plot = p1, dpi = 300, width = 12, height = 8, bg = "white")
ggsave(filename = "./texOut/figures/figure_mainEffect.png", plot = p1, dpi = 300, width = 12, height = 8, bg = "white")

p12 <- plot_df %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel (Est)")) %>% 
  filter(dif2PL %in% c(0)) %>% 
  ggplot(aes(x = x_axis_test, y = rmse, group = method, color = method, linetype = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          
                                "Omega" = "#404040",          
                                "GRM Rel" = "#808080",        
                                "Hurdle Rel (Est)" = "#C0C0C0")) +  
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel (Est)" = "solid"))+ 
  labs(title = NULL,
       x = NULL, 
       y = "RMSE",
       color = "Measure",
       linetype = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, color = "black"),  # INCREASED
        axis.text.y = element_text(color = "black", size = 16),  # INCREASED
        axis.title = element_text(color = "black", size = 18),   # INCREASED
        plot.title = element_text(color = "black", size = 20),   # INCREASED
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),   # INCREASED
        legend.position = "bottom")


## Write a csv with all of the effect sizes
#write.csv(x = all.ef.vals.iso, file = "~/Desktop/efVals.csv", quote=F, row.names = FALSE)

## Now I need to repeat this process for the raw effect sizes
long_rmse <- all.dat.collapse %>%
  select(all_of(pred.vals), alpha, omega.tot, grmRel, estInfoRel, trueGamRelSev, seedVal) %>%
  pivot_longer(cols = c(alpha, omega.tot, grmRel, estInfoRel, trueGamRelSev),
               names_to = "method", 
               values_to = "rmse") %>%
  mutate(method = factor(method, 
                         levels = c("alpha", "omega.tot", "grmRel", "estInfoRel", "trueGamRelSev"),
                         labels = c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")))

# Define factors
panel_factors <- c("discrim2pl", "dif2PL") 
x_axis_factors <- c("nCat", "facCor", "nItem", "difGrmF", "grmDiscrim")

# Create summarized data
plot_df <- summarySE(
  data = long_rmse,
  measurevar = "rmse",
  groupvars = c(panel_factors, x_axis_factors, "method")
)

# Create panel labels
plot_df$panel_label <- paste("2PL Dis:", plot_df$discrim2pl, "| 2PL Diff:", plot_df$dif2PL)

# Function to create x-axis labels for any given factor order
create_x_axis_with_order <- function(df, factor_order) {
  do.call(paste, c(df[factor_order], sep = " | "))
}

plot_df$x_axis_test = create_x_axis_with_order(plot_df, x_axis_factors)

# Create plot with different shapes and fill/color for Ground Truth
p <- plot_df %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_test, y = rmse, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  facet_wrap(~ panel_label, ncol = 2, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom")

print(p)

p2 <- plot_df %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  filter(dif2PL %in% c(0)) %>% 
  ggplot(aes(x = x_axis_test, y = rmse, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom")

ggsave(filename = "./reports/figure_mainEffect-Raw.png", plot = p, dpi = 300, width = 12, height = 8, bg = "white")
ggsave(filename = "./texOut/figures/figure_mainEffect-Raw.png", plot = p, dpi = 300, width = 12, height = 8, bg = "white")

## NOw create a fiugure with the raw and RMSE values
all_fig = ggpubr::ggarrange(p, p1, ncol=2, labels = c("A","B"), common.legend = TRUE, legend = "bottom")
ggsave(filename = "./texOut/figures/figure_mainEffect.png", plot = all_fig, dpi = 300, width = 26, height = 12, bg = "white")
all_fig2 = ggpubr::ggarrange(p2, p12, ncol=2, labels = c("A","B"), common.legend = TRUE, legend = "bottom")
ggsave(filename = "./scripts/beamer_pres/images/iso_all_effect.png", plot = all_fig2, dpi = 300, width = 14, height = 12, bg = "white")

## NOw incorporate the latest from Hao's feedback here:

# The graphs are very hard to read. Some suggestions: (a) The labels for discrimination and difficulty param-
# eters should be a range rather than a single value because they were sampled from a uniform distribution
# over a range. (b) The low 2PL difficulty columns use only a very small range of the y axis. In this case, it is
# better to align them in a row (rather than a column) so that they share the same y-axis accustomed to them.
# The other row for high 2PL difficulty conditions can use their own, separate range of y axis. (c) The graph
# can focus on reliability but leave out the RMSEs (which can be summarized in text). In addition to 2PL
# difficulty (for two rows of panels), the number of items and factor correlation should determine 4 different
# panels (within each row) so that across panels the true reliability has different profiles. (d) Within each
# panel, have 3 category conditions on the left and 6 category conditions on the right, as it affects the Hurdle
# reliability and interact with 2PL difficulty for true reliability. (e) Within the same category conditions (half
# of a panel), separate by 2PL discrimination, as it distinguishes Hurdle and GRM.
## I am going to have to use GGarrange to put together
## figures specific to: 2pl dif, 2pl dis, n itmes, & factor correlation
## Within each permutation of these I will then have to grab the mean values for the raw
## reliability within each level of n category

## First create each isolated dataset
long_rmse <- all.dat.collapse %>%
  select(all_of(pred.vals), alpha, omega.tot, grmRel, estInfoRel, trueGamRelSev, seedVal) %>%
  pivot_longer(cols = c(alpha, omega.tot, grmRel, estInfoRel, trueGamRelSev),
               names_to = "method", 
               values_to = "Reliability") %>%
  mutate(method = factor(method, 
                         levels = c("alpha", "omega.tot", "grmRel", "estInfoRel", "trueGamRelSev"),
                         labels = c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")))

## Create some values for all of the permuation names here
# nItem = 10 | 20
# cNat = 3 | 6
# facCor = 0.4 | 0.8
# difGrmF = -2 | 0
# dif2PL = -3 | 0
# grmDiscrim = 1 | 2.3
# discrim2pl = 1 | 2.3
### --- top-left ------------
iso_vals1 = long_rmse |> filter(
  discrim2pl ==1 & dif2PL == -3 & nItem == 10 & facCor == 0.4
)
iso_vals2 = long_rmse |> filter(
  discrim2pl ==1 & dif2PL == -3 & nItem == 20 & facCor == 0.4
)
iso_vals3 = long_rmse |> filter(
  discrim2pl ==1 & dif2PL == -3 & nItem == 10 & facCor == 0.8
)
iso_vals4 = long_rmse |> filter(
  discrim2pl ==1 & dif2PL == -3 & nItem == 20 & facCor == 0.8
)
## NOw take the summary SE vlaues for each of these
## and then combine them into a single plot
plot_one = summarySE(iso_vals1, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
          measurevar = "Reliability")
plot_two = summarySE(iso_vals2, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"),  
                     measurevar = "Reliability")
plot_thr = summarySE(iso_vals3, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
plot_fou = summarySE(iso_vals4, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
## Now combine these so they can be a panel within the larger figure
plot_1a = bind_rows(plot_one)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_1a$x_axis_text = create_x_axis_with_order(plot_1a, c("nCat", "grmDiscrim", "difGrmF"))

## Now make a plot for each of these
p2a <- plot_1a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==1 & dif2PL == -3 & nItem == 10 & facCor == 0.4")

## Now combine these so they can be a panel within the larger figure
plot_2a = bind_rows(plot_two)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_2a$x_axis_text = create_x_axis_with_order(plot_2a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2b <- plot_2a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==1 & dif2PL == -3 & nItem == 20 & facCor == 0.4")

## Now combine these so they can be a panel within the larger figure
plot_3a = bind_rows(plot_thr)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_3a$x_axis_text = create_x_axis_with_order(plot_2a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2c <- plot_3a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==1 & dif2PL == -3 & nItem == 10 & facCor == 0.8")


## Now combine these so they can be a panel within the larger figure
plot_4a = bind_rows(plot_fou)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_4a$x_axis_text = create_x_axis_with_order(plot_4a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2d <- plot_4a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==1 & dif2PL == -3 & nItem == 20 & facCor == 0.8")

p2a = p2a + coord_cartesian(ylim=c(.65, 1)) + theme(legend.position = "NULL")
p2b = p2b + coord_cartesian(ylim=c(.65, 1)) + theme(legend.position = "NULL")
p2c = p2c + coord_cartesian(ylim=c(.65, 1)) + theme(legend.position = "NULL")
p2d = p2d + coord_cartesian(ylim=c(.65, 1)) + theme(legend.position = "NULL")
top_left = ggarrange(p2a, p2b, p2c, p2d)

### --- top-right ------------
iso_vals1 = long_rmse |> filter(
  discrim2pl ==2.3 & dif2PL == -3 & nItem == 10 & facCor == 0.4
)
iso_vals2 = long_rmse |> filter(
  discrim2pl ==2.3 & dif2PL == -3 & nItem == 20 & facCor == 0.4
)
iso_vals3 = long_rmse |> filter(
  discrim2pl ==2.3 & dif2PL == -3 & nItem == 10 & facCor == 0.8
)
iso_vals4 = long_rmse |> filter(
  discrim2pl ==2.3 & dif2PL == -3 & nItem == 20 & facCor == 0.8
)
## NOw take the summary SE vlaues for each of these
## and then combine them into a single plot
plot_one = summarySE(iso_vals1, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
plot_two = summarySE(iso_vals2, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"),  
                     measurevar = "Reliability")
plot_thr = summarySE(iso_vals3, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
plot_fou = summarySE(iso_vals4, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
## Now combine these so they can be a panel within the larger figure
plot_1a = bind_rows(plot_one)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_1a$x_axis_text = create_x_axis_with_order(plot_1a, c("nCat", "grmDiscrim", "difGrmF"))

## Now make a plot for each of these
p2a <- plot_1a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==2.3 & dif2PL == -3 & nItem == 10 & facCor == 0.4")

## Now combine these so they can be a panel within the larger figure
plot_2a = bind_rows(plot_two)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_2a$x_axis_text = create_x_axis_with_order(plot_2a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2b <- plot_2a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==2.3 & dif2PL == -3 & nItem == 20 & facCor == 0.4")

## Now combine these so they can be a panel within the larger figure
plot_3a = bind_rows(plot_thr)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_3a$x_axis_text = create_x_axis_with_order(plot_2a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2c <- plot_3a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==2.3 & dif2PL == -3 & nItem == 10 & facCor == 0.8")


## Now combine these so they can be a panel within the larger figure
plot_4a = bind_rows(plot_fou)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_4a$x_axis_text = create_x_axis_with_order(plot_4a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2d <- plot_4a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==2.3 & dif2PL == -3 & nItem == 20 & facCor == 0.8")

p2a = p2a + coord_cartesian(ylim=c(.4, 1)) + theme(legend.position = "NULL")
p2b = p2b + coord_cartesian(ylim=c(.4, 1)) + theme(legend.position = "NULL")
p2c = p2c + coord_cartesian(ylim=c(.4, 1)) + theme(legend.position = "NULL")
p2d = p2d + coord_cartesian(ylim=c(.4, 1)) + theme(legend.position = "NULL")
top_right = ggarrange(p2a, p2b, p2c, p2d)


### --- bottom-right ------------
iso_vals1 = long_rmse |> filter(
  discrim2pl ==2.3 & dif2PL == 0 & nItem == 10 & facCor == 0.4
)
iso_vals2 = long_rmse |> filter(
  discrim2pl ==2.3 & dif2PL == 0 & nItem == 20 & facCor == 0.4
)
iso_vals3 = long_rmse |> filter(
  discrim2pl ==2.3 & dif2PL == 0 & nItem == 10 & facCor == 0.8
)
iso_vals4 = long_rmse |> filter(
  discrim2pl ==2.3 & dif2PL == 0 & nItem == 20 & facCor == 0.8
)
## NOw take the summary SE vlaues for each of these
## and then combine them into a single plot
plot_one = summarySE(iso_vals1, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
plot_two = summarySE(iso_vals2, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"),  
                     measurevar = "Reliability")
plot_thr = summarySE(iso_vals3, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
plot_fou = summarySE(iso_vals4, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
## Now combine these so they can be a panel within the larger figure
plot_1a = bind_rows(plot_one)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_1a$x_axis_text = create_x_axis_with_order(plot_1a, c("nCat", "grmDiscrim", "difGrmF"))

## Now make a plot for each of these
p2a <- plot_1a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==2.3 & dif2PL == 0 & nItem == 10 & facCor == 0.4")

## Now combine these so they can be a panel within the larger figure
plot_2a = bind_rows(plot_two)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_2a$x_axis_text = create_x_axis_with_order(plot_2a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2b <- plot_2a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==2.3 & dif2PL == 0 & nItem == 20 & facCor == 0.4")

## Now combine these so they can be a panel within the larger figure
plot_3a = bind_rows(plot_thr)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_3a$x_axis_text = create_x_axis_with_order(plot_2a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2c <- plot_3a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==2.3 & dif2PL == 0 & nItem == 10 & facCor == 0.8")


## Now combine these so they can be a panel within the larger figure
plot_4a = bind_rows(plot_fou)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_4a$x_axis_text = create_x_axis_with_order(plot_4a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2d <- plot_4a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==2.3 & dif2PL == 0 & nItem == 20 & facCor == 0.8")

p2a = p2a + coord_cartesian(ylim=c(.2, 1)) + theme(legend.position = "NULL")
p2b = p2b + coord_cartesian(ylim=c(.2, 1)) + theme(legend.position = "NULL")
p2c = p2c + coord_cartesian(ylim=c(.2, 1)) + theme(legend.position = "NULL")
p2d = p2d + coord_cartesian(ylim=c(.2, 1)) + theme(legend.position = "NULL")
bottom_right = ggarrange(p2a, p2b, p2c, p2d)


### --- bottom-left ------------
iso_vals1 = long_rmse |> filter(
  discrim2pl ==1 & dif2PL == 0 & nItem == 10 & facCor == 0.4
)
iso_vals2 = long_rmse |> filter(
  discrim2pl ==1 & dif2PL == 0 & nItem == 20 & facCor == 0.4
)
iso_vals3 = long_rmse |> filter(
  discrim2pl ==1 & dif2PL == 0 & nItem == 10 & facCor == 0.8
)
iso_vals4 = long_rmse |> filter(
  discrim2pl ==1 & dif2PL == 0 & nItem == 20 & facCor == 0.8
)
## NOw take the summary SE vlaues for each of these
## and then combine them into a single plot
plot_one = summarySE(iso_vals1, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
plot_two = summarySE(iso_vals2, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"),  
                     measurevar = "Reliability")
plot_thr = summarySE(iso_vals3, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
plot_fou = summarySE(iso_vals4, groupvars = c("nCat", "difGrmF", "grmDiscrim", 
                                              "discrim2pl", "dif2PL","nItem","facCor","method"), 
                     measurevar = "Reliability")
## Now combine these so they can be a panel within the larger figure
plot_1a = bind_rows(plot_one)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_1a$x_axis_text = create_x_axis_with_order(plot_1a, c("nCat", "grmDiscrim", "difGrmF"))

## Now make a plot for each of these
p2a <- plot_1a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==1 & dif2PL == 0 & nItem == 10 & facCor == 0.4")

## Now combine these so they can be a panel within the larger figure
plot_2a = bind_rows(plot_two)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_2a$x_axis_text = create_x_axis_with_order(plot_2a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2b <- plot_2a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==1 & dif2PL == 0 & nItem == 20 & facCor == 0.4")

## Now combine these so they can be a panel within the larger figure
plot_3a = bind_rows(plot_thr)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_3a$x_axis_text = create_x_axis_with_order(plot_2a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2c <- plot_3a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==1 & dif2PL == 0 & nItem == 10 & facCor == 0.8")


## Now combine these so they can be a panel within the larger figure
plot_4a = bind_rows(plot_fou)
## Now prepare the x axis text
x_axis_factors = c("nCat", "difGrmF", "grmDiscrim")
plot_4a$x_axis_text = create_x_axis_with_order(plot_4a, c("nCat", "grmDiscrim", "difGrmF"))
## Now make a plot for each of these
p2d <- plot_4a %>% 
  filter(method %in% c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth")) %>% 
  ggplot(aes(x = x_axis_text, y = Reliability, group = method, color = method, linetype = method)) +
  geom_line(aes(size = method)) +
  geom_point(aes(size = method, shape = method, fill = method), stroke = 1) +
  #facet_wrap(discrim2pl~ ., ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",    # Light grey
                                "Ground Truth" = "#000000")) + # Black border for Ground Truth
  scale_fill_manual(values = c("Alpha" = "#000000",           # Black
                               "Omega" = "#404040",           # Dark grey
                               "GRM Rel" = "#808080",         # Medium grey
                               "Hurdle Rel" = "#C0C0C0",     # Light grey
                               "Ground Truth" = "#FFFFFF")) + # White fill for Ground Truth
  scale_shape_manual(values = c("Alpha" = 16,                # Solid circle
                                "Omega" = 16,                # Solid circle
                                "GRM Rel" = 16,              # Solid circle
                                "Hurdle Rel" = 16,           # Solid circle
                                "Ground Truth" = 21)) +      # Fillable circle
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line and larger points
  guides(size = "none", fill = "none") +  # Hide size and fill from legend
  labs(title = NULL,
       x = NULL, 
       y = "Reliability",
       color = "Measure",
       linetype = "Measure",
       shape = "Measure") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(color = "black", size = 16),
        axis.title = element_text(color = "black", size = 18),
        plot.title = element_text(color = "black", size = 20),
        legend.text = element_text(color = "black", size = 18),      # INCREASE THIS
        legend.title = element_text(color = "black", size = 20),     # INCREASE THIS  
        strip.text = element_text(color = "black", size = 16),
        legend.position = "bottom") +
  ggtitle("discrim2pl ==1 & dif2PL == 0 & nItem == 20 & facCor == 0.8")

p2a = p2a + coord_cartesian(ylim=c(.35, 1)) + theme(legend.position = "NULL")
p2b = p2b + coord_cartesian(ylim=c(.35, 1)) + theme(legend.position = "NULL")
p2c = p2c + coord_cartesian(ylim=c(.35, 1)) + theme(legend.position = "NULL")
p2d = p2d + coord_cartesian(ylim=c(.35, 1)) + theme(legend.position = "NULL")
bottom_left = ggarrange(p2a, p2b, p2c, p2d)

### --- combine-all ------------------------
ggarrange(top_left, top_right, bottom_left, bottom_right)

### --- subfacet-approach --------------
library(ggh4x)
library(dplyr)
library(ggplot2)

plot_dat <- long_rmse |>
  filter(
    discrim2pl %in% c(1, 2.3),
    dif2PL %in% c(-3, 0),
    nItem %in% c(10, 20),
    facCor %in% c(0.4, 0.8),
    method %in% c(
      "Alpha", "Omega", "GRM Rel",
      "Hurdle Rel", "Ground Truth"
    )
  ) |>
  summarySE(
    measurevar = "Reliability",
    groupvars = c(
      "nCat", "difGrmF", "grmDiscrim",
      "discrim2pl", "dif2PL", "nItem",
      "facCor", "method"
    )
  )

plot_dat$grmDiscrim <- factor(
  plot_dat$grmDiscrim,
  levels = c(1, 2.3),
  labels = c("U(0.5, 2.0)", "U(2.0, 3.5)")
)

plot_dat$discrim2pl <- factor(
  plot_dat$discrim2pl,
  levels = c(1, 2.3),
  labels = c("U(0.5, 2.0)", "U(2.0, 3.5)")
)

plot_dat$dif2PL <- factor(
  plot_dat$dif2PL,
  levels = c(-3, 0),
  labels = c("U(-2, 0)", "U(0, 2)")
)

plot_dat$difGrmF <- factor(
  plot_dat$difGrmF,
  levels = c(-2, 0),
  labels = c("U(-2, 0)", "U(0, 2)")
)

# Reuse your existing function for the x-axis labels
plot_dat$x_axis_text <- create_x_axis_with_order(
  plot_dat,
  c("nCat", "grmDiscrim", "difGrmF")
)

p_nested <- ggplot(
  plot_dat,
  aes(
    x = x_axis_text,
    y = Reliability,
    group = method,
    color = method,
    linetype = method
  )
) +
  geom_line(aes(linewidth = method)) +
  geom_point(
    aes(size = method, shape = method, fill = method),
    stroke = 1
  ) +
  
  # Outer rows: dif2PL
  # Inner rows: facCor
  # Outer columns: discrim2pl
  # Inner columns: nItem
  ggh4x::facet_nested(
    dif2PL + facCor ~ discrim2pl + nItem,
    labeller = label_both
  ) +
  
  scale_color_manual(values = c(
    "Alpha" = "#000000",
    "Omega" = "#404040",
    "GRM Rel" = "#808080",
    "Hurdle Rel" = "#C0C0C0",
    "Ground Truth" = "#000000"
  )) +
  
  scale_fill_manual(values = c(
    "Alpha" = "#000000",
    "Omega" = "#404040",
    "GRM Rel" = "#808080",
    "Hurdle Rel" = "#C0C0C0",
    "Ground Truth" = "#FFFFFF"
  )) +
  
  scale_shape_manual(values = c(
    "Alpha" = 16,
    "Omega" = 16,
    "GRM Rel" = 16,
    "Hurdle Rel" = 16,
    "Ground Truth" = 21
  )) +
  
  scale_linetype_manual(values = c(
    "Alpha" = "solid",
    "Omega" = "dashed",
    "GRM Rel" = "dotted",
    "Hurdle Rel" = "solid",
    "Ground Truth" = "solid"
  )) +
  
  scale_linewidth_manual(values = c(
    "Alpha" = 0.5,
    "Omega" = 0.5,
    "GRM Rel" = 0.5,
    "Hurdle Rel" = 0.5,
    "Ground Truth" = 1
  )) +
  
  scale_size_manual(values = c(
    "Alpha" = 1.5,
    "Omega" = 1.5,
    "GRM Rel" = 1.5,
    "Hurdle Rel" = 1.5,
    "Ground Truth" = 2.5
  )) +
  
  guides(
    size = "none",
    linewidth = "none",
    fill = "none"
  ) +
  
  labs(
    x = NULL,
    y = "Reliability",
    color = "Measure",
    linetype = "Measure",
    shape = "Measure"
  ) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 9, color = "black"
    ),
    axis.text.y = element_text(
      size = 10, color = "black"
    ),
    axis.title = element_text(size = 12),
    strip.text = element_text(
      size = 10, color = "black"
    ),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    panel.spacing = grid::unit(0.8, "lines")
  )

p_nested


### --- fat-trim ------------------------
## NOw do a two-way interaction here:
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Calculate interaction effects for Factor Correlation × 2PL Difficulty
interaction_data <- all.dat.collapse %>%
  group_by(facCor, dif2PL) %>%
  dplyr::summarise(
    estInfoRel_mean = mean(estInfoRel, na.rm = TRUE),
    estInfoRel_se = sd(estInfoRel, na.rm = TRUE) / sqrt(n()),
    trueGamRelSev_mean = mean(trueGamRelSev, na.rm = TRUE),
    trueGamRelSev_se = sd(trueGamRelSev, na.rm = TRUE) / sqrt(n()),
    rseHurdRel_mean = mean(rseHurdRel, na.rm = TRUE),
    rseHurdRel_se = sd(rseHurdRel, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  ) %>%
  # Convert to long format
  pivot_longer(cols = ends_with(c("_mean", "_se")), 
               names_to = c("Measure", "Stat"), 
               names_pattern = "(.*)_(mean|se)$",
               values_to = "Value") %>%
  # Pivot wider to get separate mean and se columns
  pivot_wider(names_from = Stat, values_from = Value) %>%
  # Clean up measure names and create groupings
  mutate(
    Measure_Type = case_when(
      Measure %in% c("estInfoRel", "trueGamRelSev") ~ "Reliability",
      Measure == "rseHurdRel" ~ "RMSE",
      TRUE ~ "Other"
    ),
    Measure_Clean = case_when(
      Measure == "estInfoRel" ~ "Hurdle Rel",
      Measure == "trueGamRelSev" ~ "Ground Truth",
      Measure == "rseHurdRel" ~ "RMSE Hurdle Rel",
      TRUE ~ Measure
    ),
    # Calculate confidence intervals
    CI_lower = mean - se,
    CI_upper = mean + se,
    # Create factor labels
    facCor_label = paste("Factor Correlation =", facCor),
    dif2PL_label = paste("2PL Difficulty =", dif2PL)
  )

# Create a cleaner version with separate panels for each measure
create_individual_interaction_plot <- function(data, measure_name, y_label) {
  plot_data <- data[data$Measure_Clean == measure_name, ]
  
  ggplot(plot_data, aes(x = factor(facCor), y = mean, color = factor(dif2PL), 
                        group = factor(dif2PL))) +
    geom_line(size = 1.5) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                  width = 0.1, size = 1) +
    scale_color_manual(values = c("-2" = "#000000", "0" = "#666666"),
                       name = "2PL Difficulty") +
    labs(title = measure_name,
         x = "Factor Correlation",
         y = y_label) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    theme(
      text = element_text(size = 16),           # Overall text size
      axis.text = element_text(size = 14),      # Axis tick labels
      axis.title = element_text(size = 18),     # Axis titles (x and y labels)
      plot.title = element_text(size = 20),     # Main title
      strip.text = element_text(size = 16),     # Facet labels
      legend.text = element_text(size = 14),    # Legend text
      legend.title = element_text(size = 16),    # Legend title
      panel.spacing = unit(1.5, "cm")          # Increase spacing between facets
    ) +
    theme(text = element_text(size = 20, face = "bold"))
}

# Create individual plots
plot_hurdle <- create_individual_interaction_plot(interaction_data, "Hurdle Rel", "Reliability")
plot_ground <- create_individual_interaction_plot(interaction_data, "Ground Truth", "Reliability")  
plot_rmse <- create_individual_interaction_plot(interaction_data, "RMSE Hurdle Rel", "RMSE")

# Combine into composite figure
composite_interaction <- ggarrange(
  plot_hurdle, plot_ground, plot_rmse,
  ncol = 3, nrow = 1,
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

# Add overall title - FIXED: Use face = "italic" instead of style = "italic"
composite_interaction_titled <- annotate_figure(
  composite_interaction
  # top = text_grob("Factor Correlation × 2PL Difficulty Interaction Effects", 
  #                 face = "bold", size = 16),
)

print(composite_interaction_titled)

# Save the composite plot
ggsave("interaction_facCor_dif2PL_composite.png", composite_interaction_titled, 
       width = 15, height = 5, dpi = 300, bg = "white")
ggsave("./texOut/figures/interaction_facCor_dif2PL_composite.png", composite_interaction_titled, 
       width = 15, height = 5, dpi = 300, bg = "white")


## Now run an anova specific to the hurdle est
# Define outcomes and predictors
library(car)
library(effectsize)
library(xtable)
library(dplyr)
library(tidyr)

# Define outcomes and predictors
rmse.out <- c("estInfoRel", "trueGamRelSev", "rseHurdRel")

## Run all ANOVA models and extract effect sizes
all.anova.mods <- list()
all.ef.vals <- list()

for(i in 1:length(rmse.out)){
  # Create the model formula
  model.tmp <- as.formula(paste(rmse.out[i], "~(", paste(pred.vals, collapse = "+"), ")^7"))
  
  # Fit the ANOVA model
  tmp.anova <- lm(model.tmp, data = all.dat.collapse)
  all.anova.mods[[i]] <- tmp.anova
  
  # Get F statistics
  f.vals <- car::Anova(tmp.anova)
  
  # Calculate partial omega squared effect sizes
  eff.size.vals <- effectsize::omega_squared(f.vals, partial = TRUE)
  all.ef.vals[[i]] <- eff.size.vals
}

# Extract omega squared values and combine
all.ef.vals.iso <- lapply(all.ef.vals, function(x) x$Omega2)
all.ef.vals.combined <- dplyr::bind_cols(all.ef.vals[[1]]$Parameter, all.ef.vals.iso)
colnames(all.ef.vals.combined) <- c("Parameter", rmse.out)

# Function to clean up parameter names including interactions
clean_parameter_names <- function(param_name) {
  # Create mapping of old to new names
  name_mapping <- c(
    "nItem" = "Number of Items",
    "nCat" = "Number of Categories", 
    "facCor" = "Factor Correlation",
    "difGrmF" = "GRM Difficulty",
    "grmDiscrim" = "GRM Discrimination",
    "discrim2pl" = "2PL Discrimination",
    "dif2PL" = "2PL Difficulty"
  )
  
  # Split interaction terms by ":"
  components <- strsplit(param_name, ":")[[1]]
  
  # Replace each component with clean name
  clean_components <- sapply(components, function(x) {
    if(x %in% names(name_mapping)) {
      return(name_mapping[x])
    } else {
      return(x)
    }
  })
  
  # Rejoin with " × " for interactions
  return(paste(clean_components, collapse = " × "))
}

# Function to determine interaction level
get_interaction_level <- function(param_name) {
  return(length(strsplit(param_name, ":")[[1]]))
}

# Function to format values with bold for >= 0.10
format_omega_value <- function(omega_val, rmse_val = NULL) {
  if (is.na(omega_val) || omega_val <= 0.01) {
    if (!is.null(rmse_val) && !is.na(rmse_val) && rmse_val > 0.01) {
      if (rmse_val >= 0.10) {
        return(paste0("-- (\\textbf{", sprintf("%.2f", rmse_val), "})"))
      } else {
        return(paste0("-- (", sprintf("%.2f", rmse_val), ")"))
      }
    } else {
      return("--")
    }
  } else {
    main_val <- if (omega_val >= 0.10) {
      paste0("\\textbf{", sprintf("%.2f", omega_val), "}")
    } else {
      sprintf("%.2f", omega_val)
    }
    
    if (!is.null(rmse_val) && !is.na(rmse_val) && rmse_val > 0.01) {
      rmse_formatted <- if (rmse_val >= 0.10) {
        paste0("\\textbf{", sprintf("%.2f", rmse_val), "}")
      } else {
        sprintf("%.2f", rmse_val)
      }
      return(paste0(main_val, " (", rmse_formatted, ")"))
    } else if (!is.null(rmse_val)) {
      return(paste0(main_val, " (--)"))
    } else {
      return(main_val)
    }
  }
}

# Process the data
processed_data <- all.ef.vals.combined %>%
  mutate(
    Parameter_Clean = sapply(Parameter, clean_parameter_names),
    Interaction_Level = sapply(Parameter, get_interaction_level),
    # Format the columns - fix the mapply issue
    Hurdle_Col = mapply(format_omega_value, estInfoRel, rseHurdRel, USE.NAMES = FALSE),
    Ground_Truth_Col = sapply(trueGamRelSev, function(x) format_omega_value(x, NULL))
  ) %>%
  # Keep only rows where at least one value > 0.01
  filter(estInfoRel > 0.01 | trueGamRelSev > 0.01 | rseHurdRel > 0.01) %>%
  select(Parameter_Clean, Interaction_Level, Hurdle_Col, Ground_Truth_Col) %>%
  arrange(Interaction_Level, Parameter_Clean)

# Create the manual LaTeX table
create_latex_table <- function(data) {
  # Split by interaction levels
  main_effects <- data[data$Interaction_Level == 1, ]
  two_way <- data[data$Interaction_Level == 2, ]
  three_way <- data[data$Interaction_Level == 3, ]
  four_way <- data[data$Interaction_Level == 4, ]
  
  # Start building the table
  latex_lines <- c(
    "\\begin{table}[h]",
    "\t%\\centering",
    "\t\\begin{tabular}{l|cc}",
    "\t\t\\hline",
    "\t\t\\textbf{Effect} & \\textbf{$\\pi_H$} & \\textbf{Ground Truth} \\\\",
    "\t\t\\hline"
  )
  
  # Main Effects
  if (nrow(main_effects) > 0) {
    latex_lines <- c(latex_lines,
                     "\t\t\\multicolumn{3}{c}{\\textbf{Main Effects}} \\\\",
                     "\t\t\\hline")
    for (i in 1:nrow(main_effects)) {
      latex_lines <- c(latex_lines,
                       paste0("\t\t", main_effects$Parameter_Clean[i], " & ",
                              main_effects$Hurdle_Col[i], " & ",
                              main_effects$Ground_Truth_Col[i], " \\\\"))
    }
    latex_lines <- c(latex_lines, "\t\t\\hline")
  }
  
  # Two-Way Interactions
  if (nrow(two_way) > 0) {
    latex_lines <- c(latex_lines,
                     "\t\t\\multicolumn{3}{c}{\\textbf{Two-Way Interactions}} \\\\",
                     "\t\t\\hline")
    for (i in 1:nrow(two_way)) {
      latex_lines <- c(latex_lines,
                       paste0("\t\t", two_way$Parameter_Clean[i], " & ",
                              two_way$Hurdle_Col[i], " & ",
                              two_way$Ground_Truth_Col[i], " \\\\"))
    }
    latex_lines <- c(latex_lines, "\t\t\\hline")
  }
  
  # Three-Way Interactions
  if (nrow(three_way) > 0) {
    latex_lines <- c(latex_lines,
                     "\t\t\\multicolumn{3}{c}{\\textbf{Three-Way Interactions}} \\\\",
                     "\t\t\\hline")
    for (i in 1:nrow(three_way)) {
      latex_lines <- c(latex_lines,
                       paste0("\t\t", three_way$Parameter_Clean[i], " & ",
                              three_way$Hurdle_Col[i], " & ",
                              three_way$Ground_Truth_Col[i], " \\\\"))
    }
    if (nrow(four_way) > 0) {
      latex_lines <- c(latex_lines, "\t\t\\hline")
    }
  }
  
  # Four-Way and Higher Interactions
  if (nrow(four_way) > 0) {
    latex_lines <- c(latex_lines,
                     "\t\t\\multicolumn{3}{c}{\\textbf{Higher-Order Interactions}} \\\\",
                     "\t\t\\hline")
    for (i in 1:nrow(four_way)) {
      latex_lines <- c(latex_lines,
                       paste0("\t\t", four_way$Parameter_Clean[i], " & ",
                              four_way$Hurdle_Col[i], " & ",
                              four_way$Ground_Truth_Col[i], " \\\\"))
    }
  }
  
  # Close the table
  latex_lines <- c(latex_lines,
                   "\t\\end{tabular}",
                   paste0("\t\\caption{Partial $\\omega^2$ effect sizes for all effects exceeding 0.01. ",
                          "Main effects, two-way interactions, three-way interactions, and higher-order interactions are shown. ",
                          "For $\\pi_H$, RMSE values are shown in parentheses. ",
                          "Dashes (--) indicate effect sizes $\\leq$ 0.01. ",
                          "Bold values indicate effect sizes $\\geq$ 0.10.}"),
                   "\t\\label{tab:hurdle_effect_sizes}",
                   "\\end{table}")
  
  return(latex_lines)
}

# Generate and print the LaTeX table
latex_table <- create_latex_table(processed_data)
cat(paste(latex_table, collapse = "\n"))

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Define the outcome variables and predictors
reliability_vars <- c("estInfoRel", "trueGamRelSev")
rmse_vars <- c("rseHurdRel")

# Function to create a plot for a single predictor
create_single_predictor_plot <- function(data, predictor) {
  # Calculate means and standard errors for this specific predictor
  pred_means <- data %>%
    group_by(.data[[predictor]]) %>%
    dplyr::summarise(
      estInfoRel_mean = mean(estInfoRel, na.rm = TRUE),
      estInfoRel_se = sd(estInfoRel, na.rm = TRUE) / sqrt(n()),
      trueGamRelSev_mean = mean(trueGamRelSev, na.rm = TRUE),
      trueGamRelSev_se = sd(trueGamRelSev, na.rm = TRUE) / sqrt(n()),
      rseHurdRel_mean = mean(rseHurdRel, na.rm = TRUE),
      rseHurdRel_se = sd(rseHurdRel, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Fix column naming
  colnames(pred_means)[1] <- "Predictor_Value"
  
  # Convert to long format for means
  means_long <- pred_means %>%
    select(Predictor_Value, ends_with("_mean")) %>%
    pivot_longer(cols = ends_with("_mean"), 
                 names_to = "Measure", values_to = "Mean_Value") %>%
    mutate(Measure = gsub("_mean", "", Measure))
  
  # Convert to long format for standard errors
  se_long <- pred_means %>%
    select(Predictor_Value, ends_with("_se")) %>%
    pivot_longer(cols = ends_with("_se"), 
                 names_to = "Measure", values_to = "SE_Value") %>%
    mutate(Measure = gsub("_se", "", Measure))
  
  # Combine means and SEs
  plot_data <- means_long %>%
    left_join(se_long, by = c("Predictor_Value", "Measure")) %>%
    mutate(
      Measure_Type = case_when(
        Measure %in% c("estInfoRel", "trueGamRelSev") ~ "Reliability",
        Measure == "rseHurdRel" ~ "RMSE",
        TRUE ~ "Other"
      ),
      Measure_Clean = case_when(
        Measure == "estInfoRel" ~ "Hurdle Rel",
        Measure == "trueGamRelSev" ~ "Ground Truth",
        Measure == "rseHurdRel" ~ "RMSE Hurdle Rel",
        TRUE ~ Measure
      ),
      # Calculate confidence intervals
      CI_lower = Mean_Value - SE_Value,
      CI_upper = Mean_Value + SE_Value
    )
  
  # Clean predictor name for title
  predictor_clean <- case_when(
    predictor == "nItem" ~ "Number of items",
    predictor == "nCat" ~ "Number of categories", 
    predictor == "facCor" ~ "Factor correlation",
    predictor == "difGrmF" ~ "GRM difficulty",
    predictor == "grmDiscrim" ~ "GRM discrimination",
    predictor == "discrim2pl" ~ "2PL discrimination",
    predictor == "dif2PL" ~ "2PL difficulty",
    TRUE ~ predictor
  )
  
  # Create the plot
  p <- plot_data %>%
    ggplot(aes(x = factor(Predictor_Value), y = Mean_Value, fill = Measure_Clean)) +
    geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.3) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                  position = position_dodge(width = 0.9), 
                  width = 0.4, size = 0.5, color = "black") +
    facet_wrap(~ Measure_Type, scales = "free", ncol = 2) +
    scale_fill_manual(values = c("Hurdle Rel" = "#C0C0C0",
                                 "Ground Truth" = "#FFFFFF",
                                 "RMSE Hurdle Rel" = "#808080")) +
    labs(title = predictor_clean,
         x = NULL,
         y = NULL,
         fill = "Measure") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10, face = "bold"),
      legend.position = "bottom",
      legend.key.size = unit(0.4, "cm"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0.5, "lines")
    ) +
    theme(
      text = element_text(size = 16),           # Overall text size
      axis.text = element_text(size = 20),      # Axis tick labels
      axis.title = element_text(size = 18),     # Axis titles (x and y labels)
      plot.title = element_text(size = 20),     # Main title
      strip.text = element_text(size = 16),     # Facet labels
      legend.text = element_text(size = 14),    # Legend text
      legend.title = element_text(size = 16),    # Legend title
      panel.spacing = unit(1.5, "cm")          # Increase spacing between facets
    ) +
    theme(text = element_text(size = 20, face = "bold")) +
    guides(fill = guide_legend(nrow = 1))
  
  return(p)
}

# Create individual plots for each predictor
individual_plots <- list()

for(i in 1:length(pred.vals)) {
  individual_plots[[i]] <- create_single_predictor_plot(all.dat.collapse, pred.vals[i])
}

# Name the plots for reference
names(individual_plots) <- pred.vals

# Create composite figure using ggarrange
composite_plot <- ggpubr::ggarrange(
  plotlist = individual_plots,
  ncol = 3, nrow = 3,  # Adjust based on number of predictors
  common.legend = TRUE,
  legend = "bottom",
  align = "hv"
)

# Add overall title
composite_plot_titled <- ggpubr::annotate_figure(
  composite_plot
  # top = ggpubr::text_grob("Main Effects: Reliability and RMSE Across All Simulation Parameters", 
  #                         face = "bold", size = 16),
)

# Display the composite plot
print(composite_plot_titled)

# Save the composite figure
ggsave("main_effects_composite.png", composite_plot_titled, 
       width = 18, height = 15, dpi = 300, bg = "white")

ggsave("./texOut/figures/main_effects_composite.png", composite_plot_titled, 
       width = 20, height = 12, dpi = 300, bg = "white")

# # Also save individual plots if needed
# for(i in 1:length(pred.vals)) {
#   predictor_name <- pred.vals[i]
#   filename <- paste0("main_effect_", predictor_name, ".png")
#   ggsave(filename, individual_plots[[i]], 
#          width = 8, height = 4, dpi = 300, bg = "white")
# }


## Now run some quick interaction plots
# Your existing ANCOVA setup
all.ancova.mods <- list()
all.out = c("alpha", "omega.tot", "grmRel", "estInfoRel")
all.out = c("rseAlpha", "rseOmega", "rseGRMRel", "rseHurdRel")

# Storage for results
all.ancova.rsquared <- list()
all.ancova.ef.vals <- list()

for(i in 1:length(all.out)){
  ## Fit the model
  model.tmp <- as.formula(paste(all.out[i], "~(skewVal * dif2PL)"))
  tmp.anova <- lm(model.tmp, data = all.dat.collapse)
  all.ancova.mods[[i]] <- tmp.anova
  
  ## Extract R-squared
  rsq <- summary(tmp.anova)$r.squared
  adj.rsq <- summary(tmp.anova)$adj.r.squared
  all.ancova.rsquared[[i]] <- data.frame(
    metric = all.out[i],
    r.squared = rsq,
    adj.r.squared = adj.rsq
  )
  
  ## Extract effect sizes
  f.vals <- car::Anova(tmp.anova)
  eff.size.vals <- effectsize::omega_squared(f.vals, partial = TRUE)
  eff.size.vals$metric <- all.out[i]
  all.ancova.ef.vals[[i]] <- eff.size.vals
}

# Combine results
ancova.rsquared.df <- do.call(rbind, all.ancova.rsquared)
ancova.effects.df <- do.call(rbind, all.ancova.ef.vals)

# Create effect size summary table
ancova.summary <- ancova.effects.df %>%
  select(Parameter, Omega2_partial, metric) %>%
  pivot_wider(names_from = metric, values_from = Omega2_partial) %>%
  mutate(across(where(is.numeric), round, 3))

# Create R-squared summary table
rsq.summary <- ancova.rsquared.df %>%
  select(metric, r.squared, adj.r.squared) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  pivot_longer(cols = c(r.squared, adj.r.squared), 
               names_to = "Parameter", 
               values_to = "value") %>%
  pivot_wider(names_from = metric, values_from = value)

# Combine effect sizes and R-squared into one table
combined.table <- bind_rows(ancova.summary, rsq.summary)

# Reorder rows to put R-squared at the bottom
row.order <- c("skewVal", "dif2PL", "skewVal:dif2PL", "r.squared", "adj.r.squared")
combined.table <- combined.table[match(row.order, combined.table$Parameter, nomatch = 0),]

# Clean up parameter names for LaTeX
combined.table$Parameter <- case_when(
  combined.table$Parameter == "skewVal" ~ "Skewness",
  combined.table$Parameter == "dif2PL" ~ "2PL Difficulty", 
  combined.table$Parameter == "skewVal:dif2PL" ~ "Skewness $\\times$ 2PL Difficulty",
  combined.table$Parameter == "r.squared" ~ "$R^2$",
  combined.table$Parameter == "adj.r.squared" ~ "Adjusted $R^2$",
  TRUE ~ combined.table$Parameter
)

# Clean up column names to match your existing naming convention
colnames(combined.table) <- c("Parameter", "Alpha", "Omega", "GRM Rel", "Hurdle Rel")

print("Combined ANCOVA Results Table:")
print(combined.table)

# Create LaTeX table
library(kableExtra)
latex.table <- kable(combined.table, format = "latex", 
                     caption = "ANCOVA Results: Partial $\\omega^2$ Effect Sizes and Model $R^2$ Values",
                     booktabs = TRUE,
                     escape = FALSE) %>%  # Don't escape LaTeX symbols
  kable_styling() %>%
  pack_rows("Effect Sizes", 1, 3) %>%
  pack_rows("Model Fit", 4, 5)

# Write to file
cat(latex.table, file = "./reports/ancova_results_table.tex")
print("LaTeX table written to ./reports/ancova_results_table.tex")
## Now try to do this all in a single figure
all.out <- c("rseHurdRel", "rseGRMRel", "rseOmega", "rseAlpha") 
dat.tmp <- all.dat.collapse[,c("skewVal", "dif2PL", "iter", all.out)]
dat.tmp <- reshape2::melt(dat.tmp, id.vars=c("skewVal", "dif2PL", "iter"))
dat.tmp$variable <- plyr::revalue(dat.tmp$variable, c("rseHurdRel" = "Hurdle Rel", "rseGRMRel" = "GRM Rel", "rseOmega" = "Omega", "rseAlpha" = "Alpha"))
dat.tmp$variable <- factor(dat.tmp$variable, levels = c("Alpha", "Omega", "GRM Rel", "Hurdle Rel"))

tmp2 <- ggplot(data = dat.tmp, aes(x=skewVal, y=value, 
                                   group = variable, 
                                   color = variable,
                                   linetype = variable,
                                   size = variable)) +
  geom_smooth(method="lm", se = FALSE) +
  geom_point(alpha = 0.05, aes(size = variable)) +
  theme_minimal() +
  scale_color_manual(values = c("Alpha" = "#000000",          
                                "Omega" = "#404040",          
                                "GRM Rel" = "#808080",        
                                "Hurdle Rel" = "#C0C0C0")) +  
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1)) +
  facet_grid(. ~ dif2PL, scales = "free_x") + 
  xlab("Skew") +
  ylab("RMSE") +
  coord_cartesian(ylim = c(0, .6)) +
  guides(color = guide_legend(title="Reliability Metric"),
         linetype = guide_legend(title="Reliability Metric"),
         size = guide_legend(title="Reliability Metric")) +
  # --- THEME MODIFICATIONS START HERE ---
  theme(
    text = element_text(color="black", size=16), 
    axis.text = element_text(color="black"), 
    legend.key.width = unit(1, "in"),
    
    # 1. This controls the gap between panels
    panel.spacing = unit(2, "lines"), 
    
    # 2. Optional: Add spacing outside the panels
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10) 
  )
all.out <- c("alpha", "omega.tot", "grmRel", "estInfoRel", "trueGamRelSev") 
dat.tmp <- all.dat.collapse[,c("skewVal", "dif2PL", "iter", all.out)]
dat.tmp <- reshape2::melt(dat.tmp, id.vars=c("skewVal", "dif2PL", "iter"))
dat.tmp$variable <- plyr::revalue(dat.tmp$variable, c("estInfoRel" = "Hurdle Rel", "grmRel" = "GRM Rel", "omega.tot" = "Omega", "alpha" = "Alpha", "trueRel" = "True", "trueGamRelSev" = "Ground Truth"))
dat.tmp$variable <- factor(dat.tmp$variable, levels = c("Alpha", "Omega", "GRM Rel", "Hurdle Rel", "Ground Truth"))

tmp1 <- ggplot(data = dat.tmp, aes(x=skewVal, y=value, 
                                   group = variable, 
                                   color = variable,
                                   linetype = variable,
                                   size = variable)) +
  geom_smooth(method="lm", se = FALSE) +
  geom_point(alpha = 0.01, aes(size = variable)) +
  theme_minimal() +
  scale_color_manual(values = c("Alpha" = "#000000",          # Black
                                "Omega" = "#404040",          # Dark grey
                                "GRM Rel" = "#808080",        # Medium grey
                                "Hurdle Rel" = "#C0C0C0",     # Light grey
                                "Ground Truth" = "#000000")) + # Black
  scale_linetype_manual(values = c("Alpha" = "solid", 
                                   "Omega" = "dashed", 
                                   "GRM Rel" = "dotted", 
                                   "Hurdle Rel" = "solid",
                                   "Ground Truth" = "dotted")) +
  scale_size_manual(values = c("Alpha" = 1, 
                               "Omega" = 1, 
                               "GRM Rel" = 1, 
                               "Hurdle Rel" = 1,
                               "Ground Truth" = 2)) +  # Thicker line for Ground Truth
  facet_grid(. ~ dif2PL, scales = "free_x") +  # Vertical facets by dif2PL with free x-axis
  xlab("Skew") +
  ylab("Reliability") +
  coord_cartesian(ylim = c(0.5,1)) +
  guides(color = guide_legend(title="Reliability Metric"),
         linetype = guide_legend(title="Reliability Metric"),
         size = guide_legend(title="Reliability Metric")) +
  # --- THEME MODIFICATIONS START HERE ---
  theme(
    text = element_text(color="black", size=16), 
    axis.text = element_text(color="black"), 
    legend.key.width = unit(1, "in"),
    
    # 1. This controls the gap between panels
    panel.spacing = unit(2, "lines"), 
    
    # 2. Optional: Add spacing outside the panels
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10) 
  )
tmp <- ggpubr::ggarrange(tmp1, tmp2, labels = "AUTO", common.legend = TRUE, legend = "bottom")
ggsave(filename = "./reports/figureANCOR_twoWay1.png", plot=tmp2, width = 14, height = 7, dpi= 400, units="in", bg="white")
ggsave(filename = "./texOut/figures/figureANCOR_twoWay1.png", plot=tmp, width = 14, height = 7, dpi= 300, units="in", bg="white")
