## Attach real vals
sim.param.nitems <- c(8)
sim.param.ncates <- c(4)
sim.param.discri <- c(2.4,2.4)
sim.param.2pl.spread <- c(3)
sim.param.sample <- c(10000)
sim.param.faccor <- c(.3,.8)
sim.param.difgrmF <- c(-1.25,-.5)
sim.param.difgrmC <- c(2)
sim.param.discri2 <- c(1.6, 2.4)
sim.iter <- 1:50
all.sim.vals <- expand.grid(sim.param.nitems, sim.param.ncates, sim.param.discri, 
                            sim.param.2pl.spread, sim.param.sample, sim.param.faccor, 
                            sim.param.difgrmF, sim.param.difgrmC, sim.param.discri2, sim.iter)
colnames(all.sim.vals) <- c("nItem", "nCat", "grmDiscrim", "diffSpread", "n", "facCor",
                            "difGrmF", "difGrmC", "discrim", "iter")
all.sim.vals$seed <- 1:nrow(all.sim.vals)

source("./scripts/rCode/hurdleFunctions-cleanComment.r")
suppressPackageStartupMessages(library("tidyverse"))
library("numDeriv")
library("mirt")
library(ggpubr)
library(dplyr)
library(ggplot2)

seedVal <- 1
set.seed(16)

a <- rep(all.sim.vals$discrim[seedVal], all.sim.vals$nItem[seedVal])
b <- genDiffGRM(
  num_items       = all.sim.vals$nItem[seedVal],
  num_categories  = all.sim.vals$nCat[seedVal],
  min             = all.sim.vals$difGrmF[seedVal],
  max             = all.sim.vals$difGrmC[seedVal],
  min_gap         = .3
)
a_z <- rep(all.sim.vals$grmDiscrim[seedVal], all.sim.vals$nItem[seedVal])

## Item / simulation setup
n_items     <- 10
n_cat       <- 3
rho         <- 0.2
N           <- all.sim.vals$n[seedVal]
add.val.grm <- 1.5

# --------------------------------------------------------------------------
# GRM (continuous severity) item parameters — shared across all conditions
# --------------------------------------------------------------------------
a_grm1 <- rep(2, n_items)
a_grm2 <- rep(2, n_items)
a_grm3 <- rep(2, n_items)

b_grm <- genDiffGRM(
  num_items      = n_items,
  num_categories = n_cat,
  min_gap        = .75,
  min            = -1.5,
  max            = 7
)

# --------------------------------------------------------------------------
# 2PL (hurdle / zero-inflation) discrimination parameters
# --------------------------------------------------------------------------
a_2pl <- rep(1.5, n_items)

# --------------------------------------------------------------------------
# 2PL difficulty parameters for the hurdle component
#
# To guarantee a monotonically DECREASING sum-score endorsement pattern in
# the zero-inflated conditions (reps2 & reps3), difficulties must be
# sufficiently high (positive) so that P(non-zero) < 0.5 for most
# respondents.  The more positive (and the larger the spread), the stronger
# the zero-inflation and the steeper the decreasing tail.
#
#   Condition 1 – No zero inflation   : difficulties very negative → nearly
#                                        everyone clears the hurdle.
#   Condition 2 – Moderate ZI        : intermediate positive difficulties.
#   Condition 3 – Severe ZI          : high positive difficulties.
# --------------------------------------------------------------------------
b_21pl <- rep(-3.0, n_items)          # No ZI    – hurdle almost never active
b_21pl <- runif(n_items, min = -2, max = -0.5)          # No ZI    – hurdle almost never active
b_22pl <- seq(-0.5, 1.5, length.out = n_items)   # Moderate ZI
b_23pl <- seq(1, 2.5, length.out = n_items)   # Severe ZI

# --------------------------------------------------------------------------
# Convert to mirt parameterization
# --------------------------------------------------------------------------
mirt_pars <- to_mirt(a_grm1, b_grm, a_2pl, b_21pl, n_cat)
a_m1  <- mirt_pars$a_m;  d_m1  <- mirt_pars$d_m
a_z_m1 <- mirt_pars$a_z_m; d_z_m1 <- mirt_pars$d_z_m

mirt_pars <- to_mirt(a_grm2, b_grm, a_2pl, b_22pl, n_cat)
a_m2  <- mirt_pars$a_m;  d_m2  <- mirt_pars$d_m
a_z_m2 <- mirt_pars$a_z_m; d_z_m2 <- mirt_pars$d_z_m

mirt_pars <- to_mirt(a_grm3, b_grm, a_2pl, b_23pl, n_cat)
a_m3  <- mirt_pars$a_m;  d_m3  <- mirt_pars$d_m
a_z_m3 <- mirt_pars$a_z_m; d_z_m3 <- mirt_pars$d_z_m

# --------------------------------------------------------------------------
# Draw true theta from N(0, Sigma_pop)
# --------------------------------------------------------------------------
varCovMat <- matrix(c(1, rho, rho, 1), ncol = 2)
N <- 100000
theta <- MASS::mvrnorm(n = N, mu = c(0, 0), Sigma = varCovMat)
theta <- theta*1.5
# --------------------------------------------------------------------------
# Simulate hurdle responses
# --------------------------------------------------------------------------
library(dplyr)
reps1 <- suppressWarnings(sim_mirt_hurdle(a_m1, d_m1, a_z_m1, d_z_m1, rho = rho, THETA = theta, N = N, estModel = FALSE))
reps2 <- suppressWarnings(sim_mirt_hurdle(a_m2, d_m2, a_z_m2, d_z_m2, rho = rho, THETA = theta, N = N, estModel = FALSE))
reps3 <- suppressWarnings(sim_mirt_hurdle(a_m3, d_m3, a_z_m3, d_z_m3, rho = rho, THETA = theta, N = N, estModel = FALSE))

reps1$obs_responses <- data.frame(reps1$obs_responses)
reps2$obs_responses <- data.frame(reps2$obs_responses)
reps3$obs_responses <- data.frame(reps3$obs_responses)
reps1$responses <- data.frame(reps1$responses)
reps2$responses <- data.frame(reps2$responses)
reps3$responses <- data.frame(reps3$responses)

# --------------------------------------------------------------------------
# Verify monotonically decreasing sum-score endorsement for ZI conditions
# (optional diagnostic — can be removed in production)
# --------------------------------------------------------------------------
check_monotone_decrease <- function(reps, label) {
  ss   <- rowSums(reps$obs_responses)
  freq <- table(factor(ss, levels = 0:max(ss)))
  diffs <- diff(as.numeric(freq))
  if (all(diffs <= 0)) {
    message(label, ": PASS — sum-score counts are monotonically non-increasing.")
  } else {
    warning(label, ": FAIL — sum-score counts are NOT monotonically non-increasing.")
  }
}
check_monotone_decrease(reps2, "Moderate ZI")
check_monotone_decrease(reps3, "Severe ZI")

# --------------------------------------------------------------------------
# Polychoric correlations and variance-covariance matrices (first 4 items)
# --------------------------------------------------------------------------
p1  <- reshape2::melt(psych::polychoric(reps1$obs_responses[,])$rho)
p2  <- reshape2::melt(psych::polychoric(reps2$obs_responses[,])$rho)
p3  <- reshape2::melt(psych::polychoric(reps3$obs_responses[,])$rho)

p1 = p1[as.numeric(p1$Var1) > as.numeric(p1$Var2), ]
p2 = p2[as.numeric(p2$Var1) > as.numeric(p2$Var2), ]
p3 = p3[as.numeric(p3$Var1) > as.numeric(p3$Var2), ]

p12 <- reshape2::melt(var(reps1$obs_responses[,]))
p22 <- reshape2::melt(var(reps2$obs_responses[,]))
p32 <- reshape2::melt(var(reps3$obs_responses[,]))

p12 <- p12[as.numeric(p12$Var1) >= as.numeric(p12$Var2), ]
p22 <- p22[as.numeric(p22$Var1) >= as.numeric(p22$Var2), ]
p32 <- p32[as.numeric(p32$Var1) >= as.numeric(p32$Var2), ]

# --------------------------------------------------------------------------
# Row sums for sum-score histograms
# --------------------------------------------------------------------------
reps1$responses$rowSums <- rowSums(reps1$obs_responses)
reps2$responses$rowSums <- rowSums(reps2$obs_responses)
reps3$responses$rowSums <- rowSums(reps3$obs_responses)

max.val <- (n_items * n_cat) + 1

# ---- shared theme snippet -------------------------------------------------
base_hist_theme <- theme(
  axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
  panel.grid   = element_blank(), panel.border  = element_blank(),
  axis.text    = element_text(color = "black"), text = element_text(size = 16)
)

# --------------------------------------------------------------------------
# Sum-score histograms
# --------------------------------------------------------------------------
## Remove extra bolus of 0 response
# grab 2k of these
rm.index = which(reps1$responses$rowSums == 0)
rm.index = sample(rm.index, size = 2500, replace = FALSE)

d1 <- ggplot(reps1$responses[-rm.index,], aes(x = rowSums)) +
  geom_histogram(bins = max.val, color = "white", fill = "black") +
  theme_bw() + base_hist_theme +
  xlab("Total Endorsement") + ylab("n") +
  coord_cartesian(ylim = c(0, 40000), xlim = c(0, max.val)) +
  ggtitle("None")

d2 <- ggplot(reps2$responses, aes(x = rowSums)) +
  geom_histogram(bins = max.val, color = "white", fill = "black") +
  theme_bw() + base_hist_theme +
  xlab("Total Endorsement") + ylab(" ") +
  theme(axis.text.y = element_text(color = "white")) +
  coord_cartesian(ylim = c(0, 40000), xlim = c(0, max.val)) +
  ggtitle("Moderate")

d3 <- ggplot(reps3$responses, aes(x = rowSums)) +
  geom_histogram(bins = max.val, color = "white", fill = "black") +
  theme_bw() + base_hist_theme +
  xlab("Total Endorsement") + ylab(" ") +
  theme(axis.text.y = element_text(color = "white")) +
  coord_cartesian(ylim = c(0, 40000), xlim = c(0, max.val)) +
  ggtitle("Severe")

# ---- shared tile theme snippet --------------------------------------------
base_tile_theme <- theme(
  legend.position = "NULL",
  axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
  panel.grid   = element_blank(), panel.border  = element_blank(),
  axis.text    = element_text(color = "black"), text = element_text(size = 16)
)

item_labels <- c(
  "Item_11" = "Item 1", "Item_12" = "Item 2",
  "Item_13" = "Item 3", "Item_14" = "Item 4"
)

# --------------------------------------------------------------------------
# Polychoric correlation heat-maps
# --------------------------------------------------------------------------
make_cor_tile <- function(dat, hide_y = FALSE) {
  p <- ggplot(dat, aes(x = Var1, y = Var2, fill = value)) +
    theme_bw() +
    geom_tile(color = "white", lwd = 1.5, linetype = 1) +
    geom_text(aes(label = round(value, 2)), color = "white", size = 4) +
    coord_fixed() + base_tile_theme +
    scale_fill_gradientn(colours = c('black', '#FFFFFF'), limits = c(0, 1.2)) +
    xlab(NULL) + ylab(NULL) +
    scale_x_discrete(labels = item_labels) +
    scale_y_discrete(labels = item_labels)
  if (hide_y) p <- p + theme(axis.text.y = element_text(color = "white"))
  p
}

c1 <- make_cor_tile(p1, hide_y = FALSE)
c2 <- make_cor_tile(p2, hide_y = TRUE)
c3 <- make_cor_tile(p3, hide_y = TRUE)

# --------------------------------------------------------------------------
# Variance-covariance heat-maps
# --------------------------------------------------------------------------
make_vcov_tile <- function(dat, hide_y = FALSE) {
  p <- ggplot(dat, aes(x = Var1, y = Var2, fill = value)) +
    theme_bw() +
    geom_tile(color = "white", lwd = 1.5, linetype = 1) +
    geom_text(aes(label = round(value, 2)), color = "white", size = 4) +
    coord_fixed() + base_tile_theme +
    scale_fill_gradientn(colours = c('black', '#FFFFFF'), limits = c(0, 4)) +
    xlab(NULL) + ylab(NULL) +
    scale_x_discrete(labels = item_labels) +
    scale_y_discrete(labels = item_labels) +
    theme(axis.text.x = element_text(color = "white"))
  if (hide_y) p <- p + theme(axis.text.y = element_text(color = "white"))
  p
}

c12 <- make_vcov_tile(p12, hide_y = FALSE)
c22 <- make_vcov_tile(p22, hide_y = TRUE)
c32 <- make_vcov_tile(p32, hide_y = TRUE)

## Change this so it creates a barplot of the average varince and covaraince
## across all items
make_vcov_tile <- function(dat, hide_y = FALSE) {
  ## Make the varaince & covar column
  dat$Value = "Covariance"
  dat$Value[which(dat$Var1==dat$Var2)] = "Variance"
  summary_df <- dat %>%
    group_by(Value) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      se_val = sd(value, na.rm = TRUE) / sqrt(n())
    )
  summary_df$Value = factor(summary_df$Value, levels = c("Variance", "Covariance"))
  p <- ggplot(summary_df, aes(x = Value, y = mean_val, group = Value)) +
    theme_bw() +
    #geom_tile(color = "white", lwd = 1.5, linetype = 1) +
    #geom_text(aes(label = round(value, 2)), color = "white", size = 4) +
    #coord_fixed() + base_tile_theme +
    geom_bar(stat = "identity", color ="black", fill = "black") +
    #scale_fill_gradientn(colours = c('black', '#FFFFFF'), limits = c(0, 2)) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), 
                  width = 0.2) +
    xlab(NULL) + ylab(NULL) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid   = element_blank(), panel.border  = element_blank(),
          axis.text    = element_text(color = "black"), text = element_text(size = 16)) +
    coord_cartesian(ylim= c(0.2, 1.3))
    #scale_x_discrete(labels = c("Variance", "Covaraince")) +
    #scale_y_discrete(labels = item_labels) +
    #theme(axis.text.x = element_text(color = "white"))
  #if (hide_y) p <- p + theme(axis.text.y = element_text(color = "white"))
  p
}

c12 <- make_vcov_tile(p12, hide_y = FALSE)
c22 <- make_vcov_tile(p22, hide_y = TRUE)
c32 <- make_vcov_tile(p32, hide_y = TRUE)

# --------------------------------------------------------------------------
# Combine Variance Covariance and Plycor values
# --------------------------------------------------------------------------
## First create Variance, Covariance, and cor identifiers
# dat$Value = "Covariance"
# dat$Value[which(dat$Var1==dat$Var2)] = "Variance"
p12$Value = "Covariance"
p22$Value = "Covariance"
p32$Value = "Covariance"
## Now make the varaince values
p12$Value[which(p12$Var1==p12$Var2)] = "Variance"
p22$Value[which(p22$Var1==p22$Var2)] = "Variance"
p32$Value[which(p32$Var1==p32$Var2)] = "Variance"
## Now make a source value
p12$Source = "None"
p22$Source = "Moderate"
p32$Source = "Severe"
## Now do the correlation values
p1$Value = "Correlation"
p2$Value = "Correlation"
p3$Value = "Correlation"
p1$Source = "None"
p2$Source = "Moderate"
p3$Source = "Severe"
## Now combine all of these for a plot
for_plot = bind_rows(p1, p2, p3, p12, p22, p32)

## Now make the plot for this
summary_df <- for_plot %>%
  group_by(Value, Source) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    se_val = sd(value, na.rm = TRUE)
  )
summary_df$Value = factor(summary_df$Value, levels = c("Variance", "Covariance", "Correlation"))
summary_df$Source = factor(summary_df$Source, levels = c("None", "Moderate", "Severe"))

out_plotValues = ggplot(summary_df, aes(x = Source, y = mean_val, group = Value,
                                        fill = Source)) +
  theme_minimal() +
  #geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  #geom_text(aes(label = round(value, 2)), color = "white", size = 4) +
  #coord_fixed() + base_tile_theme +
  geom_bar(stat = "identity") +#, color ="black", fill = "black") +
  #scale_fill_gradientn(colours = c('black', '#FFFFFF'), limits = c(0, 2)) +
  geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val), 
                width = 0.2) +
  xlab(NULL) + ylab(NULL) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid   = element_blank(), panel.border  = element_blank(),
        axis.text    = element_text(color = "black"), 
        text = element_text(size = 16),
        legend.position = "NULL") +
  coord_cartesian(ylim= c(0.2, 1.5)) +
  facet_grid(. ~ Value) +
  scale_fill_grey()

# --------------------------------------------------------------------------
# Combine and save Figure 1
# --------------------------------------------------------------------------
out.plot1 <- ggpubr::ggarrange(d1, d2, d3, nrow = 1, ncol = 3)
out.plot1 <- annotate_figure(out.plot1, left = "Sum Scores")

out.plot2 <- ggpubr::ggarrange(c12, c22, c32, nrow = 1, ncol = 3)
out.plot2 <- annotate_figure(out.plot2, left = "Variance-Covariance Values")

out.plot3 <- ggpubr::ggarrange(c1, c2, c3, nrow = 1, ncol = 3)
out.plot3 <- annotate_figure(out.plot3, left = "Correlation Matrix")

out.plot <- ggpubr::ggarrange(out.plot1, out.plot2, out.plot3, nrow = 3)
ggsave("./reports/figure1ExampleCovar.png",          out.plot, dpi = 300, width = 12, height = 13, units = "in", bg = "white")
ggsave("./texOut/figures/figure1ExampleCovar.png",   out.plot, dpi = 300, width = 12, height = 13, units = "in", bg = "white")
ggsave("./scripts/beamer_pres/images/resp_dist.png", out.plot1, dpi = 300, width = 12, height = 8,  units = "in", bg = "white")
ggsave("./scripts/beamer_pres/images/varcov_dist.png", out.plot2, dpi = 300, width = 12, height = 8, units = "in", bg = "white")
ggsave("./scripts/beamer_pres/images/cor_dist.png",  out.plot3, dpi = 300, width = 12, height = 8,  units = "in", bg = "white")

## Now produce a reduced figure for ouhsc
out.plot.hsc <- ggpubr::ggarrange(out.plot1, out.plot3, nrow = 2)
ggsave("./scripts/beamer_pres/images/simulate_hurdle_data.png",   out.plot.hsc, dpi = 300, width = 10, height = 5, units = "in", bg = "white")

## Now make the reduced version of this figure
out.plot_rd = ggarrange(out.plot1, out_plotValues, nrow = 2)
ggsave("./reports/figure1ExampleCovarUpdate.png", out.plot_rd, 
       dpi = 300, width = 12, height = 8, units = "in", bg = "white")
ggsave("./texOut/figures/figure1ExampleCovarUpdate.png", out.plot_rd, 
       dpi = 300, width = 12, height = 8, units = "in", bg = "white")

# --------------------------------------------------------------------------
# Factor analysis (polychoric and Pearson)
# --------------------------------------------------------------------------
library(psych)
fa.1 <- fa(reps1$obs_responses[,], cor = "poly", nfactors = 1)
fa.2 <- fa(reps2$obs_responses[,], cor = "poly", nfactors = 1)
fa.3 <- fa(reps3$obs_responses[,], cor = "poly", nfactors = 1)

fa_load <- data.frame(
  items      = 1:10,
  None = as.numeric(loadings(fa.1)),  # Convert to numeric
  Moderate = as.numeric(loadings(fa.2)),
  Severe = as.numeric(loadings(fa.3))
)


fa_load_long <- reshape2::melt(fa_load, id.vars = "items")

# ---- legend extractor -----------------------------------------------------
g_legend <- function(a.gplot) {
  tmp    <- ggplot_gtable(ggplot_build(a.gplot))
  leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

# Calculate means and standard errors
fa_load_summary <- fa_load_long %>%
  group_by(variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    se = sd(value, na.rm = TRUE),
    .groups = 'drop'
  ) %>% 
  mutate(
    source = "Loading"
  )

# Bar graph for loadings
fa_vals = ggplot(fa_load_summary, aes(x = variable, y = mean, fill = variable)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, linewidth = 1) +
  labs(#title = "Factor Loadings", 
       x = "Zero Inflation", 
       y = "Mean Loading") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid   = element_blank(), panel.border  = element_blank(),
        axis.text    = element_text(color = "black"), text = element_text(size = 16)) +
  coord_cartesian(ylim= c(0.6, .8)) +
  scale_fill_grey()+
  guides(color = guide_legend(title = "Zero Inflation "), fill = guide_legend(title = "Zero Inflation "))


# --------------------------------------------------------------------------
# IRT models and test information functions
# --------------------------------------------------------------------------
mi.1 <- mirt(reps1$obs_responses, 1, itemtype = "graded")
mi.2 <- mirt(reps2$obs_responses, 1, itemtype = "graded")
mi.3 <- mirt(reps3$obs_responses, 1, itemtype = "graded")

theta_seq <- seq(-7, 7, .1)
test.info.plot <- data.frame(
  theta  = theta_seq,
  None     = mirt::testinfo(mi.1, theta_seq),
  Moderate = mirt::testinfo(mi.2, theta_seq),
  Severe   = mirt::testinfo(mi.3, theta_seq)
)
ti_long <- reshape2::melt(test.info.plot, id.vars = "theta")
ti_long$Zero_Inflation <- ti_long$variable

ti1 <- ggplot(ti_long, aes(x = theta, y = value, group = Zero_Inflation,
                           color = Zero_Inflation)) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(color = "black"))+#, panel.grid = element_blank()) +
  guides(color = guide_legend(title = "Zero Inflation "), 
         fill = guide_legend(title = "Zero Inflation ")) +
  theme_bw() +
  theme(#panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    #panel.grid   = element_blank(), panel.border  = element_blank(),
    axis.text    = element_text(color = "black"), text = element_text(size = 16)) +
  scale_color_grey() +
  geom_line(size = 2) +
  coord_cartesian(ylim = c(0, 15)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Test Information") + xlab("Theta Value")# +
  #ggtitle("Test Information Functions")


# legend from a version with legend shown
ti1p <- ggplot(ti_long, aes(x = theta, y = value, group = Zero_Inflation,
                            color = Zero_Inflation)) +
  geom_line(size = 4) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_bw() +
  ylab("Test Information") + xlab("Indicator") +
  scale_color_grey() +
  guides(color = guide_legend(title = "Extent of Zero Inflation"))

legend.val <- g_legend(ti1p)

out.plot2b <- ggarrange(fa_vals, ti1, nrow = 1, ncol = 2, common.legend = TRUE, legend = "bottom")
ggsave("./reports/figure2ExampleFAandInfo.png",        out.plot2b, dpi = 300, width = 12, height = 7, units = "in", bg = "white")
ggsave("./texOut/figures/figure2ExampleFAandInfo.png", out.plot2b, dpi = 300, width = 11, height = 6, units = "in", bg = "white")
ggsave("./scripts/beamer_pres/images/factor_loading_simulate.png", out.plot2b, dpi = 300, width = 9, height =5, units="in")

# --------------------------------------------------------------------------
# Reliability estimates
# --------------------------------------------------------------------------
omega.1 <- psych::omega(reps1$obs_responses, poly = TRUE,  plot = FALSE, nfactors = 1)
omega.2 <- psych::omega(reps2$obs_responses, poly = TRUE,  plot = FALSE, nfactors = 1)
omega.3 <- psych::omega(reps3$obs_responses, poly = TRUE,  plot = FALSE, nfactors = 1)

alpha.1 <- psych::alpha(psych::polychoric(reps1$obs_responses)$rho)
alpha.2 <- psych::alpha(psych::polychoric(reps2$obs_responses)$rho)
alpha.3 <- psych::alpha(psych::polychoric(reps3$obs_responses)$rho)

grm.1 <- grm_reliability_emp(mi.1, method = "EAP")
grm.2 <- grm_reliability_emp(mi.2, method = "EAP")
grm.3 <- grm_reliability_emp(mi.3, method = "EAP")

estHurdleRel1 <- severity_reliability_GH(sv2 = reps1$tru_mod)$rho
estHurdleRel2 <- severity_reliability_GH(sv2 = reps2$tru_mod)$rho
estHurdleRel3 <- severity_reliability_GH(sv2 = reps3$tru_mod)$rho

# --------------------------------------------------------------------------
# Figure 3 – reliability estimates
# --------------------------------------------------------------------------
to.plot <- data.frame(
  Zero_Inflation = 1:3,
  Alpha  = c(alpha.1$total$std.alpha, alpha.2$total$std.alpha, alpha.3$total$std.alpha),
  Omega  = c(omega.1$omega.tot,       omega.2$omega.tot,       omega.3$omega.tot),
  GRM    = c(grm.1, grm.2, grm.3),
  Hurdle = c(estHurdleRel1, estHurdleRel2, estHurdleRel3)
)

tp3 <- reshape2::melt(to.plot, id.vars = "Zero_Inflation")
tp3$Zero_Inflation <- factor(tp3$Zero_Inflation)
tp3$Zero_Inflation <- plyr::revalue(
  tp3$Zero_Inflation,
  replace = c("1" = "None", "2" = "Moderate", "3" = "Severe")
)

f1a <- ggplot(tp3, aes(x = variable, y = value,
                       group = Zero_Inflation, fill = Zero_Inflation)) +
  theme_bw() +
  geom_bar(stat = 'identity', position = 'dodge') +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Estimated Reliability") + xlab("Reliability Metric") +
  scale_fill_grey() +
  theme(axis.ticks = element_blank(),
        axis.text  = element_text(color = "black")) +
  labs(fill = 'Zero Inflation')

f1a
ggsave("./reports/figure3ExampleRelEsts.png",        f1a, dpi = 300, width = 8, height = 6, units = "in", bg = "white")
ggsave("./texOut/figures/figure3ExampleRelEsts.png", f1a, dpi = 300, width = 8, height = 6, units = "in", bg = "white")
ggsave("./scripts/beamer_pres/images/reliability_simulate.png", f1a, dpi = 300, width = 9, height =5, units="in")

