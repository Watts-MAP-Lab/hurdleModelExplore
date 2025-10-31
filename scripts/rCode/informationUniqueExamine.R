##### --load-library-------------
source("./scripts/rCode/hurdleFunctions.r")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(mirt))
suppressPackageStartupMessages(library(psych))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(MASS))   # mvrnorm

##### --declare-sim-params-------
source("./scripts/rCode/simParam.r")
seedVal <- as.integer(commandArgs(TRUE))
if (is.na(seedVal)) seedVal <- 1
set.seed(seedVal)

### --check-run-status----------
out.file <- paste0("./data/hurdleCollapse/seedVal_", seedVal, "_take3.RDS")
if (file.exists(out.file)) {
  message("Already done for seed ", seedVal)
  quit(save = "no")
}

##### --qol-funs-and-times-------
# Simple timestamped logger
log_msg <- function(..., level = "INFO") {
  cat(sprintf("[%s] %-5s | %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              level,
              paste(..., collapse = " ")))
}

# Timing wrapper
time_it <- function(label, expr) {
  log_msg("START:", label)
  t <- system.time(res <- eval.parent(substitute(expr)))
  log_msg("DONE:", label, "| elapsed", sprintf("%.2f s", t["elapsed"]))
  invisible(res)
}

# Format hh:mm:ss.s
format_elapsed <- function(start_time, end_time) {
  secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
  sprintf("%02d:%02d:%05.2f", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
}

start_time <- Sys.time()
log_msg("SCRIPT START for seed", seedVal)
log_msg("Output:", out.file)

##### --simulate-one-scenario----
# Pull parameters from simParam.r for this seed
n_items   <- all.sim.vals$nItem[seedVal]
n_cat     <- all.sim.vals$nCat[seedVal]
rho       <- all.sim.vals$facCor[seedVal]
N         <- all.sim.vals$n[seedVal]
add.val.2pl <- 1.5
add.val.grm <- 1.5

# Draw item parameters (traditional)
a_grm <- runif(n_items, min = all.sim.vals$grmDiscrim[seedVal],
               max = all.sim.vals$grmDiscrim[seedVal] + add.val.grm)
b_grm <- genDiffGRM(num_items = n_items, num_categories = n_cat,
                    min = all.sim.vals$difGrmF[seedVal],
                    max = all.sim.vals$difGrmF[seedVal] + 2.5,
                    rnorm_var = 0.3)
a_2pl <- runif(n_items, min = all.sim.vals$discrim2pl[seedVal],
               max = all.sim.vals$discrim2pl[seedVal] + add.val.2pl)
b_2pl <- runif(n_items, min = all.sim.vals$dif2PL[seedVal],
               max = all.sim.vals$dif2PL[seedVal] + 2.0)

# Convert to mirt parameterization
mirt_pars <- to_mirt(a_grm, b_grm, a_2pl, b_2pl, n_cat)
a_m   <- mirt_pars$a_m
d_m   <- mirt_pars$d_m
a_z_m <- mirt_pars$a_z_m
d_z_m <- mirt_pars$d_z_m

# Generate true theta from N(0, Sigma_pop)
Sigma_pop <- matrix(c(1, rho, rho, 1), ncol = 2)
theta     <- MASS::mvrnorm(n = N, mu = c(0, 0), Sigma = Sigma_pop)

# Simulate with your hurdle function (fixed-parameter true model)
reps1 <- suppressWarnings(sim_mirt_hurdle(a_m, d_m, a_z_m, d_z_m, rho = rho,
                                          THETA = theta, N = N))

# Fit an estimated hurdle model (sv1) to the NA-ed responses (this uses estimation;
# if you want scoring-only with true parameters, you can skip this and use reps1$tru_mod)
n_2pl <- n_items
n_item_total <- n_items * 2
model_str <- make_model(n_2pl, n_item_total)
item.type.rep <- c(rep("2PL", n_items), rep("graded", n_items))
sv1 <- mirt(as.data.frame(reps1$responses_na), model = model_str, itemtype = item.type.rep)

##### --reliability metrics------
vals_loop <- data.frame(seed = seedVal)

# Classic reliability metrics (alpha, omega, etc.) on full responses (no hurdle NA)
test_dat <- as.data.frame(reps1$responses)
cor.mat  <- psych::polychoric(test_dat)
rel.all  <- psych::reliability(cor.mat$rho, n.obs = N, nfactors = 3, plot = FALSE)

# Extract reliability stats (vectorized assignment)
stats <- rel.all$result.df
vals_loop <- cbind(vals_loop, t(stats[, c("omega_h","alpha","omega.tot","Uni","tau",
                                          "cong","CFI","ECV","Beta","EVR","MAP")]))
vals_loop$skewVal <- as.numeric(psych::describe(rowSums(reps1$responses))["skew"])

# Split-half (Guttman)
split.half.rel <- suppressWarnings(psych::guttman(r = cor.mat$rho))
vals_loop$lambda1Rel <- split.half.rel$lambda.1
vals_loop$lambda2Rel <- split.half.rel$lambda.2
vals_loop$lambda3Rel <- split.half.rel$lambda.3
vals_loop$lambda4Rel <- split.half.rel$lambda.4
vals_loop$lambda5Rel <- split.half.rel$lambda.5
vals_loop$lambda6Rel <- split.half.rel$lambda.6

# Single-factor omega via SEM (optional; can be slow)
omega.sem.val <- try(omegaSem(m = cor.mat$rho, nfactors = 1, plot = FALSE, n.obs = N), silent = TRUE)
if (!inherits(omega.sem.val, "try-error")) {
  vals_loop$singleFactorOmegaT <- omega.sem.val$omegaSem$omega.tot
} else {
  vals_loop$singleFactorOmegaT <- NA_real_
}

# MAP-based hurdle reliability via GH (true model with fixed parameters)
vals_loop$trueInfoRel <- severity_reliability_GH(sv2 = reps1$tru_mod)$rho

# Smoother-based empirical reliability (uses true thetas and MAP scores)
tmp.dat <- reps1$theta_vals 
rel.sev <- mgcv::gam(truP_sevMAP ~ s(trueSev) + s(trueSus), data = tmp.dat, REML = TRUE)
rel.sus <- mgcv::gam(truP_susMAP ~ s(trueSev) + s(trueSus), data = tmp.dat, REML = TRUE)
vals_loop$trueGamRelSev <- summary(rel.sev)$r.sq  # approx Var(g)/Var(est)
vals_loop$trueGamRelSus <- summary(rel.sus)$r.sq

## Variance from true score and true model EAP scores
true.score.var <- var(reps1$theta_vals$trueSev)
error.var <- var(reps1$theta_vals$truP_sevMAP - reps1$theta_vals$trueSev)
vals_loop$trueVarRel <- true.score.var / (true.score.var + error.var)

# GRM models for comparison (all items; and excluding zeros)
mod_grm <- mirt::mirt(as.data.frame(reps1$responses), 1, itemtype = "graded")
vals_loop$grmRel <- grm_reliability_emp(mod_grm, method = "EAP")

# Estimated hurdle model
vals_loop$estInfoRel <- severity_reliability_GH(sv1)$rho

# Remove structural zeros by selecting graded-only columns (>0)
iso.col <- (n_items + 1):ncol(reps1$responses_na)
mod_rm  <- mirt::mirt(as.data.frame(reps1$responses_na[, iso.col, drop = FALSE]), 1, itemtype = "graded")
vals_loop$grmRel_rmZeroOption <- grm_reliability_emp(mod_rm, method = "EAP")

##### --collect-vals----------
# Augment theta table with quick summaries (optional; avoid heavy objects if memory is tight)
theta_tab <- reps1$theta_vals
theta_tab$rowSums       <- rowSums(reps1$responses)
theta_tab$grmFacScore   <- fscores(mod_grm)
theta_tab$grmFacScoreRM <- fscores(mod_rm)
fs_sv1 <- fscores(sv1, method = "MAP")
est_df <- data.frame(sus = fs_sv1[, 1], sev = fs_sv1[, 2])
theta_tab$estSus        <- est_df$sus
theta_tab$estSev        <- est_df$sev

##### --summary print-------------
log_msg(sprintf("Reliabilities summary (seed=%d):", seedVal))
log_msg(sprintf("  Plugin MAP rho=%.4f | Var(MAP)=%.4f | E_err=%.4f", vals_loop$pluginRel, vals_loop$plugin_VarHat, vals_loop$plugin_Eerr))
log_msg(sprintf("  GAM sev r^2≈rho=%.4f; GAM sus r^2≈rho=%.4f", vals_loop$trueRelSevGam, vals_loop$trueRelSusGam))
log_msg(sprintf("  GRM all rho=%.4f; GRM nonzero rho=%.4f", vals_loop$grmRel, vals_loop$grmRel_rmZeroOption))

##### --save + DONE---------------
out.list <- list(
  allParams = vals_loop,
  thetaVals = theta_tab,
  est_mod   = sv1
)

time_it("saveRDS", saveRDS(out.list, file = out.file, compress = TRUE))

end_time <- Sys.time()
elapsed  <- format_elapsed(start_time, end_time)

log_msg("SCRIPT DONE for seed", seedVal, "| wrote:", out.file)
log_msg("TOTAL ELAPSED:", elapsed)


# Save compact output
out.list <- list(
  allParams = vals_loop,
  thetaVals = theta_tab,
  est_mod   = sv1  # if memory is tight, you can remove this or save only parameter tables
)

saveRDS(out.list, file = out.file, compress = TRUE)