# =============================================================
# Step 12: Parametric calibration experiments (referee comments 1-4)
#
# Script 11 showed the excess FWER is NOT numerical (16% at 400, 1000,
# 4000 samples and with the Gaussian method). This script asks WHY, by
# simulating data from known models and fitting the same SPDE model.
#
# All simulations use hyperparameters estimated from the REAL annual data
# (fit once at the start, standardised scale, coarse 670-node mesh), and
# simulate the fields on the same mesh, so the generating model is exactly
# the fitted model class.
#
# Scenarios (each answers a different question - do not interpret them
# interchangeably):
#   P0b  zero trend, data = w_t + noise on the standardised scale,
#        fitted with s KNOWN (= correctly specified model).
#        -> frequentist null FWER of the joint sets when the model is right.
#   P0a  as P0b but multiplied by the real cell SDs s_i and fitted with the
#        usual v2 procedure (s_i ESTIMATED from the data).
#        -> adds the effect of the standardisation step.
#   P1b  trend FIELD drawn each replicate from the fitted trend GRF
#        (b0 = 0), s known, correctly specified.
#        -> calibration averaged over the (hierarchical) prior of the trend
#           field: P(any wrong-sign cell in the joint set). NOT the same
#           question as P0.
#   P2a  zero trend, heavy-tailed years: each year's anomaly map is scaled
#        by sqrt(4 / chi2_4) (t4-type years); v2 procedure.
#   P2b  zero trend, year-to-year dependence: w_t = 0.3 w_{t-1} + innovation
#        (AR(1), variance kept); v2 procedure.
#        -> misspecified models: which real-data feature could cause excess?
#
# Per replicate and scenario: ONE fit, 4000 posterior samples, stored per
# cell: truth, posterior mean, 95% interval, excursion functions F (neg, pos,
# up to alpha = 0.10) and the alpha = 0.025 sets. FWER at any alpha <= 0.10
# can be recomputed from F -> alpha tuning can be validated out-of-sample.
#
# Order: replicate r runs all scenarios before r+1, so partial results are
# balanced across scenarios. P0a/P0b/P1b go to 100 replicates, P2a/P2b to 50.
# Runtime: ~2.5-3 min per fit -> ~12 min per replicate while all 5 run
# (first 50 reps ~10 h), then ~8 min (reps 51-100 ~7 h). Stop with Esc any
# time; finished fits are skipped on restart.
# Output: results/param/rep_XXX_<scenario>.rds, combined at the end into
#         results/param/12_cells.csv and 12_hyper.csv
# =============================================================

N_MAX   <- c(P0b = 100, P0a = 100, P1b = 100, P2a = 50, P2b = 50)
SCEN    <- names(N_MAX)            # run order within each replicate
N_SAMP  <- 4000
F_LIM   <- 0.10
SMOKE   <- TRUE                    # TRUE = 1 replicate of every scenario (~15 min); then set FALSE

suppressPackageStartupMessages({
  library(INLA); library(fmesher); library(excursions); library(Matrix)
  library(dplyr); library(tidyr); library(readr)
})
inla.setOption(num.threads = "2:1")
out <- file.path("results", "param"); dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---------------- data and structures (identical to scripts 08/11) ----------------
agg <- read_csv("results/rev/A0_seasonal_means_dayweighted.csv", show_col_types = FALSE)
cells <- read_csv("results/rev/cells.csv", show_col_types = FALSE)
cells <- cells |> mutate(X_km = (LON - 82) * 111.32 * cos(22 * pi / 180), Y_km = (LAT - 22) * 110.57)
nc <- nrow(cells); loc <- as.matrix(cells[, c("X_km", "Y_km")])
yrs <- 2001:2023
Yann <- agg |> filter(SEASON == "Annual", YEAR %in% yrs) |>
  left_join(cells |> select(LON, LAT, cell), by = c("LON", "LAT")) |>
  arrange(cell, YEAR) |> select(cell, YEAR, GHI) |>
  pivot_wider(names_from = YEAR, values_from = GHI) |> arrange(cell) |>
  select(-cell) |> as.matrix()
mesh <- fm_mesh_2d_inla(loc = loc, max.edge = c(200, 600), cutoff = 80, offset = c(200, 700))
cat("Coarse mesh nodes:", mesh$n, "(must be 670)\n")
if (mesh$n != 670) stop("Mesh differs from scripts 08/11. Stop and check the setup.")
T <- length(yrs); tc <- (yrs - mean(range(yrs))) / 10
dd <- expand.grid(cell = 1:nc, yr = 1:T); obs_loc <- loc[dd$cell, ]
A_b <- inla.spde.make.A(mesh, loc = obs_loc, weights = tc[dd$yr])
A_w <- inla.spde.make.A(mesh, loc = obs_loc, repl = dd$yr, n.repl = T)
A_pred <- inla.spde.make.A(mesh, loc = loc)
pr <- list(br = c(300, 0.05), bs = c(2, 0.05), wr = c(200, 0.05), ws = c(3, 0.05))  # v2 priors
spde_b <- inla.spde2.pcmatern(mesh, prior.range = pr$br, prior.sigma = pr$bs)
spde_w <- inla.spde2.pcmatern(mesh, prior.range = pr$wr, prior.sigma = pr$ws)
s_real <- apply(Yann, 1, function(y) sd(resid(lm(y ~ tc))))

# Fit the SPDE trend model. std = TRUE: estimate s_i from the data (usual v2);
# std = FALSE: data are already on the standardised scale (s known = 1).
fit_model <- function(Y, std) {
  s <- if (std) apply(Y, 1, function(y) sd(resid(lm(y ~ tc)))) else rep(1, nc)
  Z <- (Y - rowMeans(Y)) / s
  idx_b <- inla.spde.make.index("beta", spde_b$n.spde)
  idx_w <- inla.spde.make.index("w", spde_w$n.spde, n.repl = T)
  stk <- inla.stack(data = list(y = Z[cbind(dd$cell, dd$yr)]), A = list(A_b, A_w, 1),
                    effects = list(idx_b, idx_w, data.frame(b0 = tc[dd$yr])))
  form <- y ~ -1 + b0 + f(beta, model = spde_b) + f(w, model = spde_w, replicate = w.repl)
  fit <- inla(form, data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.compute = list(config = TRUE),
              control.inla = list(int.strategy = "eb"))
  list(fit = fit, s = s)
}
draw_trends <- function(m, n, seed) {
  smp <- inla.posterior.sample(n, m$fit, seed = seed,
                               selection = list(beta = 1:spde_b$n.spde, b0 = 1))
  rn <- rownames(smp[[1]]$latent); ib <- grep("^beta:", rn); i0 <- grep("^b0:", rn)
  stopifnot(length(ib) == spde_b$n.spde, length(i0) == 1)
  Bm <- sapply(smp, function(z) z$latent[ib, 1]); b0 <- sapply(smp, function(z) z$latent[i0, 1])
  (as.matrix(A_pred %*% Bm) + matrix(b0, nc, n, byrow = TRUE)) * m$s
}
hyp <- function(fit) {
  h <- fit$summary.hyperpar
  tibble(prec_eps = h["Precision for the Gaussian observations", "mode"],
         range_b = h["Range for beta", "mode"], sd_b = h["Stdev for beta", "mode"],
         range_w = h["Range for w", "mode"], sd_w = h["Stdev for w", "mode"])
}

# ---------------- generating hyperparameters from the real annual data ----------------
f_true <- file.path(out, "12_true_hyper.csv")
if (!file.exists(f_true)) {
  cat("Fitting real annual data once to get generating hyperparameters...\n")
  H <- hyp(fit_model(Yann, std = TRUE)$fit)
  write_csv(H, f_true)
}
H <- read_csv(f_true, show_col_types = FALSE); print(as.data.frame(H))
Q_w <- inla.spde2.precision(spde_w, theta = c(log(H$range_w), log(H$sd_w)))   # pcmatern: log range, log sigma
Q_b <- inla.spde2.precision(spde_b, theta = c(log(H$range_b), log(H$sd_b)))
sd_eps <- 1 / sqrt(H$prec_eps)
# Safety check: simulated fields must have roughly the estimated SD at the cells
chk_w <- sd(as.matrix(A_pred %*% inla.qsample(20, Q_w, seed = 99)))
chk_b <- sd(as.matrix(A_pred %*% inla.qsample(20, Q_b, seed = 98)))
cat(sprintf("Check: simulated SD w = %.3f (target %.3f) | beta = %.3f (target %.3f)\n", chk_w, H$sd_w, chk_b, H$sd_b))
if (chk_w / H$sd_w > 2 || chk_w / H$sd_w < 0.5 || chk_b / H$sd_b > 2 || chk_b / H$sd_b < 0.5)
  stop("Simulated field SD does not match the target - parametrisation problem. Stop and check the setup.")

# Simulate anomaly maps on the standardised scale (nc x T); zero mean
sim_w <- function(seed, heavy = FALSE, ar = 0) {
  X <- inla.qsample(T, Q_w, seed = seed)                         # mesh x T, independent years
  if (ar > 0) for (t in 2:T) X[, t] <- ar * X[, t - 1] + sqrt(1 - ar^2) * X[, t]
  W <- as.matrix(A_pred %*% X)
  if (heavy) W <- sweep(W, 2, sqrt(4 / rchisq(T, 4)), `*`)
  W + matrix(rnorm(nc * T, 0, sd_eps), nc, T)
}
sim_scenario <- function(sc, r) {
  sid <- match(sc, SCEN); set.seed(20000 + 1000 * sid + r); qs <- 300000 + 1000 * sid + r
  truth_z <- rep(0, nc)
  if (sc == "P0b") Z <- sim_w(qs)
  if (sc == "P0a") Z <- sim_w(qs)
  if (sc == "P1b") { truth_z <- as.numeric(A_pred %*% inla.qsample(1, Q_b, seed = qs + 500000))
                     Z <- sim_w(qs) + outer(truth_z, tc) }
  if (sc == "P2a") Z <- sim_w(qs, heavy = TRUE)
  if (sc == "P2b") Z <- sim_w(qs, ar = 0.3)
  std <- sc %in% c("P0a", "P2a", "P2b")
  Y <- if (std) Z * s_real else Z                                # v2 procedure works on the original scale
  list(Y = Y, std = std, truth = if (std) truth_z * s_real else truth_z)
}

rep_file <- function(r, sc) file.path(out, sprintf("rep_%03d_%s.rds", r, sc))
reps <- if (SMOKE) 1 else 1:max(N_MAX)
t_all <- Sys.time(); k_done <- 0
cat("Started:", format(t_all), if (SMOKE) "(SMOKE TEST)" else "", "\n")

for (r in reps) for (sc in SCEN) {
  if (r > N_MAX[[sc]] || file.exists(rep_file(r, sc))) next
  t0 <- Sys.time()
  d <- sim_scenario(sc, r)
  m <- tryCatch(fit_model(d$Y, d$std), error = function(e) { cat("fit error", sc, r, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) next
  tr <- draw_trends(m, N_SAMP, seed = 400000 + 1000 * match(sc, SCEN) + r)
  exn <- excursions.mc(tr, alpha = F_LIM, u = 0, type = "<")
  exp_ <- excursions.mc(tr, alpha = F_LIM, u = 0, type = ">")
  q <- t(apply(tr, 1, quantile, c(0.025, 0.975)))
  cellres <- tibble(rep = r, scenario = sc, cell = 1:nc, truth = d$truth, est = rowMeans(tr),
                    q025 = q[, 1], q975 = q[, 2], Fneg = as.numeric(exn$F), Fpos = as.numeric(exp_$F),
                    jneg025 = excursions.mc(tr, 0.025, 0, "<")$E == 1,
                    jpos025 = excursions.mc(tr, 0.025, 0, ">")$E == 1)
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  res <- list(cells = cellres, hyper = hyp(m$fit) |> mutate(rep = r, scenario = sc, minutes = mins, .before = 1))
  tmp <- paste0(rep_file(r, sc), ".tmp"); saveRDS(res, tmp); file.rename(tmp, rep_file(r, sc))
  k_done <- k_done + 1
  todo <- sum(sapply(SCEN, function(s2) sum(!file.exists(sapply(reps[reps <= N_MAX[[s2]]], rep_file, sc = s2)))))
  per <- as.numeric(difftime(Sys.time(), t_all, units = "mins")) / k_done
  cat(sprintf("rep %d %s done in %.1f min | %d fits left | ETA %.1f h\n", r, sc, mins, todo, per * todo / 60))
}

fs <- list.files(out, pattern = "^rep_\\d{3}_P.*\\.rds$", full.names = TRUE)
if (length(fs)) {
  all <- lapply(fs, readRDS)
  write_csv(bind_rows(lapply(all, `[[`, "cells")), file.path(out, "12_cells.csv"))
  write_csv(bind_rows(lapply(all, `[[`, "hyper")), file.path(out, "12_hyper.csv"))
}
writeLines(capture.output(sessionInfo()), file.path(out, "12_sessionInfo.txt"))
cat("\nFinished:", format(Sys.time()), "| fits saved:", length(fs), "\n")
if (SMOKE) cat("Smoke test OK -> set SMOKE <- FALSE and source again.\n") else
  cat("Run finished.\n")
