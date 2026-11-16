# =============================================================
# Step 11: Posterior-sample (Monte Carlo) convergence of joint excursion sets
#
# Question: is part of the excess FWER of the joint excursion sets
# (14-17.5% at nominal 5%, task B of script 08) caused by computing the
# sets from too few posterior samples (400)? This is a HYPOTHESIS: more
# samples may stabilise the sets without restoring nominal FWER.
#
# Design (per null replicate r = same permutation as task B, seed 5000 + r):
#   * ONE INLA fit (v2, coarse mesh, EB) - so sampling effects are isolated
#     from model-fitting variation.
#   * Two INDEPENDENT batches of posterior samples, A and B, 4000 each.
#   * Excursion sets from nested subsets of each batch: n = 400, 1000, 4000,
#     at alpha = 0.005, 0.01, 0.025, 0.05 per sign, plus the excursion
#     function F (stored per cell) so FWER at any alpha <= 0.10 can be
#     recomputed later without refitting.
#   * A third independent batch C (4000) is used ONLY to estimate a Gaussian
#     approximation (mean + covariance) of the 274 cell trends. Sets from the
#     Gaussian method (excursions(), "EB") use it. This is an approximation based
#     on posterior samples, not an exact reference, but it is independent of A, B.
#   * For every selected set: its joint probability
#       - in-sample (same samples that selected it; >= 1 - alpha by construction)
#       - out-of-sample from the OTHER batch (independent 4000 samples), with its
#         Monte Carlo standard error
#       - by Gaussian integration (gaussint) under the batch-C approximation,
#         with gaussint's own error estimate
#     A drop from in-sample to out-of-sample, larger than its Monte Carlo error,
#     indicates selection bias.
#
# Output (results/mc/): one rep_XXX_v2.rds per finished replicate (written
# atomically, so stopping and restarting never duplicates rows), combined at
# the end into 11_sets.csv, 11_F.csv, 11_timing.csv.
# Runtime: about 3-5 min per replicate -> 50 replicates ~ 3-4 h.
# Laptop: threads limited to 2:1; you can stop with Esc at any time and
# source the script again later - finished replicates are skipped.
# =============================================================

N_REP   <- 50                      # first 50 permutations of task B
SIZES   <- c(400, 1000, 4000)      # nested subsets of each 4000-sample batch
N_BATCH <- 4000
ALPHAS  <- c(0.005, 0.01, 0.025, 0.05)
F_LIM   <- 0.10                    # excursion function stored up to alpha = 0.10
MODELS  <- c("v2")                 # add "v1" later if time allows
GAUSS   <- TRUE                    # sets from Gaussian approximation (independent batch C)
SMOKE   <- FALSE                   # TRUE = quick 2-replicate test (done, passed 18 Oct)

suppressPackageStartupMessages({
  library(INLA); library(fmesher); library(excursions); library(Matrix)
  library(dplyr); library(tidyr); library(readr)
})
inla.setOption(num.threads = "2:1")
out <- file.path("results", "mc"); dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---------------- data and structures (identical to script 08) ----------------
agg <- read_csv("results/rev/A0_seasonal_means_dayweighted.csv", show_col_types = FALSE)
cells <- read_csv("results/rev/cells.csv", show_col_types = FALSE)
# recompute coordinates exactly as script 08 did (the CSV is rounded -> mesh would differ by 1 node)
cells <- cells |> mutate(X_km = (LON - 82) * 111.32 * cos(22 * pi / 180), Y_km = (LAT - 22) * 110.57)
nc <- nrow(cells); loc <- as.matrix(cells[, c("X_km", "Y_km")])
yrs <- 2001:2023
Yann <- agg |> filter(SEASON == "Annual", YEAR %in% yrs) |>
  left_join(cells |> select(LON, LAT, cell), by = c("LON", "LAT")) |>
  arrange(cell, YEAR) |> select(cell, YEAR, GHI) |>
  pivot_wider(names_from = YEAR, values_from = GHI) |> arrange(cell) |>
  select(-cell) |> as.matrix()
mesh <- fm_mesh_2d_inla(loc = loc, max.edge = c(200, 600), cutoff = 80, offset = c(200, 700))
cat("Coarse mesh nodes:", mesh$n, "(must be 670 = same as task B)\n")
if (mesh$n != 670) stop("Mesh differs from task B (", mesh$n, " nodes). Stop and check the setup.")
T <- length(yrs); tc <- (yrs - mean(range(yrs))) / 10
dd <- expand.grid(cell = 1:nc, yr = 1:T); obs_loc <- loc[dd$cell, ]
A_b <- inla.spde.make.A(mesh, loc = obs_loc, weights = tc[dd$yr])
A_w <- inla.spde.make.A(mesh, loc = obs_loc, repl = dd$yr, n.repl = T)
A_pred <- inla.spde.make.A(mesh, loc = loc)
R <- t(apply(Yann, 1, function(y) resid(lm(y ~ tc))))     # detrended residual maps

priors <- list(v2 = list(br = c(300, 0.05), bs = c(2, 0.05),   wr = c(200, 0.05), ws = c(3, 0.05)),
               v1 = list(br = c(300, 0.05), bs = c(0.3, 0.05), wr = c(200, 0.05), ws = c(0.5, 0.05)))

fit_model <- function(Y, model) {
  std <- model == "v2"
  s <- if (std) apply(Y, 1, function(y) sd(resid(lm(y ~ tc)))) else rep(1, nc)
  Z <- (Y - rowMeans(Y)) / s
  pr <- priors[[model]]
  spde_b <- inla.spde2.pcmatern(mesh, prior.range = pr$br, prior.sigma = pr$bs)
  spde_w <- inla.spde2.pcmatern(mesh, prior.range = pr$wr, prior.sigma = pr$ws)
  idx_b <- inla.spde.make.index("beta", spde_b$n.spde)
  idx_w <- inla.spde.make.index("w", spde_w$n.spde, n.repl = T)
  stk <- inla.stack(data = list(y = Z[cbind(dd$cell, dd$yr)]), A = list(A_b, A_w, 1),
                    effects = list(idx_b, idx_w, data.frame(b0 = tc[dd$yr])))
  form <- y ~ -1 + b0 + f(beta, model = spde_b) + f(w, model = spde_w, replicate = w.repl)
  fit <- inla(form, data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.compute = list(config = TRUE),
              control.inla = list(int.strategy = "eb"))
  list(fit = fit, s = s, n_spde = spde_b$n.spde)
}

# Draw n joint posterior samples of the 274 cell trends (kWh m-2 day-1 decade-1)
draw_trends <- function(m, n, seed) {
  smp <- inla.posterior.sample(n, m$fit, seed = seed,
                               selection = list(beta = 1:m$n_spde, b0 = 1))
  rn <- rownames(smp[[1]]$latent); ib <- grep("^beta:", rn); i0 <- grep("^b0:", rn)
  stopifnot(length(ib) == m$n_spde, length(i0) == 1)
  Bm <- sapply(smp, function(z) z$latent[ib, 1]); b0 <- sapply(smp, function(z) z$latent[i0, 1])
  tz <- as.matrix(A_pred %*% Bm) + matrix(b0, nc, n, byrow = TRUE)
  tz * m$s
}

set_str <- function(E) paste(which(E == 1), collapse = ";")
claim_ok <- function(x, type) if (type == "<") colSums(x >= 0) == 0 else colSums(x <= 0) == 0
joint_in <- function(tr, E, type) {           # share of samples with the whole set on the claimed side
  idx <- which(E == 1); if (!length(idx)) return(NA_real_)
  mean(claim_ok(tr[idx, , drop = FALSE], type))
}
# Gaussian approximation of the joint probability of a set. mu and Sig are estimated
# from an INDEPENDENT batch (C) of posterior samples, so this is an approximation,
# not an exact reference. Returns probability and gaussint's error estimate.
joint_gauss <- function(mu, Sig, E, type) {
  idx <- which(E == 1); if (!length(idx)) return(c(NA_real_, NA_real_))
  m <- mu[idx]; S <- Sig[idx, idx, drop = FALSE]
  if (length(idx) == 1) {
    sd1 <- sqrt(S[1, 1]); p <- if (type == "<") pnorm(0, m, sd1) else 1 - pnorm(0, m, sd1)
    return(c(as.numeric(p), 0))
  }
  Q <- Matrix(solve(S), sparse = TRUE)
  a <- if (type == "<") rep(-Inf, length(idx)) else rep(0, length(idx))
  b <- if (type == "<") rep(0, length(idx)) else rep(Inf, length(idx))
  g <- tryCatch(gaussint(mu = m, Q = Q, a = a, b = b, n.iter = 20000), error = function(e) NULL)
  if (is.null(g)) c(NA_real_, NA_real_) else c(as.numeric(g$P), as.numeric(g$E))
}
set_row <- function(r, model, bn, size, a, type, E, p_in, oos_tr, gs) {
  p_oos <- joint_in(oos_tr, E, type); n_oos <- ncol(oos_tr)
  tibble(rep = r, model = model, batch = bn, size = size, alpha = a,
         sign = ifelse(type == "<", "neg", "pos"), n_cells = sum(E), any = sum(E) > 0,
         p_in = p_in, p_oos_mc = p_oos, n_oos = n_oos,
         se_oos_mc = if (is.na(p_oos)) NA_real_ else sqrt(max(p_oos * (1 - p_oos), 1 / n_oos) / n_oos),
         p_gauss = gs[1], err_gauss = gs[2], cells = set_str(E))
}

# One .rds file per finished replicate (written atomically) -> restart-safe, no duplicates
rep_file <- function(r, model) file.path(out, sprintf("rep_%03d_%s_m670.rds", r, model))  # smoke-test files (671-node mesh) are ignored
t_all <- Sys.time(); k_done <- 0
cat("Started:", format(t_all), if (SMOKE) "(SMOKE TEST: 2 replicates)" else "", "\n")
reps <- if (SMOKE) 1:2 else 1:N_REP

for (r in reps) for (model in MODELS) {
  if (file.exists(rep_file(r, model))) next
  t0 <- Sys.time()
  set.seed(5000 + r); Ys <- R[, sample(ncol(R))]          # same permutation as task B
  m <- tryCatch(fit_model(Ys, model), error = function(e) { cat("fit error rep", r, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(m)) next
  trA <- draw_trends(m, N_BATCH, seed = 110000 + 10 * r + 1)
  trB <- draw_trends(m, N_BATCH, seed = 110000 + 10 * r + 2)
  trC <- draw_trends(m, N_BATCH, seed = 110000 + 10 * r + 3)  # only for the Gaussian approximation
  muC <- rowMeans(trC); SigC <- cov(t(trC))
  batches <- list(A = trA, B = trB); other <- list(A = trB, B = trA)

  rows <- list(); Frows <- list()
  for (bn in names(batches)) for (n in SIZES) {
    tr <- batches[[bn]][, 1:n, drop = FALSE]
    for (type in c("<", ">")) {
      ex_lim <- excursions.mc(tr, alpha = F_LIM, u = 0, type = type)
      Frows[[length(Frows) + 1]] <- tibble(rep = r, model = model, batch = bn, size = as.character(n),
                                           sign = ifelse(type == "<", "neg", "pos"), cell = 1:nc,
                                           F = as.numeric(ex_lim$F))
      for (a in ALPHAS) {
        E <- excursions.mc(tr, alpha = a, u = 0, type = type)$E; E[is.na(E)] <- 0
        rows[[length(rows) + 1]] <- set_row(r, model, bn, as.character(n), a, type, E,
                                            p_in = joint_in(tr, E, type), oos_tr = other[[bn]],
                                            gs = joint_gauss(muC, SigC, E, type))
      }
    }
  }
  if (GAUSS) {   # sets from the Gaussian approximation (batch C); evaluated on A and B (independent)
    Qc <- Matrix(solve(SigC + diag(1e-10 * mean(diag(SigC)), nc)), sparse = TRUE)
    for (type in c("<", ">")) for (a in ALPHAS) {
      E <- tryCatch(excursions(alpha = a, u = 0, mu = muC, Q = Qc, type = type, method = "EB")$E,
                    error = function(e) rep(NA, nc))
      ok <- !all(is.na(E)); E[is.na(E)] <- 0
      rows[[length(rows) + 1]] <- set_row(r, model, "C", if (ok) "gauss_approx" else "gauss_failed",
                                          a, type, E, p_in = NA_real_, oos_tr = cbind(trA, trB),
                                          gs = joint_gauss(muC, SigC, E, type))
    }
  }
  mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  res <- list(sets = bind_rows(rows), F = bind_rows(Frows),
              timing = tibble(rep = r, model = model, minutes = mins,
                              range_w = m$fit$summary.hyperpar["Range for w", "mean"],
                              range_b = m$fit$summary.hyperpar["Range for beta", "mean"]))
  tmp <- paste0(rep_file(r, model), ".tmp"); saveRDS(res, tmp); file.rename(tmp, rep_file(r, model))
  k_done <- k_done + 1
  per <- as.numeric(difftime(Sys.time(), t_all, units = "mins")) / k_done
  left <- sum(!file.exists(mapply(rep_file, rep(reps, each = length(MODELS)), MODELS)))
  cat(sprintf("rep %d (%s) done in %.1f min | %d left | ETA %.1f h\n", r, model, mins, left, per * left / 60))
}

# Combine all finished replicates into CSVs (re-created from the .rds files every time)
fs <- list.files(out, pattern = "^rep_\\d{3}_.*_m670\\.rds$", full.names = TRUE)
if (length(fs)) {
  all <- lapply(fs, readRDS)
  write_csv(bind_rows(lapply(all, `[[`, "sets")),   file.path(out, "11_sets.csv"))
  write_csv(bind_rows(lapply(all, `[[`, "F")),      file.path(out, "11_F.csv"))
  write_csv(bind_rows(lapply(all, `[[`, "timing")), file.path(out, "11_timing.csv"))
}
writeLines(capture.output(sessionInfo()), file.path(out, "11_sessionInfo.txt"))
cat("\nFinished:", format(Sys.time()), "| replicates saved:", length(fs), "\n")
if (SMOKE) cat("Smoke test OK -> set SMOKE <- FALSE at the top and source again for the full run.\n") else
  cat("Run finished.\n")
