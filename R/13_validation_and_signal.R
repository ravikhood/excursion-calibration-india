# =============================================================
# Step 13: (A) independent permutation replicates for validating a tuned
#              alpha; (B) more replicates for the mixed (S2) and sparse (S4)
#              synthetic trend surfaces, with paired MK vs SPDE comparisons.
#
# Run AFTER script 12 (needs no output of 12, but uses the laptop).
#
# Part A: null permutations r = 51..150 (same permutation rule as task B:
#         seed 5000 + r), v2 on the 670-node mesh, 4000 posterior samples,
#         excursion functions F stored per cell. Together with script 11
#         (reps 1-50, F stored) this gives 150 null replicates with F, so a
#         tuned alpha can be CHOSEN on one subset and EVALUATED on another.
# Part B: S2 and S4 (definitions copied exactly from script 08), replicates
#         r = 41..100 (seed 7000 + r, as task C), MK and v2 on the same
#         datasets (paired), 4000 samples, F stored. With task C reps 1-40
#         this gives 100 paired replicates for S2 and S4.
#
# Runtime: ~3-4 min per v2 fit. Part A ~6 h, Part B ~7 h. Stop with Esc any
# time; finished fits are skipped on restart.
# Output: results/val/A_rep_XXX.rds, results/val/B_<scen>_rep_XXX.rds,
#         combined into results/val/13_null_cells.csv, 13_signal_cells.csv
# =============================================================

PART_A_REPS <- 51:150
PART_B_REPS <- 41:100
PART_B_SCEN <- c("S2_mixed", "S4_sparse")
N_SAMP <- 4000
F_LIM  <- 0.10
SMOKE  <- TRUE          # TRUE = 1 replicate of A and 1 of each B scenario (~12 min); then FALSE

suppressPackageStartupMessages({
  library(INLA); library(fmesher); library(excursions)
  library(dplyr); library(tidyr); library(readr); library(trend)
})
inla.setOption(num.threads = "2:1")
out <- file.path("results", "val"); dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---------------- data and structures (identical to scripts 08/11/12) ----------------
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
if (mesh$n != 670) stop("Mesh differs from scripts 08/11/12. Stop and check the setup.")
T <- length(yrs); tc <- (yrs - mean(range(yrs))) / 10
dd <- expand.grid(cell = 1:nc, yr = 1:T); obs_loc <- loc[dd$cell, ]
A_b <- inla.spde.make.A(mesh, loc = obs_loc, weights = tc[dd$yr])
A_w <- inla.spde.make.A(mesh, loc = obs_loc, repl = dd$yr, n.repl = T)
A_pred <- inla.spde.make.A(mesh, loc = loc)
spde_b <- inla.spde2.pcmatern(mesh, prior.range = c(300, 0.05), prior.sigma = c(2, 0.05))
spde_w <- inla.spde2.pcmatern(mesh, prior.range = c(200, 0.05), prior.sigma = c(3, 0.05))
R <- t(apply(Yann, 1, function(y) resid(lm(y ~ tc))))      # detrended residual maps
s_ann <- apply(R, 1, sd)

# synthetic truths: copied exactly from script 08 (task C)
truths <- cells |> transmute(
  cell,
  S2_mixed  = pmax(pmin(0.12 * (X_km - median(X_km)) / sd(X_km), 0.2), -0.2),
  S4_sparse = { b <- -0.25 * (exp(-((X_km - (77 - 82) * 103.2)^2 + (Y_km - (18 - 22) * 110.57)^2) / (2 * 250^2)) +
                              exp(-((X_km - (85 - 82) * 103.2)^2 + (Y_km - (25 - 22) * 110.57)^2) / (2 * 250^2)));
                ifelse(abs(b) < 0.01, 0, b) })
chk <- read_csv("results/rev/C_truths.csv", show_col_types = FALSE)
if (max(abs(chk$S2_mixed - truths$S2_mixed)) > 1e-6 || max(abs(chk$S4_sparse - truths$S4_sparse)) > 1e-6)
  stop("Synthetic truths differ from task C. Stop and check the setup.")
cat("Synthetic truths identical to task C: OK\n")

fit_v2 <- function(Y) {
  s <- apply(Y, 1, function(y) sd(resid(lm(y ~ tc))))
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
v2_cells <- function(Y, seed) {
  tr <- draw_trends(fit_v2(Y), N_SAMP, seed)
  q <- t(apply(tr, 1, quantile, c(0.025, 0.975)))
  tibble(cell = 1:nc, est = rowMeans(tr), q025 = q[, 1], q975 = q[, 2],
         Fneg = as.numeric(excursions.mc(tr, F_LIM, 0, "<")$F),
         Fpos = as.numeric(excursions.mc(tr, F_LIM, 0, ">")$F),
         jneg025 = excursions.mc(tr, 0.025, 0, "<")$E == 1,
         jpos025 = excursions.mc(tr, 0.025, 0, ">")$E == 1)
}
mk_cells <- function(Y) tibble(cell = 1:nc,
  mk_est = apply(Y, 1, function(y) unname(sens.slope(y)$estimates)) * 10,
  mk_p   = apply(Y, 1, function(y) mk.test(y)$p.value))
save_rds <- function(x, f) { tmp <- paste0(f, ".tmp"); saveRDS(x, tmp); file.rename(tmp, f) }

a_reps <- if (SMOKE) PART_A_REPS[1] else PART_A_REPS
b_reps <- if (SMOKE) PART_B_REPS[1] else PART_B_REPS
units <- c(lapply(a_reps, function(r) list(part = "A", r = r, sc = "null")),
           unlist(lapply(b_reps, function(r) lapply(PART_B_SCEN, function(sc) list(part = "B", r = r, sc = sc))),
                  recursive = FALSE))
ufile <- function(u) if (u$part == "A") file.path(out, sprintf("A_rep_%03d.rds", u$r)) else
                     file.path(out, sprintf("B_%s_rep_%03d.rds", u$sc, u$r))
t_all <- Sys.time(); k_done <- 0
cat("Started:", format(t_all), if (SMOKE) "(SMOKE TEST)" else "", "| units:", length(units), "\n")

for (u in units) {
  f <- ufile(u); if (file.exists(f)) next
  t0 <- Sys.time()
  res <- tryCatch({
    if (u$part == "A") {
      set.seed(5000 + u$r); Ys <- R[, sample(ncol(R))]                    # same rule as task B
      v2_cells(Ys, seed = 500000 + u$r) |> mutate(rep = u$r, .before = 1)
    } else {
      set.seed(7000 + u$r); Ys <- R[, sample(ncol(R))] + outer(truths[[u$sc]], tc)   # same rule as task C
      v2_cells(Ys, seed = 600000 + 1000 * match(u$sc, PART_B_SCEN) + u$r) |>
        left_join(mk_cells(Ys), by = "cell") |>
        mutate(rep = u$r, scenario = u$sc, truth = truths[[u$sc]], .before = 1)
    }
  }, error = function(e) { cat("error", u$part, u$sc, u$r, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next
  save_rds(res, f)
  k_done <- k_done + 1
  left <- sum(!file.exists(sapply(units, ufile)))
  per <- as.numeric(difftime(Sys.time(), t_all, units = "mins")) / k_done
  cat(sprintf("%s %s rep %d done in %.1f min | %d left | ETA %.1f h\n", u$part, u$sc, u$r,
              as.numeric(difftime(Sys.time(), t0, units = "mins")), left, per * left / 60))
}

fa <- list.files(out, pattern = "^A_rep_\\d{3}\\.rds$", full.names = TRUE)
fb <- list.files(out, pattern = "^B_.*_rep_\\d{3}\\.rds$", full.names = TRUE)
if (length(fa)) write_csv(bind_rows(lapply(fa, readRDS)), file.path(out, "13_null_cells.csv"))
if (length(fb)) write_csv(bind_rows(lapply(fb, readRDS)), file.path(out, "13_signal_cells.csv"))
writeLines(capture.output(sessionInfo()), file.path(out, "13_sessionInfo.txt"))
cat("\nFinished:", format(Sys.time()), "| null reps:", length(fa), "| signal fits:", length(fb), "\n")
if (SMOKE) cat("Smoke test OK -> set SMOKE <- FALSE and source again.\n") else
  cat("Run finished.\n")
