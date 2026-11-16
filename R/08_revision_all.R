# =============================================================
# Step 8: ALL revision analyses in ONE checkpointed run
# (answers the referee report of 17 Oct; see PROJECT_LOG.md sec. 6)
#
# Rule: run once, save RAW outputs. all metrics are computed
# later from results/rev/*.csv - no need to rerun R for new tables.
#
# Tasks (set TASKS below to run a subset; order = priority):
#   A  Real data (~2 h): day-weighted seasons; v2 base (all seasons);
#      v2 + cell intercepts; v2 sum-to-zero trend field; v2 on
#      2001-2016; v2 with full CCD integration; v2 prior sensitivity;
#      mesh coarse/fine; classical MK on new data; ACF lags 1-3.
#   B  Null simulation, 200 PAIRED reps: MK, v1, v2 (coarse mesh),
#      joint sets stored at alpha 0.025 AND 0.05 per sign (~10 h)
#   C  Signal simulation: 5 synthetic truths NOT taken from any
#      model, 40 paired reps each: MK and v2 (~5 h)
#   D  Null simulation with full CCD integration, 50 reps, v2 (~2.5 h)
#   E  Propagate uncertainty of the per-cell SDs s_i (annual, 20 refits)
#
# Checkpointing: every finished unit is appended to a CSV in
# results/rev/. If R stops, just source the script again - finished
# units are skipped. Partial results are usable at any time.
# Seeds: data perms 5000+r (null), 7000+r (signal); R seed 2026.
# =============================================================

TASKS        <- c("A", "B", "C", "D", "E")
N_NULL       <- 200
N_SIG        <- 40
N_CCD        <- 50
N_SBOOT      <- 20
SAMP_REAL    <- 1000
SAMP_SIM     <- 400

suppressPackageStartupMessages({
  library(INLA); library(fmesher); library(excursions)
  library(dplyr); library(tidyr); library(readr); library(trend)
})
set.seed(2026)
out <- file.path("results", "rev"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)",
                   "Monsoon (JJAS)", "Post-monsoon (OND)")
append_csv <- function(df, f) write_csv(df, f, append = file.exists(f))
done_keys <- function(f, cols) {
  if (!file.exists(f)) return(character(0))
  x <- read_csv(f, col_select = all_of(cols), show_col_types = FALSE) |> distinct()
  do.call(paste, c(x, sep = "|"))
}
cat("Started:", format(Sys.time()), "\n")

# ---------------- A0: day-weighted seasonal means ----------------
d <- read_csv("data/power_ghi_india_monthly.csv", show_col_types = FALSE)
ndays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
d <- d |> mutate(wd = ndays[MONTH],
                 SEASON = case_when(MONTH %in% 1:2 ~ "Winter (JF)",
                                    MONTH %in% 3:5 ~ "Pre-monsoon (MAM)",
                                    MONTH %in% 6:9 ~ "Monsoon (JJAS)",
                                    TRUE ~ "Post-monsoon (OND)"))
agg <- bind_rows(
  d |> group_by(LON, LAT, YEAR) |> summarise(GHI = weighted.mean(GHI, wd), .groups = "drop") |>
    mutate(SEASON = "Annual"),
  d |> group_by(LON, LAT, YEAR, SEASON) |> summarise(GHI = weighted.mean(GHI, wd), .groups = "drop"))
write_csv(agg, file.path(out, "A0_seasonal_means_dayweighted.csv"))

cells <- agg |> distinct(LON, LAT) |> arrange(LON, LAT) |>
  mutate(cell = row_number(),
         X_km = (LON - 82) * 111.32 * cos(22 * pi / 180),
         Y_km = (LAT - 22) * 110.57,
         w_area = cos(LAT * pi / 180) / sum(cos(LAT * pi / 180)))
nc <- nrow(cells); loc <- as.matrix(cells[, c("X_km", "Y_km")])
write_csv(cells, file.path(out, "cells.csv"))

season_matrix <- function(sname, years) {
  agg |> filter(SEASON == sname, YEAR %in% years) |>
    left_join(cells |> select(LON, LAT, cell), by = c("LON", "LAT")) |>
    arrange(cell, YEAR) |> select(cell, YEAR, GHI) |>
    pivot_wider(names_from = YEAR, values_from = GHI) |> arrange(cell) |>
    select(-cell) |> as.matrix()
}

# ---------------- meshes and model structures ----------------
meshes <- list(
  base   = fm_mesh_2d_inla(loc = loc, max.edge = c(120, 500), cutoff = 40, offset = c(150, 600)),
  coarse = fm_mesh_2d_inla(loc = loc, max.edge = c(200, 600), cutoff = 80, offset = c(200, 700)),
  fine   = fm_mesh_2d_inla(loc = loc, max.edge = c(80, 400),  cutoff = 25, offset = c(150, 600)))
cat("Mesh nodes:", sapply(meshes, function(m) m$n), "\n")

make_struct <- function(mesh, years) {
  T <- length(years); tc <- (years - mean(range(years))) / 10
  dd <- expand.grid(cell = 1:nc, yr = 1:T); obs_loc <- loc[dd$cell, ]
  list(mesh = mesh, T = T, tc = tc, dd = dd,
       A_b = inla.spde.make.A(mesh, loc = obs_loc, weights = tc[dd$yr]),
       A_w = inla.spde.make.A(mesh, loc = obs_loc, repl = dd$yr, n.repl = T),
       A_pred = inla.spde.make.A(mesh, loc = loc))
}
yrs_main <- 2001:2023; yrs_pre <- 2001:2016
ST <- list(base = make_struct(meshes$base, yrs_main),
           coarse = make_struct(meshes$coarse, yrs_main),
           fine = make_struct(meshes$fine, yrs_main),
           base_pre = make_struct(meshes$base, yrs_pre))

prior_sets <- list(   # standardised (std = TRUE) units; v1 uses raw units
  base   = list(br = c(300, 0.05),  bs = c(2, 0.05),   wr = c(200, 0.05), ws = c(3, 0.05)),
  rough  = list(br = c(100, 0.05),  bs = c(6, 0.05),   wr = c(100, 0.05), ws = c(6, 0.05)),
  smooth = list(br = c(1000, 0.5),  bs = c(0.7, 0.05), wr = c(500, 0.5),  ws = c(2, 0.05)),
  v1     = list(br = c(300, 0.05),  bs = c(0.3, 0.05), wr = c(200, 0.05), ws = c(0.5, 0.05)))

# ---------------- core model fit ----------------
# Y: nc x T matrix of raw values. Returns trend samples in kWh m-2 day-1 decade-1.
fit_spde <- function(Y, S, std = TRUE, int = "eb", sumzero = FALSE, intercept = FALSE,
                     prior = if (std) "base" else "v1", n_samp = SAMP_REAL, s_override = NULL) {
  tc <- S$tc; dd <- S$dd
  s <- if (!is.null(s_override)) s_override else
       if (std) apply(Y, 1, function(y) sd(resid(lm(y ~ tc)))) else rep(1, nc)
  Yc <- if (intercept) Y else Y - rowMeans(Y)
  Z <- Yc / s
  pr <- prior_sets[[prior]]
  spde_b <- inla.spde2.pcmatern(S$mesh, prior.range = pr$br, prior.sigma = pr$bs)
  spde_w <- inla.spde2.pcmatern(S$mesh, prior.range = pr$wr, prior.sigma = pr$ws)
  idx_b <- inla.spde.make.index("beta", spde_b$n.spde)
  idx_w <- inla.spde.make.index("w", spde_w$n.spde, n.repl = S$T)
  A_con <- matrix(as.numeric(crossprod(cells$w_area, as.matrix(S$A_pred))), nrow = 1)
  stk <- inla.stack(data = list(y = Z[cbind(dd$cell, dd$yr)]), A = list(S$A_b, S$A_w, 1),
                    effects = list(idx_b, idx_w, data.frame(b0 = tc[dd$yr], cell = dd$cell)))
  f_str <- "y ~ -1 + b0 + f(beta, model = spde_b"
  if (sumzero) f_str <- paste0(f_str, ", extraconstr = list(A = A_con, e = 0)")
  f_str <- paste0(f_str, ") + f(w, model = spde_w, replicate = w.repl)")
  if (intercept) f_str <- paste0(f_str,
    " + f(cell, model = 'iid', hyper = list(prec = list(initial = -10, fixed = TRUE)))")
  form <- as.formula(f_str, env = environment())
  fit <- inla(form, data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.compute = list(config = TRUE),
              control.inla = list(int.strategy = int))
  smp <- inla.posterior.sample(n_samp, fit)
  rn <- rownames(smp[[1]]$latent); ib <- grep("^beta:", rn); i0 <- grep("^b0:", rn)
  tz <- as.matrix(S$A_pred %*% sapply(smp, function(z) z$latent[ib, 1])) +
        matrix(sapply(smp, function(z) z$latent[i0, 1]), nc, n_samp, byrow = TRUE)
  list(tr = tz * s, s = s, fit = fit)
}

cell_summary <- function(tr) {
  q <- t(apply(tr, 1, quantile, c(0.025, 0.05, 0.95, 0.975)))
  tibble(cell = 1:nc, est = rowMeans(tr), sd = apply(tr, 1, sd),
         q025 = q[, 1], q05 = q[, 2], q95 = q[, 3], q975 = q[, 4],
         jneg025 = excursions.mc(tr, 0.025, 0, "<")$E == 1,
         jpos025 = excursions.mc(tr, 0.025, 0, ">")$E == 1,
         jneg05  = excursions.mc(tr, 0.05,  0, "<")$E == 1,
         jpos05  = excursions.mc(tr, 0.05,  0, ">")$E == 1)
}
hyp <- function(fit, name) {
  h <- fit$summary.hyperpar
  if (!name %in% rownames(h)) return(c(NA, NA, NA))
  unlist(h[name, c("mean", "0.025quant", "0.975quant")])
}
mk_cells <- function(Y, years) {
  tibble(cell = 1:nc,
         est = apply(Y, 1, function(y) unname(sens.slope(y)$estimates)) * 10,
         p   = apply(Y, 1, function(y) mk.test(y)$p.value))
}
timer <- function(t0) round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

# ================= TASK A: real data =================
if ("A" %in% TASKS) {
  fA_cells <- file.path(out, "A_cells.csv"); fA_sum <- file.path(out, "A_summary.csv")
  doneA <- done_keys(fA_sum, c("variant", "SEASON"))
  run_variant <- function(variant, sname, years, S, ...) {
    if (paste(variant, sname, sep = "|") %in% doneA) return(invisible())
    t0 <- Sys.time()
    r <- fit_spde(season_matrix(sname, years), S, ...)
    cs <- cell_summary(r$tr) |> mutate(variant = variant, SEASON = sname, s = r$s, .before = 1)
    nat <- as.vector(crossprod(cells$w_area, r$tr))
    rb <- hyp(r$fit, "Range for beta"); rw <- hyp(r$fit, "Range for w")
    sb <- hyp(r$fit, "Stdev for beta"); sw <- hyp(r$fit, "Stdev for w")
    sm <- tibble(variant = variant, SEASON = sname,
                 india = mean(nat), india_lo = quantile(nat, 0.025), india_hi = quantile(nat, 0.975),
                 pct_jneg025 = 100 * mean(cs$jneg025), pct_jpos025 = 100 * mean(cs$jpos025),
                 pct_jneg05 = 100 * mean(cs$jneg05), pct_jpos05 = 100 * mean(cs$jpos05),
                 range_b = rb[1], range_b_lo = rb[2], range_b_hi = rb[3],
                 range_w = rw[1], range_w_lo = rw[2], range_w_hi = rw[3],
                 sd_b = sb[1], sd_w = sw[1], mesh_nodes = S$mesh$n, minutes = timer(t0))
    append_csv(cs, fA_cells); append_csv(sm, fA_sum)
    cat("A |", variant, "|", sname, "|", sm$minutes, "min\n")
  }
  # A1 main model, day-weighted, all seasons
  for (sn in season_levels) run_variant("v2_base", sn, yrs_main, ST$base)
  # A2 cell intercepts instead of pre-centering; A3 sum-to-zero trend field
  run_variant("v2_intercept", "Annual", yrs_main, ST$base, intercept = TRUE)
  run_variant("v2_sumzero",   "Annual", yrs_main, ST$base, sumzero = TRUE)
  # A4 pre-2017 refit, all seasons
  for (sn in season_levels) run_variant("v2_2001_2016", sn, yrs_pre, ST$base_pre)
  # A5 full integration over hyperparameters
  for (sn in c("Annual", "Winter (JF)")) run_variant("v2_ccd", sn, yrs_main, ST$base, int = "ccd")
  # A6 prior sensitivity of v2 (most prior-sensitive seasons + annual)
  for (pr in c("rough", "smooth")) for (sn in c("Annual", "Winter (JF)", "Pre-monsoon (MAM)"))
    run_variant(paste0("v2_prior_", pr), sn, yrs_main, ST$base, prior = pr)
  # A7 mesh sensitivity
  run_variant("v2_mesh_coarse", "Annual", yrs_main, ST$coarse)
  run_variant("v2_mesh_fine",   "Annual", yrs_main, ST$fine)
  # A8 v1 on day-weighted data (annual), for comparison at both alphas
  run_variant("v1_base", "Annual", yrs_main, ST$base, std = FALSE)

  # A9 classical MK on day-weighted data + ACF lags 1-3
  fA_cls <- file.path(out, "A_classical.csv")
  if (!file.exists(fA_cls)) {
    cls <- bind_rows(lapply(season_levels, function(sn) {
      Y <- season_matrix(sn, yrs_main); tc <- ST$base$tc
      ac <- t(apply(Y, 1, function(y) acf(resid(lm(y ~ tc)), lag.max = 3, plot = FALSE)$acf[2:4]))
      mk_cells(Y, yrs_main) |> mutate(SEASON = sn, acf1 = ac[, 1], acf2 = ac[, 2], acf3 = ac[, 3])
    }))
    write_csv(cls, fA_cls); cat("A | classical + ACF saved\n")
  }
}

# ---------------- simulation ingredients ----------------
Yann <- season_matrix("Annual", yrs_main)
tcm  <- ST$base$tc
R    <- t(apply(Yann, 1, function(y) resid(lm(y ~ tcm))))
s_ann <- apply(R, 1, sd)

SIM_COLS <- c("scenario", "rep", "method", "cell", "est", "sd", "q025", "q05", "q95", "q975",
              "jneg025", "jpos025", "jneg05", "jpos05", "p")
sim_rows <- function(scen, rep, method, tb) {
  tb <- tb |> mutate(scenario = scen, rep = rep, method = method)
  for (cc in setdiff(SIM_COLS, names(tb))) tb[[cc]] <- NA      # same columns in every row
  tb[, SIM_COLS]
}
run_sim_rep <- function(Ys, scen, rep, methods, f, S = ST$coarse, int = "eb") {
  if (length(methods) == 0) return(invisible())
  rows <- list()
  if ("MK" %in% methods) rows$MK <- sim_rows(scen, rep, "MK", mk_cells(Ys, yrs_main))
  if ("v1" %in% methods) rows$v1 <- sim_rows(scen, rep, "v1",
      cell_summary(fit_spde(Ys, S, std = FALSE, int = int, n_samp = SAMP_SIM)$tr))
  if ("v2" %in% methods) rows$v2 <- sim_rows(scen, rep, if (int == "ccd") "v2_ccd" else "v2",
      cell_summary(fit_spde(Ys, S, std = TRUE, int = int, n_samp = SAMP_SIM)$tr))
  append_csv(bind_rows(rows), f)
}
eta <- function(t_start, k, n) {
  m <- as.numeric(difftime(Sys.time(), t_start, units = "mins")) / k
  sprintf("%.1f min/rep, ETA %.1f h", m, m * (n - k) / 60)
}

# ================= TASK B: null, paired =================
if ("B" %in% TASKS) {
  fB <- file.path(out, "B_null_cells.csv")
  doneB <- done_keys(fB, c("rep", "method"))
  todo <- setdiff(1:N_NULL, as.integer(sub("\\|v2$", "", grep("\\|v2$", doneB, value = TRUE))))
  t_start <- Sys.time(); k <- 0
  for (r in todo) {
    set.seed(5000 + r); Ys <- R[, sample(ncol(R))]
    m <- c("MK", "v1", "v2")[!paste(r, c("MK", "v1", "v2"), sep = "|") %in% doneB]
    tryCatch(run_sim_rep(Ys, "null", r, m, fB), error = function(e) cat("B rep", r, "error:", conditionMessage(e), "\n"))
    k <- k + 1; cat("B | null rep", r, "|", eta(t_start, k, length(todo)), "\n")
  }
}

# ================= TASK C: signal, synthetic truths =================
truths <- cells |> transmute(
  cell,
  S1_smooth   = -0.15 + 0.04 * (Y_km - mean(Y_km)) / sd(Y_km),
  S2_mixed    = pmax(pmin(0.12 * (X_km - median(X_km)) / sd(X_km), 0.2), -0.2),
  S3_sharp    = ifelse(LAT < 22, -0.18, 0),
  S4_sparse   = { b <- -0.25 * (exp(-((X_km - (77 - 82) * 103.2)^2 + (Y_km - (18 - 22) * 110.57)^2) / (2 * 250^2)) +
                                exp(-((X_km - (85 - 82) * 103.2)^2 + (Y_km - (25 - 22) * 110.57)^2) / (2 * 250^2)));
                  ifelse(abs(b) < 0.01, 0, b) },
  S5_hetero   = -1.2 * s_ann)
if ("C" %in% TASKS) {
  write_csv(truths, file.path(out, "C_truths.csv"))
  fC <- file.path(out, "C_signal_cells.csv")
  doneC <- done_keys(fC, c("scenario", "rep", "method"))
  scen_names <- setdiff(names(truths), "cell")
  units <- expand.grid(rep = 1:N_SIG, scen = scen_names, stringsAsFactors = FALSE)
  units <- units[!paste(units$scen, units$rep, "v2", sep = "|") %in% doneC, ]
  t_start <- Sys.time(); k <- 0
  for (u in seq_len(nrow(units))) {
    r <- units$rep[u]; sc <- units$scen[u]
    set.seed(7000 + r); Ys <- R[, sample(ncol(R))] + outer(truths[[sc]], tcm)
    m <- c("MK", "v2")[!paste(sc, r, c("MK", "v2"), sep = "|") %in% doneC]
    tryCatch(run_sim_rep(Ys, sc, r, m, fC), error = function(e) cat("C", sc, r, "error:", conditionMessage(e), "\n"))
    k <- k + 1; cat("C |", sc, "rep", r, "|", eta(t_start, k, nrow(units)), "\n")
  }
}

# ================= TASK D: null with CCD integration =================
if ("D" %in% TASKS) {
  fD <- file.path(out, "D_null_ccd_cells.csv")
  doneD <- done_keys(fD, c("rep"))
  todo <- setdiff(1:N_CCD, as.integer(doneD))
  t_start <- Sys.time(); k <- 0
  for (r in todo) {
    set.seed(5000 + r); Ys <- R[, sample(ncol(R))]       # same perms as task B
    tryCatch(run_sim_rep(Ys, "null", r, "v2", fD, int = "ccd"),
             error = function(e) cat("D rep", r, "error:", conditionMessage(e), "\n"))
    k <- k + 1; cat("D | ccd null rep", r, "|", eta(t_start, k, length(todo)), "\n")
  }
}

# ================= TASK E: uncertainty in s_i (real data, annual) =================
if ("E" %in% TASKS) {
  fE <- file.path(out, "E_sboot_cells.csv")
  if (!file.exists(fE)) {
    df <- ST$base$T - 2; pooled <- NULL; t0 <- Sys.time()
    for (b in 1:N_SBOOT) {
      set.seed(9000 + b)
      s_star <- s_ann * sqrt(df / rchisq(nc, df))         # draw plausible SDs
      tr <- fit_spde(Yann, ST$base, std = TRUE, s_override = s_star, n_samp = 250)$tr
      pooled <- cbind(pooled, tr)
      cat("E | refit", b, "of", N_SBOOT, "|", timer(t0), "min\n")
    }
    cs <- cell_summary(pooled) |> mutate(variant = "v2_sboot", SEASON = "Annual", .before = 1)
    nat <- as.vector(crossprod(cells$w_area, pooled))
    write_csv(cs, fE)
    write_csv(tibble(variant = "v2_sboot", india = mean(nat), india_lo = quantile(nat, 0.025),
                     india_hi = quantile(nat, 0.975), pct_jneg025 = 100 * mean(cs$jneg025)),
              file.path(out, "E_sboot_summary.csv"))
  }
}

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
cat("\nFinished:", format(Sys.time()), "\nAll raw outputs are in", out,
    "- run finished.\n")
