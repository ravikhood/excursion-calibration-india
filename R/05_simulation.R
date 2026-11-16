# =============================================================
# Step 5: Simulation study (the "statistics" core of the paper)
#
# Question: when interannual anomalies are spatially correlated (as in
# the real data, range ~1000 km), how often do cell-wise Mann-Kendall
# tests (raw / FDR) and the Bayesian SPDE model (marginal / joint)
# flag trends that are not there, and how well do they recover real ones?
#
# Design (uses the REAL data's spatial dependence, no parametric guess):
#   1. Take the annual GHI series (274 cells x 23 years), remove each
#      cell's OLS trend -> residual maps R_t.
#   2. NULL scenario  : permute whole yearly maps R_t (keeps spatial
#      correlation, destroys any trend). True trend = 0 everywhere.
#   3. SIGNAL scenario: permuted maps + known trend field beta_true(s)
#      = posterior mean annual trend from Step 3.
#   Classical methods: 500 replicates per scenario (fast).
#   Bayesian model   : N_BAYES replicates per scenario (about 5 min each).
#
# Checkpointing: Bayesian results are saved after every replicate to
# results/sim/; if R stops, just run the script again and it resumes.
# Run overnight. Output: results/sim/sim_summary.csv
# =============================================================

N_CLASSIC <- 500
N_BAYES   <- 25          # per scenario; 25 x 2 x ~5 min = ~4 hours

library(INLA); library(fmesher); library(excursions)
library(dplyr); library(tidyr); library(readr); library(trend)

set.seed(42)
tag <- "2001_2023"; sim_dir <- file.path("results", "sim")
dir.create(sim_dir, recursive = TRUE, showWarnings = FALSE)

agg <- read_csv(file.path("data", paste0("ghi_seasonal_means_", tag, ".csv")),
                show_col_types = FALSE) |> filter(SEASON == "Annual")
years <- sort(unique(agg$YEAR)); n_year <- length(years)
t_mid <- mean(range(years)); tc <- (years - t_mid) / 10

cells <- agg |> distinct(LON, LAT) |> arrange(LON, LAT) |>
  mutate(cell = row_number(),
         X_km = (LON - 82) * 111.32 * cos(22 * pi / 180),
         Y_km = (LAT - 22) * 110.57)
nc  <- nrow(cells)
loc <- as.matrix(cells[, c("X_km", "Y_km")])

Y <- agg |> left_join(cells, by = c("LON", "LAT")) |> arrange(cell, YEAR) |>
  select(cell, YEAR, GHI) |> pivot_wider(names_from = YEAR, values_from = GHI) |>
  arrange(cell) |> select(-cell) |> as.matrix()
R <- t(apply(Y, 1, function(y) resid(lm(y ~ tc))))          # detrended residuals

beta_true <- read_csv(file.path("results", tag, "bayes_cell_trends.csv"),
                      show_col_types = FALSE) |>
  filter(SEASON == "Annual") |> arrange(LON, LAT) |> pull(post_mean)

make_sim <- function(scen) {
  Ys <- R[, sample(n_year)]
  if (scen == "signal") Ys <- Ys + outer(beta_true, tc)
  Ys - rowMeans(Ys)
}

# ---------- classical ----------
classic <- function(Ys) {
  p   <- apply(Ys, 1, function(y) mk.test(y)$p.value)
  sen <- apply(Ys, 1, function(y) unname(sens.slope(y)$estimates)) * 10
  pf  <- p.adjust(p, "BH")
  list(est = sen,
       neg = cbind(raw = p < 0.05 & sen < 0, fdr = pf < 0.10 & sen < 0),
       pos = cbind(raw = p < 0.05 & sen > 0, fdr = pf < 0.10 & sen > 0))
}

# ---------- Bayesian (same model as Step 3, base priors) ----------
mesh <- fm_mesh_2d_inla(loc = loc, max.edge = c(120, 500), cutoff = 40, offset = c(150, 600))
spde_b <- inla.spde2.pcmatern(mesh, prior.range = c(300, 0.05), prior.sigma = c(0.3, 0.05))
spde_w <- inla.spde2.pcmatern(mesh, prior.range = c(200, 0.05), prior.sigma = c(0.5, 0.05))
A_pred <- inla.spde.make.A(mesh, loc = loc)
dd <- expand.grid(cell = 1:nc, yr = 1:n_year)
obs_loc <- loc[dd$cell, ]
A_b <- inla.spde.make.A(mesh, loc = obs_loc, weights = tc[dd$yr])
A_w <- inla.spde.make.A(mesh, loc = obs_loc, repl = dd$yr, n.repl = n_year)
idx_b <- inla.spde.make.index("beta", spde_b$n.spde)
idx_w <- inla.spde.make.index("w", spde_w$n.spde, n.repl = n_year)

bayes <- function(Ys, n_samp = 300) {
  stk <- inla.stack(data = list(y = Ys[cbind(dd$cell, dd$yr)]), A = list(A_b, A_w, 1),
                    effects = list(idx_b, idx_w, data.frame(b0 = tc[dd$yr])))
  fit <- inla(y ~ -1 + b0 + f(beta, model = spde_b) + f(w, model = spde_w, replicate = w.repl),
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.compute = list(config = TRUE),
              control.inla = list(int.strategy = "eb"))
  smp <- inla.posterior.sample(n_samp, fit)
  rn <- rownames(smp[[1]]$latent); ib <- grep("^beta:", rn); i0 <- grep("^b0:", rn)
  tr <- as.matrix(A_pred %*% sapply(smp, function(z) z$latent[ib, 1])) +
        matrix(sapply(smp, function(z) z$latent[i0, 1]), nc, n_samp, byrow = TRUE)
  q025 <- apply(tr, 1, quantile, 0.025); q975 <- apply(tr, 1, quantile, 0.975)
  list(est = rowMeans(tr),
       neg = cbind(marg = q975 < 0, joint = excursions.mc(tr, 0.05, 0, "<")$E == 1),
       pos = cbind(marg = q025 > 0, joint = excursions.mc(tr, 0.05, 0, ">")$E == 1))
}

# ---------- metrics for one replicate ----------
metrics <- function(res, scen, rep) {
  truth_neg <- if (scen == "null") rep(FALSE, nc) else beta_true < 0
  truth_pos <- if (scen == "null") rep(FALSE, nc) else beta_true > 0
  bind_rows(lapply(colnames(res$neg), function(m) {
    flag_neg <- res$neg[, m]; flag_pos <- res$pos[, m]
    false_flag <- (flag_neg & !truth_neg) | (flag_pos & !truth_pos)
    tibble(scenario = scen, rep = rep, method = m,
           pct_flagged    = 100 * mean(flag_neg | flag_pos),
           any_false      = any(false_flag),
           pct_false      = 100 * mean(false_flag),
           power_neg      = if (scen == "null") NA else 100 * mean(flag_neg[truth_neg]),
           rmse           = sqrt(mean((res$est - if (scen == "null") 0 else beta_true)^2)))
  }))
}

# ---------- run classical ----------
cat("Classical replicates...\n")
cls <- bind_rows(lapply(c("null", "signal"), function(s)
  bind_rows(lapply(1:N_CLASSIC, function(r) metrics(classic(make_sim(s)), s, r)))))
cls$method <- paste0("MK_", cls$method)
write_csv(cls, file.path(sim_dir, "sim_classical_reps.csv"))

# ---------- run Bayesian with checkpoints ----------
ck <- file.path(sim_dir, "sim_bayes_reps.rds")
bay <- if (file.exists(ck)) readRDS(ck) else tibble()
for (s in c("null", "signal")) for (r in 1:N_BAYES) {
  if (nrow(bay) && any(bay$scenario == s & bay$rep == r)) next
  set.seed(1000 * (s == "signal") + r)
  t0 <- Sys.time()
  m <- tryCatch(metrics(bayes(make_sim(s)), s, r), error = function(e) NULL)
  if (!is.null(m)) { m$method <- paste0("Bayes_", m$method); bay <- bind_rows(bay, m); saveRDS(bay, ck) }
  cat(s, "rep", r, "of", N_BAYES, "-",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}

# ---------- summary table (goes into the paper) ----------
summ <- bind_rows(cls, bay) |>
  group_by(scenario, method) |>
  summarise(reps = n(),
            mean_pct_flagged = round(mean(pct_flagged), 1),
            FWER_any_false   = round(100 * mean(any_false), 1),
            mean_pct_false   = round(mean(pct_false), 2),
            power_neg        = round(mean(power_neg), 1),
            RMSE             = round(mean(rmse), 4), .groups = "drop") |>
  arrange(scenario, method)
write_csv(summ, file.path(sim_dir, "sim_summary.csv"))
options(width = 200); print(as.data.frame(summ))
