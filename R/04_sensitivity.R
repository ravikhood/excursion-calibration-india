# =============================================================
# Step 4: (a) India-wide area-weighted mean trend with 95% CrI
#         (b) prior sensitivity of the Bayesian SPDE trend model
#
# Why (a): beta0 is NOT the India-mean trend, because the SPDE field
# beta(s) is not constrained to sum to zero. The correct national trend
# is the area-weighted (cos latitude) average of beta0 + beta(s) over
# all cells, computed from posterior samples.
#
# Prior sets (PC priors, P(range < r0) = a, P(sd > s0) = a):
#   base   : as in Step 3
#   rough  : allows short-range, larger-amplitude trend variation
#   smooth : favours long-range, small-amplitude trend variation
#
# 15 fits (3 priors x 5 seasons), about 1.5 hours in total.
# Output: results/2001_2023/sensitivity.csv
# =============================================================

library(INLA); library(fmesher); library(excursions)
library(dplyr); library(tidyr); library(readr)

set.seed(2026)
tag <- "2001_2023"; res_dir <- file.path("results", tag)
season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)",
                   "Monsoon (JJAS)", "Post-monsoon (OND)")
agg <- read_csv(file.path("data", paste0("ghi_seasonal_means_", tag, ".csv")),
                show_col_types = FALSE)
n_year <- length(unique(agg$YEAR)); yr0 <- min(agg$YEAR) - 1
t_mid  <- mean(range(agg$YEAR))

cells <- agg |> distinct(LON, LAT) |> arrange(LON, LAT) |>
  mutate(cell = row_number(),
         X_km = (LON - 82) * 111.32 * cos(22 * pi / 180),
         Y_km = (LAT - 22) * 110.57,
         w_area = cos(LAT * pi / 180) / sum(cos(LAT * pi / 180)))
loc  <- as.matrix(cells[, c("X_km", "Y_km")])
mesh <- fm_mesh_2d_inla(loc = loc, max.edge = c(120, 500), cutoff = 40, offset = c(150, 600))
A_pred <- inla.spde.make.A(mesh, loc = loc)

priors <- list(
  base   = list(b_r = c(300, 0.05),  b_s = c(0.3, 0.05),  w_r = c(200, 0.05), w_s = c(0.5, 0.05)),
  rough  = list(b_r = c(100, 0.05),  b_s = c(1.0, 0.05),  w_r = c(100, 0.05), w_s = c(1.0, 0.05)),
  smooth = list(b_r = c(1000, 0.5),  b_s = c(0.1, 0.05),  w_r = c(500, 0.5),  w_s = c(0.3, 0.05)))

fit_one <- function(sname, pr, pname, n_samp = 400) {
  t0 <- Sys.time()
  spde_b <- inla.spde2.pcmatern(mesh, prior.range = pr$b_r, prior.sigma = pr$b_s)
  spde_w <- inla.spde2.pcmatern(mesh, prior.range = pr$w_r, prior.sigma = pr$w_s)
  dd <- agg |> filter(SEASON == sname) |> left_join(cells, by = c("LON", "LAT")) |>
    group_by(cell) |> mutate(y = GHI - mean(GHI)) |> ungroup() |>
    mutate(tc = (YEAR - t_mid) / 10, yr = YEAR - yr0)
  obs_loc <- as.matrix(dd[, c("X_km", "Y_km")])
  A_b <- inla.spde.make.A(mesh, loc = obs_loc, weights = dd$tc)
  A_w <- inla.spde.make.A(mesh, loc = obs_loc, repl = dd$yr, n.repl = n_year)
  stk <- inla.stack(data = list(y = dd$y), A = list(A_b, A_w, 1),
                    effects = list(inla.spde.make.index("beta", spde_b$n.spde),
                                   inla.spde.make.index("w", spde_w$n.spde, n.repl = n_year),
                                   data.frame(b0 = dd$tc)))
  fit <- inla(y ~ -1 + b0 + f(beta, model = spde_b) + f(w, model = spde_w, replicate = w.repl),
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.compute = list(config = TRUE, waic = TRUE),
              control.inla = list(int.strategy = "eb"))
  smp <- inla.posterior.sample(n_samp, fit)
  rn <- rownames(smp[[1]]$latent)
  ib <- grep("^beta:", rn); i0 <- grep("^b0:", rn)
  trend <- as.matrix(A_pred %*% sapply(smp, function(z) z$latent[ib, 1])) +
           matrix(sapply(smp, function(z) z$latent[i0, 1]), nrow(loc), n_samp, byrow = TRUE)
  nat <- as.vector(crossprod(cells$w_area, trend))          # India-mean trend samples
  q025 <- apply(trend, 1, quantile, 0.025); q975 <- apply(trend, 1, quantile, 0.975)
  ex <- excursions.mc(trend, alpha = 0.05, u = 0, type = "<")
  hp <- fit$summary.hyperpar
  cat(pname, "|", sname, "done in",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  tibble(prior = pname, SEASON = sname,
         india_trend = mean(nat), india_q025 = quantile(nat, 0.025), india_q975 = quantile(nat, 0.975),
         pct_marg_neg = 100 * mean(q975 < 0), pct_joint_neg = 100 * mean(ex$E == 1),
         range_beta_km = hp["Range for beta", "mean"], sd_beta = hp["Stdev for beta", "mean"],
         range_w_km = hp["Range for w", "mean"], WAIC = fit$waic$waic)
}

grid <- expand.grid(SEASON = season_levels, prior = names(priors), stringsAsFactors = FALSE)
sens <- bind_rows(lapply(seq_len(nrow(grid)), function(k)
  fit_one(grid$SEASON[k], priors[[grid$prior[k]]], grid$prior[k])))
sens <- sens |> mutate(SEASON = factor(SEASON, levels = season_levels)) |> arrange(SEASON, prior)
write_csv(sens, file.path(res_dir, "sensitivity.csv"))
options(width = 200); print(as.data.frame(sens), digits = 3)
