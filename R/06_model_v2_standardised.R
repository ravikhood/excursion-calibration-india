# =============================================================
# Step 6: Model v2 - locally standardised SPDE trend model
#
# WHY: the simulation (Step 5) showed the v1 joint excursion sets were
# anti-conservative under the null (FWER 40%). Diagnosis: interannual
# variance differs 2-3x between cells (Himalayan cells ~1.5-1.9x the
# plains), but v1 assumes ONE stationary variance for the yearly field.
# High-variance cells then get credible intervals that are too narrow.
#
# FIX (as in Haug et al. 2020): standardise each cell's series by its
# own residual SD s_i, fit the same SPDE model to z_it = y_it / s_i,
# and back-transform trend samples: trend_i = s_i * trend_z_i.
# Also: joint excursion sets now use alpha = 0.025 per sign, so the
# two-sided family-wise level is 5%; 1000 posterior samples.
#
# Part A: real data, all seasons   (~30 min)
# Part B: simulation, null + signal, N_BAYES reps each (~3-4 h,
#         checkpointed, resumes if restarted)
# =============================================================

RUN_SIM <- TRUE
N_BAYES <- 25
ALPHA_SIDE <- 0.025
N_SAMP <- 1000

library(INLA); library(fmesher); library(excursions)
library(dplyr); library(tidyr); library(readr); library(ggplot2)
library(sf); library(rnaturalearth)

set.seed(2026)
tag <- "2001_2023"
res_dir <- file.path("results", tag); fig_dir <- file.path("figures", tag)
sim_dir <- file.path("results", "sim")
season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)",
                   "Monsoon (JJAS)", "Post-monsoon (OND)")
agg <- read_csv(file.path("data", paste0("ghi_seasonal_means_", tag, ".csv")),
                show_col_types = FALSE)
years <- sort(unique(agg$YEAR)); n_year <- length(years)
t_mid <- mean(range(years)); tc <- (years - t_mid) / 10

cells <- agg |> distinct(LON, LAT) |> arrange(LON, LAT) |>
  mutate(cell = row_number(),
         X_km = (LON - 82) * 111.32 * cos(22 * pi / 180),
         Y_km = (LAT - 22) * 110.57,
         w_area = cos(LAT * pi / 180) / sum(cos(LAT * pi / 180)))
nc <- nrow(cells); loc <- as.matrix(cells[, c("X_km", "Y_km")])

mesh <- fm_mesh_2d_inla(loc = loc, max.edge = c(120, 500), cutoff = 40, offset = c(150, 600))
spde_b <- inla.spde2.pcmatern(mesh, prior.range = c(300, 0.05), prior.sigma = c(2, 0.05))  # z-units (SD per decade)
spde_w <- inla.spde2.pcmatern(mesh, prior.range = c(200, 0.05), prior.sigma = c(3, 0.05))
A_pred <- inla.spde.make.A(mesh, loc = loc)
dd <- expand.grid(cell = 1:nc, yr = 1:n_year)
obs_loc <- loc[dd$cell, ]
A_b <- inla.spde.make.A(mesh, loc = obs_loc, weights = tc[dd$yr])
A_w <- inla.spde.make.A(mesh, loc = obs_loc, repl = dd$yr, n.repl = n_year)
idx_b <- inla.spde.make.index("beta", spde_b$n.spde)
idx_w <- inla.spde.make.index("w", spde_w$n.spde, n.repl = n_year)

# Ys: nc x n_year matrix (rows = cells in 'cells' order, cols = years)
fit_v2 <- function(Ys, n_samp = N_SAMP) {
  Ys <- Ys - rowMeans(Ys)
  s  <- apply(Ys, 1, function(y) sd(resid(lm(y ~ tc))))
  Z  <- Ys / s
  stk <- inla.stack(data = list(y = Z[cbind(dd$cell, dd$yr)]), A = list(A_b, A_w, 1),
                    effects = list(idx_b, idx_w, data.frame(b0 = tc[dd$yr])))
  fit <- inla(y ~ -1 + b0 + f(beta, model = spde_b) + f(w, model = spde_w, replicate = w.repl),
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.compute = list(config = TRUE, waic = TRUE),
              control.inla = list(int.strategy = "eb"))
  smp <- inla.posterior.sample(n_samp, fit)
  rn <- rownames(smp[[1]]$latent); ib <- grep("^beta:", rn); i0 <- grep("^b0:", rn)
  tz <- as.matrix(A_pred %*% sapply(smp, function(z) z$latent[ib, 1])) +
        matrix(sapply(smp, function(z) z$latent[i0, 1]), nc, n_samp, byrow = TRUE)
  tr <- tz * s                                   # back to kWh m-2 day-1 decade-1
  q_lo <- apply(tr, 1, quantile, ALPHA_SIDE); q_hi <- apply(tr, 1, quantile, 1 - ALPHA_SIDE)
  list(tr = tr, fit = fit, s = s,
       est = rowMeans(tr), sd = apply(tr, 1, sd), q_lo = q_lo, q_hi = q_hi,
       neg = cbind(marg = q_hi < 0, joint = excursions.mc(tr, ALPHA_SIDE, 0, "<")$E == 1),
       pos = cbind(marg = q_lo > 0, joint = excursions.mc(tr, ALPHA_SIDE, 0, ">")$E == 1))
}

season_matrix <- function(sname) {
  agg |> filter(SEASON == sname) |> left_join(cells, by = c("LON", "LAT")) |>
    arrange(cell, YEAR) |> select(cell, YEAR, GHI) |>
    pivot_wider(names_from = YEAR, values_from = GHI) |> arrange(cell) |>
    select(-cell) |> as.matrix()
}

# ================= Part A: real data =================
cellres <- list(); summ <- list()
for (sname in season_levels) {
  t0 <- Sys.time()
  r <- fit_v2(season_matrix(sname))
  nat <- as.vector(crossprod(cells$w_area, r$tr))
  hp <- r$fit$summary.hyperpar
  cellres[[sname]] <- cells |> mutate(SEASON = sname, post_mean = r$est, post_sd = r$sd,
    q_lo = r$q_lo, q_hi = r$q_hi, resid_sd = r$s,
    marg_neg = r$neg[, "marg"], joint_neg = r$neg[, "joint"],
    marg_pos = r$pos[, "marg"], joint_pos = r$pos[, "joint"])
  summ[[sname]] <- tibble(SEASON = sname,
    india_trend = mean(nat), india_lo = quantile(nat, 0.025), india_hi = quantile(nat, 0.975),
    pct_marg_neg = 100 * mean(r$neg[, "marg"]), pct_joint_neg = 100 * mean(r$neg[, "joint"]),
    pct_marg_pos = 100 * mean(r$pos[, "marg"]), pct_joint_pos = 100 * mean(r$pos[, "joint"]),
    range_beta_km = hp["Range for beta", "mean"], range_w_km = hp["Range for w", "mean"],
    WAIC = r$fit$waic$waic)
  cat(sname, "done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
}
cellres <- bind_rows(cellres) |> mutate(SEASON = factor(SEASON, levels = season_levels))
summ <- bind_rows(summ)
write_csv(cellres, file.path(res_dir, "v2_cell_trends.csv"))
write_csv(summ, file.path(res_dir, "v2_summary.csv"))
options(width = 200); cat("\n===== v2 real-data summary =====\n"); print(as.data.frame(summ), digits = 3)

india <- ne_countries(country = "India", scale = "medium", returnclass = "sf")
lim <- max(abs(cellres$post_mean))
p <- ggplot(cellres) +
  geom_tile(aes(LON, LAT, fill = post_mean)) +
  geom_point(data = filter(cellres, joint_neg), aes(LON, LAT), size = 0.5) +
  geom_point(data = filter(cellres, joint_pos), aes(LON, LAT), size = 0.7, shape = 3) +
  geom_sf(data = india, fill = NA, colour = "grey30", linewidth = 0.3) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-lim, lim),
                       name = expression(paste("kWh ", m^-2, " ", day^-1, " ", decade^-1))) +
  facet_wrap(~ SEASON, ncol = 3) + coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, title = "Posterior mean trend, standardised SPDE model (v2), 2001-2023",
       caption = "Dots/crosses: joint negative/positive excursion sets, family-wise level 5%") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig_v2_trend.png"), p, width = 9, height = 7.5, dpi = 300)

# ================= Part B: simulation =================
if (RUN_SIM) {
  Y <- season_matrix("Annual")
  R <- t(apply(Y, 1, function(y) resid(lm(y ~ tc))))
  beta_true <- read_csv(file.path(res_dir, "bayes_cell_trends.csv"), show_col_types = FALSE) |>
    filter(SEASON == "Annual") |> arrange(LON, LAT) |> pull(post_mean)

  metrics <- function(res, scen, rep) {
    truth_neg <- if (scen == "null") rep(FALSE, nc) else beta_true < 0
    truth_pos <- if (scen == "null") rep(FALSE, nc) else beta_true > 0
    bind_rows(lapply(colnames(res$neg), function(m) {
      fn <- res$neg[, m]; fp <- res$pos[, m]
      ff <- (fn & !truth_neg) | (fp & !truth_pos)
      tibble(scenario = scen, rep = rep, method = paste0("Bayes2_", m),
             pct_flagged = 100 * mean(fn | fp), any_false = any(ff), pct_false = 100 * mean(ff),
             power_neg = if (scen == "null") NA else 100 * mean(fn[truth_neg]),
             rmse = sqrt(mean((res$est - if (scen == "null") 0 else beta_true)^2)))
    }))
  }
  ck <- file.path(sim_dir, "sim_bayes2_reps.rds")
  bay2 <- if (file.exists(ck)) readRDS(ck) else tibble()
  for (s in c("null", "signal")) for (r in 1:N_BAYES) {
    if (nrow(bay2) && any(bay2$scenario == s & bay2$rep == r)) next
    set.seed(1000 * (s == "signal") + r)       # same permutations as Step 5
    Ys <- R[, sample(n_year)]
    if (s == "signal") Ys <- Ys + outer(beta_true, tc)
    t0 <- Sys.time()
    m <- tryCatch(metrics(fit_v2(Ys, n_samp = 500), s, r), error = function(e) NULL)
    if (!is.null(m)) { bay2 <- bind_rows(bay2, m); saveRDS(bay2, ck) }
    cat(s, "rep", r, "-", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  }
  old <- bind_rows(read_csv(file.path(sim_dir, "sim_classical_reps.csv"), show_col_types = FALSE),
                   readRDS(file.path(sim_dir, "sim_bayes_reps.rds")))
  summ_sim <- bind_rows(old, bay2) |> group_by(scenario, method) |>
    summarise(reps = n(), mean_pct_flagged = round(mean(pct_flagged), 1),
              FWER_any_false = round(100 * mean(any_false), 1),
              mean_pct_false = round(mean(pct_false), 2),
              power_neg = round(mean(power_neg), 1), RMSE = round(mean(rmse), 4),
              .groups = "drop") |> arrange(scenario, method)
  write_csv(summ_sim, file.path(sim_dir, "sim_summary_v2.csv"))
  cat("\n===== Simulation summary incl. v2 =====\n"); print(as.data.frame(summ_sim))
}
