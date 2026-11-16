# =============================================================
# Step 3: Bayesian spatially varying trend model (INLA-SPDE)
#
# For each season, with y_it = GHI at cell i, year t, minus cell mean:
#
#   y_it = (beta0 + beta(s_i)) * t_c  +  w_t(s_i)  +  eps_it
#
#   t_c       : (year - 2012)/10, so trends are per decade
#   beta0     : India-wide mean trend (fixed effect)
#   beta(s)   : spatially varying trend deviation, Matern GRF (SPDE)
#   w_t(s)    : year-specific spatial anomaly field, Matern GRF,
#               independent replicate for each year (captures the
#               within-year spatial correlation that cell-wise tests ignore)
#   eps_it    : iid Gaussian nugget
#
# Priors: penalised-complexity (PC) priors on range and sd of both fields.
# Inference: INLA (empirical Bayes integration), posterior samples ->
#   per-cell trend summaries, marginal 95% credible intervals, and
#   JOINT excursion sets (excursions package, alpha = 0.05).
#
# Input : data/ghi_seasonal_means_2001_2023.csv   (from Step 2)
# Output: results/2001_2023/bayes_*.csv, figures/2001_2023/fig_bayes_*.png
# Runtime: roughly 2-6 minutes per season on a laptop.
# =============================================================

library(INLA); library(fmesher); library(excursions)
library(dplyr); library(tidyr); library(readr); library(ggplot2)
library(sf); library(rnaturalearth)

set.seed(2026)
tag <- "2001_2023"
res_dir <- file.path("results", tag); fig_dir <- file.path("figures", tag)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)",
                   "Monsoon (JJAS)", "Post-monsoon (OND)")
agg <- read_csv(file.path("data", paste0("ghi_seasonal_means_", tag, ".csv")),
                show_col_types = FALSE)
n_year <- length(unique(agg$YEAR))          # 23
yr0    <- min(agg$YEAR) - 1                 # 2000
t_mid  <- mean(range(agg$YEAR))             # 2012

# ---- 1. Grid cells in km (simple equirectangular projection) ----
cells <- agg |> distinct(LON, LAT) |> arrange(LON, LAT) |>
  mutate(cell = row_number(),
         X_km = (LON - 82) * 111.32 * cos(22 * pi / 180),
         Y_km = (LAT - 22) * 110.57)
loc <- as.matrix(cells[, c("X_km", "Y_km")])

# ---- 2. Mesh and SPDE models ----
mesh <- fm_mesh_2d_inla(loc = loc, max.edge = c(120, 500),
                        cutoff = 40, offset = c(150, 600))
cat("Mesh nodes:", mesh$n, "\n")
png(file.path(fig_dir, "fig_mesh.png"), width = 1600, height = 1600, res = 250)
plot(mesh, main = ""); points(loc, col = "red", pch = 16, cex = 0.4); dev.off()

spde_b <- inla.spde2.pcmatern(mesh, prior.range = c(300, 0.05),   # P(range < 300 km) = 0.05
                                    prior.sigma = c(0.3, 0.05))   # P(sd > 0.3 per decade) = 0.05
spde_w <- inla.spde2.pcmatern(mesh, prior.range = c(200, 0.05),
                                    prior.sigma = c(0.5, 0.05))
A_pred <- inla.spde.make.A(mesh, loc = loc)

# ---- 3. Fit one season ----
fit_season <- function(sname, n_samp = 500) {
  t_start <- Sys.time()
  dd <- agg |> filter(SEASON == sname) |>
    left_join(cells, by = c("LON", "LAT")) |>
    group_by(cell) |> mutate(y = GHI - mean(GHI)) |> ungroup() |>
    mutate(tc = (YEAR - t_mid) / 10, yr = YEAR - yr0)

  obs_loc <- as.matrix(dd[, c("X_km", "Y_km")])
  A_b <- inla.spde.make.A(mesh, loc = obs_loc, weights = dd$tc)
  A_w <- inla.spde.make.A(mesh, loc = obs_loc, repl = dd$yr, n.repl = n_year)
  idx_b <- inla.spde.make.index("beta", n.spde = spde_b$n.spde)
  idx_w <- inla.spde.make.index("w",    n.spde = spde_w$n.spde, n.repl = n_year)

  stk <- inla.stack(data = list(y = dd$y),
                    A = list(A_b, A_w, 1),
                    effects = list(idx_b, idx_w, data.frame(b0 = dd$tc)),
                    tag = "est")
  form <- y ~ -1 + b0 + f(beta, model = spde_b) +
                        f(w, model = spde_w, replicate = w.repl)

  fit <- inla(form, data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.compute   = list(config = TRUE, dic = TRUE, waic = TRUE),
              control.inla      = list(int.strategy = "eb"))

  # posterior samples of the total trend beta0 + beta(s) at every cell
  smp <- inla.posterior.sample(n_samp, fit)
  rn  <- rownames(smp[[1]]$latent)
  ib  <- grep("^beta:", rn); i0 <- grep("^b0:", rn)
  Bm  <- sapply(smp, function(z) z$latent[ib, 1])
  b0v <- sapply(smp, function(z) z$latent[i0, 1])
  trend <- as.matrix(A_pred %*% Bm) +
           matrix(b0v, nrow = nrow(loc), ncol = n_samp, byrow = TRUE)

  ex_neg <- excursions.mc(trend, alpha = 0.05, u = 0, type = "<")
  ex_pos <- excursions.mc(trend, alpha = 0.05, u = 0, type = ">")

  cellres <- cells |>
    mutate(SEASON    = sname,
           post_mean = rowMeans(trend),
           post_sd   = apply(trend, 1, sd),
           q025      = apply(trend, 1, quantile, 0.025),
           q975      = apply(trend, 1, quantile, 0.975),
           p_neg     = rowMeans(trend < 0),
           marg_neg  = q975 < 0,
           marg_pos  = q025 > 0,
           joint_neg = ex_neg$E == 1,
           joint_pos = ex_pos$E == 1)

  hp <- fit$summary.hyperpar
  fx <- fit$summary.fixed["b0", ]
  summ <- tibble(
    SEASON        = sname,
    beta0_mean    = fx$mean, beta0_q025 = fx$`0.025quant`, beta0_q975 = fx$`0.975quant`,
    range_beta_km = hp["Range for beta", "mean"],
    sd_beta       = hp["Stdev for beta", "mean"],
    range_w_km    = hp["Range for w", "mean"],
    sd_w          = hp["Stdev for w", "mean"],
    sd_nugget     = 1 / sqrt(hp["Precision for the Gaussian observations", "mean"]),
    pct_marg_neg  = 100 * mean(cellres$marg_neg),
    pct_joint_neg = 100 * mean(cellres$joint_neg),
    pct_marg_pos  = 100 * mean(cellres$marg_pos),
    pct_joint_pos = 100 * mean(cellres$joint_pos),
    DIC = fit$dic$dic, WAIC = fit$waic$waic,
    minutes = round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1))
  cat(sname, "done in", summ$minutes, "min\n")
  list(cellres = cellres, summ = summ)
}

# ---- 4. Run all seasons ----
out <- lapply(season_levels, fit_season)
cellres <- bind_rows(lapply(out, `[[`, "cellres")) |>
  mutate(SEASON = factor(SEASON, levels = season_levels))
summ    <- bind_rows(lapply(out, `[[`, "summ"))
write_csv(cellres, file.path(res_dir, "bayes_cell_trends.csv"))
write_csv(summ,    file.path(res_dir, "bayes_summary.csv"))

# ---- 5. Compare with classical results ----
cls <- read_csv(file.path(res_dir, "classical_trends.csv"), show_col_types = FALSE) |>
  mutate(SEASON = factor(SEASON, levels = season_levels))
cmp <- cls |> group_by(SEASON) |>
  summarise(pct_raw05_neg = 100 * mean(sig_raw05 & sen < 0),
            pct_FDR10_neg = 100 * mean(sig_fdr10 & sen < 0)) |>
  left_join(summ |> mutate(SEASON = factor(SEASON, levels = season_levels)) |>
              select(SEASON, pct_marg_neg, pct_joint_neg), by = "SEASON")
write_csv(cmp, file.path(res_dir, "comparison_classical_vs_bayes.csv"))

options(width = 200)
cat("\n===== Bayesian model summary =====\n");      print(as.data.frame(summ), digits = 3)
cat("\n===== Classical vs Bayesian (% cells with significant NEGATIVE trend) =====\n")
print(as.data.frame(cmp), digits = 3)

# ---- 6. Figures ----
india <- ne_countries(country = "India", scale = "medium", returnclass = "sf")
lim <- max(abs(cellres$post_mean))
unit_lab <- expression(paste("kWh ", m^-2, " ", day^-1, " ", decade^-1))

p1 <- ggplot(cellres) +
  geom_tile(aes(LON, LAT, fill = post_mean)) +
  geom_point(data = filter(cellres, joint_neg), aes(LON, LAT), size = 0.5) +
  geom_point(data = filter(cellres, joint_pos), aes(LON, LAT), size = 0.7, shape = 3) +
  geom_sf(data = india, fill = NA, colour = "grey30", linewidth = 0.3) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-lim, lim), name = unit_lab) +
  facet_wrap(~ SEASON, ncol = 3) + coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL, title = "Posterior mean trend (INLA-SPDE), 2001-2023",
       caption = "Dots: inside joint 95% negative excursion set; crosses: joint positive") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig_bayes_trend.png"), p1, width = 9, height = 7.5, dpi = 300)

p2 <- ggplot(cellres) +
  geom_tile(aes(LON, LAT, fill = post_sd)) +
  geom_sf(data = india, fill = NA, colour = "grey30", linewidth = 0.3) +
  scale_fill_viridis_c(name = "Posterior sd") +
  facet_wrap(~ SEASON, ncol = 3) + coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) + theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig_bayes_sd.png"), p2, width = 9, height = 7.5, dpi = 300)

sc <- cellres |> select(SEASON, LON, LAT, post_mean) |>
  left_join(cls |> select(SEASON, LON, LAT, sen), by = c("SEASON", "LON", "LAT"))
p3 <- ggplot(sc, aes(sen, post_mean)) +
  geom_abline(linetype = 2) + geom_point(size = 0.6, alpha = 0.6) +
  facet_wrap(~ SEASON, ncol = 3) +
  labs(x = "Sen slope (cell-wise)", y = "Posterior mean trend (INLA-SPDE)") +
  theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "fig_compare_scatter.png"), p3, width = 8, height = 5.5, dpi = 300)

cat("\nDone. Send me the two printed tables and the PNGs fig_bayes_trend, fig_bayes_sd,",
    "fig_compare_scatter and fig_mesh from", fig_dir, "\n")
