# =============================================================
# Step 15 (Part A only): apply the validated tuned threshold alpha* = 0.005 per sign
# to the REAL NASA POWER data, store the excursion functions F, and draw the maps.
#
# What alpha* = 0.005 was calibrated and validated on (PROJECT_LOG 4k):
#   v2 model, COARSE mesh (670 nodes), ANNUAL means, 2001-2023, 4000 posterior samples,
#   dependence-preserving permutation null (chosen on reps 1-50, validated on 51-150: FWER 5%).
# So the run below has two tiers:
#   PRIMARY   : Annual, 2001-2023, coarse mesh -> exactly the validated setting.
#   SECONDARY : Annual on the base mesh (1383 nodes, the paper's main maps) and the four
#               seasons on the base mesh. alpha* was NOT calibrated for these; they are
#               reported as descriptive/sensitivity results only.
# Runtime: ~3 min (primary) + ~6 x 5-7 min (secondary) = about 40 min. Restart-safe.
# Output: results/tuned/15_cells.csv, 15_summary.csv, 15_sessionInfo.txt
#         figures/paper/fig9_tuned_annual.(pdf|png), figS4_tuned_seasons.(pdf|png)
# =============================================================

N_SAMP <- 4000
A_STAR <- 0.005          # validated tuned alpha per sign
A_STD  <- 0.025          # standard alpha per sign (for comparison)

suppressPackageStartupMessages({
  library(INLA); library(fmesher); library(excursions)
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(sf); library(rnaturalearth)
})
inla.setOption(num.threads = "2:1")
out <- file.path("results", "tuned"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
fig <- file.path("figures", "paper"); dir.create(fig, recursive = TRUE, showWarnings = FALSE)
season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)", "Monsoon (JJAS)", "Post-monsoon (OND)")

agg <- read_csv("results/rev/A0_seasonal_means_dayweighted.csv", show_col_types = FALSE)
cells <- read_csv("results/rev/cells.csv", show_col_types = FALSE) |>
  mutate(X_km = (LON - 82) * 111.32 * cos(22 * pi / 180), Y_km = (LAT - 22) * 110.57)  # exact, not rounded
nc <- nrow(cells); loc <- as.matrix(cells[, c("X_km", "Y_km")])
yrs <- 2001:2023; T <- length(yrs); tc <- (yrs - mean(range(yrs))) / 10
meshes <- list(coarse = fm_mesh_2d_inla(loc = loc, max.edge = c(200, 600), cutoff = 80, offset = c(200, 700)),
               base   = fm_mesh_2d_inla(loc = loc, max.edge = c(120, 500), cutoff = 40, offset = c(150, 600)))
cat("Mesh nodes: coarse", meshes$coarse$n, "(must be 670) | base", meshes$base$n, "(must be 1383)\n")
if (meshes$coarse$n != 670 || meshes$base$n != 1383) stop("Mesh differs from scripts 08/11-14. Stop and check the setup.")

season_matrix <- function(sname) {
  agg |> filter(SEASON == sname, YEAR %in% yrs) |>
    left_join(cells |> select(LON, LAT, cell), by = c("LON", "LAT")) |>
    arrange(cell, YEAR) |> select(cell, YEAR, GHI) |>
    pivot_wider(names_from = YEAR, values_from = GHI) |> arrange(cell) |> select(-cell) |> as.matrix()
}
fit_v2 <- function(Y, mesh, seed) {         # identical model to scripts 08/11/13 (v2, base priors, EB)
  dd <- expand.grid(cell = 1:nc, yr = 1:T); ol <- loc[dd$cell, ]
  A_b <- inla.spde.make.A(mesh, loc = ol, weights = tc[dd$yr])
  A_w <- inla.spde.make.A(mesh, loc = ol, repl = dd$yr, n.repl = T)
  A_pred <- inla.spde.make.A(mesh, loc = loc)
  spde_b <- inla.spde2.pcmatern(mesh, prior.range = c(300, 0.05), prior.sigma = c(2, 0.05))
  spde_w <- inla.spde2.pcmatern(mesh, prior.range = c(200, 0.05), prior.sigma = c(3, 0.05))
  s <- apply(Y, 1, function(y) sd(resid(lm(y ~ tc)))); Z <- (Y - rowMeans(Y)) / s
  stk <- inla.stack(data = list(y = Z[cbind(dd$cell, dd$yr)]), A = list(A_b, A_w, 1),
                    effects = list(inla.spde.make.index("beta", spde_b$n.spde),
                                   inla.spde.make.index("w", spde_w$n.spde, n.repl = T),
                                   data.frame(b0 = tc[dd$yr])))
  fit <- inla(y ~ -1 + b0 + f(beta, model = spde_b) + f(w, model = spde_w, replicate = w.repl),
              data = inla.stack.data(stk), control.predictor = list(A = inla.stack.A(stk)),
              control.compute = list(config = TRUE), control.inla = list(int.strategy = "eb"))
  smp <- inla.posterior.sample(N_SAMP, fit, seed = seed, selection = list(beta = 1:spde_b$n.spde, b0 = 1))
  rn <- rownames(smp[[1]]$latent); ib <- grep("^beta:", rn); i0 <- grep("^b0:", rn)
  stopifnot(length(ib) == spde_b$n.spde, length(i0) == 1)
  tr <- (as.matrix(A_pred %*% sapply(smp, function(z) z$latent[ib, 1])) +
           matrix(sapply(smp, function(z) z$latent[i0, 1]), nc, N_SAMP, byrow = TRUE)) * s
  list(tr = tr, s = s)
}

units <- c(list(list(tier = "primary", mesh = "coarse", SEASON = "Annual")),
           lapply(season_levels, function(sn) list(tier = "secondary", mesh = "base", SEASON = sn)))
fc <- file.path(out, "15_cells.csv"); fs <- file.path(out, "15_summary.csv")
done <- if (file.exists(fs)) with(read_csv(fs, show_col_types = FALSE), paste(mesh, SEASON)) else character(0)
cat("Started:", format(Sys.time()), "\n")
for (k in seq_along(units)) {
  u <- units[[k]]; if (paste(u$mesh, u$SEASON) %in% done) next
  t0 <- Sys.time()
  r <- fit_v2(season_matrix(u$SEASON), meshes[[u$mesh]], seed = 800000 + 100 * k)
  Fneg <- as.numeric(excursions.mc(r$tr, 0.10, 0, "<")$F); Fpos <- as.numeric(excursions.mc(r$tr, 0.10, 0, ">")$F)
  cs <- tibble(tier = u$tier, mesh = u$mesh, SEASON = u$SEASON, cell = 1:nc, s = r$s,
               est = rowMeans(r$tr), sd = apply(r$tr, 1, sd),
               q025 = apply(r$tr, 1, quantile, 0.025), q975 = apply(r$tr, 1, quantile, 0.975),
               Fneg = Fneg, Fpos = Fpos,
               neg_std  = Fneg > 1 - A_STD,  pos_std  = Fpos > 1 - A_STD,     # E = {F > 1 - alpha} (strict)
               neg_tuned = Fneg > 1 - A_STAR, pos_tuned = Fpos > 1 - A_STAR,
               neg_direct025 = excursions.mc(r$tr, A_STD, 0, "<")$E == 1)      # cross-check of the F convention
  if (any(cs$neg_std != cs$neg_direct025)) warning("F-derived and direct 0.025 sets differ in ", u$mesh, " ", u$SEASON)
  nat <- as.vector(crossprod(cells$w_area, r$tr))
  sm <- tibble(tier = u$tier, mesh = u$mesh, SEASON = u$SEASON,
               india = mean(nat), india_lo = quantile(nat, 0.025), india_hi = quantile(nat, 0.975),
               n_neg_std = sum(cs$neg_std), n_pos_std = sum(cs$pos_std),
               n_neg_tuned = sum(cs$neg_tuned), n_pos_tuned = sum(cs$pos_tuned),
               minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  write_csv(cs, fc, append = file.exists(fc)); write_csv(sm, fs, append = file.exists(fs))
  cat(sprintf("%-9s %-6s %-19s India %.3f | cells neg/pos at 0.025: %d/%d | at 0.005: %d/%d | %.1f min\n",
              u$tier, u$mesh, u$SEASON, sm$india, sm$n_neg_std, sm$n_pos_std, sm$n_neg_tuned, sm$n_pos_tuned, sm$minutes))
}
writeLines(capture.output(sessionInfo()), file.path(out, "15_sessionInfo.txt"))

# ---------------- maps ----------------
C <- read_csv(fc, show_col_types = FALSE) |> left_join(cells |> select(cell, LON, LAT), by = "cell") |>
  mutate(SEASON = factor(SEASON, levels = season_levels))
india <- ne_countries(country = "India", scale = "medium", returnclass = "sf")
lab_tr <- expression(paste("kWh ", m^-2, " ", day^-1, " ", decade^-1))
div <- scale_fill_gradient2(name = lab_tr, low = "#2166ac", mid = "white", high = "#b2182b",
                            midpoint = 0, limits = c(-0.3, 0.3), oob = scales::squish)
base_map <- function(d) ggplot(d) + geom_tile(aes(LON, LAT, fill = est)) +
  geom_sf(data = india, fill = NA, colour = "grey30", linewidth = 0.3) + div + coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) + theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.key.width = grid::unit(2.5, "cm"))
pts <- function(d, neg, pos) list(
  geom_point(data = filter(d, {{ neg }}), aes(LON, LAT), size = 0.35, colour = "black"),
  geom_point(data = filter(d, {{ pos }}), aes(LON, LAT), size = 0.9, shape = 4, colour = "black"))
save2 <- function(p, name, w = 9, h = 5) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
# Fig 9: the validated setting (annual, coarse mesh), standard vs tuned threshold
P <- C |> filter(tier == "primary")
P2 <- bind_rows(P |> mutate(panel = "(a) alpha = 0.025 per sign (standard)",
                            neg = neg_std, pos = pos_std),
                P |> mutate(panel = "(b) alpha* = 0.005 per sign (permutation-calibrated)",
                            neg = neg_tuned, pos = pos_tuned))
save2(base_map(P2) + pts(P2, neg, pos) + facet_wrap(~ panel, ncol = 2), "fig9_tuned_annual")
# Fig S4: seasons on the base mesh at alpha* (descriptive: alpha* not calibrated for seasons)
S <- C |> filter(tier == "secondary")
save2(base_map(S) + pts(S, neg_tuned, pos_tuned) + facet_wrap(~ SEASON, ncol = 3), "figS4_tuned_seasons", h = 7.5)

cat("\nFinished:", format(Sys.time()), "| summary:\n"); print(as.data.frame(read_csv(fs, show_col_types = FALSE)), digits = 3)
cat("Run finished.\n")
