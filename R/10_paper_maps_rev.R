# =============================================================
# Step 10: Redraw the map figures from the saved revision results
# (day-weighted seasons, alpha = 0.025 per sign). NO model fitting:
# reads CSVs only, runs in under a minute.
# Input : results/rev/A_cells.csv, results/rev/A_classical.csv, results/rev/cells.csv
# Output: figures/paper/fig3_sen_fdr, fig4_v2_trend, fig5_v2_sd, figS2_v1_trend (.pdf/.png)
# Then copy the four PNGs into paper/figures/ (same file names).
# =============================================================
library(dplyr); library(readr); library(ggplot2); library(sf); library(rnaturalearth)

out <- file.path("figures", "paper"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)", "Monsoon (JJAS)", "Post-monsoon (OND)")
india <- ne_countries(country = "India", scale = "medium", returnclass = "sf")
cells <- read_csv("results/rev/cells.csv", show_col_types = FALSE) |> select(cell, LON, LAT)
lab_tr <- expression(paste("kWh ", m^-2, " ", day^-1, " ", decade^-1))

map_base <- function(d, fill, fill_scale) {
  ggplot(d) + geom_tile(aes(LON, LAT, fill = {{ fill }})) +
    geom_sf(data = india, fill = NA, colour = "grey30", linewidth = 0.3) +
    fill_scale + facet_wrap(~ SEASON, ncol = 3) + coord_sf(expand = FALSE) +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 10) + theme(legend.position = "bottom")
}
div <- scale_fill_gradient2(name = lab_tr, low = "#2166ac", mid = "white", high = "#b2182b",
                            midpoint = 0, limits = c(-0.3, 0.3), oob = scales::squish)
save2 <- function(p, name, h = 7.5) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = 9, height = h)
  ggsave(file.path(out, paste0(name, ".png")), p, width = 9, height = h, dpi = 300)
}

# ---- Fig 3: cell-wise Sen slope + BH-FDR (q = 0.10) ----
cls <- read_csv("results/rev/A_classical.csv", show_col_types = FALSE) |>
  group_by(SEASON) |> mutate(fdr = p.adjust(p, "BH") < 0.10) |> ungroup() |>
  left_join(cells, by = "cell") |> mutate(SEASON = factor(SEASON, levels = season_levels))
p3 <- map_base(cls, est, div) +
  geom_point(data = filter(cls, fdr), aes(LON, LAT), size = 0.35, colour = "black")
save2(p3, "fig3_sen_fdr")

# ---- Fig 4, 5, S2: SPDE models ----
A <- read_csv("results/rev/A_cells.csv", show_col_types = FALSE) |>
  left_join(cells, by = "cell") |> mutate(SEASON = factor(SEASON, levels = season_levels))
v2 <- filter(A, variant == "v2_base")
p4 <- map_base(v2, est, div) +
  geom_point(data = filter(v2, jneg025), aes(LON, LAT), size = 0.35, colour = "black") +
  geom_point(data = filter(v2, jpos025), aes(LON, LAT), size = 0.9, shape = 4, colour = "black")
save2(p4, "fig4_v2_trend")
p5 <- map_base(v2, sd, scale_fill_viridis_c(name = lab_tr))
save2(p5, "fig5_v2_sd")
v1 <- filter(A, variant == "v1_base")          # annual only
pS2 <- map_base(v1, est, div) +
  geom_point(data = filter(v1, jneg025), aes(LON, LAT), size = 0.35, colour = "black") +
  geom_point(data = filter(v1, jpos025), aes(LON, LAT), size = 0.9, shape = 4, colour = "black")
save2(pS2, "figS2_v1_trend", h = 4.5)
cat("Maps written to", out, "- copy the 4 PNGs into paper/figures/.\n")
