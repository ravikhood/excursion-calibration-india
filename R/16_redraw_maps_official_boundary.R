# =============================================================
# Step 16: redraw every map figure with India's official boundary
# (Natural Earth 10m "Admin 0 - Countries, point of view: India", which includes
# Pakistan-occupied Kashmir, Gilgit-Baltistan and Aksai Chin).
# NO model fitting: reads saved results only. Runtime ~2-3 min.
#
# The ANALYSIS DOMAIN IS UNCHANGED (274 cells from script 01). Areas inside the
# official boundary without analysis cells are shaded light grey and explained in
# the captions.
#
# To use a Survey of India boundary instead, set BOUNDARY_FILE to its shapefile/geojson.
#
# Inputs : results/rev/{A0_seasonal_means_dayweighted,cells,A_classical,A_cells,C_truths}.csv
#          results/tuned/15_cells.csv
# Output : figures/paper_official/  (PNG + PDF), plus copies named as in Overleaf:
#          figures/paper_official/overleaf/Figure_1.png, Figure_6.png, Figure_7.png,
#          Figure_S2.png ... Figure_S6.png   -> upload these to Overleaf "figures/"
# =============================================================

BOUNDARY_FILE <- NULL      # e.g. "data/boundary/survey_of_india_outline.shp"; NULL = Natural Earth India view

pk <- c("sf", "ggplot2", "dplyr", "readr", "tidyr", "scales")
new <- pk[!pk %in% installed.packages()[, "Package"]]; if (length(new)) install.packages(new)
suppressPackageStartupMessages({ library(sf); library(ggplot2); library(dplyr); library(readr); library(tidyr) })
sf_use_s2(FALSE)

# ---------------- boundary ----------------
bdir <- file.path("data", "boundary"); dir.create(bdir, recursive = TRUE, showWarnings = FALSE)
if (is.null(BOUNDARY_FILE)) {
  base_url <- "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/10m_cultural/ne_10m_admin_0_countries_ind"
  for (e in c("shp", "shx", "dbf", "prj")) {
    f <- file.path(bdir, paste0("ne_10m_admin_0_countries_ind.", e))
    if (!file.exists(f)) download.file(paste0(base_url, ".", e), f, mode = "wb", quiet = TRUE)
  }
  w <- st_read(file.path(bdir, "ne_10m_admin_0_countries_ind.shp"), quiet = TRUE)
  india <- w[w$ADM0_A3 == "IND", ]
  src <- "Natural Earth 10m Admin 0 (India point of view)"
} else {
  india <- st_read(BOUNDARY_FILE, quiet = TRUE); src <- BOUNDARY_FILE
}
india <- st_make_valid(st_union(st_transform(india, 4326)))
chk <- st_contains(india, st_sfc(st_point(c(73.47, 34.37)), st_point(c(79.5, 35.2)), crs = 4326), sparse = FALSE)
cat("Boundary:", src, "| contains Muzaffarabad (PoK):", chk[1], "| contains Aksai Chin:", chk[2], "\n")
if (!all(chk)) stop("Boundary does not include PoK / Aksai Chin. Stop and check the setup.")
india <- st_sf(geometry = india)

# ---------------- data ----------------
season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)", "Monsoon (JJAS)", "Post-monsoon (OND)")
cells <- read_csv("results/rev/cells.csv", show_col_types = FALSE) |> select(cell, LON, LAT)
out <- file.path("figures", "paper_official"); ovl <- file.path(out, "overleaf")
dir.create(ovl, recursive = TRUE, showWarnings = FALSE)
lab_tr <- expression(paste("kWh ", m^-2, " ", day^-1, " ", decade^-1))
lab_ghi <- expression(paste("kWh ", m^-2, " ", day^-1))
div <- scale_fill_gradient2(name = lab_tr, low = "#2166ac", mid = "white", high = "#b2182b",
                            midpoint = 0, limits = c(-0.3, 0.3), oob = scales::squish)
base_map <- function(d, fill_col, fill_scale) {
  ggplot(d) +
    geom_sf(data = india, fill = "grey88", colour = NA) +                    # official territory, no-data shade
    geom_tile(aes(LON, LAT, fill = .data[[fill_col]]), width = 1, height = 1) +
    geom_sf(data = india, fill = NA, colour = "grey20", linewidth = 0.35) +
    fill_scale + coord_sf(xlim = c(68, 97.5), ylim = c(6.5, 37.3), expand = FALSE) +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 10) +
    theme(legend.position = "bottom", legend.key.width = grid::unit(2.2, "cm"),
          panel.grid = element_line(colour = "grey93", linewidth = 0.3))
}
pts <- function(d, neg, pos) list(
  geom_point(data = d[d[[neg]], ], aes(LON, LAT), size = 0.35, colour = "black"),
  geom_point(data = d[d[[pos]], ], aes(LON, LAT), size = 0.9, shape = 4, colour = "black"))
save2 <- function(p, name, overleaf, w = 9, h = 7.5) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  file.copy(file.path(out, paste0(name, ".png")), file.path(ovl, paste0(overleaf, ".png")), overwrite = TRUE)
  cat("  ", overleaf, "<-", name, "\n")
}
fac <- function(d) mutate(d, SEASON = factor(SEASON, levels = season_levels))

# Figure 1: mean annual GHI 2001-2023
agg <- read_csv("results/rev/A0_seasonal_means_dayweighted.csv", show_col_types = FALSE)
m1 <- agg |> filter(SEASON == "Annual", YEAR %in% 2001:2023) |> group_by(LON, LAT) |>
  summarise(GHI = mean(GHI), .groups = "drop") |> inner_join(cells, by = c("LON", "LAT"))
stopifnot(nrow(m1) == 274)
save2(base_map(m1, "GHI", scale_fill_viridis_c(name = lab_ghi, option = "inferno")), "fig1_mean_ghi", "Figure_1", w = 6, h = 6.5)

# Script-15 results (4000 samples): Figure 6, Figure 7, Figure S4, Figure S6
C15 <- read_csv("results/tuned/15_cells.csv", show_col_types = FALSE) |> inner_join(cells, by = "cell") |> fac()
P <- filter(C15, tier == "primary")
P2 <- bind_rows(P |> mutate(panel = "(a) alpha = 0.025 per sign (standard)", neg = neg_std, pos = pos_std),
                P |> mutate(panel = "(b) alpha* = 0.005 per sign (permutation-calibrated)", neg = neg_tuned, pos = pos_tuned))
save2(base_map(P2, "est", div) + pts(P2, "neg", "pos") + facet_wrap(~ panel, ncol = 2), "fig9_tuned_annual", "Figure_6", h = 5.2)
B <- filter(C15, tier == "secondary")
save2(base_map(B, "sd", scale_fill_viridis_c(name = lab_tr)) + facet_wrap(~ SEASON, ncol = 3), "fig5_v2_sd", "Figure_7")
save2(base_map(B, "est", div) + pts(B, "neg_tuned", "pos_tuned") + facet_wrap(~ SEASON, ncol = 3), "figS4_tuned_seasons", "Figure_S4")
save2(base_map(B, "est", div) + pts(B, "neg_std", "pos_std") + facet_wrap(~ SEASON, ncol = 3), "fig4_v2_trend", "Figure_S6")

# Figure S2: v1 annual (script 08)
A <- read_csv("results/rev/A_cells.csv", show_col_types = FALSE) |> inner_join(cells, by = "cell") |> fac()
v1 <- A |> filter(variant == "v1_base") |> mutate(jneg025 = as.logical(jneg025), jpos025 = as.logical(jpos025))
save2(base_map(v1, "est", div) + pts(v1, "jneg025", "jpos025") + facet_wrap(~ SEASON), "figS2_v1_trend", "Figure_S2", w = 6, h = 6.5)

# Figure S3: synthetic truths
tr <- read_csv("results/rev/C_truths.csv", show_col_types = FALSE) |>
  pivot_longer(-cell, names_to = "truth", values_to = "beta") |> inner_join(cells, by = "cell") |>
  mutate(truth = factor(truth, levels = c("S1_smooth", "S2_mixed", "S3_sharp", "S4_sparse", "S5_hetero"),
                        labels = c("S1 smooth", "S2 mixed", "S3 sharp", "S4 sparse", "S5 hetero")))
divS3 <- scale_fill_gradient2(name = "True trend", low = "#2166ac", mid = "white", high = "#b2182b",
                              midpoint = 0, limits = c(-0.3, 0.3), oob = scales::squish)
save2(base_map(tr, "beta", divS3) + facet_wrap(~ truth, nrow = 1) + theme(axis.text = element_text(size = 6)),
      "figS3_truths", "Figure_S3", w = 12, h = 3.8)

# Figure S5: Sen slopes + BH (q = 0.10)
cls <- read_csv("results/rev/A_classical.csv", show_col_types = FALSE) |>
  group_by(SEASON) |> mutate(fdr = p.adjust(p, "BH") < 0.10) |> ungroup() |> inner_join(cells, by = "cell") |> fac()
save2(base_map(cls, "est", div) + geom_point(data = filter(cls, fdr), aes(LON, LAT), size = 0.35) +
        facet_wrap(~ SEASON, ncol = 3), "fig3_sen_fdr", "Figure_S5")

writeLines(c(paste("Boundary source:", src), capture.output(sessionInfo())), file.path(out, "16_sessionInfo.txt"))
cat("\nDone. Upload the 8 PNGs in", ovl, "to Overleaf 'figures/' (replace existing files).\n")
