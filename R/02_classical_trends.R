# =============================================================
# Step 2: Classical (cell-by-cell) trend analysis = BASELINE
#   - annual and seasonal mean GHI per grid cell per year
#   - Mann-Kendall test + Sen's slope at every cell
#   - lag-1 autocorrelation of detrended series
#   - Benjamini-Hochberg FDR control (Wilks 2016: alpha_FDR = 0.10)
#
# IMPORTANT: NASA POWER solar data come from GEWEX SRB (1984-2000) and
# CERES SYN1deg (2001 onward). NASA advises against trend analysis
# across that source change. So the MAIN analysis uses 2001-2023 only.
# Set START_YEAR <- 1984 to reproduce the full-record (contaminated) run.
#
# Input : data/power_ghi_india_monthly.csv   (from Step 1)
# Output: results/<tag>/..., figures/<tag>/...
# =============================================================

START_YEAR <- 2001
END_YEAR   <- 2023
tag <- paste0(START_YEAR, "_", END_YEAR)

library(dplyr); library(tidyr); library(readr); library(ggplot2)
library(trend); library(sf); library(rnaturalearth)

res_dir <- file.path("results", tag); fig_dir <- file.path("figures", tag)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

d <- read_csv("data/power_ghi_india_monthly.csv", show_col_types = FALSE) |>
  filter(YEAR >= START_YEAR, YEAR <= END_YEAR)

# ---- 1. Seasons (IMD convention) and yearly means ----
season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)",
                   "Monsoon (JJAS)", "Post-monsoon (OND)")
d <- d |>
  mutate(SEASON = case_when(MONTH %in% 1:2  ~ "Winter (JF)",
                            MONTH %in% 3:5  ~ "Pre-monsoon (MAM)",
                            MONTH %in% 6:9  ~ "Monsoon (JJAS)",
                            TRUE            ~ "Post-monsoon (OND)"))

seas <- d |> group_by(LON, LAT, YEAR, SEASON) |>
  summarise(GHI = mean(GHI), .groups = "drop")
ann  <- d |> group_by(LON, LAT, YEAR) |>
  summarise(GHI = mean(GHI), .groups = "drop") |>
  mutate(SEASON = "Annual")

agg <- bind_rows(ann, seas) |>
  mutate(SEASON = factor(SEASON, levels = season_levels)) |>
  arrange(SEASON, LON, LAT, YEAR)
write_csv(agg, file.path("data", paste0("ghi_seasonal_means_", tag, ".csv")))

# ---- 2. Mann-Kendall + Sen's slope at every cell ----
cell_trend <- function(y, yr) {
  r <- resid(lm(y ~ yr))
  c(mean_ghi = mean(y),
    sen      = unname(sens.slope(y)$estimates) * 10,   # per decade
    p_mk     = mk.test(y)$p.value,
    rho1     = acf(r, lag.max = 1, plot = FALSE)$acf[2])
}

res <- agg |>
  group_by(SEASON, LON, LAT) |>
  summarise(out = list(cell_trend(GHI, YEAR)), .groups = "drop") |>
  unnest_wider(out) |>
  group_by(SEASON) |>
  mutate(p_fdr      = p.adjust(p_mk, method = "BH"),
         sig_raw05  = p_mk  < 0.05,
         sig_raw10  = p_mk  < 0.10,
         sig_fdr10  = p_fdr < 0.10) |>
  ungroup()
write_csv(res, file.path(res_dir, "classical_trends.csv"))

# ---- 3. Summary table (goes into the paper) ----
summ <- res |>
  group_by(SEASON) |>
  summarise(cells          = n(),
            mean_GHI       = round(mean(mean_ghi), 2),
            median_slope   = round(median(sen), 3),
            pct_negative   = round(100 * mean(sen < 0), 1),
            pct_raw05      = round(100 * mean(sig_raw05), 1),
            pct_raw10      = round(100 * mean(sig_raw10), 1),
            pct_FDR10      = round(100 * mean(sig_fdr10), 1),
            median_rho1    = round(median(rho1), 2))
write_csv(summ, file.path(res_dir, "classical_summary.csv"))
print(summ, width = Inf)

# ---- 4. Maps ----
india <- ne_countries(country = "India", scale = "medium", returnclass = "sf")
lim <- max(abs(res$sen))

p1 <- ggplot(res) +
  geom_tile(aes(LON, LAT, fill = sen)) +
  geom_point(data = filter(res, sig_fdr10), aes(LON, LAT),
             size = 0.5, colour = "black") +
  geom_sf(data = india, fill = NA, colour = "grey30", linewidth = 0.3) +
  scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-lim, lim),
                       name = expression(paste("Sen slope (kWh ", m^-2, " ", day^-1, " ", decade^-1, ")"))) +
  facet_wrap(~ SEASON, ncol = 3) +
  coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL,
       title = paste0("Cell-wise Sen slopes, ", START_YEAR, "-", END_YEAR),
       caption = "Dots: significant after FDR control (Benjamini-Hochberg, alpha = 0.10)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "fig_classical_sen_slope.png"), p1, width = 9, height = 7.5, dpi = 300)

p2 <- ggplot(filter(res, SEASON == "Annual")) +
  geom_tile(aes(LON, LAT, fill = mean_ghi)) +
  geom_sf(data = india, fill = NA, colour = "grey30", linewidth = 0.3) +
  scale_fill_viridis_c(option = "inferno",
                       name = expression(paste("Mean GHI (kWh ", m^-2, " ", day^-1, ")"))) +
  coord_sf(expand = FALSE) + labs(x = NULL, y = NULL) + theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "fig_mean_ghi_annual.png"), p2, width = 5.5, height = 5.5, dpi = 300)

p3 <- ggplot(res, aes(rho1)) +
  geom_histogram(bins = 30, fill = "grey60", colour = "white") +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~ SEASON, ncol = 3) +
  labs(x = "Lag-1 autocorrelation of detrended series", y = "Number of cells") +
  theme_bw(base_size = 10)
ggsave(file.path(fig_dir, "fig_rho1_hist.png"), p3, width = 8, height = 5, dpi = 300)

cat("\nDone:", tag, "- send me the printed table and the PNGs in", fig_dir, "\n")
