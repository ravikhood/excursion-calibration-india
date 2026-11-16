# =============================================================
# Step 9: Independent check with ERA5 reanalysis (reviewer comments 2, 3, 13)
#
# ERA5 is produced independently of NASA POWER/CERES (ECMWF model +
# data assimilation). If ERA5 shows similar dimming over India, the
# trend is not just a CERES/POWER product artefact.
#
# What it does:
#   1. Downloads ERA5 monthly surface solar radiation downwards (SSRD),
#      2001-2023, India box, from the Copernicus CDS (about 50 MB).
#   2. Converts J m-2 (daily accumulation in monthly means) to kWh m-2 day-1.
#   3. Averages the 0.25 deg ERA5 grid onto the same 274 one-degree cells.
#   4. Day-weighted seasonal means (same as script 08, task A0).
#   5. Classical trends (India-mean OLS + Sen; per-cell Sen + MK + BH-FDR)
#      for 2001-2023 and 2001-2016, and the same for NASA POWER.
#
# Before running: free CDS account + accept the licence of
# "ERA5 monthly averaged data on single levels from 1940 to present"
# (Download tab, bottom of page). The script asks for your API token once;
# never paste the token into a chat or into this file.
#
# Output: data/era5/*.nc, data/era5_ghi_india_monthly.csv,
#         results/era5/*.csv   (analysed by the Python scripts in analysis/)
# Runtime: download 5-30 min (CDS queue), processing ~2 min.
# =============================================================

pk <- c("ecmwfr", "terra", "dplyr", "tidyr", "readr", "trend")
new <- pk[!pk %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)
suppressPackageStartupMessages({
  library(ecmwfr); library(terra); library(dplyr); library(tidyr); library(readr); library(trend)
})
dir.create("data/era5", recursive = TRUE, showWarnings = FALSE)
dir.create("results/era5", recursive = TRUE, showWarnings = FALSE)
nc_file <- file.path("data/era5", "era5_ssrd_india_2001_2023.nc")

# ---- 1. API key (stored securely in your system keyring, asked once) ----
has_key <- tryCatch({ nchar(wf_get_key()) > 0 }, error = function(e) FALSE)
if (!has_key) {
  tok <- readline("Paste your CDS API token and press Enter: ")
  wf_set_key(key = trimws(tok))
  rm(tok)
}

# ---- 2. Download (skipped if the file already exists) ----
if (!file.exists(nc_file)) {
  request <- list(
    dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
    product_type    = "monthly_averaged_reanalysis",
    variable        = "surface_solar_radiation_downwards",
    year            = as.character(2001:2023),
    month           = sprintf("%02d", 1:12),
    time            = "00:00",
    area            = c(38, 67, 6, 98),          # North, West, South, East
    data_format     = "netcdf",
    download_format = "unarchived",
    target          = basename(nc_file))
  cat("Submitting CDS request (may wait in the queue)...\n")
  wf_request(request, path = "data/era5", transfer = TRUE)
}
if (!file.exists(nc_file)) stop("Download did not produce ", nc_file,
                                ". Check the CDS website (licence accepted? request finished?).")

# ---- 3. Read, convert units, regrid to the 274 one-degree cells ----
r <- rast(nc_file)
cat("ERA5 layers:", nlyr(r), "\n")
# The new CDS netCDF stores time as 'valid_time' (seconds since 1970), which terra
# does not always pick up. Read it directly with ncdf4; fall back to the request order.
tt <- tryCatch(as.Date(time(r)), error = function(e) NULL)
if (is.null(tt) || all(is.na(tt))) {
  if (!requireNamespace("ncdf4", quietly = TRUE)) install.packages("ncdf4")
  nc <- ncdf4::nc_open(nc_file)
  tname <- intersect(c("valid_time", "time"), names(nc$var) |> c(names(nc$dim)))[1]
  if (!is.na(tname)) {
    tv <- ncdf4::ncvar_get(nc, tname)
    un <- ncdf4::ncatt_get(nc, tname, "units")$value
    cat("Time variable:", tname, "| units:", un, "\n")
    origin <- as.POSIXct(sub(".*since ", "", un), tz = "UTC")
    mult <- if (grepl("^hours", un)) 3600 else if (grepl("^days", un)) 86400 else 1
    tt <- as.Date(origin + tv * mult)
  }
  ncdf4::nc_close(nc)
}
if (is.null(tt) || all(is.na(tt)) || length(tt) != nlyr(r)) {
  stopifnot(nlyr(r) == 23 * 12)                      # 2001-2023 x 12 months, in request order
  tt <- seq(as.Date("2001-01-01"), by = "month", length.out = nlyr(r))
  cat("Time axis rebuilt from request order.\n")
}
cat("Time range:", format(min(tt)), "to", format(max(tt)), "| unique months:",
    length(unique(format(tt, "%Y-%m"))), "\n")
stopifnot(length(unique(format(tt, "%Y-%m"))) == nlyr(r))  # no duplicated months
df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
names(df) <- c("x", "y", format(tt, "%Y-%m"))
cells <- read_csv("results/rev/cells.csv", show_col_types = FALSE) |> select(LON, LAT, cell, w_area)

era <- df |>
  mutate(LON = floor(x) + 0.5, LAT = floor(y) + 0.5) |>
  select(-x, -y) |>
  pivot_longer(-c(LON, LAT), names_to = "ym", values_to = "ssrd") |>
  group_by(LON, LAT, ym) |> summarise(ssrd = mean(ssrd), n_pts = n(), .groups = "drop") |>
  inner_join(cells, by = c("LON", "LAT")) |>
  mutate(YEAR = as.integer(substr(ym, 1, 4)), MONTH = as.integer(substr(ym, 6, 7)),
         GHI = ssrd / 3.6e6) |>                     # J m-2 per day -> kWh m-2 day-1
  select(LON, LAT, YEAR, MONTH, GHI, n_pts) |> arrange(LON, LAT, YEAR, MONTH)
cat("Cells matched:", n_distinct(paste(era$LON, era$LAT)), "of", nrow(cells),
    "| ERA5 points per cell (median):", median(era$n_pts), "\n")
write_csv(era, "data/era5_ghi_india_monthly.csv")

# ---- 4. Day-weighted seasonal means for both datasets ----
ndays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
seasonal <- function(d) {
  d <- d |> mutate(wd = ndays[MONTH],
                   SEASON = case_when(MONTH %in% 1:2 ~ "Winter (JF)",
                                      MONTH %in% 3:5 ~ "Pre-monsoon (MAM)",
                                      MONTH %in% 6:9 ~ "Monsoon (JJAS)",
                                      TRUE ~ "Post-monsoon (OND)"))
  bind_rows(
    d |> group_by(LON, LAT, YEAR) |> summarise(GHI = weighted.mean(GHI, wd), .groups = "drop") |>
      mutate(SEASON = "Annual"),
    d |> group_by(LON, LAT, YEAR, SEASON) |> summarise(GHI = weighted.mean(GHI, wd), .groups = "drop"))
}
era_s  <- seasonal(era) |> mutate(dataset = "ERA5")
nasa_s <- read_csv("results/rev/A0_seasonal_means_dayweighted.csv", show_col_types = FALSE) |>
  filter(YEAR >= 2001, YEAR <= 2023) |> mutate(dataset = "NASA_POWER")
both <- bind_rows(era_s, nasa_s) |> inner_join(cells, by = c("LON", "LAT"))
write_csv(both, "results/era5/seasonal_means_both.csv")

# ---- 5. Trends ----
bhq <- function(p, q = 0.10) p.adjust(p, "BH") < q
periods <- list(`2001-2023` = 2001:2023, `2001-2016` = 2001:2016)

nat <- both |> group_by(dataset, SEASON, YEAR) |>
  summarise(GHI = sum(GHI * w_area) / sum(w_area), .groups = "drop")
write_csv(nat, "results/era5/national_series.csv")

nat_tr <- bind_rows(lapply(names(periods), function(pn) {
  nat |> filter(YEAR %in% periods[[pn]]) |> group_by(dataset, SEASON) |>
    summarise(period = pn,
              mean_GHI = mean(GHI),
              ols = 10 * coef(lm(GHI ~ YEAR))[2],
              ols_lo = 10 * confint(lm(GHI ~ YEAR))[2, 1],
              ols_hi = 10 * confint(lm(GHI ~ YEAR))[2, 2],
              sen = 10 * unname(sens.slope(GHI)$estimates),
              mk_p = mk.test(GHI)$p.value, .groups = "drop")
}))
write_csv(nat_tr, "results/era5/national_trends.csv")

cell_tr <- bind_rows(lapply(names(periods), function(pn) {
  both |> filter(YEAR %in% periods[[pn]]) |> arrange(YEAR) |>
    group_by(dataset, SEASON, LON, LAT, cell) |>
    summarise(period = pn, sen = 10 * unname(sens.slope(GHI)$estimates),
              p = mk.test(GHI)$p.value, .groups = "drop")
})) |> group_by(dataset, SEASON, period) |> mutate(sig_fdr10 = bhq(p)) |> ungroup()
write_csv(cell_tr, "results/era5/cell_trends.csv")

# Correlation of yearly India-mean anomalies between datasets (agreement in variability)
cors <- nat |> pivot_wider(names_from = dataset, values_from = GHI) |>
  group_by(SEASON) |> summarise(r_interannual = cor(ERA5, NASA_POWER), .groups = "drop")
write_csv(cors, "results/era5/interannual_correlation.csv")

options(width = 200)
cat("\n===== India-mean trends (kWh m-2 day-1 per decade) =====\n")
print(as.data.frame(nat_tr |> filter(SEASON == "Annual") |> select(-SEASON)), digits = 3)
cat("\nDone. Run finished.\n")
