# =============================================================
# Step 1: Download NASA POWER monthly surface solar radiation
# (ALLSKY_SFC_SW_DWN, kWh/m^2/day) over India, 1984-2023
# Output: data/power_ghi_india_monthly.csv
# =============================================================

# ---- 0. Install packages (run once) ----
pkgs <- c("nasapower", "sf", "rnaturalearth", "rnaturalearthdata",
          "dplyr", "tidyr", "readr", "ggplot2", "trend")
new <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new)) install.packages(new)

# INLA is not on CRAN (needed from Step 3 onward)
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages("INLA",
    repos = c(getOption("repos"),
              INLA = "https://inla.r-inla-download.org/R/stable"),
    type = "binary", dep = TRUE)   # binary avoids the Windows 'cp not found' error
}
if (!requireNamespace("excursions", quietly = TRUE)) install.packages("excursions")

library(nasapower); library(dplyr); library(tidyr); library(readr)
library(sf); library(rnaturalearth)

dir.create("data", showWarnings = FALSE)

# ---- 1. Split India's bounding box into tiles ----
# Solar parameters are on a 1-degree grid, so the POWER regional API
# accepts boxes of at most 5 x 5 degrees (42 tiles in total)
lon_breaks <- seq(68, 98, by = 5)
lat_breaks <- seq(6, 41, by = 5)
tiles <- expand.grid(i = 1:(length(lon_breaks) - 1),
                     j = 1:(length(lat_breaks) - 1))

out <- list()
for (k in seq_len(nrow(tiles))) {
  bb <- c(lon_breaks[tiles$i[k]],     lat_breaks[tiles$j[k]],
          lon_breaks[tiles$i[k] + 1], lat_breaks[tiles$j[k] + 1])
  cat("Tile", k, "of", nrow(tiles), ":", bb, "\n")
  res <- tryCatch(
    get_power(community    = "re",
              lonlat       = bb,
              pars         = "ALLSKY_SFC_SW_DWN",
              dates        = c("1984-01-01", "2023-12-31"),
              temporal_api = "monthly"),
    error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
  if (!is.null(res)) out[[k]] <- res
  Sys.sleep(2)   # be polite to the server
}
raw <- bind_rows(out) |> distinct()

# ---- 2. Reshape to long format: one row per cell-year-month ----
mon <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
long <- raw |>
  select(LON, LAT, YEAR, any_of(mon)) |>
  pivot_longer(any_of(mon), names_to = "MONTH", values_to = "GHI") |>
  mutate(MONTH = match(MONTH, mon),
         GHI   = ifelse(GHI < 0, NA, GHI))   # -999 = missing

# ---- 3. Keep only grid cells inside India ----
india <- ne_countries(country = "India", scale = "medium", returnclass = "sf")
cells <- long |> distinct(LON, LAT) |>
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, remove = FALSE)
inside <- cells[lengths(st_intersects(cells, india)) > 0, ] |> st_drop_geometry()
long <- semi_join(long, inside, by = c("LON", "LAT"))

write_csv(long, "data/power_ghi_india_monthly.csv")

# ---- 4. Quick check: report these numbers back ----
cat("\nGrid cells inside India:", nrow(inside), "\n")
cat("Years:", range(long$YEAR), "\n")
cat("Missing values:", sum(is.na(long$GHI)), "\n")
print(summary(long$GHI))
