# =============================================================
# Step 14: Compare surface solar radiation products (referee comments 7-8)
#
# Products (all monthly, regridded to the same 274 one-degree cells):
#   NASA POWER (CERES SYN1deg Ed4.1)  - already in data/  (not for trend studies per CERES)
#   ERA5                              - already in data/  (script 09)
#   CERES EBAF-Surface Ed4.2.1        - YOU download (see below); recommended by CERES for trends
#   MERRA-2 (M2TMNXRAD, SWGDN)        - YOU download (see below); assimilates aerosols
# NOTE: EBAF-Surface is built from SYN1deg inputs (constrained to EBAF-TOA), so
# EBAF and POWER are NOT independent. ERA5 and MERRA-2 are reanalyses.
#
# DOWNLOADS (free NASA Earthdata login; never paste passwords into chat):
#  1. EBAF-Surface: https://ceres.larc.nasa.gov/data/  -> EBAF -> "Order data"
#     Product: EBAF Ed4.2.1 (latest; Ed4.2 also works), Surface, monthly, parameter "Surface Shortwave Flux Down,
#     All-Sky" (sfc_sw_down_all_mon), region lat 6-38 N, lon 67-98 E, 2001-01 to 2023-12,
#     format NetCDF. Save the .nc file into  data/ebaf/
#  2. MERRA-2: https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary -> "Subset / Get Data"
#     Variable SWGDN, region 67,6,98,38 (W,S,E,N), 2001-01 to 2023-12, NetCDF.
#     Download all monthly files into  data/merra2/   (one file per month is fine)
# The script skips a product whose folder is empty, so you can run it with one product.
#
# Output: results/products/
#   products_seasonal_means.csv, products_national_series.csv,
#   products_national_trends.csv, products_cell_trends.csv,
#   products_difference_tests.csv, products_v2_india.csv
# Runtime: ~5 min without the SPDE fits; ~45 min with them (4000 samples, base mesh).
# =============================================================

RUN_V2 <- TRUE        # SPDE v2 India-mean trend for EBAF and MERRA-2 (2001-2023 and 2001-2016)
N_SAMP <- 4000

pk <- c("ncdf4", "dplyr", "tidyr", "readr", "trend")
new <- pk[!pk %in% installed.packages()[, "Package"]]; if (length(new)) install.packages(new)
suppressPackageStartupMessages({
  library(ncdf4); library(dplyr); library(tidyr); library(readr); library(trend)
  if (RUN_V2) { library(INLA); library(fmesher) }
})
out <- file.path("results", "products"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
cells <- read_csv("results/rev/cells.csv", show_col_types = FALSE) |>
  mutate(X_km = (LON - 82) * 111.32 * cos(22 * pi / 180), Y_km = (LAT - 22) * 110.57)
ndays <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
W_TO_KWH <- 24 / 1000                                    # W m-2 -> kWh m-2 day-1

# ---- read a gridded monthly variable from one or more netCDF files ----
read_nc_monthly <- function(files, var_pattern) {
  bind_rows(lapply(files, function(f) {
    nc <- nc_open(f); on.exit(nc_close(nc))
    v <- grep(var_pattern, names(nc$var), value = TRUE)[1]
    if (is.na(v)) stop("No variable matching '", var_pattern, "' in ", basename(f),
                       ". Variables: ", paste(names(nc$var), collapse = ", "))
    lon <- ncvar_get(nc, grep("^lon", names(nc$dim), value = TRUE)[1])
    lat <- ncvar_get(nc, grep("^lat", names(nc$dim), value = TRUE)[1])
    tn <- grep("^time", names(nc$dim), value = TRUE)[1]
    tv <- ncvar_get(nc, tn); un <- ncatt_get(nc, tn, "units")$value
    origin <- as.POSIXct(sub(".*since ", "", un), tz = "UTC")
    mult <- if (grepl("^hours", un)) 3600 else if (grepl("^days", un)) 86400 else if (grepl("^minutes", un)) 60 else 1
    tt <- as.Date(origin + tv * mult)
    x <- ncvar_get(nc, v)                                   # lon x lat (x time)
    if (length(dim(x)) == 2) dim(x) <- c(dim(x), 1)
    expand.grid(lon = as.numeric(lon), lat = as.numeric(lat), k = seq_along(tt)) |>
      mutate(value = as.numeric(x), date = tt[k]) |> select(-k)
  }))
}
to_cells <- function(g, product) {      # area-overlap (conservative) average of grid boxes into our 1-degree cells
  # Works for any regular grid: EBAF (1 deg, same centres -> weight 1 on one box) and
  # MERRA-2 (0.625 x 0.5 deg, boxes straddle cell edges -> fractional weights).
  g <- g |> mutate(lon = ifelse(lon > 180, lon - 360, lon)) |>
    filter(as.integer(format(date, "%Y")) >= 2001, as.integer(format(date, "%Y")) <= 2023)
  lons <- sort(unique(g$lon)); lats <- sort(unique(g$lat))
  dlon <- min(diff(lons)); dlat <- min(diff(lats))
  ov <- function(p, d, c) pmax(0, pmin(p + d / 2, c + 0.5) - pmax(p - d / 2, c - 0.5))   # overlap length (deg)
  wx <- crossing(cells |> select(cell, LON), lon = lons) |> mutate(wx = ov(lon, dlon, LON)) |> filter(wx > 1e-9)
  wy <- crossing(cells |> select(cell, LAT), lat = lats) |> mutate(wy = ov(lat, dlat, LAT)) |> filter(wy > 1e-9)
  W <- inner_join(wx, wy, by = "cell", relationship = "many-to-many") |> mutate(w = wx * wy) |> select(cell, lon, lat, w)
  g |> inner_join(W, by = c("lon", "lat"), relationship = "many-to-many") |>
    group_by(cell, date) |>
    summarise(GHI = sum(w * value, na.rm = TRUE) / sum(w * !is.na(value)) * W_TO_KWH, n_pts = n(), .groups = "drop") |>
    mutate(YEAR = as.integer(format(date, "%Y")), MONTH = as.integer(format(date, "%m"))) |>
    inner_join(cells |> select(LON, LAT, cell), by = "cell") |>
    select(LON, LAT, YEAR, MONTH, GHI, n_pts, cell) |> mutate(product = product)
}
seasonal <- function(d) {
  d <- d |> mutate(wd = ndays[MONTH], SEASON = case_when(MONTH %in% 1:2 ~ "Winter (JF)", MONTH %in% 3:5 ~ "Pre-monsoon (MAM)",
                                                         MONTH %in% 6:9 ~ "Monsoon (JJAS)", TRUE ~ "Post-monsoon (OND)"))
  bind_rows(d |> group_by(product, LON, LAT, cell, YEAR) |> summarise(GHI = weighted.mean(GHI, wd), .groups = "drop") |> mutate(SEASON = "Annual"),
            d |> group_by(product, LON, LAT, cell, YEAR, SEASON) |> summarise(GHI = weighted.mean(GHI, wd), .groups = "drop"))
}

# ---- assemble products ----
monthly <- list()
pw <- read_csv("data/power_ghi_india_monthly.csv", show_col_types = FALSE)
monthly$NASA_POWER <- pw |> filter(YEAR >= 2001, YEAR <= 2023) |>
  inner_join(cells |> select(LON, LAT, cell), by = c("LON", "LAT")) |> mutate(product = "NASA_POWER") |>
  select(LON, LAT, cell, YEAR, MONTH, GHI, product)
if (file.exists("data/era5_ghi_india_monthly.csv"))
  monthly$ERA5 <- read_csv("data/era5_ghi_india_monthly.csv", show_col_types = FALSE) |>
    inner_join(cells |> select(LON, LAT, cell), by = c("LON", "LAT")) |> mutate(product = "ERA5") |>
    select(LON, LAT, cell, YEAR, MONTH, GHI, product)
fe <- list.files("data/ebaf", pattern = "\\.nc4?$", full.names = TRUE)
if (length(fe)) monthly$EBAF <- to_cells(read_nc_monthly(fe, "sfc_sw_down_all_mon"), "EBAF")
fm <- list.files("data/merra2", pattern = "\\.nc4?$", full.names = TRUE)
if (length(fm)) monthly$MERRA2 <- to_cells(read_nc_monthly(fm, "^SWGDN$"), "MERRA2")
cat("Products found:", paste(names(monthly), collapse = ", "), "\n")
for (p in names(monthly)) {
  m <- monthly[[p]]
  cat(sprintf("  %-10s cells %3d | months %3d | mean %.2f kWh/m2/day\n", p, n_distinct(m$cell),
              n_distinct(paste(m$YEAR, m$MONTH)), mean(m$GHI, na.rm = TRUE)))
  if (n_distinct(paste(m$YEAR, m$MONTH)) != 276) warning(p, ": expected 276 months (2001-2023)")
  if (n_distinct(m$cell) != nrow(cells)) warning(p, ": not all 274 cells matched")
}
both <- seasonal(bind_rows(monthly)) |> inner_join(cells |> select(cell, w_area), by = "cell")
write_csv(both, file.path(out, "products_seasonal_means.csv"))

# ---- national series and trends ----
nat <- both |> group_by(product, SEASON, YEAR) |> summarise(GHI = sum(GHI * w_area) / sum(w_area), .groups = "drop")
write_csv(nat, file.path(out, "products_national_series.csv"))
periods <- list(`2001-2023` = 2001:2023, `2001-2016` = 2001:2016)
nat_tr <- bind_rows(lapply(names(periods), function(pn) nat |> filter(YEAR %in% periods[[pn]]) |>
  group_by(product, SEASON) |> summarise(period = pn, mean_GHI = mean(GHI),
    ols = 10 * coef(lm(GHI ~ YEAR))[2], ols_lo = 10 * confint(lm(GHI ~ YEAR))[2, 1],
    ols_hi = 10 * confint(lm(GHI ~ YEAR))[2, 2], sen = 10 * unname(sens.slope(GHI)$estimates),
    mk_p = mk.test(GHI)$p.value, .groups = "drop")))
write_csv(nat_tr, file.path(out, "products_national_trends.csv"))

# ---- difference series vs POWER (step after 2016?) ----
ann <- nat |> filter(SEASON == "Annual") |> pivot_wider(names_from = product, values_from = GHI)
others <- setdiff(names(monthly), "NASA_POWER")
diffs <- bind_rows(lapply(others, function(p) {
  d <- ann$NASA_POWER - ann[[p]]; pre <- ann$YEAR <= 2016
  pt <- pettitt.test(d)
  tibble(product = p, mean_diff_2001_16 = mean(d[pre]), mean_diff_2017_23 = mean(d[!pre]),
         welch_p = t.test(d[pre], d[!pre])$p.value,
         drift_pre_per_dec = 10 * coef(lm(d[pre] ~ ann$YEAR[pre]))[2],
         pettitt_year = ann$YEAR[pt$estimate[1]], pettitt_p = pt$p.value,
         r_detrended_2001_16 = cor(resid(lm(ann$NASA_POWER[pre] ~ ann$YEAR[pre])), resid(lm(ann[[p]][pre] ~ ann$YEAR[pre]))))
}))
write_csv(diffs, file.path(out, "products_difference_tests.csv"))

# ---- cell-wise Sen / MK / BH-FDR ----
cell_tr <- bind_rows(lapply(names(periods), function(pn) both |> filter(YEAR %in% periods[[pn]]) |>
  arrange(YEAR) |> group_by(product, SEASON, cell) |>
  summarise(period = pn, sen = 10 * unname(sens.slope(GHI)$estimates), p = mk.test(GHI)$p.value, .groups = "drop"))) |>
  group_by(product, SEASON, period) |> mutate(sig_fdr10 = p.adjust(p, "BH") < 0.10) |> ungroup()
write_csv(cell_tr, file.path(out, "products_cell_trends.csv"))

# ---- SPDE v2 India-mean trend (annual) for EBAF and MERRA-2 ----
if (RUN_V2 && any(c("EBAF", "MERRA2") %in% names(monthly))) {
  inla.setOption(num.threads = "2:1")
  loc <- as.matrix(cells[, c("X_km", "Y_km")]); nc <- nrow(cells)
  mesh <- fm_mesh_2d_inla(loc = loc, max.edge = c(120, 500), cutoff = 40, offset = c(150, 600))
  cat("Base mesh nodes:", mesh$n, "(should be 1383, as in the real-data analysis)\n")
  spde_b <- inla.spde2.pcmatern(mesh, prior.range = c(300, 0.05), prior.sigma = c(2, 0.05))
  spde_w <- inla.spde2.pcmatern(mesh, prior.range = c(200, 0.05), prior.sigma = c(3, 0.05))
  A_pred <- inla.spde.make.A(mesh, loc = loc)
  v2_india <- function(Y, years, seed) {
    T <- length(years); tc <- (years - mean(range(years))) / 10
    dd <- expand.grid(cell = 1:nc, yr = 1:T); ol <- loc[dd$cell, ]
    A_b <- inla.spde.make.A(mesh, loc = ol, weights = tc[dd$yr])
    A_w <- inla.spde.make.A(mesh, loc = ol, repl = dd$yr, n.repl = T)
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
    tr <- (as.matrix(A_pred %*% sapply(smp, function(z) z$latent[ib, 1])) +
             matrix(sapply(smp, function(z) z$latent[i0, 1]), nc, N_SAMP, byrow = TRUE)) * s
    natv <- as.vector(crossprod(cells$w_area, tr))
    c(india = mean(natv), lo = unname(quantile(natv, 0.025)), hi = unname(quantile(natv, 0.975)))
  }
  res <- list()
  for (p in intersect(c("EBAF", "MERRA2", "NASA_POWER"), names(monthly))) for (pn in names(periods)) {
    yy <- periods[[pn]]
    Y <- both |> filter(product == p, SEASON == "Annual", YEAR %in% yy) |> arrange(cell, YEAR) |>
      select(cell, YEAR, GHI) |> pivot_wider(names_from = YEAR, values_from = GHI) |> arrange(cell) |>
      select(-cell) |> as.matrix()
    t0 <- Sys.time(); r <- v2_india(Y, yy, seed = 700000 + nchar(p) * 10 + length(yy))
    res[[length(res) + 1]] <- tibble(product = p, period = pn, india = r["india"], lo = r["lo"], hi = r["hi"])
    cat(sprintf("v2 %s %s: %.3f (%.3f, %.3f) [%.1f min]\n", p, pn, r["india"], r["lo"], r["hi"],
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  write_csv(bind_rows(res), file.path(out, "products_v2_india.csv"))
}

options(width = 200)
cat("\n===== Annual India-mean OLS trends (kWh m-2 day-1 per decade) =====\n")
print(as.data.frame(nat_tr |> filter(SEASON == "Annual") |> select(product, period, ols, ols_lo, ols_hi, sen)), digits = 3)
cat("\n===== POWER minus other products =====\n"); print(as.data.frame(diffs), digits = 3)
cat("\nDone. Run finished.\n")
