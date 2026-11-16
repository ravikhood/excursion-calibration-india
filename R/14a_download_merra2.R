# =============================================================
# Step 14a: download the MERRA-2 SWGDN monthly subsets listed in a GES DISC links file.
#
# Before running:
#  1. In your Earthdata profile (urs.earthdata.nasa.gov -> Applications -> Authorized Apps)
#     approve "NASA GESDISC DATA ARCHIVE" (one-time; without it every download returns a login page).
#  2. On https://disc.gsfc.nasa.gov/datasets/M2TMNXRAD_5.12.4/summary click "Subset / Get Data":
#       Download method: Get File Subsets (OPeNDAP or subsetter, either is fine)
#       Dates 2001-01-01 to 2023-12-31; Region 67, 6, 98, 38 (W, S, E, N)
#       Variable: SWGDN only; Format: netCDF
#     Click "Get Data", then "Download links list". Save the .txt file into  data/merra2/
#  3. Run:  source("14a_download_merra2.R")
#     A box asks for your Earthdata username and password (typed into RStudio, never stored,
#     never printed; do not paste them into the chat).
# Safe to rerun: files already downloaded and valid are skipped.
# =============================================================

if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
library(httr)
dir.create("data/merra2", recursive = TRUE, showWarnings = FALSE)
lf <- list.files("data/merra2", pattern = "\\.txt$", full.names = TRUE)
if (!length(lf)) stop("No links .txt file in data/merra2. See step 2 at the top of this script.")
urls <- unique(trimws(unlist(lapply(lf, readLines, warn = FALSE))))
urls <- urls[grepl("^https?://", urls) & grepl("tavgM_2d_rad|M2TMNXRAD", urls)]   # drop README/doc links
cat("Data links found:", length(urls), "(expect 276 = 23 years x 12 months)\n")

fname <- function(u) {
  lab <- regmatches(u, regexpr("LABEL=[^&]+", u))
  f <- if (length(lab)) sub("LABEL=", "", lab) else basename(sub("\\?.*$", "", u))
  f <- URLdecode(f); if (!grepl("\\.nc4?$", f)) f <- paste0(f, ".nc4"); file.path("data/merra2", f)
}
is_netcdf <- function(f) {                      # HDF5 (netCDF-4) or classic netCDF signature
  if (!file.exists(f) || file.size(f) < 1000) return(FALSE)
  b <- readBin(f, "raw", 4); identical(b[2:4], charToRaw("HDF")) || identical(b[1:3], charToRaw("CDF"))
}

user <- rstudioapi::showPrompt("Earthdata", "Earthdata username:")
pass <- rstudioapi::askForPassword("Earthdata password:")
cj <- tempfile("edl_cookies")
ok <- 0; bad <- character()
for (i in seq_along(urls)) {
  f <- fname(urls[i]); if (is_netcdf(f)) { ok <- ok + 1; next }
  tmp <- paste0(f, ".part")
  r <- tryCatch(GET(urls[i], authenticate(user, pass), write_disk(tmp, overwrite = TRUE), timeout(300),
                    config(followlocation = 1L, unrestricted_auth = 1L, cookiefile = cj, cookiejar = cj)),
                error = function(e) NULL)
  if (!is.null(r) && status_code(r) == 200 && is_netcdf(tmp)) { file.rename(tmp, f); ok <- ok + 1 }
  else { unlink(tmp); bad <- c(bad, basename(f)) }
  if (i %% 12 == 0 || i == length(urls)) cat(sprintf("%d / %d checked | ok %d | failed %d\n", i, length(urls), ok, length(bad)))
  if (i == 3 && ok == 0) { rm(pass); stop("First 3 downloads failed: check the username/password and step 1 (Authorized Apps).") }
}
rm(pass, user)
cat("\nValid MERRA-2 files in data/merra2:", sum(sapply(list.files("data/merra2", "\\.nc4?$", full.names = TRUE), is_netcdf)), "\n")
if (length(bad)) cat("Failed (rerun the script to retry):", head(bad, 10), if (length(bad) > 10) "...", "\n") else
  cat("All downloaded. Next: put the EBAF .nc file in data/ebaf/, then source(\"14_radiation_products.R\")\n")
