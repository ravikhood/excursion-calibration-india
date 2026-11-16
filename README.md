# Frequentist calibration of Bayesian joint excursion sets for gridded trend detection

Code and result tables for the paper

> Ravi Kumar and Sumanta Pasari. *Frequentist calibration of Bayesian joint excursion sets
> for gridded trend detection: evidence from surface solar radiation over India.*

Department of Mathematics, BITS Pilani, Pilani Campus, India.
Contact: Ravi Kumar (p20250014@pilani.bits-pilani.ac.in)

## What is here

| Folder | Content |
|---|---|
| `R/` | R scripts, numbered in the order they were run |
| `analysis/` | Python scripts that compute the error rates, tests and figures from the saved outputs |
| `tables/` | Result tables behind the numbers in the paper (CSV) |

Input data and large simulation outputs are not stored here. The scripts download the
input data and recreate the outputs (see below).

## Scripts

| Script | Purpose | Paper section |
|---|---|---|
| `01_download_power.R` | Download NASA POWER monthly GHI; keep the 274 cells over India | 2.1 |
| `02_classical_trends.R`, `02b_breakpoint_check.R` | Sen slope, Mann-Kendall, BH adjustment; Pettitt check | 3.1 |
| `03`-`07` | Earlier model versions and figures (kept for the record; not used for the final numbers) | - |
| `08_revision_all.R` | Real-data fits and robustness; 200 permutation-null replicates; synthetic signals (40 replicates); CCD check; scale-perturbation sensitivity | 4.1, 4.2, 4.6, 4.7 |
| `09_era5_validation.R` | Download ERA5 and compare with NASA POWER | 4.8 |
| `10_paper_maps_rev.R` | Map figures from script 08 outputs | - |
| `11_mc_convergence.R` | Monte Carlo convergence of excursion sets (50 replicates) | 4.3 |
| `12_parametric_calibration.R` | Parametric scenarios P0a, P0b, P1b, P2a, P2b | 4.4 |
| `13_validation_and_signal.R` | Held-out null replicates 51-150; S2/S4 signals, replicates 41-100 | 4.5, 4.6 |
| `14a_download_merra2.R`, `14_radiation_products.R` | CERES EBAF, MERRA-2, ERA5 and POWER comparison | 4.8 |
| `15_tuned_alpha_real_data.R` | Real data at the calibrated threshold | 4.7 |
| `16_redraw_maps_official_boundary.R` | Final map figures with the boundary of India (Natural Earth, India point of view) | Figures |

Random seeds are fixed inside each script (permutation null: 5000 + r; synthetic signals: 7000 + r).
Long scripts save one file per replicate and continue from where they stopped if run again.

## Software

R 4.5.1 with INLA 25.10.19, fmesher 0.8.0, excursions 2.5.11, trend, sf, dplyr, tidyr, readr, ggplot2, ncdf4.
Python 3 with numpy, pandas, scipy, statsmodels, matplotlib and rdata.

INLA is installed from its own repository:
`install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)`

## Data sources

| Product | Source | Access |
|---|---|---|
| NASA POWER GHI | NASA POWER Project | downloaded by `01_download_power.R` |
| ERA5 monthly means, surface solar radiation downwards | Copernicus Climate Data Store, DOI 10.24381/cds.f17050d7 | `09_era5_validation.R` (free CDS account and API key) |
| CERES EBAF Edition 4.2.1, `sfc_sw_down_all_mon` | NASA LaRC ASDC, DOI 10.5067/TERRA-AQUA-NOAA20/CERES/EBAF_L3B004.2.1 | manual order from ceres.larc.nasa.gov, region 6-38 N, 67-98 E, 2001-2023; save in `data/ebaf/` |
| MERRA-2 M2TMNXRAD v5.12.4, `SWGDN` | NASA GES DISC, DOI 10.5067/OU3HJDS973O0 | subset links list saved in `data/merra2/`, then `14a_download_merra2.R` (free Earthdata account) |
| Boundary of India | Natural Earth 1:10m Admin 0, India point of view (public domain) | downloaded by `16_redraw_maps_official_boundary.R` |

Never put passwords or API keys in a script; the scripts ask for them when they run.

## How to reproduce

1. Set the working directory to the repository folder and run the R scripts in number order.
   Scripts 08 and 11-13 are long (several hours each on a laptop).
2. Run the Python scripts from the repository folder, for example
   `python analysis/analyse_val_A.py results tables`.
3. The tables in `tables/` are the outputs we obtained; rerunning should reproduce them.

## Use of AI tools

The analysis code was developed and checked with the help of Claude (Anthropic), and ChatGPT (OpenAI)
assisted in drafting the manuscript. The authors reviewed the code, the results and the text and are
responsible for them.

## Licence

Code: MIT licence (see `LICENSE`). The input datasets remain under the terms of their providers.
