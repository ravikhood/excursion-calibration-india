# =============================================================
# Step 2b: Is there a jump at the 2000/2001 data-source change?
#   NASA POWER solar = GEWEX SRB (1984-2000) + CERES SYN1deg (2001-)
#   1) Pettitt change-point test on each cell's annual series
#   2) India-mean annual GHI time series with the switch marked
# Input : data/power_ghi_india_monthly.csv
# Output: results/breakpoint_check.csv, figures/fig_breakpoint.png
# =============================================================

library(dplyr); library(readr); library(ggplot2); library(trend)
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)

ann <- read_csv("data/power_ghi_india_monthly.csv", show_col_types = FALSE) |>
  group_by(LON, LAT, YEAR) |> summarise(GHI = mean(GHI), .groups = "drop") |>
  arrange(LON, LAT, YEAR)

# ---- 1. Pettitt test per cell ----
pt <- ann |> group_by(LON, LAT) |>
  summarise(cp_year = YEAR[pettitt.test(GHI)$estimate[1]],
            p_pett  = pettitt.test(GHI)$p.value,
            jump    = mean(GHI[YEAR >= 2001 & YEAR <= 2005]) -
                      mean(GHI[YEAR >= 1996 & YEAR <= 2000]),
            .groups = "drop")
write_csv(pt, "results/breakpoint_check.csv")

cat("Cells with significant change point (p<0.05):",
    round(100 * mean(pt$p_pett < 0.05), 1), "%\n")
cat("Of those, change point in 1999-2002:",
    round(100 * mean(pt$cp_year[pt$p_pett < 0.05] %in% 1999:2002), 1), "%\n")
cat("Most common change-point years:\n")
print(head(sort(table(pt$cp_year[pt$p_pett < 0.05]), decreasing = TRUE), 5))
cat("Median jump (mean 2001-05 minus mean 1996-2000):",
    round(median(pt$jump), 3), "kWh/m2/day\n")

# ---- 2. India-mean series ----
nat <- ann |> group_by(YEAR) |> summarise(GHI = mean(GHI))
p <- ggplot(nat, aes(YEAR, GHI)) +
  geom_line() + geom_point(size = 1) +
  geom_vline(xintercept = 2000.5, linetype = 2, colour = "red") +
  annotate("text", x = 1992, y = max(nat$GHI), label = "GEWEX SRB", size = 3) +
  annotate("text", x = 2012, y = max(nat$GHI), label = "CERES SYN1deg", size = 3) +
  labs(x = NULL, y = expression(paste("India-mean annual GHI (kWh ", m^-2, " ", day^-1, ")"))) +
  theme_bw(base_size = 10)
ggsave("figures/fig_breakpoint.png", p, width = 7, height = 3.5, dpi = 300)
