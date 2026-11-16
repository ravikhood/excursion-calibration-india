# =============================================================
# Step 7: Final figures for the manuscript (PDF + PNG)
# Output: figures/paper/
#   fig5_v2_sd.pdf     posterior SD of the v2 trend (5 seasons)
#   fig6_simulation.pdf simulation: FWER under null, power & RMSE
# (Figures 1-4 are copied from earlier steps; see list at the end.)
# =============================================================

library(dplyr); library(tidyr); library(readr); library(ggplot2)
library(sf); library(rnaturalearth)

out <- file.path("figures", "paper"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
season_levels <- c("Annual", "Winter (JF)", "Pre-monsoon (MAM)",
                   "Monsoon (JJAS)", "Post-monsoon (OND)")
india <- ne_countries(country = "India", scale = "medium", returnclass = "sf")

# ---- Fig 5: posterior SD of v2 trend ----
v2 <- read_csv("results/2001_2023/v2_cell_trends.csv", show_col_types = FALSE) |>
  mutate(SEASON = factor(SEASON, levels = season_levels))
p5 <- ggplot(v2) +
  geom_tile(aes(LON, LAT, fill = post_sd)) +
  geom_sf(data = india, fill = NA, colour = "grey30", linewidth = 0.3) +
  scale_fill_viridis_c(name = expression(paste("Posterior SD (kWh ", m^-2, " ", day^-1, " ", decade^-1, ")"))) +
  facet_wrap(~ SEASON, ncol = 3) + coord_sf(expand = FALSE) +
  labs(x = NULL, y = NULL) + theme_bw(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out, "fig5_v2_sd.pdf"), p5, width = 9, height = 7.5)
ggsave(file.path(out, "fig5_v2_sd.png"), p5, width = 9, height = 7.5, dpi = 300)

# ---- Fig 6: simulation summary ----
lab <- c(MK_raw = "MK (raw)", MK_fdr = "MK + BH-FDR",
         Bayes_marg = "SPDE v1 (marginal)", Bayes_joint = "SPDE v1 (joint)",
         Bayes2_marg = "SPDE v2 (marginal)", Bayes2_joint = "SPDE v2 (joint)")
sim <- read_csv("results/sim/sim_summary_v2.csv", show_col_types = FALSE) |>
  mutate(method = factor(lab[method], levels = lab))
long <- bind_rows(
  sim |> filter(scenario == "null")   |> transmute(method, panel = "(a) Null: FWER (%)", value = FWER_any_false),
  sim |> filter(scenario == "signal") |> transmute(method, panel = "(b) Signal: power (%)", value = power_neg),
  sim |> filter(scenario == "signal") |> transmute(method, panel = "(c) Signal: wrong-sign FWER (%)", value = FWER_any_false),
  sim |> transmute(method, panel = ifelse(scenario == "null", "(d) Null: RMSE", "(e) Signal: RMSE"), value = RMSE))
p6 <- ggplot(long, aes(method, value)) +
  geom_col(fill = "grey45", width = 0.7) +
  geom_hline(data = data.frame(panel = "(a) Null: FWER (%)", y = 5), aes(yintercept = y),
             linetype = 2, colour = "red") +
  facet_wrap(~ panel, scales = "free_x", ncol = 3) + coord_flip() +
  labs(x = NULL, y = NULL) + theme_bw(base_size = 9)
ggsave(file.path(out, "fig6_simulation.pdf"), p6, width = 9, height = 5)
ggsave(file.path(out, "fig6_simulation.png"), p6, width = 9, height = 5, dpi = 300)

# ---- Copy earlier figures under paper names ----
file.copy("figures/2001_2023/fig_mean_ghi_annual.png",     file.path(out, "fig1_mean_ghi.png"), overwrite = TRUE)
file.copy("figures/fig_breakpoint.png",                    file.path(out, "fig2_breakpoint.png"), overwrite = TRUE)
file.copy("figures/2001_2023/fig_classical_sen_slope.png", file.path(out, "fig3_sen_fdr.png"), overwrite = TRUE)
file.copy("figures/2001_2023/fig_v2_trend.png",            file.path(out, "fig4_v2_trend.png"), overwrite = TRUE)
file.copy("figures/2001_2023/fig_mesh.png",                file.path(out, "figS1_mesh.png"), overwrite = TRUE)
file.copy("figures/2001_2023/fig_bayes_trend.png",         file.path(out, "figS2_v1_trend.png"), overwrite = TRUE)
cat("Paper figures written to", out, "\n")
