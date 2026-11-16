"""Revision figures that need no map boundary (Python, run from the saved CSV files).
Usage: python3 make_rev_figures.py <results dir (contains rev/ and era5/)> <analysis out dir> <figure out dir>
Writes fig6_null_calibration, fig7_signal, fig8_era5, figS3_truths (.pdf + .png)."""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES, AOUT, FIG = sys.argv[1:4]; os.makedirs(FIG, exist_ok=True)
BLUE, ORANGE, AQUA, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#8a8984", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6, "axes.axisbelow": True})
def save(fig, name):
    for ext in ("pdf", "png"): fig.savefig(f"{FIG}/{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
def wilson(k, n, z=1.959964):
    p = k / n; d = 1 + z*z/n; c = (p + z*z/(2*n))/d; h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return 100*(c-h), 100*(c+h)

# ---- Fig 6: null calibration ----
B = pd.read_csv(f"{AOUT}/B_null_metrics.csv"); D = pd.read_csv(f"{AOUT}/D_ccd_metrics.csv")
rows = [("MK, BH-FDR q = 0.10", "MK", "MK BH q=0.10", BLUE),
        ("MK, BH-FDR q = 0.05", "MK", "MK BH q=0.05", BLUE),
        ("v1, joint (α = 0.025 per sign)", "v1", "joint a=0.025/sign", GREY),
        ("v2, joint (α = 0.025 per sign)", "v2", "joint a=0.025/sign", ORANGE),
        ("v1, marginal 95%", "v1", "marginal 95%", GREY),
        ("v2, marginal 95%", "v2", "marginal 95%", ORANGE),
        ("v1, joint (α = 0.05 per sign)", "v1", "joint a=0.05/sign", GREY),
        ("v2, joint (α = 0.05 per sign)", "v2", "joint a=0.05/sign", ORANGE)]
fig, ax = plt.subplots(figsize=(6.3, 3.6))
for i, (lab, m, r, col) in enumerate(rows):
    x = B[(B.method == m) & (B.rule == r)].iloc[0]
    ax.plot([x.FWER_lo, x.FWER_hi], [i, i], color=col, lw=2, solid_capstyle="round")
    ax.plot(x.FWER, i, "o", ms=7, color=col, mec="white", mew=1.5)
    ax.text(max(x.FWER_hi, 6) + 1, i, f"{x.FWER:.1f}", va="center", fontsize=8, color=INK2)
d = D[D.rule == "joint a=0.025/sign"].iloc[0]
i = len(rows)
ax.plot([d.lo, d.hi], [i, i], color=ORANGE, lw=2, ls=(0, (2, 1.5)))
ax.plot(d.FWER_ccd, i, "s", ms=7, color=ORANGE, mec="white", mew=1.5)
ax.text(d.hi + 1, i, f"{d.FWER_ccd:.1f} (50 reps)", va="center", fontsize=8, color=INK2)
labels = [r[0] for r in rows] + ["v2 + CCD integration, joint (α = 0.025)"]
ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels); ax.set_ylim(len(labels) + 0.1, -0.6)
ax.axvline(5, color=INK, lw=1, ls="--"); ax.text(5.5, len(rows) + 0.45, "nominal 5%", fontsize=8, color=INK, va="center")
raw = B[(B.method == "MK") & (B.rule == "MK raw p<0.05")].iloc[0]
ax.set_xlim(0, 48); ax.set_xlabel("Family-wise error rate under the null (%) with 95% Wilson interval")
ax.set_title(f"Raw cell-wise MK (p < 0.05): {raw.FWER:.0f}% (off scale)", fontsize=8, color=INK2, loc="left")
ax.grid(axis="y", visible=False)
save(fig, "fig6_null_calibration")

# ---- Fig 7: synthetic truths ----
C = pd.read_csv(f"{AOUT}/C_signal_metrics.csv")
sc = ["S1_smooth", "S2_mixed", "S3_sharp", "S4_sparse", "S5_hetero"]
scl = ["S1", "S2", "S3", "S4", "S5"]
mk = C[(C.method == "MK") & (C.rule == "MK BH q=0.10")].set_index("scenario").loc[sc]
v2 = C[(C.method == "v2") & (C.rule == "joint a=0.025/sign")].set_index("scenario").loc[sc]
panels = [("FDR", "(a) False discovery rate (%)"), ("FWER", "(b) P(any false flag) (%)"),
          ("dirFWER", "(c) P(any wrong sign) (%)"), ("power", "(d) Power (%)"), ("RMSE", "(e) RMSE of trend")]
fig, axs = plt.subplots(1, 5, figsize=(10.5, 2.6))
xx = np.arange(len(sc))
for ax, (col, title) in zip(axs, panels):
    ax.plot(xx - 0.12, mk[col], "o", ms=7, color=BLUE, mec="white", mew=1.5, label="MK + BH-FDR (q = 0.10)")
    ax.plot(xx + 0.12, v2[col], "s", ms=7, color=ORANGE, mec="white", mew=1.5, label="v2 joint (α = 0.025/sign)")
    ax.set_xticks(xx); ax.set_xticklabels(scl, fontsize=7.5); ax.set_title(title, fontsize=8.5, loc="left")
    ax.grid(axis="x", visible=False)
    if col in ("FDR", "FWER", "dirFWER"): ax.set_ylim(bottom=0)
h, l = axs[0].get_legend_handles_labels()
fig.subplots_adjust(left=0.05, right=0.99, top=0.88, bottom=0.25, wspace=0.38)
fig.legend(h, l, loc="lower center", ncol=2, frameon=False, fontsize=8.5)
save(fig, "fig7_signal")

# ---- Fig 8: ERA5 vs NASA ----
s = pd.read_csv(f"{RES}/era5/national_series.csv")
a = s[s.SEASON == "Annual"].pivot(index="YEAR", columns="dataset", values="GHI")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.1))
ax1.plot(a.index, a.NASA_POWER, "-o", color=BLUE, lw=2, ms=4, label="NASA POWER (CERES SYN1deg)")
ax1.plot(a.index, a.ERA5, "-s", color=ORANGE, lw=2, ms=4, label="ERA5")
ax1.text(2023.4, a.NASA_POWER.iloc[-1], "NASA", color=INK2, va="center", fontsize=8)
ax1.text(2023.4, a.ERA5.iloc[-1], "ERA5", color=INK2, va="center", fontsize=8)
ax1.set_ylabel("kWh m$^{-2}$ day$^{-1}$"); ax1.set_title("(a) India-mean annual GHI", loc="left", fontsize=9)
ax1.legend(frameon=False, fontsize=8, loc="lower left"); ax1.set_xlim(2000, 2026)
dd = a.NASA_POWER - a.ERA5
ax2.bar(dd.index, dd.values, width=0.7, color=GREY)
for (y0, y1) in [(2001, 2016), (2017, 2023)]:
    m = dd.loc[y0:y1].mean(); ax2.plot([y0 - 0.4, y1 + 0.4], [m, m], color=INK, lw=1.5)
    ax2.text(y0 - 0.4, -0.235, f"{y0}–{str(y1)[2:]} mean {m:.3f}", ha="left", va="bottom", fontsize=8, color=INK)
ax2.set_ylim(-0.245, 0.005); ax2.axvline(2016.5, color=INK2, lw=1, ls="--")
ax2.set_title("(b) Difference NASA − ERA5", loc="left", fontsize=9); ax2.set_ylabel("kWh m$^{-2}$ day$^{-1}$")
ax2.grid(axis="x", visible=False)
fig.tight_layout(); save(fig, "fig8_era5")

# ---- Fig S3: synthetic truth fields (tiles on the analysis grid) ----
t = pd.read_csv(f"{RES}/rev/C_truths.csv").merge(pd.read_csv(f"{RES}/rev/cells.csv"), on="cell")
fig, axs = plt.subplots(1, 5, figsize=(11, 2.9))
cmap = plt.get_cmap("RdBu_r"); vmax = 0.3
for ax, c, l in zip(axs, sc, ["S1 smooth", "S2 mixed", "S3 sharp", "S4 sparse", "S5 hetero"]):
    im = ax.scatter(t.LON, t.LAT, c=t[c], cmap=cmap, vmin=-vmax, vmax=vmax, marker="s", s=17, linewidths=0)
    ax.set_aspect("equal"); ax.set_title(l, fontsize=9); ax.set_xticks([70, 80, 90]); ax.set_yticks([10, 20, 30])
    ax.grid(False); ax.tick_params(labelsize=7)
cb = fig.colorbar(im, ax=axs, orientation="horizontal", fraction=0.06, pad=0.12, shrink=0.5)
cb.set_label("True trend (kWh m$^{-2}$ day$^{-1}$ decade$^{-1}$)")
save(fig, "figS3_truths")
print("figures written to", FIG)
