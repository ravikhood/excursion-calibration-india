# Final summary figures from saved CSVs (no model fitting).
#   fig10_parametric      : error rate by parametric scenario (script 12)
#   fig11_alpha_validation: FWER vs alpha, selection (reps 1-50) vs validation (reps 51-150) (scripts 11/13)
#   fig12_products        : India-mean annual GHI for four products + POWER-minus-product differences (script 14)
# Usage: python3 make_final_figures.py <results_root> <analysis_out_dir> <fig_dir>
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = sys.argv[1] if len(sys.argv) > 1 else "results"
AOUT = sys.argv[2] if len(sys.argv) > 2 else "tables"
FIG = sys.argv[3] if len(sys.argv) > 3 else "figures"
import os; os.makedirs(FIG, exist_ok=True)
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
S = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # validated categorical order (light)
MK = ["o", "s", "^", "D"]
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "savefig.dpi": 300})
def save(f, name):
    for ext in ("png", "pdf"): f.savefig(f"{FIG}/{name}.{ext}", bbox_inches="tight")
    plt.close(f)
def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)

# ---- fig10: parametric scenarios (from 12_cells.csv, strict F > 1 - alpha) ----
C = pd.read_csv(f"{ROOT}/param/12_cells.csv")
order = [("P1b", "P1: trend from fitted\nGRF (prior-averaged)"), ("P0b", "P0: zero trend,\ns known"),
         ("P0a", "P0: zero trend,\ns estimated (v2)"), ("P2a", "P2: heavy-tailed\nyears (v2)"),
         ("P2b", "P2: AR(1)=0.3\nyears (v2)")]
f, ax = plt.subplots(figsize=(7.2, 3.4))
for j, (a, lab) in enumerate([(0.025, "nominal 5% (alpha 0.025/sign)"), (0.05, "nominal 10% (alpha 0.05/sign)")]):
    for i, (sc, _) in enumerate(order):
        g = C[C.scenario == sc]
        r = (((g.Fneg > 1 - a) & (g.truth >= 0)) | ((g.Fpos > 1 - a) & (g.truth <= 0))).groupby(g.rep).any()
        k, n = int(r.sum()), len(r); lo, hi = wilson(k, n); x = i + (j - 0.5) * 0.28
        ax.errorbar(x, 100 * k / n, yerr=[[100 * k / n - lo], [hi - 100 * k / n]], fmt=MK[j], color=S[j], ms=6,
                    capsize=0, elinewidth=1.4, label=lab if i == 0 else None)
    ax.axhline(100 * 2 * a, color=S[j], lw=1, ls="--")
ax.set_xticks(range(len(order))); ax.set_xticklabels([o[1] for o in order])
ax.set_ylabel("Replicates with any wrong flag (%)"); ax.set_ylim(0, 103)
ax.legend(frameon=False, loc="upper left"); ax.text(4.45, 6, "dashed: nominal", color=INK2, ha="right", fontsize=8)
save(f, "fig10_parametric")

# ---- fig11: alpha validation curve ----
T = pd.read_csv(f"{AOUT}/val_A_fwer_curve.csv")
f, ax = plt.subplots(figsize=(5.4, 3.6))
sets = [("select_reps1-50_batchA", "Selection: reps 1-50 (batch A)"), ("select_reps1-50_batchB", "Selection: reps 1-50 (batch B)"),
        ("validate_reps51-150", "Validation: reps 51-150")]
for i, (s, lab) in enumerate(sets):
    g = T[T.set == s].sort_values("alpha_per_sign")
    ax.plot(g.alpha_per_sign, g.FWER, marker=MK[i], color=S[i], lw=2, ms=5, label=lab)
    if s.startswith("validate"):
        ax.fill_between(g.alpha_per_sign, g.lo, g.hi, color=S[i], alpha=0.12, lw=0)
ax.axhline(5, color=INK2, lw=1, ls="--"); ax.text(0.0505, 5.8, "5%", color=INK2, ha="right", fontsize=8)
ax.axvline(0.005, color=INK2, lw=1, ls=":"); ax.text(0.0056, 33.5, "alpha* = 0.005", color=INK2, fontsize=8)
ax.axvline(0.025, color=INK2, lw=1, ls=":"); ax.text(0.0258, 33.5, "standard 0.025", color=INK2, fontsize=8)
ax.set_xlabel("Per-sign alpha"); ax.set_ylabel("Null FWER (%)"); ax.set_ylim(0, 36)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.9), fontsize=8)
save(f, "fig11_alpha_validation")

# ---- fig12: four products ----
N = pd.read_csv(f"{ROOT}/products/products_national_series.csv")
A = N[N.SEASON == "Annual"].pivot(index="YEAR", columns="product", values="GHI")
prods = [("NASA_POWER", "NASA POWER"), ("EBAF", "EBAF Ed4.2.1"), ("ERA5", "ERA5"), ("MERRA2", "MERRA-2")]
f, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4), gridspec_kw={"width_ratios": [1.25, 1], "wspace": 0.38})
for i, (p, lab) in enumerate(prods):
    a1.plot(A.index, A[p], marker=MK[i], color=S[i], lw=2, ms=4)
    a1.annotate(lab, (A.index[-1], A[p].iloc[-1]), xytext=(4, {"NASA_POWER": 5, "EBAF": -5}.get(p, 0)), textcoords="offset points", va="center",
                fontsize=8, color=INK)
a1.set_xlim(2000.5, 2027.5); a1.set_ylabel("India-mean annual GHI (kWh m$^{-2}$ day$^{-1}$)")
a1.set_title("(a) Four products", loc="left", fontsize=9, color=INK)
for i, (p, lab) in enumerate(prods[1:], start=1):
    d = A["NASA_POWER"] - A[p]; d = d - d[d.index <= 2016].mean()
    a2.plot(A.index, d, marker=MK[i], color=S[i], lw=2, ms=4, label=f"POWER minus {lab}")
a2.axvline(2016.5, color=INK2, lw=1, ls=":"); a2.axhline(0, color=INK2, lw=0.8)
a2.set_ylabel("Difference minus its 2001-2016 mean"); a2.set_title("(b) POWER minus other products", loc="left", fontsize=9, color=INK)
a2.legend(frameon=False, fontsize=7.5, loc="lower left")
save(f, "fig12_products")
print("written:", FIG)
