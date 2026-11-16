"""Analysis of revision runs (script 08) - computed from saved raw CSVs.
Run: python3 analyse_rev.py <folder with results/rev files> <output folder>
Writes: B_null_metrics.csv, B_paired.csv, B_false_location.csv,
        C_signal_metrics.csv, A_acf.csv, D_ccd_metrics.csv (if present)
"""
import sys, os
import numpy as np, pandas as pd

IN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
cells = pd.read_csv(f"{IN}/cells.csv")

def bh(p, q):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    thr = q * np.arange(1, n + 1) / n
    ok = p[o] <= thr
    k = np.max(np.where(ok)[0]) + 1 if ok.any() else 0
    rej = np.zeros(n, bool); rej[o[:k]] = True
    return rej

def flags(df):
    """Return dict rule -> (neg flags, pos flags) for one rep/method (rows ordered by cell)."""
    out = {}
    if df.method.iloc[0] == "MK":
        p, s = df.p.values, df.est.values
        for name, rej in [("MK raw p<0.05", p < 0.05), ("MK BH q=0.05", bh(p, 0.05)),
                          ("MK BH q=0.10", bh(p, 0.10))]:
            out[name] = (rej & (s < 0), rej & (s > 0))
    else:
        out["marginal 95%"] = (df.q975.values < 0, df.q025.values > 0)
        out["marginal 90%"] = (df.q95.values < 0, df.q05.values > 0)
        out["joint a=0.025/sign"] = (df.jneg025.values.astype(bool), df.jpos025.values.astype(bool))
        out["joint a=0.05/sign"] = (df.jneg05.values.astype(bool), df.jpos05.values.astype(bool))
    return out

def ci(k, n, z=1.959964):
    """Wilson 95% interval for a proportion, in percent."""
    ph = k / n; den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den; h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return 100 * (c - h), 100 * (c + h)

# ---------------- B: null ----------------
b = pd.read_csv(f"{IN}/B_null_cells.csv").sort_values(["method", "rep", "cell"])
rows, anyfalse, loc_count = [], {}, {}
for (m, r), g in b.groupby(["method", "rep"]):
    for rule, (fn, fp) in flags(g).items():
        f = fn | fp
        key = (m, rule)
        anyfalse.setdefault(key, {})[r] = f.any()
        loc_count.setdefault(key, np.zeros(len(cells)))
        loc_count[key] += f
        rows.append(dict(method=m, rule=rule, rep=r, pct_flagged=100 * f.mean(), any=f.any(),
                         rmse=np.sqrt(np.mean(g.est.values ** 2))))
B = pd.DataFrame(rows)
summ = []
for (m, rule), g in B.groupby(["method", "rule"]):
    k, n = int(g["any"].sum()), len(g); lo, hi = ci(k, n)
    summ.append(dict(method=m, rule=rule, reps=n, FWER=100 * k / n, FWER_lo=lo, FWER_hi=hi,
                     mean_pct_flagged=g.pct_flagged.mean(), RMSE=g.rmse.mean()))
Bs = pd.DataFrame(summ)
Bs.round(3).to_csv(f"{OUT}/B_null_metrics.csv", index=False)

# paired comparison v1 vs v2 (same permutations)
pairs = []
for rule in ["joint a=0.025/sign", "joint a=0.05/sign", "marginal 95%"]:
    a1, a2 = anyfalse[("v1", rule)], anyfalse[("v2", rule)]
    common = sorted(set(a1) & set(a2))
    x = np.array([a1[r] for r in common]); y = np.array([a2[r] for r in common])
    n10, n01 = int((x & ~y).sum()), int((~x & y).sum())
    from scipy.stats import binomtest
    pv = binomtest(n10, n10 + n01, 0.5).pvalue if n10 + n01 > 0 else 1.0
    pairs.append(dict(rule=rule, reps=len(common), v1_only=n10, v2_only=n01, both=int((x & y).sum()),
                      diff_FWER=100 * (x.mean() - y.mean()), mcnemar_exact_p=pv))
pd.DataFrame(pairs).round(4).to_csv(f"{OUT}/B_paired.csv", index=False)

# where false flags occur (v2 joint 0.025) vs local variability / latitude
loc = cells.copy()
for key in [("v2", "joint a=0.025/sign"), ("v1", "joint a=0.05/sign"), ("MK", "MK BH q=0.10")]:
    loc[f"{key[0]}|{key[1]}"] = loc_count[key] / len(anyfalse[key])
loc.to_csv(f"{OUT}/B_false_location.csv", index=False)

# ---------------- C: signal ----------------
c = pd.read_csv(f"{IN}/C_signal_cells.csv").sort_values(["scenario", "method", "rep", "cell"])
truths = pd.read_csv(f"{IN}/C_truths.csv").sort_values("cell")
rows = []
for (sc, m, r), g in c.groupby(["scenario", "method", "rep"]):
    t = truths[sc].values
    tneg, tpos, tzero = t < 0, t > 0, t == 0
    for rule, (fn, fp) in flags(g).items():
        disc = fn | fp
        false = (fn & ~tneg) | (fp & ~tpos)            # zero-truth or wrong sign
        wrong_sign = (fn & tpos) | (fp & tneg)
        nd = disc.sum()
        d = dict(scenario=sc, method=m, rule=rule, rep=r,
                 FDP=false.sum() / nd if nd else 0.0,
                 dirFDP=wrong_sign.sum() / nd if nd else 0.0,
                 any_false=false.any(), any_wrong_sign=wrong_sign.any(),
                 power=((fn & tneg) | (fp & tpos)).sum() / max((~tzero).sum(), 1),
                 pct_flagged=100 * disc.mean(),
                 rmse=np.sqrt(np.mean((g.est.values - t) ** 2)),
                 bias=np.mean(g.est.values - t),
                 n_zero_truth=int(tzero.sum()))
        if m != "MK":
            d["coverage95"] = np.mean((g.q025.values <= t) & (t <= g.q975.values))
        rows.append(d)
C = pd.DataFrame(rows)
Cs = C.groupby(["scenario", "method", "rule"]).agg(
    reps=("rep", "nunique"), FDR=("FDP", "mean"), dirFDR=("dirFDP", "mean"),
    FWER=("any_false", "mean"), dirFWER=("any_wrong_sign", "mean"),
    power=("power", "mean"), pct_flagged=("pct_flagged", "mean"),
    RMSE=("rmse", "mean"), bias=("bias", "mean"), coverage95=("coverage95", "mean"),
    n_zero_truth=("n_zero_truth", "first")).reset_index()
for col in ["FDR", "dirFDR", "FWER", "dirFWER", "power", "coverage95"]:
    Cs[col] = 100 * Cs[col]
Cs.round(3).to_csv(f"{OUT}/C_signal_metrics.csv", index=False)

# ---------------- A: autocorrelation ----------------
a = pd.read_csv(f"{IN}/A_classical.csv")
acf = a.groupby("SEASON")[["acf1", "acf2", "acf3"]].describe(percentiles=[0.05, 0.5, 0.95])
acf.round(3).to_csv(f"{OUT}/A_acf.csv")

# ---------------- D: CCD null (if present) ----------------
if os.path.exists(f"{IN}/D_null_ccd_cells.csv"):
    d = pd.read_csv(f"{IN}/D_null_ccd_cells.csv").sort_values(["rep", "cell"])
    rows = []
    for r, g in d.groupby("rep"):
        for rule, (fn, fp) in flags(g).items():
            rows.append(dict(rule=rule, rep=r, any=(fn | fp).any()))
    D = pd.DataFrame(rows)
    ds = []
    for rule, g in D.groupby("rule"):
        k, n = int(g["any"].sum()), len(g); lo, hi = ci(k, n)
        # same reps under EB (task B) for a paired view
        eb = B[(B.method == "v2") & (B.rule == rule) & (B.rep.isin(g.rep))]
        ds.append(dict(rule=rule, reps=n, FWER_ccd=100 * k / n, lo=lo, hi=hi,
                       FWER_eb_same_reps=100 * eb["any"].mean()))
    pd.DataFrame(ds).round(2).to_csv(f"{OUT}/D_ccd_metrics.csv", index=False)
print("done")
