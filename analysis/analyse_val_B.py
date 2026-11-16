# Script 13 Part B + task C: MK/BH-FDR vs v2 joint excursion sets on S2_mixed and S4_sparse, 100 paired reps.
# Convention: E_alpha = {F > 1-alpha} (strict), matching excursions.mc sets exactly (ties at F = 1-alpha occur with 4000 samples).
# reps 1-40 from task C (results/rev/C_signal_cells.csv), reps 41-100 from script 13 (results/val/13_signal_cells.csv).
# Tuned alpha* = 0.005/sign (chosen on null reps 1-50, validated on 51-150, see PROJECT_LOG 4k) needs F -> reps 41-100 only.
import sys, numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon
ROOT = sys.argv[1] if len(sys.argv) > 1 else "results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "tables"
SC = ["S2_mixed", "S4_sparse"]
tr = pd.read_csv(f"{ROOT}/rev/C_truths.csv").sort_values("cell")
def bh(p, q):
    p = np.asarray(p); n = len(p); o = np.argsort(p); ok = p[o] <= q*np.arange(1, n+1)/n
    k = np.max(np.where(ok)[0]) + 1 if ok.any() else 0; r = np.zeros(n, bool); r[o[:k]] = True; return r
def wilson(k, n, z=1.96):
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d; h = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return 100*(c-h), 100*(c+h)
old = pd.read_csv(f"{ROOT}/rev/C_signal_cells.csv"); old = old[old.scenario.isin(SC)]
new = pd.read_csv(f"{ROOT}/val/13_signal_cells.csv")
recs = []
def add(sc, r, t, rules):
    tneg, tpos, tz = t < 0, t > 0, t == 0
    for rule, (fn, fp, est) in rules.items():
        disc = fn | fp; false = (fn & ~tneg) | (fp & ~tpos); ws = (fn & tpos) | (fp & tneg); nd = disc.sum()
        recs.append(dict(scenario=sc, rep=r, rule=rule, any_false=false.any(), any_wrong_sign=ws.any(),
                         n_false=false.sum(), FDP=false.sum()/nd if nd else 0.0,
                         power=((fn & tneg) | (fp & tpos)).sum()/max((~tz).sum(), 1),
                         pct_flagged=100*disc.mean(), rmse=np.sqrt(np.mean((est-t)**2))))
def mk_rules(p, s, est):
    return {f"MK BH q={q:.2f}": (bh(p, q) & (s < 0), bh(p, q) & (s > 0), est) for q in (0.05, 0.10)} | \
           {"MK raw p<0.05": ((p < .05) & (s < 0), (p < .05) & (s > 0), est)}
for (sc, r), g in old.groupby(["scenario", "rep"]):
    t = tr[sc].values; m = g[g.method == "MK"].sort_values("cell"); v = g[g.method == "v2"].sort_values("cell")
    add(sc, r, t, mk_rules(m.p.values, m.est.values, m.est.values) | {
        "v2 joint 0.025/sign": (v.jneg025.values.astype(bool), v.jpos025.values.astype(bool), v.est.values),
        "v2 joint 0.05/sign": (v.jneg05.values.astype(bool), v.jpos05.values.astype(bool), v.est.values),
        "v2 marginal 95%": (v.q975.values < 0, v.q025.values > 0, v.est.values)})
for (sc, r), g in new.groupby(["scenario", "rep"]):
    g = g.sort_values("cell"); t = tr[sc].values; assert np.allclose(g.truth.values, t)
    add(sc, r, t, mk_rules(g.mk_p.values, g.mk_est.values, g.mk_est.values) | {
        "v2 joint 0.025/sign": (g.Fneg.values > .975, g.Fpos.values > .975, g.est.values),
        "v2 joint 0.05/sign": (g.Fneg.values > .95, g.Fpos.values > .95, g.est.values),
        "v2 joint 0.005/sign (tuned)": (g.Fneg.values > .995, g.Fpos.values > .995, g.est.values),
        "v2 marginal 95%": (g.q975.values < 0, g.q025.values > 0, g.est.values)})
    assert ((g.Fneg > .975) == g.jneg025).all() and ((g.Fpos > .975) == g.jpos025).all()
R = pd.DataFrame(recs); R["block"] = np.where(R.rep <= 40, "1-40", "41-100")
R.to_csv(f"{OUT}/val_B_per_rep.csv", index=False)
def summ(D):
    rows = []
    for (sc, rule), g in D.groupby(["scenario", "rule"]):
        k, n = int(g.any_false.sum()), len(g); lo, hi = wilson(k, n)
        rows.append(dict(scenario=sc, rule=rule, reps=n, FWER=100*k/n, FWER_lo=lo, FWER_hi=hi,
                         dirFWER=100*g.any_wrong_sign.mean(), FDR=100*g.FDP.mean(), mean_false_cells=g.n_false.mean(),
                         power=100*g.power.mean(), pct_flagged=g.pct_flagged.mean(), RMSE=g.rmse.mean()))
    return pd.DataFrame(rows)
S = summ(R); S.round(3).to_csv(f"{OUT}/val_B_metrics_100.csv", index=False)
S60 = summ(R[R.rep > 40]); S60.round(3).to_csv(f"{OUT}/val_B_metrics_41_100.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20)
print("=== all reps (tuned rule: reps 41-100 only) ==="); print(S.round(2).to_string(index=False))
print("=== block consistency: FWER / power by block ===")
print(R.groupby(["scenario", "rule", "block"]).agg(FWER=("any_false", "mean"), power=("power", "mean")).unstack().mul(100).round(1))
# paired tests
pr = []
def paired(sc, a, b, D):
    A = D[(D.scenario == sc) & (D.rule == a)].set_index("rep"); B = D[(D.scenario == sc) & (D.rule == b)].set_index("rep")
    c = A.index.intersection(B.index); x, y = A.loc[c, "any_false"].values, B.loc[c, "any_false"].values
    n10, n01 = int((x & ~y).sum()), int((~x & y).sum())
    p_m = binomtest(n10, n10+n01, .5).pvalue if n10+n01 else 1.0
    dp = (A.loc[c, "power"] - B.loc[c, "power"]).values
    p_w = wilcoxon(dp).pvalue if np.any(dp != 0) else 1.0
    pr.append(dict(scenario=sc, A=a, B=b, reps=len(c), FWER_A=100*x.mean(), FWER_B=100*y.mean(), A_only=n10, B_only=n01,
                   mcnemar_p=p_m, power_A=100*A.loc[c, "power"].mean(), power_B=100*B.loc[c, "power"].mean(), wilcoxon_p=p_w))
for sc in SC:
    paired(sc, "MK BH q=0.05", "v2 joint 0.025/sign", R)
    paired(sc, "MK BH q=0.10", "v2 joint 0.05/sign", R)
    paired(sc, "MK BH q=0.05", "v2 joint 0.005/sign (tuned)", R)
    paired(sc, "v2 joint 0.025/sign", "v2 joint 0.005/sign (tuned)", R)
P = pd.DataFrame(pr); P.round(4).to_csv(f"{OUT}/val_B_paired.csv", index=False)
print("=== paired ==="); print(P.round(4).to_string(index=False))
