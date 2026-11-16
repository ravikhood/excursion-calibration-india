"""Analysis of script 11 (MC convergence). Run: python3 analyse_mc.py <results dir> <out dir>
Writes mc_fwer.csv, mc_selection_bias.csv, mc_ab_agreement.csv, mc_fwer_curve_F.csv"""
import sys, os, numpy as np, pandas as pd
RES, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
s = pd.read_csv(f"{RES}/mc/11_sets.csv", keep_default_na=False, na_values=["NA"]); s["size"] = s["size"].astype(str)
F = pd.read_csv(f"{RES}/mc/11_F.csv"); F["size"] = F["size"].astype(str)
def wil(k, n, z=1.959964):
    p = k / n; d = 1 + z*z/n; c = (p + z*z/(2*n))/d; h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return 100*(c-h), 100*(c+h)
# FWER by batch/size/alpha (any flag, either sign)
g = s.groupby(["batch", "size", "alpha", "rep"])["any"].any().groupby(["batch", "size", "alpha"]).agg(["sum", "count"]).reset_index()
g["FWER"] = 100*g["sum"]/g["count"]; g[["lo", "hi"]] = [wil(k, n) for k, n in zip(g["sum"], g["count"])]
g.round(2).to_csv(f"{OUT}/mc_fwer.csv", index=False)
# selection bias: in-sample vs out-of-sample joint probability of non-empty sets
rows = []
for (sz, a), q in s[s["any"]].groupby(["size", "alpha"]):
    d = q.p_in - q.p_oos_mc
    rows.append(dict(size=sz, alpha=a, n_sets=len(q), target=1-a, p_in=q.p_in.mean(), p_oos_mc=q.p_oos_mc.mean(),
                     p_gauss=q.p_gauss.mean(), mean_drop=d.mean(), mean_se_oos=q.se_oos_mc.mean(),
                     mean_z=(d/q.se_oos_mc).mean(), share_oos_below_target=(q.p_oos_mc < 1-a).mean()))
pd.DataFrame(rows).round(4).to_csv(f"{OUT}/mc_selection_bias.csv", index=False)
# A vs B agreement (same fit, independent samples)
cs = lambda x: set() if x == "" or pd.isna(x) else set(map(int, str(x).split(";")))
s["cs"] = s.cells.apply(cs); rows = []
for (sz, a), q in s[s.batch.isin(["A", "B"])].groupby(["size", "alpha"]):
    A = q[q.batch == "A"].set_index(["rep", "sign"]).cs; B = q[q.batch == "B"].set_index(["rep", "sign"]).cs
    J = [len(A[k] & B[k]) / len(A[k] | B[k]) for k in A.index if A[k] or B[k]]
    D = [len(A[k] ^ B[k]) for k in A.index if A[k] or B[k]]
    ra = q[q.batch == "A"].groupby("rep")["any"].any(); rb = q[q.batch == "B"].groupby("rep")["any"].any()
    rows.append(dict(size=sz, alpha=a, cases=len(J), identical=sum(j == 1 for j in J), mean_jaccard=np.mean(J) if J else np.nan,
                     mean_cells_differing=np.mean(D) if D else np.nan, rep_flag_disagree=int((ra != rb).sum())))
pd.DataFrame(rows).round(3).to_csv(f"{OUT}/mc_ab_agreement.csv", index=False)
# FWER curve from excursion function (max F over cells, both signs)
mx = F.groupby(["batch", "size", "rep"]).F.max().fillna(0).reset_index(); rows = []
for (bt, sz), q in mx.groupby(["batch", "size"]):
    for a in [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.015, 0.02, 0.025, 0.05]:
        k = int((q.F >= 1 - a).sum()); lo, hi = wil(k, len(q))
        rows.append(dict(batch=bt, size=sz, alpha=a, nominal_two_sided=2*a, FWER=100*k/len(q), lo=lo, hi=hi))
pd.DataFrame(rows).round(2).to_csv(f"{OUT}/mc_fwer_curve_F.csv", index=False)
print("done")
