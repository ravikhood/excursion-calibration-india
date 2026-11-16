# Script 13 Part A: out-of-sample validation of a tuned alpha for v2 joint excursion sets.
# Convention: E_alpha = {F > 1-alpha} (strict), matching excursions.mc sets exactly (ties at F = 1-alpha occur with 4000 samples).
# Selection set: script 11 null reps 1-50 (v2, 4000 samples, batch A; batch B as sensitivity).
# Validation set: script 13 null reps 51-150 (v2, 4000 samples, independent permutations).
# Rule (fixed before looking at 51-150): alpha* = largest grid alpha per sign with selection FWER <= 5%.
import glob, sys, numpy as np, pandas as pd, rdata
ROOT = sys.argv[1] if len(sys.argv) > 1 else "results"
OUT = sys.argv[2] if len(sys.argv) > 2 else "tables"
GRID = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.025, 0.05]
def wilson(k, n, z=1.96):
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d; h = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return 100*(c-h), 100*(c+h)
def fwer(any_by_rep):
    k, n = int(any_by_rep.sum()), len(any_by_rep); lo, hi = wilson(k, n); return k, n, 100*k/n, lo, hi
# selection set
F = pd.read_csv(f"{ROOT}/mc/11_F.csv")
F = F[(F["model"] == "v2") & (F["size"] == 4000)]
W = F.pivot_table(index=["batch", "rep", "cell"], columns="sign", values="F").reset_index()
# validation set
V = pd.concat([pd.DataFrame(rdata.read_rds(f)) for f in sorted(glob.glob(f"{ROOT}/val/A_rep_*.rds"))])
assert V.rep.nunique() == 100 and V.groupby("rep").size().eq(274).all()
rows = []
for a in GRID:
    for b in ["A", "B"]:
        g = W[W.batch == b]; r = ((g.neg > 1-a) | (g.pos > 1-a)).groupby(g.rep).any()
        rows.append(("select_reps1-50_batch"+b, a, *fwer(r)))
    r = ((V.Fneg > 1-a) | (V.Fpos > 1-a)).groupby(V.rep).any()
    rows.append(("validate_reps51-150", a, *fwer(r)))
T = pd.DataFrame(rows, columns=["set", "alpha_per_sign", "k", "n", "FWER", "lo", "hi"])
T["nominal_two_sided"] = 2*T.alpha_per_sign
T.to_csv(f"{OUT}/val_A_fwer_curve.csv", index=False)
print(T.pivot(index="alpha_per_sign", columns="set", values="FWER"))
for b in ["A", "B"]:
    s = T[(T.set == "select_reps1-50_batch"+b) & (T.FWER <= 5)]
    a_star = s.alpha_per_sign.max()
    v = T[(T.set == "validate_reps51-150") & (T.alpha_per_sign == a_star)].iloc[0]
    print(f"batch {b}: alpha* = {a_star} -> validation FWER {v.FWER:.0f}% ({v.k}/{v.n}), 95% CI {v.lo:.1f}-{v.hi:.1f}")
# direct check at 0.025 and wrong-flag size
d = (V.jneg025 | V.jpos025).groupby(V.rep)
print("direct E(0.025) FWER:", 100*d.any().mean(), "| mean flagged cells:", d.sum().mean(), "| median when any:", d.sum()[d.sum()>0].median())
