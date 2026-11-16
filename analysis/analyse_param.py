"""Analysis of script 12 (parametric calibration). Run: python3 analyse_param.py <param dir> <out dir>"""
import sys, os, glob, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore"); import rdata
D, OUT = sys.argv[1], sys.argv[2]; os.makedirs(OUT, exist_ok=True)
cells, hyp = [], []
for f in sorted(glob.glob(f"{D}/rep_*_P*.rds")):
    x = rdata.read_rds(f); cells.append(pd.DataFrame(x["cells"])); hyp.append(pd.DataFrame(x["hyper"]))
C = pd.concat(cells, ignore_index=True); H = pd.concat(hyp, ignore_index=True)
C["rep"] = C["rep"].astype(int); H["rep"] = H["rep"].astype(int)
def wil(k, n, z=1.959964):
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d; h = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return 100*(c-h), 100*(c+h)
rows = []
for sc, g in C.groupby("scenario"):
    for a in [0.005, 0.01, 0.025, 0.05]:
        wrongneg = (g.Fneg >= 1-a) & (g.truth >= 0); wrongpos = (g.Fpos >= 1-a) & (g.truth <= 0)
        r = (wrongneg | wrongpos).groupby(g.rep).any(); k, n = int(r.sum()), len(r); lo, hi = wil(k, n)
        rows.append(dict(scenario=sc, alpha=a, nominal_two_sided=2*a, reps=n, error_rate=100*k/n, lo=lo, hi=hi))
    d = (((g.jneg025.astype(bool)) & (g.truth >= 0)) | ((g.jpos025.astype(bool)) & (g.truth <= 0))).groupby(g.rep).any()
    rows.append(dict(scenario=sc, alpha="0.025 direct E", nominal_two_sided=0.05, reps=len(d), error_rate=100*d.mean(), lo=np.nan, hi=np.nan))
R = pd.DataFrame(rows); R.round(2).to_csv(f"{OUT}/param_error_rates.csv", index=False)
cov = C.assign(cov=(C.q025 <= C.truth) & (C.truth <= C.q975)).groupby("scenario")["cov"].mean().mul(100)
H.to_csv(f"{OUT}/param_hyper.csv", index=False)
pd.set_option("display.width", 200)
print(R.round(1).to_string()); print("\nmarginal 95% coverage (%):"); print(cov.round(1))
print("\nhyperparameter estimates (median over reps):"); print(H.groupby("scenario")[["prec_eps","range_b","sd_b","range_w","sd_w"]].median().round(3))
