"""Global sensitivity of model v2 (Monte Carlo + Spearman rank correlation).

Which uncertain inputs change the answer most? → what to measure first at the
next pilot. Usage: python3 python/sensitivity_v2.py
"""
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from model_v2 import params, simulate, at, O_LIMIT

rng = np.random.default_rng(42)
N = 400
# name: (low, high, log-scale?, label, how to measure)
RANGES = {
    "Q":         (0.5, 8.0, True,  "Flow through the filter (m³/h)", "bucket-and-stopwatch or float method, 3× per visit"),
    "V3":        (0.5, 4.0, True,  "Water volume of the biochar stage (m³)", "measure bed length × width × water depth × porosity"),
    "k_bio_max": (0.02, 0.5, True, "Biodegradation rate of mature biofilm (1/h)", "weekly in/out KMnO4 or COD after week 4"),
    "t_bio":     (7*24, 60*24, True, "Biofilm maturation time (h)", "weekly in/out samples from day 0 of a new bed"),
    "f_ref":     (0.1, 0.6, False, "Refractory share of organics (–)", "lab: COD vs BOD5 ratio of inflow water"),
    "q_max":     (5.0, 50.0, True, "Biochar adsorption capacity (mg/g)", "jar test: biochar + canal water, 3 doses, 24 h"),
    "K_L":       (0.02, 1.0, True, "Biochar adsorption affinity (L/mg)", "same jar test"),
    "k_ads":     (0.001, 0.05, True, "Adsorption rate (1/h)", "jar test sampled at 1, 4, 24 h"),
    "leach0":    (0.0, 0.04, False, "Leachable organics in fresh biochar (g/g)", "rinse test: soak biochar, measure KMnO4 of rinse water"),
    "M":         (1e5, 2e6, True,  "Biochar mass in the stage (g)", "weigh or estimate volume × bulk density"),
}

def sample(lo, hi, log):
    return float(np.exp(rng.uniform(np.log(lo), np.log(hi)))) if log else float(rng.uniform(lo, hi))

X, y180, y20 = [], [], []
for i in range(N):
    draw = {k: sample(lo, hi, lg) for k, (lo, hi, lg, _, _) in RANGES.items()}
    p = params(**draw)
    d, t, o, s = simulate(p, 200)
    X.append([draw[k] for k in RANGES]); y180.append(at(d, o, 180)); y20.append(at(d, o, 20))
X = np.array(X); y180 = np.array(y180); y20 = np.array(y20)

def rank(y):
    rows = []
    for j, k in enumerate(RANGES):
        rho = spearmanr(X[:, j], y).correlation
        rows.append((k, rho))
    return sorted(rows, key=lambda r: -abs(r[1]))

r180, r20 = rank(y180), rank(y20)
print(f"{N} runs. Organics day 180: median {np.median(y180):.1f}, 90% range {np.percentile(y180,5):.1f}–{np.percentile(y180,95):.1f} mg/L; "
      f"below 10 in {100*np.mean(y180<O_LIMIT):.0f}% of runs")
print("\nLong term (day 180) — Spearman rho:")
for k, r in r180: print(f"  {k:10s} {r:+.2f}  {RANGES[k][3]}")
print("\nStart-up (day 20) — Spearman rho:")
for k, r in r20: print(f"  {k:10s} {r:+.2f}  {RANGES[k][3]}")

out = Path(__file__).resolve().parent.parent / "figures"
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)
for a, rows, title in ((ax[0], r180, "Long term (day 180)"), (ax[1], r20, "Start-up (day 20)")):
    names = [RANGES[k][3] for k, _ in rows][::-1]; vals = [r for _, r in rows][::-1]
    a.barh(names, vals, color=["#c0392b" if v > 0 else "#2471a3" for v in vals])
    a.axvline(0, color="k", lw=0.8); a.set_xlim(-1, 1); a.set_title(title); a.set_xlabel("effect on organics after filter (Spearman ρ)")
    a.tick_params(axis="y", labelsize=8)
plt.suptitle("What drives the result? Red = raises organics, blue = lowers them")
plt.tight_layout(); plt.savefig(out / "v2_sensitivity.png", dpi=130); plt.close()

import json
json.dump({"N": N, "day180": r180, "day20": r20,
           "median180": float(np.median(y180)), "p5": float(np.percentile(y180,5)), "p95": float(np.percentile(y180,95)),
           "share_below_limit": float(np.mean(y180<O_LIMIT))},
          open(out.parent / "python" / "sensitivity_v2.json", "w"), indent=1)
