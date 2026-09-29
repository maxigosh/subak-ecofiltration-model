"""Numbers and figures for model v2. Usage: python3 python/run_v2.py"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from model_v2 import *

out = Path(__file__).resolve().parent.parent / "figures"
out.mkdir(exist_ok=True)
kT = params()["k_T"]
print(f"Calibrated particle capture k_T = {kT:.3f} 1/h → turbidity {T_IN} → {T_OUT} NTU")

rows = []
plt.figure(figsize=(8, 4.5))
for name, kw in SCENARIOS.items():
    p = params(k_T=kT, **kw)
    d, t, o, s = simulate(p, 365)
    plt.plot(d, o, lw=2, label=name)
    rows.append((name, at(d, o, 20), at(d, o, 180), at(d, t, 180), s[-1]))
plt.axhline(O_LIMIT, color="k", ls="--", lw=1, label="Limit 10 mg/L (Permenkes 32/2017)")
plt.scatter([20, 20], [O_IN, O_OUT], color="red", zorder=5, label="Udayana test: in / out")
plt.xlabel("day"); plt.ylabel("organic matter (KMnO4), mg/L")
plt.title("Dissolved organics after the filter — scenarios")
plt.legend(fontsize=8); plt.tight_layout(); plt.savefig(out / "v2_organics_scenarios.png", dpi=130); plt.close()

V3s = np.array([1, 2, 4, 6, 8, 12, 16, 20, 24, 32])
vals = []
for V3 in V3s:
    d, t, o, s = simulate(params(k_T=kT, leach0=0.0, V3=float(V3)), 200)
    vals.append(at(d, o, 180))
taus = V3s / params()["Q"]
plt.figure(figsize=(7, 4))
plt.plot(taus, vals, "o-", lw=2)
plt.axhline(O_LIMIT, color="k", ls="--", lw=1)
plt.xlabel("contact time in the biochar stage, hours"); plt.ylabel("organics after filter, mg/L (day 180)")
plt.title("How long must water stay in the biochar stage?")
plt.tight_layout(); plt.savefig(out / "v2_contact_time.png", dpi=130); plt.close()

print("\nscenario | organics day 20 | organics day 180 | turbidity day 180 | biochar saturation day 365")
for r in rows:
    print(f"{r[0]} | {r[1]:.1f} | {r[2]:.1f} | {r[3]:.2f} | {100*r[4]:.0f}%")
need = next((tau for tau, v in zip(taus, vals) if v < O_LIMIT), None)
print("\ncontact time h → organics day 180:", ", ".join(f"{tau:g}h:{v:.1f}" for tau, v in zip(taus, vals)))
print("contact time needed for < 10 mg/L:", f"~{need:g} h" if need else "> range")
