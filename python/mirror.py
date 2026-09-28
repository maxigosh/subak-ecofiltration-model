"""Reference mirror of src/SubakFilter.jl in Python (SciPy), used to check the
Julia model and render figures where Julia is not installed. Same equations,
same default parameters."""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

def params(**kw):
    p = dict(Q_base=2.0, V1=4.0, V2=2.0, V3=1.0,
             storm_times_h=np.array([20, 55, 90, 125, 160]) * 24.0, storm_width_h=6.0,
             storm_flow_mult=4.0, storm_tss_mult=6.0, storm_c_mult=3.0,
             TSS_base=50.0, C_base=2.0, k_settle=0.30, k_filter=0.80,
             k_bio1=0.005, k_bio2=0.02, k_bio3=0.01,
             M_biochar=500_000.0, q_max=40.0, K_L=1.0, k_ads=0.05)
    p.update(kw); return p

def storm(t, p): return float(np.sum(np.exp(-((t - p["storm_times_h"]) / p["storm_width_h"]) ** 2)))

def inflow(t, p):
    s = storm(t, p)
    return (p["Q_base"] * (1 + (p["storm_flow_mult"] - 1) * s),
            p["TSS_base"] * (1 + (p["storm_tss_mult"] - 1) * s),
            p["C_base"] * (1 + (p["storm_c_mult"] - 1) * s))

def rhs(t, u, p):
    S1, C1, S2, C2, C3, q = u
    Q, TSS, Cin = inflow(t, p)
    d1, d2, d3 = Q / p["V1"], Q / p["V2"], Q / p["V3"]
    c = max(C3, 0.0)
    qeq = p["q_max"] * p["K_L"] * c / (1 + p["K_L"] * c)
    r = p["k_ads"] * (qeq - q)
    return [d1 * (TSS - S1) - p["k_settle"] * S1,
            d1 * (Cin - C1) - p["k_bio1"] * C1,
            d2 * (S1 - S2) - p["k_filter"] * S2,
            d2 * (C1 - C2) - p["k_bio2"] * C2,
            d3 * (C2 - C3) - p["M_biochar"] / (p["V3"] * 1000) * r - p["k_bio3"] * C3,
            r]

def simulate(p, days=180):
    t = np.arange(0, days * 24 + 1, 1.0)
    u0 = [p["TSS_base"], p["C_base"], 0, 0, 0, 0]
    return solve_ivp(rhs, (0, days * 24), u0, args=(p,), method="Radau", t_eval=t, rtol=1e-6, atol=1e-8, max_step=2.0)

def breakthrough(sol, p, limit=None):
    limit = 0.5 * p["C_base"] if limit is None else limit
    for i, t in enumerate(sol.t):
        if storm(t, p) < 0.01 and sol.y[4, i] > limit: return t / 24
    return None

if __name__ == "__main__":
    out = Path(__file__).resolve().parent.parent / "figures"; out.mkdir(exist_ok=True)
    p = params(); s = simulate(p); d = s.t / 24
    cin = np.array([inflow(t, p)[2] for t in s.t])
    tss_in = np.array([inflow(t, p)[1] for t in s.t])

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(d, cin, label="inflow", lw=1.6); ax.plot(d, s.y[4], label="after 3 stages", lw=1.8)
    ax.set_xlabel("day"); ax.set_ylabel("dissolved pollutant, mg/L"); ax.set_title("Outflow vs inflow (spikes = storms)")
    ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(out / "outflow.png", dpi=150)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(d, 100 * s.y[5] / p["q_max"], lw=2, color="tab:green")
    bt = breakthrough(s, p)
    if bt: ax.axvline(bt, ls="--", color="tab:red"); ax.text(bt + 2, 10, f"breakthrough ≈ day {bt:.0f}", color="tab:red")
    ax.set_ylim(0, 100); ax.set_xlabel("day"); ax.set_ylabel("% of Langmuir q_max")
    ax.set_title("Biochar saturation — when does it need replacing?"); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(out / "saturation.png", dpi=150)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.plot(d, tss_in, label="inflow", lw=1.4); ax.plot(d, s.y[0], label="after settling", lw=1.4); ax.plot(d, s.y[2], label="after sand filter", lw=1.6)
    ax.set_xlabel("day"); ax.set_ylabel("suspended sediment, mg/L"); ax.set_title("Sediment through the first two stages")
    ax.legend(); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(out / "sediment.png", dpi=150)

    rows = []
    for m in [100, 250, 500, 1000, 2000]:
        pm = params(M_biochar=m * 1000.0); sm = simulate(pm)
        cinm = np.array([inflow(t, pm)[2] for t in sm.t])
        eff = 100 * np.mean(1 - sm.y[4] / cinm); b = breakthrough(sm, pm)
        rows.append((m, eff, b))
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar([str(r[0]) for r in rows], [r[1] for r in rows], color="tab:blue")
    for i, r in enumerate(rows):
        ax.text(i, r[1] + 1, f"{r[1]:.0f}%\n{'>180 d' if r[2] is None else f'{r[2]:.0f} d'}", ha="center")
    ax.set_ylim(0, 110); ax.set_xlabel("biochar in stage 3, kg"); ax.set_ylabel("mean removal over 180 days, %")
    ax.set_title("Configuration comparison: removal and time to replacement"); fig.tight_layout(); fig.savefig(out / "configurations.png", dpi=150)
    print("biochar_kg mean_removal_% breakthrough_day")
    for r in rows: print(r[0], round(r[1], 1), r[2] if r[2] is None else round(r[2]))
