"""Model v2 — calibrated on the Udayana water tests (16 Aug 2026).

Python mirror of src/SubakFilterV2.jl (same equations, same parameters), used
to produce the numbers and figures while Julia is being set up.

What the data say (irrigation water, before → after the pilot filter):
  turbidity         3.77 → 1.78 NTU   (−53%)   particles are removed well
  organic (KMnO4)   16.42 → 17.68 mg/L (+8%)   dissolved organics are NOT removed,
                                               and are above the 10 mg/L limit
So v2 separates two processes the v1 demo lumped together:
  * particulate (turbidity): straining/settling in every stage, first order;
  * dissolved organics: adsorption on biochar (Langmuir, linear driving force)
    + biodegradation by a biofilm that matures over weeks
    + leaching of organics from fresh (un-rinsed, un-charged) biochar
    + a refractory fraction nothing in the filter removes.

One before/after pair taken at the same moment cannot identify kinetics. The
"as-built" scenario is fitted to reproduce it; everything else is a scenario
to be tested with a time series of samples at the new pilot.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

T_IN, T_OUT = 3.77, 1.78       # NTU, measured
O_IN, O_OUT = 16.42, 17.68     # mg/L KMnO4, measured
O_LIMIT = 10.0                 # Permenkes 32/2017


def params(**kw):
    p = dict(
        Q=2.0,                       # m³/h through the unit (placeholder — measure on site)
        V1=4.0, V2=2.0, V3=1.0,      # water volume per stage, m³ (placeholder)
        k_T=None,                    # 1/h particle capture, calibrated below
        f_ref=0.35,                  # refractory share of KMnO4 organics
        M=500_000.0,                 # g biochar in stage 3
        q_max=15.0,                  # mg/g Langmuir capacity for dissolved organics
        K_L=0.15,                    # L/mg Langmuir affinity
        k_ads=0.004,                 # 1/h linear-driving-force rate (contact limited)
        k_bio_max=0.10,              # 1/h biodegradation in stages 2–3 when biofilm is mature
        t_bio=21 * 24.0,             # h, biofilm maturation time constant (~3 weeks)
        leach0=0.0189,               # g leachable organics per g fresh biochar
        k_leach=1 / (20 * 24.0),     # 1/h, leaching decay (~20 days)
        sample_day=20.0,             # bed age when the Udayana samples were taken (assumed)
    )
    p.update(kw)
    if p["k_T"] is None:
        p["k_T"] = calibrate_kT(p)
    return p


def calibrate_kT(p):
    """k_T such that steady-state turbidity removal equals the measured 53%."""
    taus = [p["V1"] / p["Q"], p["V2"] / p["Q"], p["V3"] / p["Q"]]
    target = T_OUT / T_IN
    f = lambda k: np.prod([1 / (1 + k * t) for t in taus]) - target
    return brentq(f, 1e-6, 100.0)


def rhs(t, u, p):
    T1, T2, T3, O1, O2, O3, q, L = u
    Q = p["Q"]
    d1, d2, d3 = Q / p["V1"], Q / p["V2"], Q / p["V3"]
    kb = p["k_bio_max"] * (1 - np.exp(-t / p["t_bio"]))
    O_ads_in = (1 - p["f_ref"]) * O_IN          # only the degradable/adsorbable part
    # particulate
    dT1 = d1 * (T_IN - T1) - p["k_T"] * T1
    dT2 = d2 * (T1 - T2) - p["k_T"] * T2
    dT3 = d3 * (T2 - T3) - p["k_T"] * T3
    # dissolved organics (adsorbable part), mg/L
    C3 = max(O3, 0.0)
    q_eq = p["q_max"] * p["K_L"] * C3 / (1 + p["K_L"] * C3)
    r_ads = p["k_ads"] * (q_eq - q)             # mg/g/h
    leach = p["k_leach"] * L                    # g/h released into stage 3
    dO1 = d1 * (O_ads_in - O1)
    dO2 = d2 * (O1 - O2) - kb * O2
    dO3 = d3 * (O2 - O3) - (p["M"] / (p["V3"] * 1000)) * r_ads - kb * O3 + leach / p["V3"]
    dL = -leach
    return [dT1, dT2, dT3, dO1, dO2, dO3, r_ads, dL]


def simulate(p, days=365):
    u0 = [T_IN, T_IN, T_IN, (1 - p["f_ref"]) * O_IN, (1 - p["f_ref"]) * O_IN, (1 - p["f_ref"]) * O_IN,
          0.0, p["leach0"] * p["M"]]
    t_eval = np.arange(0, days * 24 + 1, 6.0)
    sol = solve_ivp(rhs, (0, days * 24), u0, args=(p,), method="BDF", t_eval=t_eval, rtol=1e-6, atol=1e-8)
    days_ = sol.t / 24
    turb = sol.y[2]
    organic = sol.y[5] + p["f_ref"] * O_IN      # add back the refractory part
    sat = sol.y[6] / p["q_max"]
    return days_, turb, organic, sat


def at(days, series, d):
    return float(np.interp(d, days, series))


def first_day_below(days, series, limit, after=0.0):
    idx = np.where((series < limit) & (days >= after))[0]
    return float(days[idx[0]]) if len(idx) else None


def last_day_below(days, series, limit, after):
    """Day the organics climb back above the limit (biochar exhausted)."""
    below = (series < limit) & (days >= after)
    if not below.any():
        return None
    start = np.argmax(below)
    above = np.where(~below[start:])[0]
    return float(days[start + above[0]]) if len(above) else None


SCENARIOS = {
    "A. As built (fitted to Udayana test)": dict(),
    "B. Rinsed + charged biochar (no leaching)": dict(leach0=0.0),
    "C. B + 16 h contact time (slower flow or bigger bed)": dict(leach0=0.0, V3=32.0),
    "D. C + 2× biochar": dict(leach0=0.0, V3=32.0, M=1_000_000.0),
}
