"""
    SubakFilterV2

Model v2, calibrated on the Udayana water tests of 16 Aug 2026
(irrigation water before/after the pilot filter):

    turbidity        3.77 → 1.78 NTU   (−53%)  particles removed well
    organic (KMnO4)  16.42 → 17.68 mg/L (+8%)  dissolved organics not removed,
                                               above the 10 mg/L limit

Two process groups, three stages in series (tanks-in-series):
  * particulate (turbidity): first-order capture k_T in every stage,
    k_T calibrated so steady-state removal equals the measured 53%;
  * dissolved organics: Langmuir adsorption on biochar (linear driving force),
    biodegradation by a biofilm that matures over weeks, leaching from fresh
    (un-rinsed / un-charged) biochar, and a refractory share nothing removes.

A single before/after pair cannot identify kinetics: scenario "A" is fitted to
reproduce it, the rest are scenarios to test with a time series at the new
pilot. Python mirror: python/model_v2.py (same equations and parameters).
"""
module SubakFilterV2

using OrdinaryDiffEq
using Roots: find_zero, Bisection

export T_IN, T_OUT, O_IN, O_OUT, O_LIMIT, params, simulate, scenarios

const T_IN, T_OUT = 3.77, 1.78      # NTU, measured
const O_IN, O_OUT = 16.42, 17.68    # mg/L KMnO4, measured
const O_LIMIT = 10.0                # Permenkes 32/2017

function calibrate_kT(Q, V1, V2, V3)
    taus = (V1 / Q, V2 / Q, V3 / Q)
    f(k) = prod(1 / (1 + k * t) for t in taus) - T_OUT / T_IN
    find_zero(f, (1e-6, 100.0), Bisection())
end

"""Parameters (placeholders marked). `k_T = nothing` → calibrated from turbidity."""
function params(; Q = 2.0, V1 = 4.0, V2 = 2.0, V3 = 1.0, k_T = nothing,
                f_ref = 0.35, M = 500_000.0, q_max = 15.0, K_L = 0.15, k_ads = 0.004,
                k_bio_max = 0.10, t_bio = 21 * 24.0, leach0 = 0.0189, k_leach = 1 / (20 * 24.0))
    kT = k_T === nothing ? calibrate_kT(Q, V1, V2, V3) : k_T
    (; Q, V1, V2, V3, k_T = kT, f_ref, M, q_max, K_L, k_ads, k_bio_max, t_bio, leach0, k_leach)
end

function rhs!(du, u, p, t)
    T1, T2, T3, O1, O2, O3, q, L = u
    d1, d2, d3 = p.Q / p.V1, p.Q / p.V2, p.Q / p.V3
    kb = p.k_bio_max * (1 - exp(-t / p.t_bio))
    O_ads_in = (1 - p.f_ref) * O_IN
    du[1] = d1 * (T_IN - T1) - p.k_T * T1
    du[2] = d2 * (T1 - T2) - p.k_T * T2
    du[3] = d3 * (T2 - T3) - p.k_T * T3
    C3 = max(O3, 0.0)
    q_eq = p.q_max * p.K_L * C3 / (1 + p.K_L * C3)
    r_ads = p.k_ads * (q_eq - q)
    leach = p.k_leach * L
    du[4] = d1 * (O_ads_in - O1)
    du[5] = d2 * (O1 - O2) - kb * O2
    du[6] = d3 * (O2 - O3) - (p.M / (p.V3 * 1000)) * r_ads - kb * O3 + leach / p.V3
    du[7] = r_ads
    du[8] = -leach
    nothing
end

"""Run `days` days. Returns (days, turbidity, organics incl. refractory, biochar saturation)."""
function simulate(p; days = 365)
    Oa = (1 - p.f_ref) * O_IN
    u0 = [T_IN, T_IN, T_IN, Oa, Oa, Oa, 0.0, p.leach0 * p.M]
    sol = solve(ODEProblem(rhs!, u0, (0.0, days * 24.0), p), Rodas5P(); saveat = 6.0, reltol = 1e-6, abstol = 1e-8)
    d = sol.t ./ 24
    turb = [u[3] for u in sol.u]
    org = [u[6] + p.f_ref * O_IN for u in sol.u]
    sat = [u[7] / p.q_max for u in sol.u]
    d, turb, org, sat
end

"""Scenarios A–D (k_T fixed from the as-built calibration)."""
function scenarios()
    kT = params().k_T
    [
        "A. As built (fitted to Udayana test)" => params(k_T = kT),
        "B. Rinsed + charged biochar (no leaching)" => params(k_T = kT, leach0 = 0.0),
        "C. B + 16 h contact time (slower flow or bigger bed)" => params(k_T = kT, leach0 = 0.0, V3 = 32.0),
        "D. C + 2× biochar" => params(k_T = kT, leach0 = 0.0, V3 = 32.0, M = 1_000_000.0),
    ]
end

end # module
