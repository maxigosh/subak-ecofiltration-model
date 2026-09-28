"""
    SubakFilter

A first, deliberately simple dynamic model of a three-stage nature-based
filtration train for Subak irrigation water:

    inflow → [1] settling basin → [2] sand/gravel filter → [3] biochar bed → outflow

Each stage is treated as a well-mixed reactor (tanks-in-series). This is the
standard first step before moving to a spatially resolved (PDE) model: it
captures residence time, removal kinetics, media saturation and storm events
with few parameters, so it can be calibrated early against sparse field data.

All parameter values in `default_params()` are illustrative placeholders, not
measurements. They are meant to be replaced by lab/field data from Bali.
"""
module SubakFilter

using DifferentialEquations

export default_params, inflow, simulate, removal_efficiency, breakthrough_day

# ── Forcing: tropical flow with storm events ────────────────────────────────

"""
Storm pulses (time in hours). Each storm multiplies flow and brings a
"first flush" of sediment and dissolved pollutant.
"""
function storm_factor(t, p)
    s = 0.0
    for tc in p.storm_times_h
        s += exp(-((t - tc) / p.storm_width_h)^2)
    end
    return s
end

"""Inflow rate Q (m³/h) and inflow concentrations (mg/L) at time t (h)."""
function inflow(t, p)
    s = storm_factor(t, p)
    Q = p.Q_base * (1 + (p.storm_flow_mult - 1) * s)
    TSS_in = p.TSS_base * (1 + (p.storm_tss_mult - 1) * s)
    C_in = p.C_base * (1 + (p.storm_c_mult - 1) * s)
    return Q, TSS_in, C_in
end

# ── Model ───────────────────────────────────────────────────────────────────

"""
State vector u:
  u[1] S1  suspended solids in settling basin      (mg/L)
  u[2] C1  dissolved pollutant in settling basin   (mg/L)
  u[3] S2  suspended solids in sand filter          (mg/L)
  u[4] C2  dissolved pollutant in sand filter       (mg/L)
  u[5] C3  dissolved pollutant in biochar bed       (mg/L)
  u[6] q   pollutant loading on biochar             (mg/g)

Processes:
  settling          k_settle * S1          (= settling velocity / depth)
  filtration        k_filter * S2
  biodegradation    k_bio_i  * C_i          (biofilm, first order)
  adsorption        k_L * (q_eq(C3) - q)   (linear driving force)
  Langmuir isotherm q_eq(C) = q_max * K * C / (1 + K * C)

Units: volumes m³, flow m³/h, concentrations mg/L (= g/m³), biochar mass g.
"""
function rhs!(du, u, p, t)
    S1, C1, S2, C2, C3, q = u
    Q, TSS_in, C_in = inflow(t, p)

    d1 = Q / p.V1
    d2 = Q / p.V2
    d3 = Q / p.V3

    # Stage 1 — settling basin
    du[1] = d1 * (TSS_in - S1) - p.k_settle * S1
    du[2] = d1 * (C_in - C1) - p.k_bio1 * C1

    # Stage 2 — sand/gravel filter (particles strained, biofilm degrades)
    du[3] = d2 * (S1 - S2) - p.k_filter * S2
    du[4] = d2 * (C1 - C2) - p.k_bio2 * C2

    # Stage 3 — biochar bed (adsorption with saturation + biofilm)
    q_eq = p.q_max * p.K_L * max(C3, 0.0) / (1 + p.K_L * max(C3, 0.0))
    r_ads = p.k_ads * (q_eq - q)                 # mg/g/h onto the media
    du[5] = d3 * (C2 - C3) - (p.M_biochar / (p.V3 * 1000)) * r_ads - p.k_bio3 * C3
    du[6] = r_ads
    return nothing
end

"""Illustrative defaults. Replace with Bali lab/field values."""
default_params(; kw...) = merge((
    # hydraulics
    Q_base = 2.0,            # m³/h base irrigation flow through the unit
    V1 = 4.0, V2 = 2.0, V3 = 1.0,   # stage water volumes, m³
    # storms (hours from start) and their effect
    storm_times_h = [20, 55, 90, 125, 160] .* 24.0,
    storm_width_h = 6.0,
    storm_flow_mult = 4.0,
    storm_tss_mult = 6.0,
    storm_c_mult = 3.0,
    # inflow quality
    TSS_base = 50.0,         # mg/L suspended sediment
    C_base = 2.0,            # mg/L dissolved pollutant (e.g. phosphate)
    # kinetics
    k_settle = 0.30,         # 1/h
    k_filter = 0.80,         # 1/h
    k_bio1 = 0.005, k_bio2 = 0.02, k_bio3 = 0.01,   # 1/h
    # biochar
    M_biochar = 500_000.0,   # g (500 kg)
    q_max = 40.0,            # mg/g Langmuir capacity
    K_L = 1.0,               # L/mg Langmuir affinity
    k_ads = 0.05,            # 1/h mass-transfer rate
), (; kw...))

"""Run the model for `days` days. Returns the ODE solution (time in hours)."""
function simulate(p = default_params(); days = 180)
    u0 = [p.TSS_base, p.C_base, 0.0, 0.0, 0.0, 0.0]
    prob = ODEProblem(rhs!, u0, (0.0, days * 24.0), p)
    # Rodas5P: stiff solver — adsorption and storm pulses make the system stiff.
    solve(prob, Rodas5P(); saveat = 1.0, reltol = 1e-6, abstol = 1e-8)
end

"""Instantaneous removal of the dissolved pollutant, 0–1."""
function removal_efficiency(sol, p)
    [1 - sol.u[i][5] / inflow(sol.t[i], p)[3] for i in eachindex(sol.t)]
end

"""First day when outflow exceeds `limit` mg/L outside storms (biochar exhausted)."""
function breakthrough_day(sol, p; limit = 0.5 * p.C_base)
    for i in eachindex(sol.t)
        storm_factor(sol.t[i], p) < 0.01 && sol.u[i][5] > limit && return sol.t[i] / 24
    end
    return nothing
end

end # module
