# Model v2 for the Denpasar Barat channel (SITE_DENPASAR_BARAT.md). From the repo root:
#   julia --project=. -e 'using Pkg; Pkg.instantiate()'
#   julia --project=. scripts/site_denpasar.jl
#
# Bed: the 20 m straight reach, 0.71 m wide, water held at 0.30 m by a weir,
# porosity 0.4 → 1.7 m³ of pore water in 3 equal stages. Rinsed and charged
# biochar (no leaching), 250 kg in the last stage. Particle capture rate kept
# from the Udayana calibration. Values on day 180 (mature biofilm).
include(joinpath(@__DIR__, "..", "src", "SubakFilterV2.jl"))
using .SubakFilterV2

interp(x, y, x0) = y[searchsortedfirst(x, x0)]

const W, DEPTH, L, POROSITY = 0.71, 0.30, 20.0, 0.4
const VPORE = W * DEPTH * L * POROSITY
kT = params().k_T

cases = [
    ("whole channel, 0.1 m/s", 0.71 * 0.12 * 0.1 * 3600),
    ("side stream 3.3 L/min", 0.2),
    ("side stream 1.8 L/min", 0.106),
    ("side stream 1.0 L/min", 0.06),
]

println("pore volume = ", round(VPORE; digits = 2), " m³")
println("case | Q m³/h | contact h | organics day 180 mg/L | turbidity day 180 NTU")
for (name, Q) in cases
    v = VPORE / 3
    p = params(k_T = kT, Q = Q, V1 = v, V2 = v, V3 = v, M = 250_000.0, leach0 = 0.0)
    d, t, o, _ = simulate(p; days = 200)
    println(name, " | ", round(Q; digits = 3), " | ", round(VPORE / Q; digits = 1), " | ",
            round(interp(d, o, 180); digits = 1), " | ", round(interp(d, t, 180); digits = 2))
end
