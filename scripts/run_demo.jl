# Runs the baseline scenario and a configuration comparison, and writes
# figures to ../figures. Usage (from the repo root):
#   julia --project=. -e 'using Pkg; Pkg.add(["OrdinaryDiffEq","Plots","Roots"])'
#   julia --project=. scripts/run_demo.jl
include(joinpath(@__DIR__, "..", "src", "SubakFilter.jl"))
using .SubakFilter
using Plots

p = default_params()
sol = simulate(p; days = 180)
days = sol.t ./ 24
Cin = [inflow(t, p)[3] for t in sol.t]
Cout = [u[5] for u in sol.u]
sat = [u[6] / p.q_max for u in sol.u]

mkpath(joinpath(@__DIR__, "..", "figures"))

plt1 = plot(days, Cin; label = "inflow", lw = 2, xlabel = "day", ylabel = "dissolved pollutant, mg/L",
    title = "Outflow vs inflow (storms = spikes)")
plot!(plt1, days, Cout; label = "after 3 stages", lw = 2)
savefig(plt1, joinpath(@__DIR__, "..", "figures", "outflow.png"))

plt2 = plot(days, 100 .* sat; lw = 2, label = "biochar saturation", xlabel = "day", ylabel = "% of capacity",
    title = "When does the biochar need replacing?", ylim = (0, 100))
savefig(plt2, joinpath(@__DIR__, "..", "figures", "saturation.png"))

# Configuration comparison: biochar mass vs. mean removal and breakthrough.
masses_kg = [100, 250, 500, 1000, 2000]
println("biochar_kg  mean_removal_%  breakthrough_day")
for m in masses_kg
    pm = default_params(M_biochar = m * 1000.0)
    s = simulate(pm; days = 180)
    eff = removal_efficiency(s, pm)
    bt = breakthrough_day(s, pm)
    println(rpad(m, 12), rpad(round(100 * sum(eff) / length(eff); digits = 1), 16), bt === nothing ? ">180" : round(bt; digits = 0))
end
