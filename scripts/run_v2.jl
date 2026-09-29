# Model v2 scenarios and figures. From the repo root:
#   julia --project=. -e 'using Pkg; Pkg.instantiate()'
#   GKSwstype=100 julia --project=. scripts/run_v2.jl
include(joinpath(@__DIR__, "..", "src", "SubakFilterV2.jl"))
using .SubakFilterV2
using Plots

interp(x, y, x0) = y[searchsortedfirst(x, x0)]
figs = joinpath(@__DIR__, "..", "figures"); mkpath(figs)

plt = plot(xlabel = "day", ylabel = "organic matter (KMnO4), mg/L", title = "Dissolved organics after the filter", legend = :topright)
println("scenario | organics day 20 | organics day 180 | turbidity day 180 | saturation day 365")
for (name, p) in scenarios()
    d, t, o, s = simulate(p; days = 365)
    plot!(plt, d, o; lw = 2, label = name)
    println(name, " | ", round(interp(d, o, 20); digits = 1), " | ", round(interp(d, o, 180); digits = 1),
            " | ", round(interp(d, t, 180); digits = 2), " | ", round(Int, 100 * s[end]), "%")
end
hline!(plt, [O_LIMIT]; ls = :dash, color = :black, label = "limit 10 mg/L")
scatter!(plt, [20, 20], [O_IN, O_OUT]; color = :red, label = "Udayana test in/out")
savefig(plt, joinpath(figs, "v2_organics_scenarios_julia.png"))
