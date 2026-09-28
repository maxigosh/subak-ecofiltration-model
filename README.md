# Subak ecofiltration — first dynamic model (draft)

A first-pass computational model of a three-stage, nature-based filtration
train for Subak irrigation water in Bali:

```
inflow → [1] settling basin → [2] sand/gravel filter → [3] biochar bed → outflow
```

It is a **starting point for discussion**, built to show how the physical
system could be translated into a calibratable model. All parameter values are
**illustrative placeholders**, not measurements — the next step is to replace
them with lab and field data.

## What the model represents

Each stage is a well-mixed reactor (tanks-in-series), which keeps the number of
parameters small enough to calibrate against sparse early field data. A
spatially resolved (advection–dispersion PDE) version of the biochar bed is a
natural second step once breakthrough curves are measured.

| Process | Form |
|---|---|
| Flow & residence time | inflow Q(t) with storm pulses; residence time V/Q per stage |
| Sedimentation | first-order, `k_settle = v_s / h` |
| Filtration of particles | first-order straining, `k_filter` |
| Biological degradation | first-order biofilm decay per stage, `k_bio,i` |
| Adsorption onto biochar | linear-driving-force kinetics towards a Langmuir isotherm, `q_eq = q_max·K·C/(1+K·C)` |
| Media saturation | biochar loading `q(t)` tracked explicitly; breakthrough detected |
| Storm / high-flow events | flow ×4, sediment ×6, dissolved load ×3 (first flush), Gaussian pulses |

State variables: suspended solids and dissolved pollutant in stages 1–2,
dissolved pollutant in stage 3, and pollutant loading on the biochar.

## First results (illustrative parameters)

![Outflow vs inflow](figures/outflow.png)

![Biochar saturation](figures/saturation.png)

![Sediment removal](figures/sediment.png)

Configuration comparison — biochar mass in stage 3 vs mean removal over 180 days
and day of breakthrough (outflow above 50 % of the inflow baseline):

| Biochar, kg | Mean removal, % | Breakthrough |
|---|---|---|
| 100 | 14 | day 22 |
| 250 | 35 | day 59 |
| 500 | 64 | day 124 |
| 1000 | 85 | > 180 days |
| 2000 | 93 | > 180 days |

![Configurations](figures/configurations.png)

This is the kind of question the model is meant to answer with real data: *how
much media, in which order, replaced how often, for a target removal under
tropical flow?*

## Run it

Julia (main model, DifferentialEquations.jl with a stiff `Rodas5P` solver):

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/run_demo.jl
```

`python/mirror.py` is a line-by-line SciPy mirror of the same equations and
defaults, used to cross-check the Julia model; the figures above were rendered
with it.

## What would make it real

1. The field setup: stage volumes, media masses, flow rates, hydraulic paths.
2. Target pollutants (e.g. phosphate, nitrate, TSS, pesticides) — one state per pollutant.
3. Lab isotherms and batch kinetics for the actual biochar (q_max, K, k_ads).
4. Inflow time series and rainfall; outflow measurements for calibration.
5. Then: parameter estimation, global sensitivity analysis, and a PDE model of the biochar column.

## Next steps

- Separate pollutants (N, P, organics) and add nutrient uptake in the ecological polishing stage.
- Replace the biochar CSTR by a 1-D advection–dispersion–adsorption column (method of lines).
- Calibrate against the first field measurements; rank parameters by sensitivity.
- Compare modular configurations under a rainfall record for Bali.

— Max Igoshev
