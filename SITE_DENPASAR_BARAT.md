# Site: subak channel, Denpasar Barat

*Measured by Fen's team, 30 Sep 2026. Location: https://maps.app.goo.gl/rdKFXZH5pAL2Rhhn6*

| Measure | Value |
|---|---|
| Channel width | 0.70–0.72 m (0.71 used) |
| Channel height (wall) | 0.50 m |
| Current water level | 0.12 m |
| Straight length available | 20 m |

## What this gives the model

- Wetted cross-section: 0.71 × 0.12 = **0.085 m²**; water in the 20 m reach: **1.7 m³**.
- The flow rate is still unknown, and it decides everything. Contact time over the 20 m reach at typical channel speeds:

| Water speed | Flow | Contact time over 20 m |
|---|---|---|
| 0.05 m/s | 4.3 L/s (15 m³/h) | 6.7 min |
| 0.1 m/s | 8.5 L/s (31 m³/h) | 3.3 min |
| 0.2 m/s | 17 L/s (61 m³/h) | 1.7 min |
| 0.3 m/s | 26 L/s (92 m³/h) | 1.1 min |

Model v2 needs about **16 h** of contact to bring organics under 10 mg/L. Treating the whole channel flow in this reach gives minutes, not hours.

## Design consequence: treat a side stream

If the 20 m reach is filled with media (porosity ~0.4), the flow it can treat for 16 h is small:

| Water depth in the filter | Pore volume | Max flow for 16 h contact |
|---|---|---|
| 0.12 m (as now) | 0.68 m³ | 0.7 L/min |
| 0.30 m (small weir) | 1.7 m³ | 1.8 L/min |
| 0.40 m (weir near wall top) | 2.3 m³ | 2.4 L/min |

So the realistic pilot is a **side stream**: divert 1–2 L/min through the media bed (raised with a weir), let the rest of the channel pass.

## Model v2 run for this site

Model runs for this project are done in Julia only. The site run is `scripts/site_denpasar.jl` (same bed and assumptions as the side-stream table above: 20 m reach, 0.30 m weir, rinsed biochar, day 180). Results go here once it has run:

    julia --project=. scripts/site_denpasar.jl

| Flow through the filter | Contact time | Organics, mg/L | Turbidity, NTU |
|---|---|---|---|
| Whole channel (0.1 m/s) | — | pending Julia run | pending |
| Side stream 3.3 L/min | — | pending | pending |
| Side stream 1.8 L/min | — | pending | pending |
| Side stream 1.0 L/min | — | pending | pending |

## Still needed from the site

1. **Flow**: float method over the 20 m (time a floating object, 3 runs) → speed × 0.085 m² × 0.8.
2. Whether a weir or a bypass pipe is allowed in this channel.
3. Weekly in/out samples once the pilot runs (see FIELD_PLAN.md).
