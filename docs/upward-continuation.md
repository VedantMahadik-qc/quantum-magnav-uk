# Field at constant altitude (Phase 2b)

`scripts/upward_bgs.py`, region E 100–250 km, N 680–780 km (Isle of Mull and surroundings,
52,737 points, surveys CA63 and HG62, observation heights 305–1350 m above OD).

## Method

Equivalent sources (harmonica 0.7, `EquivalentSourcesGB`): a layer of point sources is fitted to
every data point **at its own observation height** (`aod`). The field is then predicted at any
altitude. This levels the mixed survey surfaces and continues the field upward in one step.

The fit uses 250 m block medians of the data, sources in 500 m blocks 1000 m below the data,
damping 1, and 20 km windows. These values were picked on a separate set of tuning lines (seed
100), excluding the evaluation lines below. Output grids are 250 m, masked to within 2 km of
data and at least 150 m above the ground.

## Results (pass criteria fixed before the first run)

**1. Recovering a known field: FAIL.** A synthetic harmonic field (300 point masses 0.5–5 km
below Ordnance Datum) is sampled at the real BGS positions and heights. The criterion is a
relative RMS error ≤ 5% at every altitude.

| Case | 1000 m AOD | 2000 m AOD | 4000 m AOD |
|---|---|---|---|
| Real BGS sampling, sources 0.5–5 km deep (**criterion**) | 12.8% | 6.9% | 5.0% |
| Diagnostic: real sampling, sources 2.5–5 km deep | 2.3% | 1.9% | 2.6% |
| Diagnostic: gap-free sampling, sources 0.5–5 km deep | 2.1% | 2.0% | 3.0% |

The method recovers the field to 2–3% when the field is smooth, or when the sampling has no
gaps. It fails only when shallow sources meet 2 km line spacing. So the failure is a limit
of **the survey's line spacing, not of the method**. Short-wavelength anomalies between lines
were never measured, and no processing can recover them. Real magnetic sources reach the ground
surface, so for real geology this error at 1000 m is probably larger than the 13% here.
Excluding a 10 km border changes the numbers by at most a point, so edge effects are minor.

**2. Survey-level hold-out against the cubic grid: pass.** On the lines of fold 0 of
`holdout_bgs.py` (9 lines, 4,168 points, never used for tuning), equivalent sources predict
**118.7 nT RMS** against 144.7 nT for the cubic grid. The criterion was ES ≤ 1.2× cubic. This box
is the steepest part of Britain: the field's own standard deviation there is about 280 nT.

**3. Navigation against altitude: pass.** A 60 km flight east from E 144 km, N 730 km. HG4930 INS,
2 nT sensor noise, truth and filter on the same map, 3 seeds. Prediction: the median final error
grows with altitude at every step.

| Map | RMS gradient | Final error (seeds 1–3) | Median |
|---|---|---|---|
| 305 m draped (cubic grid) | – | 6.6, 4.5, 5.8 m | 5.8 m |
| 1000 m AOD | 87 nT/km | 16.5, 29.8, 11.0 m | 16.5 m |
| 2000 m AOD | 39 nT/km | 103.5, 43.9, 74.1 m | 74.1 m |
| 4000 m AOD | 16 nT/km | 191.6, 409.3, 310.3 m | 310.3 m |

Caveats: the 305 m row uses a different method (cubic, draped surface) from the others
(equivalent sources, constant altitude). Truth and filter share one map, so these numbers
measure how much information the map holds at each altitude, not the effect of map error.
Given result 1, the 1000 m map itself carries ≳ 13% field error that the filter never sees.

## Not done yet

National grids at altitude, which Phase 3 needs. At this setting the cost scales to roughly
10 minutes per altitude for all of Britain.
