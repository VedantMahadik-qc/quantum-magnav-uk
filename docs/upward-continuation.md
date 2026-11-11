# Field at constant altitude (Phase 2b)

`scripts/upward_bgs.py`, region E 100–250 km, N 680–780 km (Isle of Mull and surroundings,
52,737 points, surveys CA63 and HG62, observation heights 305–1350 m above OD).

## Method

Equivalent sources (harmonica 0.7, `EquivalentSourcesGB`): a layer of point sources is fitted to
every data point **at its own observation height** (`aod`). The field is then predicted at any
altitude. This levels the mixed survey surfaces and continues the field upward in one step.

The fit uses 250 m block medians of the data, sources in 500 m blocks 1000 m below the data,
damping 1, and 80 km windows. Depth, damping and block sizes were picked on a separate set of
tuning lines (seed 100), excluding the evaluation lines below.

**Window size (changed 2026-11-11).** The first version used 20 km windows. Building the
national grids showed that the gradient-boosted fit had not converged at that size: two fits
of the same Mull data that differ only in `random_state` (the order windows are visited)
disagreed by **12% RMS at every altitude**, and by up to 176 nT RMS over the Mull igneous
centre. With 40 km windows the seed-to-seed difference is 2.7%, and with 80 km it is 2.0%. The
data misfit falls from 34 to 24 nT. A fit costs about 5 times more. All results below use
80 km windows; the 20 km values are given in brackets where they were published before.
Depth and damping were not retuned for the new window. Output grids are 250 m, masked to within 2 km of
data and at least 150 m above the ground.

## Results (pass criteria fixed before the first run)

**1. Recovering a known field: FAIL.** A synthetic harmonic field (300 point masses 0.5–5 km
below Ordnance Datum) is sampled at the real BGS positions and heights. The criterion is a
relative RMS error ≤ 5% at every altitude.

| Case | 1000 m AOD | 2000 m AOD | 4000 m AOD |
|---|---|---|---|
| Real BGS sampling, sources 0.5–5 km deep (**criterion**) | 12.7% (12.8) | 6.7% (6.9) | 4.7% (5.0) |
| Diagnostic: real sampling, sources 2.5–5 km deep | 2.2% | 1.7% | 2.4% |
| Diagnostic: gap-free sampling, sources 0.5–5 km deep | 1.6% | 1.8% | 3.1% |

The method recovers the field to 2–3% when the field is smooth, or when the sampling has no
gaps. It fails only when shallow sources meet 2 km line spacing. So the failure is a limit
of **the survey's line spacing, not of the method**. Short-wavelength anomalies between lines
were never measured, and no processing can recover them. Real magnetic sources reach the ground
surface, so for real geology this error at 1000 m is probably larger than the 13% here.
Excluding a 10 km border changes the numbers by at most a point, so edge effects are minor.

**2. Survey-level hold-out against the cubic grid: pass.** On the lines of fold 0 of
`holdout_bgs.py` (9 lines, 4,168 points, never used for tuning), equivalent sources predict
**117.3 nT RMS** (118.7 with 20 km windows) against 144.7 nT for the cubic grid. The criterion was ES ≤ 1.2× cubic. This box
is the steepest part of Britain: the field's own standard deviation there is about 280 nT.

**3. Navigation against altitude: pass.** A 60 km flight east from E 144 km, N 730 km. HG4930 INS,
2 nT sensor noise, truth and filter on the same map, 3 seeds. Prediction: the median final error
grows with altitude at every step.

| Map | RMS gradient | Final error (seeds 1–3) | Median |
|---|---|---|---|
| 305 m draped (cubic grid) | – | 6.6, 4.5, 5.8 m | 5.8 m |
| 1000 m AOD | 89 nT/km | 16.7, 27.3, 12.0 m | 16.7 m (16.5) |
| 2000 m AOD | 39 nT/km | 93.8, 41.1, 55.9 m | 55.9 m (74.1) |
| 4000 m AOD | 16 nT/km | 139.4, 352.4, 326.5 m | 326.5 m (310.3) |

Caveats: the 305 m row uses a different method (cubic, draped surface) from the others
(equivalent sources, constant altitude). Truth and filter share one map, so these numbers
measure how much information the map holds at each altitude, not the effect of map error.
Given result 1, the 1000 m map itself carries ≳ 13% field error that the filter never sees.

## Not done yet

National grids at altitude, which Phase 3 needs. At this setting the cost scales to roughly
10 minutes per altitude for all of Britain.
